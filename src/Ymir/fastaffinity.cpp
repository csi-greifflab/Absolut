#include "fastaffinity.h"
#include "proteins.h"
#include "receptorligand.h"
#include <fstream>
#include <sstream>
#include <cmath>

#include "../Tools/distribution.h"
#include "../Tools/md5.h"

#include "plot3d.h"

#include "pthread.h"
extern pthread_mutex_t lockAccessPrecompAffinities = PTHREAD_MUTEX_INITIALIZER;




// Philippe 2019-10-5 todo:
// 1- hash memory by start position
// 2- the only important is the occupied positions from the ligand. Change the recursion to use 'occupied positions'
// 3- make a multistruct class with sequences and structures + numbering scheme
// 4- Change the definition of compact, with position-to-position
// 5- Find the optimal sequence
// 6- find where Latfit expects a missing residue

affinityOneLigand::affinityOneLigand(string _ligandStructureSeq, string _ligandAAseq, int startPosition, int _sizeReceptors, int _minimalNInteract, int _minSelfFoldings, double _KT, vector<int> listForbiddenPositions)
{
    if(_ligandAAseq.size() != _ligandStructureSeq.size() + 1){
        cerr << "ERR: affinityOneLigand, incompatible structure " << _ligandStructureSeq << " with AA sequence (size should be +1) " << _ligandAAseq << endl;
        return;
    }
    //ligandSeq = _ligandStructureSeq;
    //ligandAAseq = _ligandAAseq;
    struct3D _ligand = struct3D(_ligandStructureSeq, UnDefined, startPosition);
    //cout << print(_ligand) << endl;
    superProtein* _ligandProt = new superProtein(_ligand);
    //cout << print(*_ligandProt) << endl;
    _ligandProt->setAAs(_ligandAAseq);
    //cout << print(*_ligandProt) << endl;
    initialize(_ligandProt, _sizeReceptors, _minimalNInteract, _minSelfFoldings, _KT, listForbiddenPositions);
}

affinityOneLigand::affinityOneLigand(superProtein* _ligand, int _sizeReceptors, int _minimalNInteract, int _minSelfFoldings, double _KT, vector<int> listForbiddenPositions){
    initialize(_ligand, _sizeReceptors, _minimalNInteract, _minSelfFoldings, _KT, listForbiddenPositions);
}

void affinityOneLigand::initialize(superProtein* _ligand, int _sizeReceptors, int _minimalNInteract, int _minSelfFoldings, double _KT, vector<int> listForbiddenPositions){

    // reading or regenerating the binding structures

    // 3 files are generated each time a new ligand has to be precomputed:
    //          - fStruct:  list of structures touching this ligand
    //          - fAll:     list of structures touching + list of interactions (AA in the ligand - position in the receptor)
    //          - fCompact: only list of interactions (AA ligand - position receptor) with the number of times they are repeated.
    // When creating a new AffinityOneLigand class,
    //          - if fCompact exist for this ligand (structure AND AA sequence), nothing to do,
    //          - if fStruct exist for this ligand (structure) but not fCompact (other AA sequence than previously computed ligand with same structure)
    //                  then reads fStruct and regenerates fCompact for this new AA sequence
    //          - if none exist, then generates all files from scratch (takes quite some time ...)
    //          Note, fAll is a readout to see concretely the structures, not needed for affinities.

    modeUltraFast = false;
    ligand = _ligand;
    sizeReceptors = _sizeReceptors;
    minimalNInteract = _minimalNInteract;
    minNrSelfInteractions = _minSelfFoldings;
    KT = _KT;

    // Generates the file names for this ligand and options
    if(!ligand->structure) {cerr << "NULL ligand structure" << endl; return;}
    ligandSeq = ligand->structure->sequence;
    ligandAAseq = ligand->getAAseq();

    cout << "====== Initialization of receptor structures ====== \n";
    cout << "   ->  Requesting all possible receptor structures for ligand:\n       Struct=" << ligandSeq << ", start pos " << ligand->structure->startingPosition
         << ",\n       AAs   =" << ligandAAseq << ", with receptors of size " << sizeReceptors << ",\n       with min " << minimalNInteract << " interactions and pooled with self-folding with at least " << _minSelfFoldings << " self-interactions,\n       and forbidden positions:" << codePos(listForbiddenPositions) << endl;

    fStruct =    fnameStructures(ligand, sizeReceptors, minimalNInteract, listForbiddenPositions);
    fAll =       fnameStructuresAndCompactForAASeqLigand(ligand, sizeReceptors, minimalNInteract, listForbiddenPositions);
    fCompact =   fileNameCompactForAASeqLigand(ligand, sizeReceptors, minimalNInteract, listForbiddenPositions);



    // ============= Part 1: try to read existing file with list of structures for this antigen structure. =========

    // Does fCompact exist ?: look for the file with all structures and interactions for this antigen with this AA sequence
    ifstream fread0(fCompact.c_str());

    // If no, does fStruct exist ?: look for the file with all structures (but no interactions) for this antigen structure, but not knowing the AA sequence of the antigen
    if(!fread0){
        ifstream fread007(fStruct.c_str());
        // if yes, regenerates the fCompact and fALL files
        if(fread007){
            fread007.close();
            cout << "   ... The structures for the same ligand structure but different sequence are already known.\n"
                 << "       Regenerating the interactions with this new AA sequence" << endl;
            reGenerateCompressedStructures(fStruct, ligandAAseq, fAll, fCompact, ligand);
        } else { // try as well in the parent folder.
            ifstream fread008(string("../") + fStruct.c_str());
            // if yes, regenerates the fCompact and fALL files
            if(fread008){
                fread008.close();
                cout << "   ... The structures for the same ligand structure but different sequence are already known (Parent Folder).\n"
                     << "       Regenerating the interactions with this new AA sequence" << endl;
                reGenerateCompressedStructures(string("../") + fStruct, ligandAAseq, fAll, fCompact, ligand);
            }
        }
    } else {
        fread0.close();
    }

    // ============= Part 2: If no file was found, regenerates all receptor structures around the ligand . =========

    // Now, if fCompact still doesn"t exist, it means none of the precompted files exist, and everything should be recomputed
    ifstream fread1(fCompact.c_str());
    if(!fread1){

        cout << "   ... Has never been computed before. Will be stored once for all inside files:\n"
             << "       " << fStruct << " for structures \n"
             << "       " << fAll << " for structures and interaction profiles\n"
             << "       " << fCompact << " for interaction profiles \n";
        set<int>* merged = generateForbidden(listForbiddenPositions);

        receptorLigand* a = new receptorLigand(*ligand->structure, sizeReceptors, minimalNInteract, *merged);
        set<int>* forbid = generateForbidden(listForbiddenPositions);
        a->setForbiddenVolume(*forbid);
        delete forbid;

        // ================================================= Here is the generation of receptors ==================================================

        // This function generates all possible structures. It stores them inside the field vector<struct*> possibleReceptors
        a->generateReceptors();

        // ========================================================================================================================================

        cout << "Receptor structures generated for L= " << sizeReceptors << "\tminI=" << minimalNInteract << "\t --> " << a->possibleReceptors.size()  << " structures " << endl;
        //stringstream fname; fname << ligand->sequence << "-" << ligandAAseq << "-" << sizeReceptors << "-" << minimalNInteract << codePos(listForbiddenPositions);

        // not necessary anymore : the struct3D already has the AA sequence
        a->putSequenceLigand(ligandAAseq);


        // ================================================= Here is the generation of receptors ==================================================

        // This function saves the result into the three files files .

        a->printToFile(fStruct);
        int nbNonredundantCodes = reGenerateCompressedStructures(fStruct, ligandAAseq, fAll, fCompact, ligand);
        cout << "Got non-redundant interactions: " <<nbNonredundantCodes << endl;

        // ========================================================================================================================================

        size_t NRs = a->possibleReceptors.size();
        vector<int> classesNbInteract(50, 0);
        for(size_t i = 0; i < NRs; ++i){
            classesNbInteract[nbTouchPoints(*(a->possibleReceptors[i]), *ligand->structure)]++;
            delete a->possibleReceptors[i];
        }
        int sum = 0;
        for(size_t i = 0; i < classesNbInteract.size(); ++i){
            if( classesNbInteract[i] > 0){
                cout << "   " << i << " Interact, \t" << classesNbInteract[i] << endl;
                sum += classesNbInteract[i];
            }
        }
        delete a;
        //delete ligand;
        delete merged;

        if(sum != NRs) cerr << "ERROR, size inconsistency " << endl;
    } else {
        fread1.close();
    }
    // now all three files should exist => The data is saved, the program can crash it's no problem
    // Note: the field 'possibleReceptors' is now deleted of all its members (they take too much memory) and instead everything is saved as strings into the text files.





    // ===================== Part 3: Loading the receptor structures and interaction profiles from the text files =========================

    // Loading the interaction codes from the compact file
    ifstream fread(fCompact.c_str());
    string useless;
    if(fread){
        interactions = new vector<string>;
        cout << "   ... Structures have been computed - reading from file " << fCompact << endl;
        fread >> nInterCodes;
        fread >> useless;
        fread >> sizeReceptors;     // idem
        interactions->resize(nInterCodes);
        nbRepeats.resize(nInterCodes);
        for(int i = 0; i < nInterCodes; ++i){
            fread >> (*interactions)[i] >> nbRepeats[i];
        }
        cout << "       loaded " << nInterCodes << " non redundant possible interaction profiles" << endl;
        fread.close();
    } else {
        cerr << "ERR: affinityOneLigand, could not read the generated file for receptor structures (" << fCompact << ") Maybe problem with folders ?" << endl;
        return;
    }
    if(nInterCodes > 1e12) cerr << "ERR: read " << nInterCodes << " Interaction codes when initializing affinityOneLigand class. Seems to be too much" << endl;



//    // ===================== Part 4: Self-foldings of the receptor: Reading from file or regenerating all self-foldings with this size. This is independent of the ligand.

    // reading or regenerating self-foldings.
    // reading or regenerating the binding structures
    nFoldingCodes = 0;
    if(minNrSelfInteractions == -1) minNrSelfInteractions = minimalNInteract;
    if(minNrSelfInteractions < 1) cerr << "ERR: forbidden to call loadOrGenerateSelfFoldings with less than one interaction" << endl;

    stringstream fnameb; fnameb << "selfFoldingsL=" << sizeReceptors << "minI=" << minimalNInteract << "Compact.txt";
    string fname = fnameb.str();
    ifstream fread2(fname.c_str());
    if(!fread2){
        cout << "   ... The self-foldings for L=" << sizeReceptors << " and minI=" << minimalNInteract << " are being generated " << endl;
        vector<struct3D*> selfFoldingsForThisSize = generateSelfFoldings(sizeReceptors,minNrSelfInteractions);
        //cout << "Generated " << selfFoldingsForThisSize.size() << " self-foldings" << endl;
        stringstream fname2; fname2 << "selfFoldingsL=" << sizeReceptors << "minI=" << minimalNInteract << ".txt"; // without compact
        cout << "       saving self-foldings in " << fname2.str() << " and associated compact file " << endl;
        int newNS = exportSelfInteractions(selfFoldingsForThisSize, fname2.str());
        for(size_t i = 0; i < selfFoldingsForThisSize.size(); ++i){
            if(selfFoldingsForThisSize[i]) delete selfFoldingsForThisSize[i];
        }
        cout << "       saved " << newNS << " non-redundant structures (sometimes there is none if minInteract is high)" << endl;
    } else {
        fread2.close();
    }
    ifstream fread3(fname.c_str());
    if(fread3){
        selfInteractions = new vector<string>;
        cout << "   ... Structures have been computed - reading from file " << fname << endl;
        fread3 >> nFoldingCodes;
        selfInteractions->resize(nFoldingCodes);
        nbSelfRepeats.resize(nFoldingCodes);
        for(int i = 0; i < nFoldingCodes; ++i){
            fread3 >> (*selfInteractions)[i] >> nbSelfRepeats[i];
        }
        cout << "       loaded " << nFoldingCodes << " non redundant possible self-folding interaction profiles" << endl;
        fread3.close();
    } else {
        cerr << "ERR: affinityOneLigand, could not read the generated file for self-foldings. Maybe problem with folders ?" << endl;
        return;
    }
    if(nFoldingCodes > 1e10) cerr << "ERR: when initializing affinityOneLigand, got " << nFoldingCodes << " folding codes. This seems too much." << endl;
    // names without 'Compact' for the full details of the structures
    //stringstream fnamec; fnamec << ligandSeq << "-" << ligandAAseq << "-" <<sizeReceptors << "-" << minimalNInteract << codePos(listForbiddenPositions) << ".txt";
    stringstream fnamed; fnamed << "selfFoldingsL=" << sizeReceptors << "minI=" << minimalNInteract << codePos(listForbiddenPositions) << ".txt";
    fileStructures = fAll;
    fileSelfFoldings = fnamed.str();


    // ============= Part 5 (optional): Generates a dictionary between interaction profile and structure to retrieve later all structures with the best interaction profile ===================

//#ifdef showBestStructures
    cout << "   ... (optional:) Creating a dictionary between compacted structures to original structure [only for debugging]" << endl;
    cout << "       Reading file " << fileStructures << endl;

    // Generates a dictionnary mapping each interaction profile to the list of structures that had this interaction code.
    profileToStructure.clear();
    ifstream f(fileStructures.c_str());
    if(!f){
        cerr << "ERR: Can not open " << fileStructures << endl;
        exit(-1);
    } else {
        string useless;
        int nForbidden = 0;
        f >> useless >> useless >> useless;
        f >> useless >> nForbidden;
        for(int i = 0; i < nForbidden; ++i){
            f >> useless;
        }
        f >> useless >> useless >> useless >> useless;
        f >> useless >> useless >> useless >> useless >> useless;

        f >> useless; // nStructures, don't need, will just do a while

        int startPos;
        string structure;
        string interactions;

        if(!interactions.compare(string("fLhLhViGaKiKcPiPdChCeYbgbidg"))) {cerr << "FOUND!" ; exit(-1);}

        while((f >> startPos)){
            f >> structure;
            f >> interactions;
            std::map<string, vector<std::pair<int,string> > >::iterator it = profileToStructure.find(interactions);

            if(it != profileToStructure.end()){
                //toAdd = *it;
                //toAdd.push_back(structure);
                //profileToStructure[interactions] = toAdd;
                it->second.push_back(std::pair<int, string>(startPos,structure));
            } else {
                vector<std::pair<int,string> > toAdd;
                toAdd.push_back(std::pair<int, string>(startPos,structure));
                profileToStructure[interactions] = toAdd;
            }
        }
        cout << "       Dictionnary size " << profileToStructure.size() << endl;
    }

    cout << "   === Structures ready! ===\n" << endl;
//#endif
}









// Predefined, for 1 ligand
//affinityOneLigand::affinityOneLigand(string fname); // loads ligand sequence and all interactions
//vector<string>* interactions; //storage from the file => put to Multivector
// map<string, double> knownAffinities;

int code (char c1, char c2){
    return 256*(c1 - 'a') + (c2 - 'a');
}


// boolean comparators for sorting structures
bool compBind(std::pair<double, std::pair<double, string> > a, std::pair<double, std::pair<double, string> > b){
    return (a.first < b.first);
}
bool compTot(std::pair<double, std::pair<double, string> > a, std::pair<double, std::pair<double, string> > b){
    return (a.second.first < b.second.first);
}
bool compNegTot(std::pair<double, std::pair<double, string> > a, std::pair<double, std::pair<double, string> > b){
    return (a.second.first > b.second.first);
}

int nrBindings(string interactionCode){
    int res = 0;
    int nSelf = 0;
    size_t ICS = interactionCode.size();
    if(((ICS / 2) * 2) != ICS) cerr << "ERR: nrBindings, non-even interaction code " << interactionCode;
    for(size_t j = 0; j < ICS; j = j + 2){
        if(interactionCode[j+1] < 'a'){
            res++;
        } else {
            nSelf++;
        }
    }
    return res;
}

// tool function: the list of optimal interaction codes will be stored inside one string, with two spaces "  " as delimiter. This function is done to convert it back into a vector of string
vector<string> slice(string toSlice, string delimiter){
    size_t pos = 0;
    string token;
    vector<string> res;
    while ((pos = toSlice.find(delimiter)) != string::npos) {
        token = toSlice.substr(0, pos);
        res.push_back(token);
        toSlice.erase(0, pos + delimiter.length());
    }
    if(toSlice.size() > 0) res.push_back(toSlice);
    return res;
}

void testSLice(){
    cout << "|" << printVector(slice("ABCD  EFGHI  JKL", "  "));
    cout << "|" <<  printVector(slice("ABCD   JKL", "  "));
    cout << "|" <<  printVector(slice("ABCD", "  "));
    cout << "|" <<  printVector(slice("", "  "));
    cout << "|" <<  printVector(slice("      ", "  ")); // in this case it performs badly, don't understand why, but not important for us
    cout << "|" <<  printVector(slice("ABC  ", "  "));
}



// ===================================    one function to rule them all !  ======================================

// This function takes a receptor sequence. It browses all possible interaction profiles between ligand and receptor, fills them
// one by one with the receptor sequence, getting a binding and total energy. The structure with best total energy is kept (best energy),
// and a boltzmann average of the binding energies of each interaction profile (weighted by total energies) returns a statistical affinity.

// Option: when showing all structures, will do with the best ones first. Alternately, can show it without sorting (put false)
#define showSorted false
std::pair<double, double> affinityOneLigand::affinity(string receptorAASeq, bool showStructures, vector<string>* returnBestStructures){

    //if(modeUltraFast) cout << "Mode Ultrafast is ON" << endl;
    if(modeUltraFast) showStructures = false;

    // if wrong size of sequence
    if(static_cast<int>(receptorAASeq.size()) != sizeReceptors+1){
        cerr << "ERR: affinityOneLigand::affinity(" << receptorAASeq << "), wrong size of AA sequence ! the interactions profiles have only been loaded for receptors of size " << sizeReceptors + 1 << " AAs (i.e. with -1 number of moves)" << endl;
        return std::pair<double, double>(NAN, NAN);
    }


    // --------- Part 1: During the lifetime of the program, it stores all affinities it has computed so far in memory, to avoid recomputing.

    // if already computed before
    pthread_mutex_lock(&lockAccessPrecompAffinities);

    std::map<string,double>::iterator it = (knownBestAffinities.find(receptorAASeq));
    std::map<string,double>::iterator itEnd = knownBestAffinities.end();

    pthread_mutex_unlock(&lockAccessPrecompAffinities);
    if (it != itEnd){ // found

        pthread_mutex_lock(&lockAccessPrecompAffinities);
        std::map<string, string>::iterator it3 = (knownBestInteractions.find(receptorAASeq));
        if(it3 == knownBestInteractions.end()) cerr << "ERR: affinity(), the best affinity for sequence " << receptorAASeq << " is known but not the best interaction codes => Will return empty interaction codes, might cause problem" << endl;
        else {
            // BAAAAAAAAAAAAAAAAAAAAD, this is a copy paste...
            // --------- Part Optional: If requests the list of optimal structures, retrives it.
            if(returnBestStructures != nullptr){
                returnBestStructures->clear();
                vector<string> tokens = slice(it3->second, "  ");
                size_t nCodes = tokens.size();
                for(size_t k = 0; k < nCodes; ++k){
                    string token = tokens[k];
                    if(showStructures) cout << token;
                    vector<std::pair<int,string> > structuresForThisProfile = profileToStructure[token];
                    size_t K = structuresForThisProfile.size();
                    for(size_t k = 0; k < K; ++k){
                        std::pair<int,string> structure = structuresForThisProfile[k];
                        stringstream s;
                        s << structure.first << "-" << structure.second; // << "\t" << it3->second;
                        if(showStructures)cout << "-" << structure.first << "-" << structure.second << "-" << it3->second;
                        returnBestStructures->push_back(s.str());
                    }
                    if(showStructures) cout << endl;
                }
            }
        }
        std::pair<double, double> res;
        if(!modeUltraFast){
            std::map<string,double>::iterator it2 = (knownStatisticalAffinities.find(receptorAASeq));
            if(it2 == knownStatisticalAffinities.end()) cerr << "ERR: bestAffinity known but statistical affinity not known for " << receptorAASeq << endl;
            //cout << "Affinity for " << receptorAASeq << " already computed: " << it->second << endl;
            res = std::pair<double, double> (it->second, it2->second);
        } else {
            res = std::pair<double, double> (it->second, NAN);
        }
        // make sure to unlok before using return.
            pthread_mutex_unlock(&lockAccessPrecompAffinities);
        return res;
    }
    // ---------  Part 2:To speed-up the filling of interaction profiles with the AA sequence of the receptor (the input), makes a dictionnary: position in the receptor => Amino Acid at this position.


    // makes a dictionnary for each possible interactions, according to the receptor sequence.
    // an interaction is 2 letters: example: aA -> affinity(position 1, A),
    // inside an interaction either 'a-z': position in the receptor, or 'A-Z': AA on the ligand side.
    // does not put it as field of the affinityOneLigand class because each thread needs to use and access its own
    map<int, double> affSingleInteractions; // affSingleInteractions.clear();
    for(int i = 0; i < sizeReceptors+1; ++i){
        // cases of interaction between receptor and ligand (always position first and ligand second): [a..a+size-1][A-Z]
        for(int j = 0; j < 26; ++j){
            if((AA_ID((char) 'A' + j) != NB_AAs) && (AA_ID((char) 'A' + j) != UndefinedYet)) {
                int singleCode = code('a' + i, 'A' + j);
                double val = AAaffinity(AA_ID(receptorAASeq[i]), AA_ID((char) 'A' + j));
                affSingleInteractions.insert(pair<int, double>(singleCode, val));
                //cout << (char) ('a' + i) << (char) ('A' + j) << "(code" << singleCode << ")->" << receptorAASeq[i] << AAname(AA_ID((char) 'A' + j)) << " aff=" << val << endl;
            }
        }
        // cases of interaction between positions inside the receptor (easier)
        for(int j = 0; j < sizeReceptors+1; ++j){
            int singleCode = code('a' + i, 'a' + j); // any order.[a..a+size-1][a..a+size-1]
            double val = AAaffinity(AA_ID(receptorAASeq[i]), AA_ID(receptorAASeq[j]));
            affSingleInteractions.insert(pair<int, double>(singleCode, val)); // now, only interaction affinities.
            //cout << (char) ('a' + i) << (char) ('a' + j) << "->" << receptorAASeq[i] << receptorAASeq[j] << " aff=" << val << endl;
        }
    }

    //  --------- Part 3: Main job: browses all structures and computes the affinity (either best affinity or statistical affinity at the same time)

    double minAffinity = 1e6;
    double minBinding = 1e6;
    double avgBinding = 1e6;
    double summedZ = 0;
    int nOptSeq = 0;
    string bestAffProfile = string(""); // Will store the equally optimal profiles separated with two spaces (to keep even number of letters)

    // When computing statistical energies, will need to sum boltzmann weights. Summing big and small values in an arbitrary order creates precision errors. Ex: 1e6 + 1e-6 = 1e6 (error...)
    // when adding double numbers of different magnitudes, make a list, sort from smaller and then to add them later in this order.
    vector<double> summed_weights; summed_weights.reserve(nInterCodes + nFoldingCodes);
    vector<double> summed_statistical_affinity; summed_statistical_affinity.reserve(nInterCodes + nFoldingCodes);

    //#ifdef showBestStructures
    vector<std::pair<double, std::pair<double, string> > > bindTotSeqList;
    bindTotSeqList.reserve(nInterCodes + nFoldingCodes);
    //#endif


    //if(showStructures) cout << "Structure         \tInteract\tFolding\tTotal\tProba\tsumZ" << endl;
    for(int i = 0; i < nInterCodes; ++i){

        // for this folding/interaction profile,
        double affFold = 0;
        double affInteract = 0;
        int nSelf = 0;

        string s = (*interactions)[i];
        int ICS = s.size();
        if(((ICS / 2) * 2) != ICS) cerr << "ERR: non-even interaction code " << s;

        // gets the affinity of this particular interaction
        for(int j = 0; j < ICS; j = j + 2){
            int singleCode = code(s[j], s[j+1]);
            std::map<int, double>::iterator it = (affSingleInteractions.find(singleCode));
            if (it != affSingleInteractions.end()){ // found
                if(s[j+1] < 'a'){
                    affInteract += it->second; // A..Z : outside interactions  // note: char to int: a:97 -> z ; A:65 ->...
                    //nLinks++;
                } else {
                    affFold += it->second; // a..z : inside interactions
                    nSelf++;
                }
            } else {
                cerr << s << endl;
                cerr << s[j] << s[j+1] << "(code " << code(s[j], s[j+1]) << "), single interaction code not found [1] Maybe a nonexisting AA" << endl;
                exit(-1);
            }
        }

        // now we have 2 informations: the affInteraction towards the ligand and the self-binding affinity
        double totalAff = affInteract + affFold;

        if(showStructures){ //#ifdef showBestStructures
            bindTotSeqList.push_back(std::pair<double, std::pair<double, string> > ( affInteract, std::pair<double, string>( totalAff ,s)));
        }//#endif

        //cout << "Mayday" << totalAff << "," << KT << "," <<  nbRepeats[i] << "," << exp(-(totalAff / KT)) << endl;
        double proba_weight = exp(-(totalAff / KT)) * nbRepeats[i];
        if(std::isnan(proba_weight) || std::isinf(proba_weight)) proba_weight = 0;
        //if(nSelf < minNrSelfInteractions) {
            summed_weights.push_back(proba_weight); // if not, these structures will also be enumerated from the self-foldings. Problem here: the same structure can be enumerated twice around the ligand but only once as self-folding. Need to think about it, or quantify it ...
            summedZ += proba_weight;
        //}
        summed_statistical_affinity.push_back(proba_weight * affInteract); // the average binding energy only counts the binding affinity, not self-folding


        // keeps track of the best profiles so far (then the structures would need to be retrieved)
        if(totalAff == minAffinity){
            bestAffProfile.append(string("  ")+s);
            minBinding = min(minBinding, affInteract);
            avgBinding = ((avgBinding * (double) nOptSeq) + affInteract) / (double) (nOptSeq + 1);
            nOptSeq++;
        }
        if(totalAff < minAffinity){
            minAffinity = totalAff;
            bestAffProfile = s;
            minBinding = affInteract;
            avgBinding = affInteract;
            nOptSeq = 1;
        }
        if(showStructures) {
            //cout << s << "   \t" << affInteract << "\t" << affFold << "\t" << totalAff << "\t" << proba_weight << "\t" << ((nSelf < minNrSelfInteractions) ? summedZ : 0) << /*"\t" << summed_statistical_affinity << */ endl;
        }
    }
    //if(fabs(minBinding - avgBindiing) > 1e-12) cerr << "Inconsistency error in affinities " << endl;


/*
    //  --------- Part 4: Checks if a self-folding would produce a better total energy
    bool bestFromSelfFold = false;
    if(showStructures) cout << "Self-foldings       \tInteract\tFolding\tTotal\tProba\tsumZ\tsumAff" << endl;
    for(size_t i = 0; i < static_cast<size_t>(nFoldingCodes); ++i){

        double affFold = 0;
        string s = (*selfInteractions)[i];
        size_t ICS = s.size();
        if(((ICS / 2) * 2) != ICS) cerr << "ERR: non-even interaction code " << s;

        // gets the affinity of this particular interaction
        for(size_t j = 0; j < ICS; j = j + 2){
            int singleCode = code(s[j], s[j+1]);
            std::map<int, double>::iterator it = (affSingleInteractions.find(singleCode));
            if (it != affSingleInteractions.end()){ // found
                affFold += it->second; // a..z : inside interactions
            } else {
                cerr << s[j] << s[j+1] << "(code " << code(s[j], s[j+1]) << "), single interaction code not found [2]" << endl;
            }
        }

        double proba_weight = exp(-(affFold / KT)) * nbSelfRepeats[i];
        if(std::isnan(proba_weight) || std::isinf(proba_weight)) proba_weight = 0;
        summed_weights.push_back(proba_weight); // if not, these structures will also be enumerated from the self-foldings. Problem here: the same structure can be enumerated twice around the ligand but only once as self-folding. Need to think about it, or quantify it ...
        summedZ += proba_weight;

        // keeps track of the best profiles so far (then the structures would need to be retrieved)
        if(affFold == minAffinity){
            bestAffProfile.append(string("  ")+s);
            // what to do in that case ? Not obvious ...
        }
        if(affFold < minAffinity){
            bestFromSelfFold = true;
            minAffinity = affFold;
            bestAffProfile = s;
        }
        if(showStructures) cout << s << "   \t  -  \t" << affFold << "\t  -  \t" << proba_weight << "\t" << summedZ << "\t  -  " << endl;
    }

    if(bestFromSelfFold){
        cerr << "WRN: " << receptorAASeq << " a self folding has better energy than any of the receptor foldings around the ligand ... " << endl;
        minAffinity = 1e-12;
    }
*/

    // ======> Finished! Best binding energy is now known.

    //cout << "bestTotalAff " << minAffinity << " and resulting bindingAff " << avgBinding << " (min " << minBinding << " from profiles: " << bestAffProfile << endl;
    // cout << receptorAASeq << "\t" << avgBinding << "\t" << bestAffProfile << endl;
        pthread_mutex_lock(&lockAccessPrecompAffinities);

    knownBestAffinities.insert(std::pair<string,double>(receptorAASeq, avgBinding));

        pthread_mutex_unlock(&lockAccessPrecompAffinities);



    // --------- Part Optional: If requests the list of optimal structures, retrives it.
    if(returnBestStructures != nullptr){
        returnBestStructures->clear();

        vector<string> tokens = slice(bestAffProfile, "  ");
        size_t nCodes = tokens.size();
        for(size_t k = 0; k < nCodes; ++k){
            string token = tokens[k];
            if(showStructures) cout << token;
            vector<std::pair<int,string> > structuresForThisProfile = profileToStructure[token];
            size_t K = structuresForThisProfile.size();
            for(size_t k = 0; k < K; ++k){
                std::pair<int,string> structure = structuresForThisProfile[k];
                stringstream s;
                s << structure.first << "-" << structure.second; // << "\t" << bestAffProfile;
                if(showStructures)cout << "\t" << structure.first << "-" << structure.second << "-" << bestAffProfile;
                returnBestStructures->push_back(s.str());
            }
            if(showStructures) cout << endl;
        }
    }

        pthread_mutex_lock(&lockAccessPrecompAffinities);

    knownBestInteractions.insert(std::pair<string, string>(receptorAASeq, bestAffProfile));

        pthread_mutex_unlock(&lockAccessPrecompAffinities);



    // In ultrafast mode, just returns best energy
    if(modeUltraFast) return std::pair<double, double> (avgBinding, NAN);




    //  --------- Part 3: Now calculating the statistical affinity over all structures.

    // Now, for the statistical affinity, the energy of the self-folding structures has to be included in the computation
    // of the total possible weight. Namely structures that fold on themselves and that would have not been included
    // in the list of structures.

    std::sort(summed_statistical_affinity.begin(), summed_statistical_affinity.end());
    // numbers are negative, so should have smaller magnitude first
    std::reverse(summed_statistical_affinity.begin(), summed_statistical_affinity.end() );
    size_t NS1 = summed_statistical_affinity.size();
    double morePreciseSumAff = 0;
    for(size_t i = 0; i < NS1; ++i){
        morePreciseSumAff += summed_statistical_affinity[i];
    }

    // this one is positive, not reverted
    std::sort(summed_weights.begin(), summed_weights.end());
    size_t NS2 = summed_weights.size();
    double morePreciseSumWei = 0;
    for(size_t i = 0; i < NS2; ++i){
        morePreciseSumWei += summed_weights[i];
    }

    //double statAff = summed_statistical_affinity / (summed_weights + 1e-12); // this one would have all the computing errors
    double statAff = morePreciseSumAff / (morePreciseSumWei + 1e-12);


        pthread_mutex_lock(&lockAccessPrecompAffinities);

    knownStatisticalAffinities.insert(std::pair<string,double>(receptorAASeq, statAff));

        pthread_mutex_unlock(&lockAccessPrecompAffinities);





    // =============== Now, the computation is finished. The remaining block is only to display info about best structures ===================
    if(showStructures){
        // Show the best structures from binding affinity
        // This sorting is computationally very expensive.
        if(showSorted) std::sort(bindTotSeqList.begin(), bindTotSeqList.end(), compBind);
        if(bindTotSeqList.size() > 0){
            cout << receptorAASeq << "\t" << avgBinding << "\t" << bindTotSeqList[0].second.second;
            std::map<string, vector<std::pair<int, string> > >::iterator it = profileToStructure.find(bindTotSeqList[0].second.second);
            if(it !=profileToStructure.end() ){
                vector<std::pair<int,string>> listStr = it->second;
                cout << "\t" << listStr.size();
                for(size_t j = 0; j < listStr.size(); ++j){
                    cout << "\t" << listStr[j].first << "-" << listStr[j].second;
                    //bestSequences.push_back(new struct3D(listStr[j].second, UnDefined, listStr[j].first));
                }
            }
            cout << endl;
        }

        #ifdef ALLOW_GRAPHICS
        if(true){
            cerr << "Now displaying 2500 / " << bindTotSeqList.size() << " found structures " << endl;
            char *c[] = {(char*)"Hello",nullptr};
            glDisplay(0,c);
            //addToDisplay(merged);
            superProtein* ligandProt = new superProtein(*ligand);
            addToDisplay(ligandProt, true);

            vector<struct3D*> subSequences;
            int NS4 = bindTotSeqList.size();
            for(int i = 0; i < min(300000, NS4); ++i){
                std::map<string, vector<std::pair<int, string> > >::iterator it = profileToStructure.find(bindTotSeqList[i].second.second);
                if(it !=profileToStructure.end() ){
                    vector<std::pair<int,string>> listStr = it->second;
                    for(int j = 0; j < (int) listStr.size(); ++j){
                        subSequences.push_back(new struct3D(listStr[j].second, UnDefined, listStr[j].first));
                    }
                } else {
                    cerr << "ERR: " << bindTotSeqList[i].second.second << " not found in dictionary [1]" << endl;
                }
            }


            int NRs = subSequences.size();
            for(int i = 0; i < min(2500,NRs); ++i){
                int rd = i;// random::uniformInteger(0, NRs-1);
                if(!subSequences[rd]) {cerr << "WTF for best sequences " << rd << endl; }
                else {
                    superProtein* receptorProt = new superProtein(*subSequences[rd]);
                    addToDisplay(receptorProt, false);
                }
            }
            // Note: the program will stop here.
            glutMainLoop();
        }
        #endif



        /* This piece of code was made to make a heatmap where the possible structures of one receptor sequence bind (according to boltzman).
         * The code works but is a bit useless because the only very similar structures get a good total energy.

        cout << " ============= The receptor seaquence is " << receptorAASeq << " , Aff best " << avgBinding << ", stat " << statAff << endl;
        cerr << "   -> Now, will show the details of best structures. Files written, program can be switched off if stuck." << endl;
        // Now for each position of the ligand, gets the probability of being bound.
        set<int> ligandSpace = ligand->occupiedPositions;
        std::map<int, int> nbTouchingReceptors;
        std::map<int, double> weightTouchingReceptors;

        vector<struct3D*> bestSequences;

        //#ifdef showBestStructures
        int NS4 = bindTotSeqList.size();
        std::sort(bindTotSeqList.begin(), bindTotSeqList.end(), compBind);
        cout << "List of structures(profiles) from the best ones, sorted by decreasing binding affinity" << endl;
        cout << "BindAff \tTotAff  \tProfile\tstructures\n";
        for(int i = 0; i < min(250, NS4); ++i){
            if(minBinding)
            cout << bindTotSeqList[i].first << "\t" << bindTotSeqList[i].second.first << "\t" << bindTotSeqList[i].second.second;
            std::map<string, vector<std::pair<int, string> > >::iterator it = profileToStructure.find(bindTotSeqList[i].second.second);
            if(it !=profileToStructure.end() ){
                vector<std::pair<int,string>> listStr = it->second;
                for(int j = 0; j < (int) listStr.size(); ++j){
                    cout << "\t" << listStr[j].first << "-" << listStr[j].second;
                    //bestSequences.push_back(new struct3D(listStr[j].second, UnDefined, listStr[j].first));
                }
            } else {
                cerr << "ERR: " << bindTotSeqList[i].second.second << " not found in dictionary [1]" << endl;
            }
            cout << endl;
        }
        std::sort(bindTotSeqList.begin(), bindTotSeqList.end(), compTot);
        cout << "List of structures(profiles) from the best ones, sorted by decreasing TOTAL affinity" << endl;
        cout << "BindAff \tTotAff  \tProfile\tstructures\n";
        for(int i = 0; i < min(250, NS4); ++i){
            cout << bindTotSeqList[i].first << "\t" << bindTotSeqList[i].second.first << "\t" << bindTotSeqList[i].second.second;
            std::map<string, vector<std::pair<int,string> > >::iterator it = profileToStructure.find(bindTotSeqList[i].second.second);
            if(it != profileToStructure.end() ){
                vector<std::pair<int,string>> listStr = it->second;
                for(int j = 0; j < (int) listStr.size(); ++j){
                    cout << "\t" << listStr[j].first << "-" << listStr[j].second;
                    bestSequences.push_back(new struct3D(listStr[j].second, UnDefined, listStr[j].first));
                }

            } else {
                cerr << "ERR: " << bindTotSeqList[i].second.second << " not found in dictionary [2]" << endl;
            }
            cout << endl;
        }



        // preparing outputs for this block
        ofstream plot1("plotBindInteractToEnergy");
        vector<double> bindEnergiesStructures;  // for histogram
        bindEnergiesStructures.reserve(1e6);
        vector<double> totEnergiesStructures;  // for histogram
        totEnergiesStructures.reserve(1e6);
        vector<double> listWeightsPerBindingEnergy(500, 0.);
        vector<double> listWeightsPerTotalEnergy(500, 0.);

        // bindTotSeqList stores, for each interaction code (.second.second), total affinity(.second.first) and binding affinity(.first).
        // note: smaller first because needs to sum
        std::sort(bindTotSeqList.begin(), bindTotSeqList.end(), compNegTot);
        for(int i = 0; i < NS4; ++i){       // for each interaction code,
            //cout << bindTotSeqList[i].first << "\t" << bindTotSeqList[i].second.first << "\t" << bindTotSeqList[i].second.second;
            string interactionCode = bindTotSeqList[i].second.second;
            std::map<string, vector<std::pair<int, string> > >::iterator it = profileToStructure.find(interactionCode);
            if(it != profileToStructure.end() ){

                // for each structure that represents this interaction code,
                vector<std::pair<int,string>> listStr = it->second;
                for(int j = 0; j < (int) listStr.size(); ++j){
                    struct3D a = struct3D(listStr[j].second, UnDefined, listStr[j].first);
                    set<int> touched = neighborPositions(a);
                    double totAff = bindTotSeqList[i].second.first;
                    double proba_weight = exp(-(totAff / KT)); // * nbRepeats[i]; no need for nb repeats because do it for all structs, a bit tedious

                    // Writing outputs
                    double bindAff = bindTotSeqList[i].first;
                    bindEnergiesStructures.push_back(bindAff);
                    totEnergiesStructures.push_back(totAff);
                    plot1 << nrBindings(interactionCode) << "\t" << bindAff << "\t" << totAff << "\t" << proba_weight << "\t" <<  proba_weight / morePreciseSumWei << "\t" << interactionCode << "\n";
                    listWeightsPerBindingEnergy[250 + int(bindAff-0.5)] += proba_weight; // hope doesnt exceed -250 or +250 .... Note, they are negative so need -0.5 and not +0.5
                    listWeightsPerTotalEnergy[250 + int(totAff-0.5)] += proba_weight;

                    // list of positions touched
                    set<int> overlap ;
                    set_intersection(ligandSpace.begin(),ligandSpace.end(),touched.begin(),touched.end(),
                                      std::inserter(overlap,overlap.begin()));
                    //cout << listStr[j].second << " -> ";
                    for(set<int>::iterator it = overlap.begin(); it != overlap.end(); ++it){
                        if(nbTouchingReceptors.find(*it) != nbTouchingReceptors.end()){
                            nbTouchingReceptors[*it] += 1;
                        } else {
                            nbTouchingReceptors[*it] = 1;
                        }
                        if(weightTouchingReceptors.find(*it) != weightTouchingReceptors.end()){
                            weightTouchingReceptors[*it] += proba_weight / max(morePreciseSumWei, 1e-12);
                        } else {
                            weightTouchingReceptors[*it] = proba_weight / max(morePreciseSumWei, 1e-12);
                        }
                        //cout << *it << "\t" << totAff << "\t" << proba_weight << "\t" << proba_weight / (morePreciseSumWei + 1e-12);
                    }
                    //cout << endl;
                    //cout << "\t" << listStr[j].first << "-" << listStr[j].second;
                }
            } else {
                cerr << "ERR: " << bindTotSeqList[i].second.second << " not found in dictionary [3]" << endl;
            }
            //cout << endl;
        }

        cout << "   -> Writing the number of bindings and energies for each structure" << endl;
        plot1.close();

        vector<double> intBoundaries;
        for(double d = -120; d < +120; ++d){
            intBoundaries.push_back(d-0.5);
        }
        cout << "   -> Distributions of binding energy among ALL possible structuresm for this receptor AA sequence" << endl;
        histogramFromDistrib h1(bindEnergiesStructures, intBoundaries);
        cout << h1.print(true) << endl;

        cout << "   -> Distributions of total energy among ALL possible structuresm for this receptor AA sequence" << endl;
        histogramFromDistrib h2(totEnergiesStructures, intBoundaries);
        cout << h2.print(true) << endl;

        cout << "   -> Average statistical weight of each class of binding energy" << endl;
        for(int i = 100; i < min(300,(int)listWeightsPerBindingEnergy.size()); ++i){
            cout << (int) - 250 + i << "\t" << listWeightsPerBindingEnergy[i] / morePreciseSumWei << endl;
        }
        cout << "   -> Average statistical weight of each class of total energy" << endl;
        for(int i = 100; i < min(300,(int)listWeightsPerTotalEnergy.size()); ++i){
            cout << (int) - 250 + i << "\t" << listWeightsPerTotalEnergy[i] / morePreciseSumWei << endl;
        }


        cout << "Details of binding per ligand position\n";
        cout << "Pos\tx\ty\tz\tNrRecept\tProba\n";
        for(set<int>::iterator it = ligandSpace.begin(); it != ligandSpace.end(); ++it){
            cout << "Pos " << *it << " ( " << printVector(lattice::positionFromID(*it)) << "\t";
            if(nbTouchingReceptors.find(*it) != nbTouchingReceptors.end()){
                cout << nbTouchingReceptors[*it];
            } else cout << "0";
            if(weightTouchingReceptors.find(*it) != weightTouchingReceptors.end()){
                cout << "\t" << weightTouchingReceptors[*it];
            } else cout << "\t0";
            cout << endl;
        }
        for(map<int, double>::iterator it = weightTouchingReceptors.begin(); it != weightTouchingReceptors.end(); ++it){
            cout << it->first << " -> " << it->second << endl;
        }
        */






    }

    return std::pair<double, double> (avgBinding, statAff);
}

// Question : does removing internal interaction changes something ???
//map<char[2], double> affSingleInteractions; // to learn once per sequence, and will be applied to all interaction codes
// saves the affinity for computed sequences
//map<string, double> knownAffinities;
//};

void affinityOneLigand::printInfos(){
    stringstream res;
    res << "affinityOneLigand: Ligand=" << ligandSeq << ", receptors expected of size sizeReceptors=L=" << sizeReceptors << " BOUNDS, i.e. " << sizeReceptors + 1 << " residues " << endl;
    res << "number of interactions stored: " << nInterCodes << endl;
    res << "Single Interaction matrix for the previous receptor sequence" << endl;
}


/*void oldtestFastAffinity(){
    cout << "Testing the FastAffinity class from a file containing a list of interactions (put inside 'TestFileAffinity.txt'" << endl;
    ofstream test("TestFileAffinity.txt");
    test << "13	AGNTGYMPARNW	6"
    "aAaNbWdWfWbg	4"
    "aAcAfReNcf	2"
    "aNeMfPdRfRbe	2"
    "bGcNcGgGadcf	2"
    "bYdYdPfPgA	2"
    "cPgPbAgRadcf	1"
    "dMdAcRbNaWbedg	2"
    "eGdNdGcYbgdg	2"
    "eTdGaYcYcPcf	1"
    "fMfAeRdN	17"
    "gGdYfYcMbPfPaAbe	1"
    "gNeMdRcNbW	3"
    "gYbMcPgP	3";
    test.close();
    affinityOneLigand aol = affinityOneLigand("TestFileAffinity.txt", 1.0);
    cout << "Affinity for receptor ARGTHC:" << aol.affinity(string("ARGTHC")).first << endl;
    cout << "Affinity for receptor ARGTHCK:" << aol.affinity(string("ARGTHCK")).first << endl;
    cout << "Affinity for receptor ARGTHCK:" << aol.affinity(string("ARGTHCK")).first << endl;
    cout << "Affinity for receptor GHXKLMK:" << aol.affinity(string("GHXKLMK")).first << endl;
    cout << "Affinity for receptor GHXKLMK:" << aol.affinity(string("GHTKLMK")).first << endl;
}*/



void testFastAffinity(){
    string simpleAccessible = string("UDRLLRLLRLLR");
    string AAsimple = string("AAAAAAAAAAAAA");
    int receptorSize = 8;
    affinityOneLigand T1 = affinityOneLigand(simpleAccessible, AAsimple, -1, receptorSize, 4, 4, 0.4);

    cout << "Computing different affinities for ligand " << simpleAccessible << " (" << AAsimple << "), receptors " << receptorSize << " minI=4" << endl;
    for(int i = 0; i < 100; ++i){
        string Px = randomProt(receptorSize+1);
        std::pair<double, double> res = T1.affinity(Px);
        cout << "affinity(" << Px << ") = \t" << res.first << "\t" << res.second << endl;
    }

    cout << "Details of the structures and affinities for " << simpleAccessible << " (" << AAsimple << "), receptors " << receptorSize << " minI=4" << endl;
    for(int i = 0; i < 1; ++i){
        string Px = randomProt(receptorSize+1);
        std::pair<double, double> res = T1.affinity(Px, true);
        cout << "affinity(" << Px << ") = \t" << res.first << "\t" << res.second << endl;
    }
}


/*
#define stepMultiV 10000
struct multiVector {
    multiVector() : _size(0), NL(0) {}
    vector< vector<string>*> contents;
    int size(){return _size;}
    int _size;
    int NL;
    string operator [](int i){
        if((i < 0) || (u >= _size)){
            return string("OUT OF BOUNDS");
        }
        int id = i / stepMultiV;
        return contents[id][i - id * stepMultiV];
    }
    void push_back(string & toAdd){
        if(((_size + 1) / stepMultiV) > contents.size()
    }
};*/

/*affinityOneLigand::affinityOneLigand(string fname, double _KT)
{
    KT = _KT;
    interactions = new vector<string>;
    ifstream fread(fname.c_str());
    if(!fread){
        cerr << "ERR: " << fname << ", file not found ";
    }
    fread >> nInterCodes;
    fread >> ligandSeq;
    fread >> sizeReceptors;
    interactions->resize(nInterCodes);
    nbRepeats.resize(nInterCodes);
    for(int i = 0; i < nInterCodes; ++i){
        fread >> (*interactions)[i] >> nbRepeats[i];
    }
    cout << "File " << fname << " loaded, with " << nInterCodes << " non redundant possible interaction profiles" << endl;
}*/


