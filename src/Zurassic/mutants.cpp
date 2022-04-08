#include <vector>
#include <string>
#include <map>
#include <set>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include "../Absolut/Absolut.h"
using namespace std;

#define allAAs string("CMFILVWYAGTSNQDEHRKP")



// Already defined within Absolut!
/*enum AA{Cys, Met, Phe, Ile, Leu, Val, Trp, Tyr, Ala, Gly, Thr, Ser, Asn, Gln, Asp, Glu, His, Arg, Lys, Pro, UndefinedYet, NB_AAs};

AA AA_ID(char AA_name){
    switch(AA_name){
    case 'C': return Cys;
    case 'M': return Met;
    case 'F': return Phe;
    case 'I': return Ile;
    case 'L': return Leu;
    case 'V': return Val;
    case 'W': return Trp;
    case 'Y': return Tyr;
    case 'A': return Ala;
    case 'G': return Gly;
    case 'T': return Thr;
    case 'S': return Ser;
    case 'N': return Asn;
    case 'Q': return Gln;
    case 'D': return Asp;
    case 'E': return Glu;
    case 'H': return His;
    case 'R': return Arg;
    case 'K': return Lys;
    case 'P': return Pro;
    case '?': return UndefinedYet;
    default: {return NB_AAs;}
    }
}

char randomAA(){
    return allAAs[random::uniformUInteger(0,19)];
}
*/



// mask=1 = mutable - can be derived from paratope or epitope encodings.
string mutateNPositions(string seq, int nMut, string mask) {
    size_t size = seq.size();
    if(mask.size() != size) {
        cerr << "ERR:mutateOnePosition, wrong mask size " << endl;
        return string("Wrong Mask");
    }
   int cpt = 0;
   int success = 0;
   while(cpt < 25){
       char test = randomAA();
       size_t position = random::uniformUInteger(0,size-1);
       if((test != seq[position]) && (mask[position] == '1')){
            seq[position] = test;
            success++;
            mask[position] = '0'; // do not mutate again there
            if(success == nMut){
                return seq; // success
            }
       }
       cpt++;
   }
   cerr << "ERR: foldedFree::mutateOnePosition, got trouble to mutate an AA" << endl;
   return string("MutateFailed"); // silent
}



void testAAs(){
    for(size_t i = 0; i < 20; ++i){
          cout << allAAs[i] << "\t" << allAAs[AA_ID(allAAs[i])] << endl;
    }
}


// Gets a different one, uniform proba
char mutateAA(char AA){
    size_t ID = AA_ID(AA);
    size_t newAAid = random::uniformUInteger(0,19-1); // boundaries included
    if(newAAid >= ID){newAAid++;}
    return allAAs[newAAid];
}




//#include <algorithm>
//std::string s = "a_b_c";
//size_t n = std::count(s.begin(), s.end(), '_');




// Strategy 1: vicinity of selected sequences, depending on their affinity
// includes the original sequence
set<string> allMutants(string seq, int nMut, string mask = ""){

    set<string> res = set<string>();

    // Stop condition, since the function is recursive
    if(nMut == 0) return set<string>({seq});

    size_t L = seq.size();
    if(mask.size() != L){mask = string(L, '0');}

    if(count(mask.begin(), mask.end(), '0') < nMut){
        cerr << "ERR: allMutants, seq=" << seq << ", nMut=" << nMut << ", mask=" << mask << ", not enough mutable positions" << endl;
    }

    // For every point mutation, will recursively call to mutate other positions
    for(size_t i = 0; i < L; ++i){

        // Will mutate only not masked positions
        if(mask[i] == '0'){

            // The AA to modify
            size_t ID = AA_ID(seq[i]);
            string newMask = mask;
            newMask[i] = '1';  // blocking this position to not be mutated further

            // replacing AA at position i by any other AA
            for(size_t mutatedID = 0; mutatedID < 20; ++mutatedID){
                if(mutatedID != ID){
                    string newSeq = seq;
                    newSeq[i] = allAAs[mutatedID];

                    set<string> generated = allMutants(newSeq, nMut-1, newMask);
                    res.insert(generated.begin(), generated.end());
                }
            }
        }
    }
    return res;
}





// From proteins.h/cpp
string randomProt(size_t size){
    string res = string(size, '?');
    for(size_t i = 0; i < size; ++i){
        size_t newAAid = random::uniformUInteger(0,19);
        res[i] = allAAs[newAAid];
    }
    return res;
}


void randomize(string content, string mask = "") {
    if(mask.size() == 0){
        mask = string(content.size(),'0');
    }
    string randomContent = randomProt(content.size());
    if(mask.size() != randomContent.size()) cerr << "ERR: FoldedFree, inconsistent size of content and conserved positions" << endl;
    for(unsigned int i = 0; i < mask.size(); ++i){
        if(!mask[i]) content[i] = randomContent.at(i);
    }
}




// Meaningful functions. returns 1: success, 0: silent, -1: fail
int mutateOnePosition(string content, string conservedPositions) {
   size_t size = content.size();
   int cpt = 0;
   while(cpt < 100){
       char test = randomAA();
       size_t position = random::uniformUInteger(0,size-1);
       if((test != content[position]) && (conservedPositions[position] == '0')){
            content[position] = test;
            return +1; // success
       }
       cpt++;
   }
   cerr << "ERR: mutateOnePosition, got trouble to mutate an AA" << endl;
   return 0; // silent
}




// Algos around one antigen, information: antigen, binders, non-binders, or sequences + affinity

struct singleAntigenBinary {
    int length; // put -1 if variable length
    vector<string> binders;
    vector<string> nonBinders;


};

struct singleAntigenRegression {
    double thresholdBinder;     // just define one threshold, we can make more later
    vector<map<string, double>> dataset;

};



/*antigenInfo getAntigenInfos(std::string ID){
    // note: the global fields are cleared inside getAntigen()

    std::pair<superProtein*, vector<int> > AG = getAntigen(ID);

    antigenInfo res;
    res.ID = ID;
    res.first = AG.first;
    res.second = AG.second;
    res.glycans = lastSearchGlycan;
    res.hotspotsCore = lastSearchHotspotsCore;
    res.hotspotsLarge = lastSearchHotspotsLarge;
    res.antibodyChains = lastAntibodyChains;
    return res;
}*/


string nDigitsNumbers(int n, int L){
    stringstream s;
    for(int k = 0; k < L; ++k){
        if(L < std::pow(10, k)) s << "0";
    }
    s << n;
    return s.str();
}



// Create a name for a mutated sequence: 1ADQ_A+P050Q+R1052P" ...
string generateMutatedAntigenCode(string nameAntigen, string mutatedSeq, superProtein* P){
    string initialSeq = P->getAAseq();
    if(initialSeq.size() != mutatedSeq.size()){cerr << "ERR: generateCode, initial and mutated Seq of different sizes" << endl;}

    if(!initialSeq.compare(mutatedSeq)) {return nameAntigen;}
    stringstream res;

    res << nameAntigen;

    size_t L = initialSeq.size();
    for(size_t i = 0; i < L; ++i){
        if(initialSeq[i] != mutatedSeq[i]){
            res << "+" << initialSeq[i] << P->points[i].IDresidue /* nDigitsNumbers(P->points[i].IDresidue, 4)*/ << mutatedSeq[i];
        }
    }

    return res.str();
}

// probably very unefficient... as long as it works!
set<string> randomSubset(set<string> bigSet, size_t maxKeep){
    size_t L = bigSet.size();
    if(L < maxKeep) return bigSet;
    vector<string> shuf = vector<string>(bigSet.begin(), bigSet.end());
    random::shuffle(shuf);
    shuf.resize(maxKeep);
    return set<string>(shuf.begin(), shuf.end());
}



// starts from an AA sequence, test the quality, then looks as N-points mutants and looks at the quality as well,
// Mask = 0100011 with 1 are positions that can be mutated.
// note: the random sequence is NOT following the mask, it's completely independent of the original one
// use NbMuts < 0 tto completely randomize each time

// InitialSeq allows to start with another sequence directly.
void createAntigenMutants(string antigenID, int muts, size_t maxEachType = 50000, string outputFile = "ListGeneratedAntigens.txt", bool mutatePerHotspot = true){

    precise_stopwatch stopwatch; // careful, doesnt work for more than 2 hours

    antigenInfo AG = getAntigenInfos(antigenID);

    // Checking the antigen
    superProtein* AGstruct = AG.first;
    vector<int> blockedPositions = AG.second;
    string AAseq = AG.first->getAAseq();
    size_t L = static_cast<size_t>(AGstruct->size());
    if(L != AAseq.size()){cerr << "ERR on library antigen, different structure and sequence sizes" << endl;}

    // Checking its hotspots
    vector<vector<int>> hotspotsLarge = AG.hotspotsLarge;
    vector<vector<int>> hotspotsCore = AG.hotspotsCore;
    size_t numberHotspots = hotspotsCore.size();
    if(numberHotspots != hotspotsLarge.size()) cerr << "ERR, antigen" << antigenID << ", different nr of core and large hotspots" << endl;

    // Now preparing the output
    ofstream ListAntigens(outputFile.c_str(), ios::app);
    set<string> generated;


    for(size_t nH = 0; nH < numberHotspots; nH++){

        set<int> HCore;
        set<int> HLargeAndCore;
        stringstream nameHotspot; nameHotspot << "H" << nH+1;
        string IDhot = nameHotspot.str();

        // if not per hotspot, merge all positions of all hotspots
        if(!mutatePerHotspot){
            nH = numberHotspots;    // will not loop, just once
            IDhot = "All";
            for(size_t i = 0; i < hotspotsCore.size(); ++i){
                HCore.insert(hotspotsCore[i].begin(),hotspotsCore[i].end());
            }
            for(size_t i = 0; i < hotspotsLarge.size(); ++i){
                HLargeAndCore.insert(hotspotsLarge[i].begin(),hotspotsLarge[i].end());
            }
        } else { // only positions of current hotspots
            HCore = set<int>(hotspotsCore[nH].begin(),hotspotsCore[nH].end());
            HLargeAndCore = set<int>(hotspotsLarge[nH].begin(),hotspotsLarge[nH].end());
        }


        // Creates mask for hotspot large and small - all positions start with '1' (blocked) and only hotspot positions will be mutable
        string maskCore = string(L, '1');
        string maskCoreAndSide = string(L, '1');
        string maskSideOnly = string(L, '1');

        for(size_t i = 0; i < L; ++i){
            int residueID = AGstruct->points[i].IDresidue;
            //int positionID = AGstruct->points[i].IDposition;
            //AA nameAA = AGstruct->points[i].TypeResidue;

            if(HCore.find(residueID) != HCore.end()){
                maskCore[i] = '0';
            }
            if(HLargeAndCore.find(residueID) != HLargeAndCore.end()){
                maskCoreAndSide[i] = '0';
            }
            if((maskCore[i] == '1') && (maskCoreAndSide[i] == '0')) {
                maskSideOnly = '0';
            }
        }

        // store the original AA sequence
        string originalSeq = AAseq;

        //ListAntigens << "# Mutants " << antigenID << "\t" << maskCore << "\t" << maskCoreAndSide << endl;
        //for(int muts = 0; muts <= NbMuts; ++muts){

            // Mutating only the core hotspot residues
            set<string> resCore = randomSubset(allMutants(originalSeq, muts, maskCore), maxEachType);

            set<string> resAll = randomSubset(allMutants(originalSeq, muts, maskCoreAndSide), maxEachType);
            resAll.insert(resCore.begin(), resCore.end());

            for(set<string>::iterator it = resAll.begin(); it != resAll.end(); ++it){
                bool onlyCore = (resCore.find(*it) != resCore.end());
                ListAntigens << generateMutatedAntigenCode(antigenID, *it, AGstruct) << "\t" << antigenID << "\t" << muts << "\t";
                ListAntigens << IDhot << "\t" << ((onlyCore)? "C": "A") << "\t" << *it << endl;
            }
        //}
    }
    ListAntigens.close();
}



 void generateLibraryMutants(){
    vector<string> availableAGs = {"1ADQ_A",
                                   "1FBI_X",
                                   "1FNS_A",
                                   "1FSK_A",
                                   "1H0D_C",
                                   "1JPS_T",
                                   "1KB5_AB",
                                   "1NCB_N",
                                   "1NSN_S",
                                   "1OAZ_A",
                                   "1OB1_C",
                                   "1OSP_O",
                                   "1PKQ_J",
                                   "1QFW_AB",
                                   "1RJL_C",
                                   "1S78_B",
                                   "1TQB_A",
                                   "1WEJ_F",
                                   "1YJD_C",
                                   "1ZTX_E",
                                   "2ARJ_RQ",
                                   "2B2X_A",
                                   "2FD6_AU",
                                   "2HFG_R",
                                   "2IH3_C",
                                   "2JEL_P",
                                   "2Q8A_A",
                                   "2R29_A",
                                   "2R4R_A",
                                   "2R56_A",
                                   "2UZI_R",
                                   "2VXQ_A",
                                   "2VXT_I",
                                   "2W9E_A",
                                   "2WUC_I",
                                   "2XQB_A",
                                   "2XWT_C",
                                   "2YC1_C",
                                   "2YPV_A",
                                   "2ZCH_P",
                                   "3BGF_S",
                                   "3BN9_A",
                                   "3CVH_ABC",
                                   "3DVG_XY",
                                   "3EFD_K",
                                   "3GI9_C",
                                   "3HI6_A",
                                   "3JBQ_B",
                                   "3KJ4_A",
                                   "3KR3_D",
                                   "3KS0_J",
                                   "3L5X_A",
                                   "3L95_X",
                                   "3MJ9_A",
                                   "3NCY_A",
                                   "3NFP_I",
                                   "3NH7_A",
                                   "3Q3G_E",
                                   "3R08_E",
                                   "3R1G_B",
                                   "3RAJ_A",
                                   "3RKD_A",
                                   "3RVV_A",
                                   "3SKJ_E",
                                   "3SQO_A",
                                   "3TT1_A",
                                   "3U9P_C",
                                   "3UBX_A",
                                   "3V6O_A",
                                   "3VG9_A",
                                   "3VRL_C",
                                   "3WD5_A",
                                   "4AEI_A",
                                   "4CAD_C",
                                   "4DKE_A",
                                   "4H88_A",
                                   "4HC1_B",
                                   "4HJ0_B",
                                   "4I18_R",
                                   "4I77_Z",
                                   "4K24_A",
                                   "4K3J_A",
                                   "4K9E_C",
                                   "4KI5_M",
                                   "4KXZ_A",
                                   "4LU5_B",
                                   "4MXV_B",
                                   "4N9G_C",
                                   "4NP4_A",
                                   "4OII_A",
                                   "4OKV_E",
                                   "4PP1_A",
                                   "4QCI_D",
                                   "4QEX_A",
                                   "4QNP_A",
                                   "4QWW_A",
                                   "4R9Y_D",
                                   "4RGM_S",
                                   "4U1G_A",
                                   "4U6V_A",
                                   "4WV1_F",
                                   "4Y5V_C",
                                   "4YPG_C",
                                   "4YUE_C",
                                   "4ZFG_A",
                                   "4ZFO_F",
                                   "4ZSO_E",
                                   "5B8C_C",
                                   "5BVP_I",
                                   "5C0N_A",
                                   "5C7X_A",
                                   "5CZV_A",
                                   "5D93_A",
                                   "5DFV_A",
                                   "5DHV_M",
                                   "5DMI_A",
                                   "5DO2_B",
                                   "5E8D_A",
                                   "5E8E_LH",
                                   "5E94_G",
                                   "5EII_G",
                                   "5EPM_C",
                                   "5EU7_A",
                                   "5EZO_A",
                                   "5F3B_C",
                                   "5FB8_C",
                                   "5H35_C",
                                   "5HDQ_A",
                                   "5HI4_B",
                                   "5IKC_M",
                                   "5J13_A",
                                   "5JW4_A",
                                   "5JZ7_A",
                                   "5KN5_C",
                                   "5KTE_A",
                                   "5L0Q_A",
                                   "5LQB_A",
                                   "5MES_A",
                                   "5T5F_A",
                                   "5TH9_A",
                                   "5TLJ_X",
                                   "5TZ2_C"};
    for(size_t i = 0; i < availableAGs.size(); ++i){
        cout << availableAGs[i] << endl;
        cout << "   -> 1 mutation" << endl;
        createAntigenMutants(availableAGs[i], 1, 500, "ListMutantAntigens2.5k.txt", true); // We always mutate per hotspot, so that we know on which position the mutation was... duplicates possible, it's fine
        //cout << "   -> 2 mutations anywhere" << endl;
        //createAntigenMutants(availableAGs[i], 2, 2000, "ListMutantAntigens2.5k.txt", false); // any two mutations wherever (control, less interesting)
        cout << "   -> 2 mutations same hotspot" << endl;
        createAntigenMutants(availableAGs[i], 2, 500, "ListMutantAntigens2.5k.txt", true); // any two mutations on the same hotspot
        cout << "   -> 3 mutations " << endl;
        createAntigenMutants(availableAGs[i], 3, 1500, "ListMutantAntigens2.5k.txt", true);  // any three mutations on the same hotspot
    }
}



pair<superProtein*, vector<int>> getFileAntigen(string antigenID, string fileName){
    pair<superProtein*, vector<int>> res(nullptr, vector<int>());

    ifstream f(fileName);
    if(!f) cerr << "ERR: getFileAntigen, file " << fileName << endl;

    string AGname = "start", AGorigin, location, sequence;
    char CoreOrAll;
    int nMuts;

    int cpt = 0;
    bool found = false;
    while((AGname.size() > 0) && (cpt < 100000000) && (!found)){
        AGname = "";
        f >> AGname >> AGorigin >> nMuts >> location >> CoreOrAll >> sequence;
        //if(cpt < 100){cout << AGname << "\t" << AGorigin << "\t" << nMuts << "\t" << location << "\t" << CoreOrAll << "\t" << sequence << "\n";}
        if(!AGname.compare(antigenID)){found = true;}
        cpt++;
    }
    f.close();
    if(!found){
        cerr << "ERR: antigen " << antigenID << ", not found in file " << fileName << endl;
        return res;
    }

    pair<superProtein*, vector<int>> AG = getAntigen(AGorigin);
    AG.first->setAAs(sequence);
    return AG;
}







#include <sstream>
#include <algorithm>
#include <vector>

// Returns the list of moves (each impacting distance one)
string LevenshteinDistance(string source, string target) {


    if (source.size() > target.size()) {
        return LevenshteinDistance(target, source);
    }
    // then it is reverse order, be careful!
    random::initialize();

    size_t min_size = source.size(), max_size = target.size();
    std::vector<size_t> lev_dist(min_size + 1);
    std::vector<string> history(min_size + 1); // for each level, all possible strings

    //cout << "j=0";
    for (size_t i = 0; i <= min_size; ++i) {
        lev_dist[i] = i;
        //history[i] = source; //.substr(0, i);
        //stringstream hst;
        history[i] = "";
        for(size_t k = 0; k < min(i, min_size); ++k){
             history[i] = history[i] + "D" + string(1, 'a' + k) + source[k] + "+";
        }
        //cout << "\t" << i << "_" << history[i];
    }
    //cout << endl;

    for (size_t  j = 1; j <= max_size; ++j) {
        //cout << "j=" << j;

        size_t  previous_diagonal = lev_dist[0];
        string previous_diagonal_history = history[0];
        size_t previous_diagonal_save;
        string previous_diag_history_save;

        history[0] = history[0] + "Ia" + string(1, target[j-1]) + "+";
        //cout << "\t" << 0 << "_" << history[0];

        ++lev_dist[0];


        for (size_t  i = 1; i <= min_size; ++i) {
            previous_diagonal_save = lev_dist[i];
            previous_diag_history_save = history[i];
            if (source[i - 1] == target[j - 1]) {
                lev_dist[i] = previous_diagonal;
                history[i] = previous_diagonal_history;
                // history[i] = history[i]; // + string(1, target[i - 1]);
            } else {
                size_t minVal = std::min(std::min(lev_dist[i - 1], lev_dist[i]), previous_diagonal) + 1;

                // will make a list of possible choices, shuffle and take the first one
                vector<string> possibles;

                // we are at position i,
                if(minVal == lev_dist[i-1]+1){ // delete one letter from the source
                    possibles.push_back(history[i-1] + "D" + string(1, 'a' + i-1) + string(1, source[i-1]) + "+"); //
                }
                if(minVal == lev_dist[i]+1){ // insert one letter from the target target[j-1] at position [source[i-1]], but then it will get longer,
                    possibles.push_back(history[i] + "I" + string(1, 'a' + i-1) + string(1, target[j-1]) + "+"); // insert letter from target
                }
                if(minVal == previous_diagonal+1){ // substitution: delete and Insert new
                    possibles.push_back(previous_diagonal_history + "S" + string(1, 'a' + i-1) + string(1, target[j-1]) + "+");
                                        //"SI" + string(1, 'a' + i-1) + string(1, target[j-1]) + "/D" + string(1, '0' + i-1) + string(1, source[i-1]) + "+");
                }

                if(possibles.size() < 1) {
                    cerr << "ERR: impossible null amount of possibles. God, please restore the possibles." << endl;
                    possibles.push_back("ERROR");
                }
                random::shuffle<string>(possibles);
                history[i] = possibles[0];
                lev_dist[i] = minVal;

            }
            //cout << "\t" << i << "_" << history[i];
            previous_diagonal = previous_diagonal_save;
            previous_diagonal_history = previous_diag_history_save;
        }
        //cout << endl;
    }

    return history[min_size];
}

// Note: there might be asymmetry in the choice of paths to reconstitute the actions in the function above, but here we always have same size so we are good!

vector<string> split_str(string str, char delim){
    vector<string> out;
    stringstream s(str);
    string s2;
    while (getline (s, s2, delim) )
    {
        out.push_back(s2); // store the string in s2
    }
    return out;
}

string pulledprint(vector<string> sliced){
    stringstream res;
    for(size_t j = 0; j < sliced.size(); ++j){
        res << sliced[j];
    }
    return res.str();
}

// Assumptions:
// 1- the insertions are decided after the deletions
// 2- there can not be deletion + insertion on the same place BUT
//      there can be substitution + insertions
// 3 - there can not be more than one deletion at a place
// 4 - if more than one insertion at a place, they can be done at random sequence BUT keep their relative order.

// we keep the first position for the original Letter, if deleted becomes _ and if substituted it is stil lthe first positions (insertions come after)
// we prepare the insertions by creating "_" as many as insertions, and they will happen at the good subposition inside these ___ in random order

enum typePath {FROM_LEFT, RANDOM, MAINTAIN_SIZE};

// Note: maintain size has only sense if the sequences are +/-1 the same size!!
vector<string> applySteps(string start, string steps, typePath typeShuffling = RANDOM){


    // transform the start as a vector of strings (since each position can be enlarged or deleted)
    vector<string> splittedStart;
    for(size_t i = 0; i < start.size(); ++i){
        splittedStart.push_back(string(1, start[i]));
    }

    // now get the orders!
    vector<string> todos = split_str(steps, '+');

    // Now we will count the insertions at each position and annotate the tasks themselves with their relative ID
    int cptperpos = 0;
    int lastpos = -1;
    for(size_t i = 0; i < todos.size(); ++i){
        if(todos[i].size() != 3) cerr << "ERR: Unknown task: " << todos[i] << endl;
        int position = todos[i][1] - 'a';
        if(position != lastpos) {
            if(cptperpos > 0){
                splittedStart[lastpos] = splittedStart[lastpos] + string(cptperpos, '_');
            }
            cptperpos = 0;
        }
        if(todos[i][0] == 'I'){
            cptperpos++;
            todos[i] = todos[i] + string(1, 'a' + cptperpos); // these will start at 1, (the position 0 is the original letter)
        }
        lastpos = position;
    }
    if(cptperpos > 0){
        splittedStart[lastpos] = splittedStart[lastpos] + string(cptperpos, '_');
    }


    vector<string> res = {pulledprint(splittedStart)};

    if(typeShuffling != FROM_LEFT){
        random::shuffle<string>(todos);
    }

    // now if we want to keep size, every time there is an insertion, we need to either substitute or delete next, and vice versa.
    // so we use todos as a pile, and only accept the first element if not generating too much imbalance. If not we reshuffle the remaining tasks.

    if(typeShuffling == MAINTAIN_SIZE){
        int imbalance = 0;
        vector<string> remainingTodos = todos;
        todos = vector<string>();

        int cpt = 0;
        while((remainingTodos.size() > 0) && (cpt < 10000)){
            bool keepTask = true;

            if((remainingTodos[0][0] == 'I') && (imbalance >= 1)) keepTask = false;
            if((remainingTodos[0][0] == 'D') && (imbalance <= - 1)) keepTask = false;
            if((remainingTodos[0][0] == 'S') && (imbalance != 0)) keepTask = false; // refuse substition is size has changed


            if(keepTask) {
                imbalance = imbalance + ((remainingTodos[0][0] == 'I') ? 1 : 0) + ((remainingTodos[0][0] == 'D')? -1 : 0);
                todos.push_back(remainingTodos[0]);
                remainingTodos.erase(remainingTodos.begin());
            } else {
                random::shuffle<string>(remainingTodos);
            }
            cpt++;
        }
    }

    for(size_t j = 0; j < todos.size(); ++j){
        string task = todos[j];
        //cout << task;
        int position = task[1] - 'a';
        //cout << "position " << position << endl;
        if((position < 0) || (position >= start.size())) cerr << "ERR: out of range position for task " << task << endl;

        // Note, we will insert AFTER and delete BEFORE, to make sure that if there is insert then delete at the same position we are still good
        if(task[0] == 'D'){
            if(splittedStart[position].size() < 1) cerr << "ERR: removing empty position in " << task << endl;
            if(splittedStart[position][0] != task[2]) cerr << "ERR: removing the wrong AA in " << task << " from " << splittedStart[position] << endl;
            //splittedStart[position].erase(0,1); // removes the first position
            splittedStart[position][0] = '_';
        } else if(task[0] == 'I'){
            if(task.size() < 4) cerr << "ERR: insertions should have been expanded with their relative position in task " << task << endl;
            int insertPos = task[3] - 'a';
            if((insertPos < 1) || (insertPos >= splittedStart[position].size())) cerr << "ERR: relative position of insertion is out of bounds in task " << task << " for " << splittedStart[position] << endl;
            if(splittedStart[position][insertPos] != '_') cerr << "ERR: performing twice an insertion at the same place in task " << task << " for " << splittedStart[position] << endl;
            splittedStart[position][insertPos] = task[2];
        } else if(task[0] == 'S'){
            // I assume if we substitute, we don't have inserted before, I am not sure this is true. There might be a problem in the order of insertions at the same place...
            if(splittedStart[position].size() < 1) cerr << "ERR: Substituting non-single position in " << task << endl;
            splittedStart[position][0] = task[2];
        } else {
            cerr << "ERR: Unknown task " << task << endl;
        }
        //cout << "Did " << task << " got " << pulledprint(splittedStart) << endl;
        res.push_back(pulledprint(splittedStart));
    }

    return res;
}

//#include <algorithm>
string unalign(string s){
    s.erase(std::remove(s.begin(), s.end(), '_'), s.end());
    return s;
}

// https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance
template<typename T>
typename T::size_type OfficialLevenshteinDistance(const T &source, const T &target) {
    if (source.size() > target.size()) {
        return OfficialLevenshteinDistance(target, source);
    }

    using TSizeType = typename T::size_type;
    const TSizeType min_size = source.size(), max_size = target.size();
    std::vector<TSizeType> lev_dist(min_size + 1);

    for (TSizeType i = 0; i <= min_size; ++i) {
        lev_dist[i] = i;
    }

    for (TSizeType j = 1; j <= max_size; ++j) {
        TSizeType previous_diagonal = lev_dist[0], previous_diagonal_save;
        ++lev_dist[0];

        for (TSizeType i = 1; i <= min_size; ++i) {
            previous_diagonal_save = lev_dist[i];
            if (source[i - 1] == target[j - 1]) {
                lev_dist[i] = previous_diagonal;
            } else {
                lev_dist[i] = std::min(std::min(lev_dist[i - 1], lev_dist[i]), previous_diagonal) + 1;
            }
            previous_diagonal = previous_diagonal_save;
        }
    }

    return lev_dist[min_size];
}

size_t LD(string s1, string s2){
    return OfficialLevenshteinDistance<string>(s1,s2);
}

void showTest(string s1, string s2, typePath t){
    if(s1.size() > s2.size()) showTest(s2, s1, t);
    string steps = LevenshteinDistance(s1, s2);
    cout << "Final Result:" << steps << endl;
    vector<string> reconstituted = applySteps(s1, steps, t);
    cout << "Initial:" << s1 << endl;
    for(size_t i = 0; i < reconstituted.size(); ++i){
        cout << LD(unalign(s1), unalign(reconstituted[i])) << "\t" << reconstituted[i] << "\t" << unalign(reconstituted[i]) << "\t" << LD(unalign(reconstituted[i]), unalign(s2)) << endl;
    }
    cout << "Ending:" << s2 << endl;
}

// Same function but does not need to add "_" before
vector<string> inBetween(string s1, string s2, typePath t){
    vector<string> res;
    if((t == MAINTAIN_SIZE) && (s1.size() != s2.size())) {return res;}
    if(s1.size() > s2.size()) return inBetween(s2, s1, t);

    string steps = LevenshteinDistance("_" + s1, "_" + s2);
    vector<string> reconstituted = applySteps("_" + s1, steps, t);
    for(size_t i = 1; i < reconstituted.size()-1; ++i){
        string newSeq = unalign(reconstituted[i]); // note: unalign removes the "_", so also the one we introduced
        //cout << newSeq << endl;
        if((t != MAINTAIN_SIZE) || (newSeq.size() == s1.size())){
            res.push_back(newSeq);
        }
    }
    return res;
}


// will do: 2x for 1 mutation, then 1x for 2 mutations, then 10% of 2-point mutated are going to get mutated 0.1x
void generateMutantCDRs(string repertoireFile, size_t nSeqsPerMut, size_t nPairs, string outputFileName){
// --- code to generate mutant
    repertoire rep = repertoire(repertoireFile, true);
    cout << "Got " << rep.nLines() << " sequences " << endl;
    if(rep.nLines() == 0) return;
    if(rep.nLines() < 1) {
        cerr << "ERR: No sequences to process, trying next task" << endl;
    }

    repertoire outputrep = repertoire();

    //  Step 1, generate mutants, names will be ID_xm_y, x = nb mutations, y is the ID of the mutated slide
    for(size_t i = 0; i < rep.nLines(); ++i){
        pair<string, string > IDandSeq = rep.getLine(static_cast<int>(i));
        outputrep.addSequence(IDandSeq.first + "_0m_0", IDandSeq.second);

        int y = 0;
        set<string> resCore1 = randomSubset(allMutants(IDandSeq.second, 1, ""), 2 * nSeqsPerMut);
        for(set<string>::iterator it = resCore1.begin(); it != resCore1.end(); it++){
            ++y;
            stringstream makeID; makeID << IDandSeq.first << "_1m_" << y;
            outputrep.addSequence(makeID.str(), *it);
        }
        cout << i << "/" << rep.nLines() << ", " << resCore1.size() << " / " << 2* nSeqsPerMut << " asked seqs generated with 1 mutation for" << IDandSeq.second << endl;

        y = 0;
        int z = 0;
        set<string> resCore2 = randomSubset(allMutants(IDandSeq.second, 2, ""), nSeqsPerMut);
        for(set<string>::iterator it = resCore2.begin(); it != resCore2.end(); it++){
            ++y;
            {
                stringstream makeID; makeID << IDandSeq.first << "_2m_" << y;
                outputrep.addSequence(makeID.str(), *it);
            }

            // Generating 3 mutations by adding 1 mutation to 2-points. only does it on 10% of sequences
            if(y % 10 == 0){
                set<string> resCore3b = randomSubset(allMutants(*it, 1, ""), 10);
                for(set<string>::iterator it = resCore3b.begin(); it != resCore3b.end(); it++){
                    stringstream makeID; makeID << IDandSeq.first << "_3m_" << z;
                    if(!outputrep.containsID(makeID.str())){
                        ++z;
                        outputrep.addSequence(makeID.str(), *it);
                    } else {cout << "^";}
                }
            }
        }
        cout << i << "/" << rep.nLines() <<  ", " << y << " / " << nSeqsPerMut << " seqs generated with 2 mutations for" << IDandSeq.second << endl;
        cout << i << "/" << rep.nLines() << ", " << z << " seqs generated with 3 mutation for" << IDandSeq.second << endl;
        /*y = 0;
        set<string> resCore3 = randomSubset(allMutants(IDandSeq.second, 3, ""), nSeqsPerMut);
        for(set<string>::iterator it = resCore3.begin(); it != resCore3.end(); it++){
            ++y;
            stringstream makeID; makeID << IDandSeq.first << "_3m_" << y;
            outputrep.addSequence(makeID.str(), *it);
        }
        cout << i << "/" << rep.nLines() << ", " << nSeqsPerMut << " seqs generated with 3 mutations for" << IDandSeq.second << endl;*/

    }

    // Step 2, generate sequence interpolations
    // we select a few sequences and take pairs between them, so statiscally we try different ways to interpolate from the same sequences
    vector< pair<string, string > > poolSeqs = rep.getRandomLines(nPairs/2);
    int ngenerated = 0;
    for(size_t i = 0; i < nPairs; ++i){

        // Pick two random sequences in the pre-selection
        size_t selected1 = random::uniformUInteger(0, poolSeqs.size()-1);
        size_t selected2 = random::uniformUInteger(0, poolSeqs.size()-1);
        int cpt = 0;
        while((selected1 == selected2) && (cpt < 100)){
            selected2 = random::uniformUInteger(0, poolSeqs.size()-1);
            cpt++;
        }

        pair<string, string> seq1 = poolSeqs[selected1];
        pair<string, string> seq2 = poolSeqs[selected2];

        //cout << seq1.second << "," << seq2.second << endl;
        vector<string> interpolated = inBetween(seq1.second, seq2.second, RANDOM);
        for(size_t j = 0; j < interpolated.size(); ++j){
            stringstream makeID; makeID << seq1.first << "_to_" << seq2.first << "_v" << j;
            if(!outputrep.containsID(makeID.str())){
                ngenerated++;
                outputrep.addSequence(makeID.str(), interpolated[j]);
                //cout << "+";
            } else cout << "!";
        }
    }
    cout << ngenerated << " interpolated sequences from a pool of " << nPairs/2 << endl;

    outputrep.write(outputFileName);
}

vector<string> generateRandomMutantsBigHD(string initialSequence, size_t nMuts, size_t nSeqs){
    size_t L = initialSequence.size();
    if(nMuts > L) {
        cerr << "ERR: generateRandomMutantsBigHD, requested more mutations than sequence size" << endl;
        return vector<string>();
    }
    vector<string> res = vector<string>();
    for(size_t i = 0; i < nSeqs; ++i){
        // Creates a vector with 1111 (nMut time) then 0000. Then shuffles it, it will be the list of positions to mutate
        vector<bool> posMut = vector<bool>(nMuts, true);
        posMut.resize(L, false);
        random::shuffle<bool>(posMut);

        string newSeq = initialSequence;
        for(size_t j = 0; j < L; ++j){
            if(posMut[j] == true){
                int cpt = 0;
                while((newSeq[j] == initialSequence[j]) && (cpt < 100)){
                    newSeq[j] = randomAA();
                    cpt++;
                }
            }
        }
        res.push_back(newSeq);
    }
    return res;
}

void generateRepertoireMasonHD(string initialSequence, string outputFileName, size_t nTimes){


    vector<string> seqsHD6 = generateRandomMutantsBigHD(initialSequence, 6, 1240 * nTimes);
    vector<string> seqsHD7 = generateRandomMutantsBigHD(initialSequence, 7, 4480 * nTimes);
    vector<string> seqsHD8 = generateRandomMutantsBigHD(initialSequence, 8, 8533 * nTimes);
    vector<string> seqsHD9 = generateRandomMutantsBigHD(initialSequence, 9, 9047 * nTimes);
    vector<string> seqsHD10 = generateRandomMutantsBigHD(initialSequence, 10, 4037 * nTimes);

    vector<string> allSeqs;
    allSeqs.insert(allSeqs.end(), seqsHD6.begin(), seqsHD6.end());
    allSeqs.insert(allSeqs.end(), seqsHD7.begin(), seqsHD7.end());
    allSeqs.insert(allSeqs.end(), seqsHD8.begin(), seqsHD8.end());
    allSeqs.insert(allSeqs.end(), seqsHD9.begin(), seqsHD9.end());
    allSeqs.insert(allSeqs.end(), seqsHD10.begin(), seqsHD10.end());

    repertoire outputrep = repertoire();
    for(size_t i = 0; i < allSeqs.size(); ++i){
        stringstream ID; ID << i;
        outputrep.addSequence("Mut" + ID.str(), allSeqs[i]);
    }
    outputrep.write(outputFileName);
}


int main(int argc, char* argv[]){

    generateRepertoireMasonHD("LLFITTVVAPF", "mutantsMasonHDsize10.txt", 3);

    //cout << printVector(generateRandomMutantsBigHD("AAAA", 1, 10)) << endl;

    return 0;
    // Task7
    generateMutantCDRs("D:/pprobert/Datasets/Task7/1ADQ_H1_TwoMascottes.txt", 70000, 0, "D:/pprobert/Datasets/Task7/3MascotteMutants_1ADQ_H1.txt");


    cout << printVector(inBetween("ABIGKOALA", "KANGOUROU", RANDOM));
    cout << printVector(inBetween("ABIGKOALA", "KANGOUROU", MAINTAIN_SIZE));
    return 0;

// Task7, full
// generateMutantCDRs("D:/pprobert/Datasets/Task7/1ADQ_H1.txt", 25000, 0, "D:/pprobert/Datasets/Task7/8Mutants_1ADQ_H1.txt");
// generateMutantCDRs("D:/pprobert/Datasets/Task7/1ADQ_H2.txt", 25000, 0, "D:/pprobert/Datasets/Task7/8Mutants_1ADQ_H2.txt");


// Task7, glimpse
// generateMutantCDRs("D:/pprobert/Datasets/Task7/1ADQ_H1.txt", 25000, 0, "D:/pprobert/Datasets/Task7/8Mutants_1ADQ_H1.txt");
// generateMutantCDRs("D:/pprobert/Datasets/Task7/1ADQ_H2.txt", 25000, 0, "D:/pprobert/Datasets/Task7/8Mutants_1ADQ_H2.txt");


// Task 6
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1ADQ_A_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1ADQ_A_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1ADQ_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1ADQ_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1ADQ_A_Classes-H3.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1ADQ_A_InterMut-H3.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1ADQ_A_Classes-H4.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1ADQ_A_InterMut-H4.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1ADQ_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1ADQ_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1FBI_X_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1FBI_X_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1FBI_X_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1FBI_X_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1FBI_X_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1FBI_X_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1FNS_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1FNS_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1FNS_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1FNS_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1FSK_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1FSK_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1FSK_A_Classes-H1+H3.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1FSK_A_InterMut-H1+H3.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1FSK_A_Classes-H3.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1FSK_A_InterMut-H3.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1FSK_A_Classes-H4.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1FSK_A_InterMut-H4.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1FSK_A_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1FSK_A_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1FSK_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1FSK_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1H0D_C_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1H0D_C_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1H0D_C_Classes-H1+H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1H0D_C_InterMut-H1+H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1H0D_C_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1H0D_C_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1H0D_C_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1H0D_C_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1JPS_T_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1JPS_T_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1JPS_T_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1JPS_T_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1KB5_AB_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1KB5_AB_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1KB5_AB_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1KB5_AB_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1KB5_AB_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1KB5_AB_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1NCB_N_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1NCB_N_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1NSN_S_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1NSN_S_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1NSN_S_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1NSN_S_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1OAZ_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1OAZ_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1OB1_C_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1OB1_C_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1OB1_C_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1OB1_C_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1OSP_O_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1OSP_O_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1OSP_O_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1OSP_O_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1OSP_O_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1OSP_O_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1PKQ_J_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1PKQ_J_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1PKQ_J_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1PKQ_J_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1PKQ_J_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1PKQ_J_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1QFW_AB_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1QFW_AB_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1QFW_AB_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1QFW_AB_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1RJL_C_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1RJL_C_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1RJL_C_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1RJL_C_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1S78_B_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1S78_B_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1TQB_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1TQB_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1TQB_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1TQB_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1WEJ_F_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1WEJ_F_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1WEJ_F_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1WEJ_F_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1WEJ_F_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1WEJ_F_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1YJD_C_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1YJD_C_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1ZTX_E_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1ZTX_E_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1ZTX_E_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1ZTX_E_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2ARJ_RQ_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2ARJ_RQ_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2ARJ_RQ_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2ARJ_RQ_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2ARJ_RQ_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2ARJ_RQ_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2B2X_A_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2B2X_A_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2B2X_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2B2X_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2B2X_A_Classes-H3.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2B2X_A_InterMut-H3.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2B2X_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2B2X_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2FD6_AU_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2FD6_AU_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2HFG_R_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2HFG_R_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2HFG_R_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2HFG_R_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2IH3_C_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2IH3_C_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2IH3_C_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2IH3_C_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2JEL_P_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2JEL_P_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2JEL_P_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2JEL_P_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2JEL_P_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2JEL_P_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2Q8A_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2Q8A_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2Q8A_A_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2Q8A_A_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2Q8A_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2Q8A_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2R29_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2R29_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2R29_A_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2R29_A_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2R29_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2R29_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2R4R_A_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2R4R_A_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2R4R_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2R4R_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2R4R_A_Classes-H1+H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2R4R_A_InterMut-H1+H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2R4R_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2R4R_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2R56_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2R56_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2R56_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2R56_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2UZI_R_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2UZI_R_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2UZI_R_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2UZI_R_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2VXQ_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2VXQ_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2VXQ_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2VXQ_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2VXT_I_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2VXT_I_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2VXT_I_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2VXT_I_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2VXT_I_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2VXT_I_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2W9E_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2W9E_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2W9E_A_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2W9E_A_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2W9E_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2W9E_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2WUC_I_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2WUC_I_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2WUC_I_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2WUC_I_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2XQB_A_Classes-H1+H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2XQB_A_InterMut-H1+H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2XQB_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2XQB_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2XQB_A_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2XQB_A_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2XQB_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2XQB_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2XWT_C_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2XWT_C_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2YC1_C_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2YC1_C_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2YC1_C_Classes-H3.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2YC1_C_InterMut-H3.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2YC1_C_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2YC1_C_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2YC1_C_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2YC1_C_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2YPV_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2YPV_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2YPV_A_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2YPV_A_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2YPV_A_Classes-H3.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2YPV_A_InterMut-H3.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2YPV_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2YPV_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_2ZCH_P_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_2ZCH_P_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3BGF_S_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3BGF_S_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3BGF_S_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3BGF_S_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3BN9_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3BN9_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3BN9_A_Classes-H3.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3BN9_A_InterMut-H3.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3BN9_A_Classes-H4.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3BN9_A_InterMut-H4.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3BN9_A_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3BN9_A_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3BN9_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3BN9_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3CVH_ABC_Classes-H1+H3.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3CVH_ABC_InterMut-H1+H3.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3CVH_ABC_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3CVH_ABC_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3CVH_ABC_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3CVH_ABC_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3CVH_ABC_Classes-H1+H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3CVH_ABC_InterMut-H1+H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3CVH_ABC_Classes-H3.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3CVH_ABC_InterMut-H3.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3CVH_ABC_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3CVH_ABC_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3DVG_XY_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3DVG_XY_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3DVG_XY_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3DVG_XY_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3DVG_XY_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3DVG_XY_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3EFD_K_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3EFD_K_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3EFD_K_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3EFD_K_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3GI9_C_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3GI9_C_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3GI9_C_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3GI9_C_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3GI9_C_Classes-H3.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3GI9_C_InterMut-H3.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3GI9_C_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3GI9_C_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3HI6_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3HI6_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3HI6_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3HI6_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3JBQ_B_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3JBQ_B_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3JBQ_B_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3JBQ_B_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3JBQ_B_Classes-H3.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3JBQ_B_InterMut-H3.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3JBQ_B_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3JBQ_B_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3KJ4_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3KJ4_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3KJ4_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3KJ4_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3KR3_D_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3KR3_D_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3KR3_D_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3KR3_D_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3KS0_J_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3KS0_J_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3KS0_J_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3KS0_J_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3L5X_A_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3L5X_A_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3L5X_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3L5X_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3L5X_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3L5X_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3L95_X_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3L95_X_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3L95_X_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3L95_X_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3L95_X_Classes-H3.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3L95_X_InterMut-H3.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3L95_X_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3L95_X_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3MJ9_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3MJ9_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3NCY_A_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3NCY_A_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3NCY_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3NCY_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3NCY_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3NCY_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3NFP_I_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3NFP_I_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3NFP_I_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3NFP_I_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3NH7_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3NH7_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3NH7_A_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3NH7_A_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3NH7_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3NH7_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3Q3G_E_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3Q3G_E_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3Q3G_E_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3Q3G_E_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3R08_E_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3R08_E_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3R08_E_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3R08_E_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3R08_E_Classes-H1+H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3R08_E_InterMut-H1+H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3R08_E_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3R08_E_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3R1G_B_Classes-H1+H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3R1G_B_InterMut-H1+H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3R1G_B_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3R1G_B_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3R1G_B_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3R1G_B_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3R1G_B_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3R1G_B_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3RAJ_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3RAJ_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3RAJ_A_Classes-H1+H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3RAJ_A_InterMut-H1+H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3RAJ_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3RAJ_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3RAJ_A_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3RAJ_A_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3RKD_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3RKD_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3RKD_A_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3RKD_A_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3RKD_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3RKD_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3RVV_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3RVV_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3RVV_A_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3RVV_A_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3RVV_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3RVV_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3SKJ_E_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3SKJ_E_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3SKJ_E_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3SKJ_E_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3SKJ_E_Classes-H3.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3SKJ_E_InterMut-H3.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3SKJ_E_Classes-H4.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3SKJ_E_InterMut-H4.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3SKJ_E_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3SKJ_E_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3SQO_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3SQO_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3SQO_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3SQO_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3TT1_A_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3TT1_A_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3TT1_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3TT1_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3TT1_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3TT1_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3U9P_C_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3U9P_C_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3UBX_A_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3UBX_A_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3UBX_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3UBX_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3UBX_A_Classes-H3.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3UBX_A_InterMut-H3.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3UBX_A_Classes-H4.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3UBX_A_InterMut-H4.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3UBX_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3UBX_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3V6O_A_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3V6O_A_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3V6O_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3V6O_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3V6O_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3V6O_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3VG9_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3VG9_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3VG9_A_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3VG9_A_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3VG9_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3VG9_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3VRL_C_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3VRL_C_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3VRL_C_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3VRL_C_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3WD5_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3WD5_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3WD5_A_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3WD5_A_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_3WD5_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_3WD5_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4AEI_A_Classes-H1+H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4AEI_A_InterMut-H1+H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4AEI_A_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4AEI_A_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4AEI_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4AEI_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4CAD_C_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4CAD_C_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4DKE_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4DKE_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4H88_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4H88_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4H88_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4H88_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4HC1_B_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4HC1_B_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4HJ0_B_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4HJ0_B_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4HJ0_B_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4HJ0_B_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4HJ0_B_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4HJ0_B_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4I18_R_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4I18_R_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4I77_Z_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4I77_Z_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4K24_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4K24_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4K24_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4K24_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4K3J_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4K3J_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4K3J_A_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4K3J_A_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4K3J_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4K3J_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4K9E_C_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4K9E_C_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4KI5_M_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4KI5_M_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4KI5_M_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4KI5_M_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4KI5_M_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4KI5_M_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4KXZ_A_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4KXZ_A_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4KXZ_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4KXZ_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4KXZ_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4KXZ_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4LU5_B_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4LU5_B_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4MXV_B_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4MXV_B_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4N9G_C_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4N9G_C_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4N9G_C_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4N9G_C_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4N9G_C_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4N9G_C_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4NP4_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4NP4_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4NP4_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4NP4_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4OII_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4OII_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4OKV_E_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4OKV_E_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4OKV_E_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4OKV_E_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4PP1_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4PP1_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4PP1_A_Classes-H6.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4PP1_A_InterMut-H6.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4PP1_A_Classes-H1+H7.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4PP1_A_InterMut-H1+H7.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4PP1_A_Classes-H4.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4PP1_A_InterMut-H4.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4PP1_A_Classes-H3.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4PP1_A_InterMut-H3.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4PP1_A_Classes-H1+H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4PP1_A_InterMut-H1+H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4PP1_A_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4PP1_A_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4PP1_A_Classes-H7.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4PP1_A_InterMut-H7.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4PP1_A_Classes-H5.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4PP1_A_InterMut-H5.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4PP1_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4PP1_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4QCI_D_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4QCI_D_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4QCI_D_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4QCI_D_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4QEX_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4QEX_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4QNP_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4QNP_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4QWW_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4QWW_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4R9Y_D_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4R9Y_D_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4R9Y_D_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4R9Y_D_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4R9Y_D_Classes-H1+H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4R9Y_D_InterMut-H1+H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4R9Y_D_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4R9Y_D_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4RGM_S_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4RGM_S_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4U1G_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4U1G_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4U1G_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4U1G_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4U6V_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4U6V_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4U6V_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4U6V_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4WV1_F_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4WV1_F_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4WV1_F_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4WV1_F_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4WV1_F_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4WV1_F_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4Y5V_C_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4Y5V_C_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4Y5V_C_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4Y5V_C_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4YPG_C_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4YPG_C_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4YPG_C_Classes-H3.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4YPG_C_InterMut-H3.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4YPG_C_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4YPG_C_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4YPG_C_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4YPG_C_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4YUE_C_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4YUE_C_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4YUE_C_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4YUE_C_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4YUE_C_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4YUE_C_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4ZFG_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4ZFG_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4ZFG_A_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4ZFG_A_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4ZFG_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4ZFG_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4ZFO_F_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4ZFO_F_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4ZFO_F_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4ZFO_F_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4ZSO_E_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4ZSO_E_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4ZSO_E_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4ZSO_E_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_4ZSO_E_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_4ZSO_E_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5B8C_C_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5B8C_C_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5B8C_C_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5B8C_C_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5BVP_I_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5BVP_I_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5BVP_I_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5BVP_I_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5BVP_I_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5BVP_I_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5C0N_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5C0N_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5C0N_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5C0N_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5C7X_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5C7X_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5C7X_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5C7X_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5CZV_A_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5CZV_A_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5CZV_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5CZV_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5CZV_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5CZV_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5D93_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5D93_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5D93_A_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5D93_A_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5D93_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5D93_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5DFV_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5DFV_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5DFV_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5DFV_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5DHV_M_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5DHV_M_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5DHV_M_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5DHV_M_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5DHV_M_Classes-H3.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5DHV_M_InterMut-H3.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5DHV_M_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5DHV_M_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5DMI_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5DMI_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5DMI_A_Classes-H1+H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5DMI_A_InterMut-H1+H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5DMI_A_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5DMI_A_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5DMI_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5DMI_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5DO2_B_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5DO2_B_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5DO2_B_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5DO2_B_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5DO2_B_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5DO2_B_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5E8D_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5E8D_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5E8D_A_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5E8D_A_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5E8D_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5E8D_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5E8E_LH_Classes-H1+H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5E8E_LH_InterMut-H1+H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5E8E_LH_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5E8E_LH_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5E8E_LH_Classes-H3.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5E8E_LH_InterMut-H3.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5E8E_LH_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5E8E_LH_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5E8E_LH_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5E8E_LH_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5E8E_LH_Classes-H4.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5E8E_LH_InterMut-H4.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5E94_G_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5E94_G_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5E94_G_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5E94_G_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5E94_G_Classes-H1+H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5E94_G_InterMut-H1+H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5E94_G_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5E94_G_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5EII_G_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5EII_G_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5EII_G_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5EII_G_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5EII_G_Classes-H3.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5EII_G_InterMut-H3.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5EII_G_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5EII_G_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5EPM_C_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5EPM_C_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5EPM_C_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5EPM_C_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5EU7_A_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5EU7_A_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5EU7_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5EU7_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5EU7_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5EU7_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5EZO_A_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5EZO_A_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5EZO_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5EZO_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5EZO_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5EZO_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5F3B_C_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5F3B_C_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5F3B_C_Classes-H1+H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5F3B_C_InterMut-H1+H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5F3B_C_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5F3B_C_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5F3B_C_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5F3B_C_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5FB8_C_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5FB8_C_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5FB8_C_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5FB8_C_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5H35_C_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5H35_C_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5H35_C_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5H35_C_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5HDQ_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5HDQ_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5HDQ_A_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5HDQ_A_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5HDQ_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5HDQ_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5HI4_B_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5HI4_B_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5HI4_B_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5HI4_B_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5HI4_B_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5HI4_B_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5IKC_M_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5IKC_M_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5IKC_M_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5IKC_M_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5IKC_M_Classes-H3.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5IKC_M_InterMut-H3.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5IKC_M_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5IKC_M_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5J13_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5J13_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5J13_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5J13_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5JW4_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5JW4_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5JZ7_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5JZ7_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5KN5_C_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5KN5_C_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5KN5_C_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5KN5_C_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5KN5_C_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5KN5_C_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5KTE_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5KTE_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5KTE_A_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5KTE_A_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5KTE_A_Classes-H3.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5KTE_A_InterMut-H3.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5KTE_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5KTE_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5L0Q_A_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5L0Q_A_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5L0Q_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5L0Q_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5L0Q_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5L0Q_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5LQB_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5LQB_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5LQB_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5LQB_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5MES_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5MES_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5MES_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5MES_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5T5F_A_Classes-H1+H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5T5F_A_InterMut-H1+H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5T5F_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5T5F_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5T5F_A_Classes-H2.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5T5F_A_InterMut-H2.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5T5F_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5T5F_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5TH9_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5TH9_A_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5TH9_A_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5TH9_A_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5TLJ_X_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5TLJ_X_InterMut-H1.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5TLJ_X_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5TLJ_X_InterMut-Unknown.txt");
//    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_5TZ2_C_Classes-Unknown.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_5TZ2_C_InterMut-Unknown.txt");

    return 0;
    generateMutantCDRs("D:/pprobert/Datasets/Task6/T6_1ADQ_A_Classes-H1.txt", 100, 20000, "D:/pprobert/Datasets/Task6/Mutated/T6_1ADQ_A_InterMut-H1.txt");

    return 0;
    cout << printVector(inBetween("ABIGKOALA", "KALANGOUROU", RANDOM));
    //cout << printVector(inBetween("ABIGKOALA", "KANGOUROU", MAINTAIN_SIZE));
    //cout << printVector(inBetween("ABIGKOALA", "KANGOUROU", MAINTAIN_SIZE));
    //return 0;

    //showTest("CARSGYYGAMDYW", "CGREGLYVGFDYW", FROM_LEFT);
    //showTest("CARSGYYGAMDYW", "CGREGLYVGFDYW", RANDOM);
    //showTest("CARSGYYGAMDYW", "CGREGLYVGFDYW", MAINTAIN_SIZE);

    // apparently, need to start with _
    // showTest("_ABIGKOALA", "_KANGOUROU", FROM_LEFT);
    // showTest("_ABIGKOALA", "_KANGOUROU", RANDOM);
    showTest("_ABIGKOALA", "_KANGOUROU", MAINTAIN_SIZE);


    //LevenshteinDistance("MARCOU", "AMARI");
    //generateLibraryMutants();
    return 0;
    //int main(void){
        //testAAs();

        /*set<string> s1({"A", "B", "C", "D", "E", "F", "G", "H", "I"});

        set<string> s2 = randomSubset(s1, 20);
        cout << print(s2) << endl;
        set<string> s3 = randomSubset(s1, 5);
        cout << print(s3) << endl;
        return 0;*/

        //set<string> res = allMutants("CARICATURE", 2);
        //for(set<string>::iterator it = res.begin(); it != res.end(); ++it){
        //    cout << *it << "\n";
        //}

    pair<superProtein*, vector<int>> AG = getFileAntigen("1FNS_A+V48H+L70T", "ListMutantAntigens.txt");
    if(AG.first != nullptr){
        cout << AG.first->getAAseq() << endl;
    }

    return 0;
}





/*
vector<string> LevenshteinDistance(string source, string target) {
    if (source.size() > target.size()) {
        return LevenshteinDistance(target, source);
    }

    size_t min_size = source.size(), max_size = target.size();
    std::vector<size_t> lev_dist(min_size + 1);
    std::vector<string> history(min_size + 1); // for each level, all possible strings

    for (size_t i = 0; i <= min_size; ++i) {
        lev_dist[i] = i;
        history[i] = source.substr(0, i);
    }

    for (size_t  j = 1; j <= max_size; ++j) {
        size_t  previous_diagonal = lev_dist[0], previous_diagonal_save;
        ++lev_dist[0];

        for (size_t  i = 1; i <= min_size; ++i) {
            previous_diagonal_save = lev_dist[i];
            if (source[i - 1] == target[j - 1]) {
                lev_dist[i] = previous_diagonal;
                history[i] = history[i] + string(1, target[i - 1]);
            } else {
                size_t minVal = std::min(std::min(lev_dist[i - 1], lev_dist[i]), previous_diagonal) + 1;
                string possible;

                // we are at position i,
                if(minVal == lev_dist[i-1]){ // deletion
                    possible = history[i] + string(1, target[i - 1]);
                }
                if(minVal == lev_dist[i]){ // doesnt change to add?
                    possible = history[i]; // deletion
                }
                if(minVal == previous_diagonal){
                    possible = history[i-1] + string(1, target[i - 1]);
                }
                lev_dist[i] = minVal;
                cout << i << "\t" << j << "\t" << history[i] << endl;
            }
            previous_diagonal = previous_diagonal_save;
        }
    }

    return history[min_size];
}



size_t LevenshteinDistance(string source, string target) {
    if (source.size() > target.size()) {
        return LevenshteinDistance(target, source);
    }


    size_t min_size = source.size(), max_size = target.size();
    std::vector<size_t> lev_dist(min_size + 1);
    std::vector<string> history(min_size + 1);


    for (size_t i = 0; i <= min_size; ++i) {
        lev_dist[i] = i;
    }

    for (size_t  j = 1; j <= max_size; ++j) {
        size_t  previous_diagonal = lev_dist[0], previous_diagonal_save;
        ++lev_dist[0];

        for (size_t  i = 1; i <= min_size; ++i) {
            previous_diagonal_save = lev_dist[i];
            if (source[i - 1] == target[j - 1]) {
                lev_dist[i] = previous_diagonal;
            } else {
                lev_dist[i] = std::min(std::min(lev_dist[i - 1], lev_dist[i]), previous_diagonal) + 1;
            }
            previous_diagonal = previous_diagonal_save;
        }
    }

    return lev_dist[min_size];
}
*/
#ifdef longlongcode


int sizeReceptors = 10;
int minInteract = 11;

        AG.first->setAAs(AAseq);
        affinityOneLigand T1 = affinityOneLigand(AG.first, sizeReceptors, minInteract, -1, 1, AG.second);

        unsigned int actual_wait_time = stopwatch.elapsed_time<unsigned int, std::chrono::microseconds>();
        cerr << actual_wait_time / 1000. << "ms Elapsed to compute structures" << endl;

        cout << "Computing 1000 affinities for ligand " << antigenID << " with sequence " << AAseq << " as reference, receptorSize " << sizeReceptors << " minI=" << minInteract << endl;

        vector<double> listAff;
        double vmin = 0;
        double avg = 0;
        int cpt = 0;
        for(int i = 0; i < nBCRseqPerAG; ++i){
            string Px = randomProt(sizeReceptors+1);
            std::pair<double, double> res = T1.affinity(Px);
            vmin = min(vmin, res.first);
            avg += res.first;
            listAff.push_back(res.first);
            //cout << "affinity(" << Px << ") = \t" << res.first << "\t" << res.second << endl;
        }
        cout << "AG:" << AAseq << ", avg binding " << avg / (double) cpt << endl;

        unsigned int actual_wait_time2 = stopwatch.elapsed_time<unsigned int, std::chrono::microseconds>();
        cerr << actual_wait_time2 / 1000.<< " ms Elapsed to compute 1000 affinities" << endl;

        /*        vector<double> intBoundaries;
        for(double d = -120; d < -10; ++d){
            intBoundaries.push_back(d-0.5);
        }
        cout << "   -> Distributions of best energies" << endl;
        histogramFromDistrib h1(listAff, intBoundaries);
        cout << h1.print(true) << endl; */
    }

    //    cout << "Best Aff(type Best) so far:" << vmin << endl;
    //    for(int i = 0; i < 1; ++i){
    //        string Px = randomProt(sizeReceptors+1);
//        cout << "Details of the structures and affinities" << endl; // for " << simpleAccessible << " (Ag=" << AAsimple << ", BCR=" << Px << "), receptors " << receptorSize << " minI=4" << endl;
//        std::pair<double, double> res = T1->affinity(Px, true);
//        cout << "affinity(" << Px << ") = \t" << res.first << "\t" << res.second << endl;
//    }

}




/*
vector<double> distribEnergies(
AG.first->setAAs(AAseq);
affinityOneLigand T1 = affinityOneLigand(AG.first, sizeReceptors, minInteract, -1, 1, AG.second);

unsigned int actual_wait_time = stopwatch.elapsed_time<unsigned int, std::chrono::microseconds>();
cerr << actual_wait_time / 1000. << "ms Elapsed to compute structures" << endl;

cout << "Computing 1000 affinities for ligand " << antigenID << " with sequence " << AAseq << " as reference, receptorSize " << sizeReceptors << " minI=" << minInteract << endl;

vector<double> listAff;
double vmin = 0;
double avg = 0;
int cpt = 0;
for(int i = 0; i < nBCRseqPerAG; ++i){
    string Px = randomProt(sizeReceptors+1);
    std::pair<double, double> res = T1.affinity(Px);
    vmin = min(vmin, res.first);
    avg += res.first;
    listAff.push_back(res.first);
    //cout << "affinity(" << Px << ") = \t" << res.first << "\t" << res.second << endl;
}
cout << "AG:" << AAseq << ", avg binding " << avg / (double) cpt << endl;

unsigned int actual_wait_time2 = stopwatch.elapsed_time<unsigned int, std::chrono::microseconds>();
cerr << actual_wait_time2 / 1000.<< " ms Elapsed to compute 1000 affinities" << endl;

*/


//string binarySequence::testAffinityFunctions(double L, double R, int maxClusters, int typeAffinityFunction) {

//    stringstream subFolder;
//    subFolder << "L" << L << "_r" << R << nameTypeAff(typeAffinityFunction) << "_C" << maxClusters;
//    string folder = string("C:/Users/Philippe/Desktop/Sequences/") + subFolder.str() + string("/");
//    createFolder(folder);

//    //#define out cout
//    stringstream out;
//    out << "Testing the properties ot the affinity function for the following parameters : \n";
//    out << "   ->     L= " << L << "\t(Size of sequences)" << endl;
//    out << "   ->     R= " << R << "\t(specificity parameter)" << endl;
//    out << "   -> maxCl= " << maxClusters << "\t(cluster size scale)" << endl;
//    switch (typeAffinityFunction) {
//    case seqAff: {out << "   -> Using standard affinity (Saham's)\n"; break;}
//    case seqAffNorm: {out << "   -> Using standard affinity normalized by maxCl^r\n"; break;}
//    case seqAffWindow: {out << "   -> Using the maximum affinity of a sliding window\n"; break;}
//    }
//    out << "==== Part 1 : enumerates all (if possible), or a lot of sequences and sort them by affinity to get the best ones : ====" << endl;

//#define resolutiondistrib 100
//#define maxSequencesToEnumeate 500000

//    // The reference antigen is 00000...
//    binarySequence * ref = new binarySequence(L);

//    // will store a large list of random sequences, with their affinity to ref
//    vector<pair<double, binarySequence*> > store;


//    out << "1 ----------- Distribution of affinities ----------------------" << endl;
//    vector<double> distribution;
//    vector<double> logDistrib;
//    distribution.resize(resolutiondistrib + 1);

//    // in case L is big, will only sample maxSequencesToEnumeate sequences.
//    int total = 0;
//    int maxim = pow(2, L);
//    if (L > 26) {
//        maxim = maxSequencesToEnumeate + 1;           // to avoid it to become negative ...
//    }

//    // if possible to enumerate all, do it one by one.
//    bool enumerateAll = (maxim < maxSequencesToEnumeate);
//    if (enumerateAll) {
//        for (int i = 0; i < maxim; ++i) {
//            binarySequence * a = new binarySequence(L, (long) i);
//            double affi = binarySequence::affinity(a, ref, R, maxClusters, typeAffinityFunction);
//            distribution[(int) (((double) resolutiondistrib) * affi)] += 1.0; // put into the histogram
//            store.push_back(pair<double, binarySequence*> (affi,a));
//            logDistrib.push_back(log10(affi+1e-6));
//            total++;
//        }
//    } else {
//        // if not, than sample randomly
//        for (int i = 0; i < maxSequencesToEnumeate; ++i) {
//            if(((i%100000) == 0) || (i == maxSequencesToEnumeate-1)) cout << i << "/" << maxSequencesToEnumeate << endl;
//            binarySequence * a = new binarySequence(L);
//            a->randomize();
//            double affi = binarySequence::affinity(a, ref, R, maxClusters, typeAffinityFunction);
//            distribution[(int) (((double) resolutiondistrib) * affi)] += 1.0;
//            store.push_back(pair<double, binarySequence*> (affi,a));
//            logDistrib.push_back(log10(affi+1e-6));
//            total++;
//        }
//    }

//    out << "Distribution of affinities\n";
//    for (int i = 0; i < (int) distribution.size(); ++i) {
//        distribution[i] /= (double) total;
//        out << i << "\t"  << double (i) * (1.0 / (double) resolutiondistrib)
//            << "\t" << double (i + 1) * (1.0 / (double) resolutiondistrib)
//            << "\t" << distribution[i] <<  endl;
//    }
//    out << "Distribution of affinities in LOG scale\n";
//    histogramFromDistrib v(logDistrib, 100);
//    out << v.print() << endl;

//    out << "\nSequences and affinity, "
//        << ((enumerateAll) ? " in the order of ID\n" : " randomly generated\n");
//    for (int i = 0; i < 200; ++i) {
//        out << i << "\t" << store[i].second->print() << "\t" << store[i].first << "\n";
//    }

//    out << "\nSequences, sorted from the best, out of the "
//        << min(maxim,(int) maxSequencesToEnumeate) << " evaluated sequences\n";
//    std::sort(store.begin(), store.end(), compSequences);
//    for (int i = 0; i < 200; ++i) {
//        out << i << "\t" << store[i].second->print() << "\t" << store[i].first << "\n";
//    }

//    if (enumerateAll) {
//        out << "\nAffinity of sequences taken randomly\n";
//        for (int i = 0; i < 100; ++i) {
//            binarySequence * seqtmp = new binarySequence(L);
//            seqtmp->randomize();
//            out << i << "\t" << binarySequence::affinity(seqtmp,ref,R,maxClusters,typeAffinityFunction) << "\t" << seqtmp->print() << "\n";
//        }
//    }

//    out << "2 ------------------------- Mutation histograms -----------------------" << endl;

//    // Cuts the distribution of affinities in blocs of 5% of the sequences
//#define sizeClasses 0.05
//    out << "\nEqual classes of affinities (best percentiles)" << endl;
//    vector<double> classes;
//    vector<double> valuesClasses;

//    std::reverse(store.begin(), store.end()); // now to put it increasing

//    // every 5% of sequences browsed (fraction), store the value of affinity
//    double fraction = 0.00;
//    for(unsigned int i = 0; i < store.size(); ++i){
//        //out << store[i].second->print() << " Aff " << store[i].first << endl;
//        if((((double) i / (double) (store.size())) >= fraction - 1e-6) || (i == store.size() - 1)){
//            //out << "]" << fraction - sizeClasses << "," << fraction << "]\t" << store[i].first << endl;
//            valuesClasses.push_back(store[i].first);
//            classes.push_back(fraction);
//            fraction += sizeClasses;
//        }
//    }
//    out << "Affinity classes: cell at percent XX has the affinity value of: (last line=log, 0->-10)" << endl;
//    out << printVec(classes) << endl;
//    out << printVec(valuesClasses) << endl;
//    for(unsigned int i = 0; i < valuesClasses.size(); ++i){
//        out << "\t" << log10(valuesClasses[i] + 1e-10);
//    } out << endl;

//    out << " Now, for each class of sequences, makes the distribution of affinities." << endl;

//    vector<double> newAffinitiesByMutation;      // collecting new affinities only inside one class
//    vector<double> foldAffinitiesByMutation;
//    vector< vector<double>> tableAbsMutations;   // to store the histograms inside each class.
//    vector< vector<double>> tableFoldMutations;
//    vector<double> TOTALnewAffinitiesByMutation; // collecting new affinities for all sequences.
//    vector<double> TOTALfoldAffinitiesByMutation;
//    // For making histograms in term of fold induction, the following groups/classes will be used to make histograms
//    vector<double> classesFoldInd = {0,0.1,0.2,0.4,0.6,0.8,0.9,0.99,1.01,1.1,1.2,1.4,2,2.5,3.25,5,7.5,10};
//    // For making histograms in term of affinity, the classes inside valueclasses will be used.

//    classes.push_back(1.0); // to avoid seg fault

//    // Now will browse again, class by class. Values of sequences inside valuesClasses[currentClass-1] and valuesClasses[currentClass]
//    fraction = sizeClasses;
//    int currentClass = 1;

//    for(unsigned int i = 0; i < store.size(); ++i){
//        binarySequence* thisSeq = store[i].second;
//        double oldAff = binarySequence::affinity(thisSeq, ref, R, maxClusters, typeAffinityFunction);
//        if((oldAff > valuesClasses[currentClass] + 1e-6) || (oldAff < valuesClasses[currentClass-1] - 1e-6) ) cerr << "Class problems, sequence outside its class" << endl;

//        // within one class, add the new possible affinities
//        for(int j = 0; j < thisSeq->size * 2; ++j){ // might get sequence of different lengths !
//            binarySequence seq(thisSeq);
//            seq.mutateOnePosition();
//            double newAff = binarySequence::affinity(&seq, ref, R, maxClusters, typeAffinityFunction);
//            newAffinitiesByMutation.push_back(newAff);
//            if(oldAff > 0) foldAffinitiesByMutation.push_back(newAff / oldAff);
//            TOTALnewAffinitiesByMutation.push_back(newAff);
//            if(oldAff > 0) TOTALfoldAffinitiesByMutation.push_back(newAff / oldAff);
//            //out << thisSeq->print() << " Aff " << oldAff << "->" << seq.print() << " Newaff " << newAff << " Fold " << newAff / oldAff << endl;
//        }
//        if((((double) i / (double) (store.size())) >= fraction - 1e-6) || (i == store.size() - 1)){
//            //out << "Histograms for mutation for the class [" << valuesClasses[currentClass-1] << "," << valuesClasses[currentClass] << "]" << endl;
//            histogramFromDistrib resAbsForThisClass(newAffinitiesByMutation,valuesClasses);
//            tableAbsMutations.push_back(resAbsForThisClass.densities);
//            //out << "New affinities" << endl;
//            //out << resAbsForThisClass.print(true);
//            histogramFromDistrib resFoldForThisClass(foldAffinitiesByMutation,classesFoldInd);
//            tableFoldMutations.push_back(resFoldForThisClass.densities);
//            //out << "Fold increase in the affinity" << endl;
//            //out << resFoldForThisClass.print(true);
//            //out << endl << endl << endl;
//            newAffinitiesByMutation.clear();
//            foldAffinitiesByMutation.clear();
//            currentClass++;
//            fraction += sizeClasses;
//        }

//    }
//    out << "Outputing the affinity changes from all sequences" << endl;
//    histogramFromDistrib TOTALresAbsForThisClass(TOTALnewAffinitiesByMutation,valuesClasses);
//    out << "New affinities" << endl;
//    out << TOTALresAbsForThisClass.print(true);
//    histogramFromDistrib TOTALresFoldForThisClass(TOTALfoldAffinitiesByMutation,classesFoldInd);
//    out << "Fold increase in the affinity" << endl;
//    out << TOTALresFoldForThisClass.print(true);

//    out << "Outputing the affinity changes (new affinity value), from each class of affinity" << endl;
//    out << "\t";
//    for(unsigned int j = 0; j < tableAbsMutations.size(); ++j){
//        out << "[" << valuesClasses[j] << "," << valuesClasses[j+1] << "]\t";
//    }
//    out << "\n";
//    for(unsigned int i = 0; i < TOTALresAbsForThisClass.densities.size(); ++i){
//        out << "[" << TOTALresAbsForThisClass.lowBoundsXs[i] << "," << TOTALresAbsForThisClass.highBoundsXs[i] << "]";
//        for(unsigned int j = 0; j < tableAbsMutations.size(); ++j){
//            out << "\t" << tableAbsMutations[j][i];
//        }
//        out << endl;
//    }

//    out << "Outputing the affinity changes (new affinity value), from each class of affinity" << endl;
//    out << "\t";
//    for(unsigned int j = 0; j < tableFoldMutations.size(); ++j){
//        out << "[" << valuesClasses[j] << "," << valuesClasses[j+1] << "]\t";
//    }
//    out << "\n";
//    for(unsigned int i = 0; i < TOTALresFoldForThisClass.densities.size(); ++i){
//        out << "[" << TOTALresFoldForThisClass.lowBoundsXs[i] << "," << TOTALresFoldForThisClass.highBoundsXs[i] << "]";
//        for(unsigned int j = 0; j < tableFoldMutations.size(); ++j){
//            out << "\t" << tableFoldMutations[j][i];
//        }
//        out << endl;
//    }




//    out << "==== Part 2 : Evaluating cross-reactivity in the system : ====" << endl;

//    for(int k = 0; k < 10; ++k){

//        binarySequence ref2(ref);
//        ref2.mutateOnePosition();
//        ref2.mutateOnePosition();

//        stringstream fname;
//        fname << folder << "DotPlot2CloseAG-L" << L << "_r" << R << nameTypeAff(typeAffinityFunction) << "_C" << maxClusters << "rep" << k << "_" << ref->print() << "to" << ref2.print() << ".txt";
//        fstream of(fname.str().c_str(), ios::out);

//        vector<double> affRef;
//        vector<double> affRef2;

//        for(int i = 0; i < 20000; ++i){
//            binarySequence test(ref2);
//            test.randomize();
//            double aff1 = binarySequence::affinity(&test, ref, R, maxClusters, typeAffinityFunction);
//            double aff2 = binarySequence::affinity(&test, &ref2, R, maxClusters, typeAffinityFunction);
//            affRef.push_back(aff1);
//            affRef2.push_back(aff2);
//            of << aff1 << "\t" << aff2 << endl;
//        }
//        of.close();
//    }

//    for(int k = 0; k < 10; ++k){

//        binarySequence ref2(ref);
//        ref2.mutateOnePosition();
//        ref2.mutateOnePosition();
//        ref2.mutateOnePosition();
//        ref2.mutateOnePosition();
//        ref2.mutateOnePosition();
//        ref2.mutateOnePosition();
//        ref2.mutateOnePosition();

//        stringstream fname;
//        fname << folder << "DotPlot6CloseAG-L" << L << "_r" << R << nameTypeAff(typeAffinityFunction) << "_C" << maxClusters << "rep" << k << "_" << ref->print() << "to" << ref2.print() << ".txt";
//        fstream of(fname.str().c_str(), ios::out);

//        vector<double> affRef;
//        vector<double> affRef2;

//        for(int i = 0; i < 20000; ++i){
//            binarySequence test(ref2);
//            test.randomize();
//            double aff1 = binarySequence::affinity(&test, ref, R, maxClusters, typeAffinityFunction);
//            double aff2 = binarySequence::affinity(&test, &ref2, R, maxClusters, typeAffinityFunction);
//            affRef.push_back(aff1);
//            affRef2.push_back(aff2);
//            of << aff1 << "\t" << aff2 << endl;
//        }
//        of.close();
//    }

//    for(int k = 0; k < 10; ++k){
//        binarySequence ref2(ref);
//        ref2.randomize();

//        stringstream fname;
//        fname << folder << "DotPlotRandCloseAG-L" << L << "_r" << R << nameTypeAff(typeAffinityFunction) << "_C" << maxClusters << "rep" << k << "_" << ref->print() << "to" << ref2.print() << ".txt";
//        fstream of(fname.str().c_str(), ios::out);


//        vector<double> affRef;
//        vector<double> affRef2;

//        for(int i = 0; i < 20000; ++i){
//            binarySequence test(ref2);
//            test.randomize();
//            double aff1 = binarySequence::affinity(&test, ref, R, maxClusters, typeAffinityFunction);
//            double aff2 = binarySequence::affinity(&test, &ref2, R, maxClusters, typeAffinityFunction);
//            affRef.push_back(aff1);
//            affRef2.push_back(aff2);
//            of << aff1 << "\t" << aff2 << endl;
//        }
//        of.close();
//    }

//    ofstream ffin(folder + string("Output.txt"));
//    ffin << out.str();
//    ffin.close();
//    return out.str();


//    int nbAntigens = 10;
//    out << "Generating randomly " << nbAntigens << " antigens " << endl;
//    vector<binarySequence*> ags;
//    for (int i = 0; i < nbAntigens; ++i) {
//        binarySequence * seq = new binarySequence(L);
//        seq->randomize();
//        ags.push_back(seq);
//        out << "\tAg nr " << i << "\t" << seq->print() << endl;
//    }
//    out << "\nNumber of antigens recognized by randomly generated sequences, based on threshold\n";

//    out << "  -> (for the first 100 sequences : ) In the case of random sequences" << endl;
//    total = 0;
//#define thresholdRecoAg 0.1
//    int nbDiscardedSeq = 0;  // sequences that don't recognize anything
//    int countprint = 0;
//    for (int k = 0; k < min(maxim, (int) maxSequencesToEnumeate); ++k) {
//        if (k == 100) {
//            out
//                    <<
//                       "  -> (for the remaining sequences) for sequences recognizing at least an antigen with affinity 0.1"
//                    << endl;
//        }
//        total++;

//        // for each sequence,
//        bool recoAtLeastOne = false;
//        vector<double> nbRecoDepThresh(10, 0.0);
//        vector<double> affinityEach(nbAntigens, 0.0);
//        binarySequence * seqtmp = new binarySequence(L);
//        seqtmp->randomize();
//        for (int j = 0; j < nbAntigens; ++j) {
//            double thisAff = binarySequence::affinity(seqtmp,
//                                                      ags[j],
//                                                      R,
//                                                      maxClusters,
//                                                      typeAffinityFunction);
//            if ((thisAff > thresholdRecoAg) || (k < 100)) {
//                recoAtLeastOne = true;
//            } else { nbDiscardedSeq++; }
//            affinityEach[j] = thisAff;
//            for (int i = 0; i <= (int) (9.99 * thisAff); ++i) {
//                if (i < 10) { nbRecoDepThresh[i]++; }
//            }
//        }
//        if (recoAtLeastOne && (countprint < 5000)) {
//            countprint++;
//            out << "RandSeq " << k << ", " << seqtmp->print() << " ";
//            out << "nbAgPerThreshold:";
//            for (int i = 0; i < 10; ++i) {
//                out << "\t" << nbRecoDepThresh[i];
//            }
//            out << "\taffPerAg:";
//            for (int i = 0; i < nbAntigens; ++i) {
//                out << "\t" << affinityEach[i];
//            }
//            out << endl;
//        }
//        delete seqtmp;
//    }
//    out << "   ... Nb of sequences analyzed: " << total << endl;
//    out << "   ... Nb of sequences discarded: " << nbDiscardedSeq
//        << "(except the 100 first ones, i.e. among the :" << total - 100 << " remaining)" << endl;

//    out << "==== Part 3 : Evaluating the effect of mutations : ====" << endl;

//    binarySequence * start = new binarySequence(L);    // starting by '0000' : the best sequence
//    out << "NbMut\tsequence\taffinity\n";
//    for (int i = 0; i < 2 * L; ++i) {
//        out << i << "\t" << start->print() << "\t" << binarySequence::affinity(start,
//                                                                               ref,
//                                                                               R,
//                                                                               maxClusters,
//                                                                               typeAffinityFunction)
//            << endl;
//        start->mutateOnePosition();
//    }

//    out << "\tReaching a good affinity\t" << endl;
//    start->randomize();
//    double prevaff = binarySequence::affinity(start, ref, R, maxClusters, typeAffinityFunction);

//    bool stop = false;
//    for (int i = 0; (i < L) && (!stop); ++i) {
//        out << "sequence : " << start->print() << "\tAff:\t" << prevaff << "\t";
//        out << "PossibleMut:";
//        vector<int> posGoodMutations;
//        for (int i = 0; i < L; ++i) {
//            binarySequence stmp = binarySequence(start);
//            stmp.content[i] = !stmp.content[i];
//            double newaff
//                    = binarySequence::affinity(&stmp, ref, R, maxClusters, typeAffinityFunction);
//            out << "\t" << newaff;
//            if (newaff > prevaff) { posGoodMutations.push_back(i); }
//        }
//        out << endl;
//        if (posGoodMutations.size() > 0) {
//            int nextmut = random::uniformInteger(0,posGoodMutations.size() - 1);
//            start->content[posGoodMutations[nextmut]] = !start->content[posGoodMutations[nextmut]];
//            prevaff = binarySequence::affinity(start, ref, R, maxClusters, typeAffinityFunction);
//        } else {
//            stop = true;
//        }
//    }


//    for (int i = 0; i < (int) store.size(); ++i) {
//        delete store[i].second;
//    }


//    out << "\tReaching a good affinity\t" << endl;
//    start->randomize();
//    /*double prevaff = binarySequence::seq_affinity(start, ref, R, maxClusters, typeAffinityFunction);

//       bool stop = false;
//       for (int i = 0; (i < L) && (!stop); ++i) {
//          out << "sequence : " << start->print() << "\tAff:\t" << prevaff << "\t";
//          out << "PossibleMut:";
//          vector<int> posGoodMutations;
//          for (int i = 0; i < L; ++i) {
//             binarySequence stmp = binarySequence(start);
//             stmp.content[i] = !stmp.content[i];
//             double newaff
//                = binarySequence::seq_affinity(&stmp, ref, R, maxClusters, typeAffinityFunction);
//             out << "\t" << newaff;
//             if (newaff > prevaff) { posGoodMutations.push_back(i); }
//          }
//          out << endl;
//          if (posGoodMutations.size() > 0) {
//             int nextmut = irandom(posGoodMutations.size() - 1);
//             start->content[posGoodMutations[nextmut]] = !start->content[posGoodMutations[nextmut]];
//             prevaff = binarySequence::seq_affinity(start, ref, R, maxClusters, typeAffinityFunction);
//          } else {
//             stop = true;
//          }
//       }*/



//    out << "2D plot of affinities to 2 antigens" << endl;
//    out << "Case 1 : two close antigens, one mutation away" << endl;

//    return out.str();
//}







//#include "quality.h"
//#include <vector>
//#include <string>
//#include <iostream>
//#include <sstream>
//#include <fstream>
//#include <antigenLib.h>
//#include "../Tools/stopwatch.h"
//using namespace std;

//void testQuality(){
//    evaluateAntigen("L10paper2AGs(V1a)", true, "101101111101011111100000000000100000000", 2);
//}

//// starts from an AA sequence, test the quality, then looks as N-points mutants and looks at the quality as well,
//// Mask = 0100011 with 1 are positions that can be mutated.
//// note: the random sequence is NOT following the mask, it's completely independent of the original one
//void evaluateAntigen(string antigenID, bool startFromRandomAAseq, string maskMut, int NbMuts){
//    cout << "Now real mask" << endl;
//    for(int KL = 0; KL < 100; ++KL){
//        precise_stopwatch stopwatch; // careful, doesnt work for more than 2 hours

//        std::pair<superProtein*, vector<int> > AG = getAntigen(antigenID);

//        string AAseq = AG.first->getAAseq();
//        size_t nAAs = AAseq.size();
//        if(maskMut.size() != nAAs){
//            cerr << "ERR: evaluateAntigen(" << antigenID << ", mask= " << maskMut << " is not same size as AA sequence from antigen " << AAseq << endl;
//        }
//        if(startFromRandomAAseq){
//            string newProt = randomProt(nAAs);
//            for(size_t i = 0; i < newProt.size(); ++i){
//                if(maskMut[i] == '1'){
//                    AAseq[i] = newProt[i];
//                }
//            }
//            AG.first->setAAs(AAseq);
//        }

//        int sizeReceptors = 8;
//        int minInteract = 4;
//        affinityOneLigand T1 = affinityOneLigand(AG.first, sizeReceptors, minInteract, -1, 1, AG.second);

//        unsigned int actual_wait_time = stopwatch.elapsed_time<unsigned int, std::chrono::microseconds>();
//        cerr << actual_wait_time / 1000. << "ms Elapsed to compute structures" << endl;

//        cout << "Computing 1000 affinities for ligand " << antigenID << " with sequence " << AAseq << " as reference, receptorSize " << sizeReceptors << " minI=" << minInteract << endl;

//        vector<double> listAff;
//        double vmin = 0;
//        double avg = 0;
//        int cpt = 0;
//        for(int i = 0; i < 200; ++i){
//            string Px = randomProt(sizeReceptors+1);
//            std::pair<double, double> res = T1.affinity(Px);
//            vmin = min(vmin, res.first);
//            avg += res.first;
//            listAff.push_back(res.first);
//            cout << "affinity(" << Px << ") = \t" << res.first << "\t" << res.second << endl;
//        }
//        cout << "AG:" << AAseq << ", avg binding " << avg / (double) cpt << endl;

//        unsigned int actual_wait_time2 = stopwatch.elapsed_time<unsigned int, std::chrono::microseconds>();
//        cerr << actual_wait_time2 / 1000.<< " ms Elapsed to compute 1000 affinities" << endl;

///*        vector<double> intBoundaries;
//        for(double d = -120; d < -10; ++d){
//            intBoundaries.push_back(d-0.5);
//        }
//        cout << "   -> Distributions of best energies" << endl;
//        histogramFromDistrib h1(listAff, intBoundaries);
//        cout << h1.print(true) << endl; */
//    }

////    cout << "Best Aff(type Best) so far:" << vmin << endl;
////    for(int i = 0; i < 1; ++i){
////        string Px = randomProt(sizeReceptors+1);
////        cout << "Details of the structures and affinities" << endl; // for " << simpleAccessible << " (Ag=" << AAsimple << ", BCR=" << Px << "), receptors " << receptorSize << " minI=4" << endl;
////        std::pair<double, double> res = T1->affinity(Px, true);
////        cout << "affinity(" << Px << ") = \t" << res.first << "\t" << res.second << endl;
////    }

//}




////string binarySequence::testAffinityFunctions(double L, double R, int maxClusters, int typeAffinityFunction) {

////    stringstream subFolder;
////    subFolder << "L" << L << "_r" << R << nameTypeAff(typeAffinityFunction) << "_C" << maxClusters;
////    string folder = string("C:/Users/Philippe/Desktop/Sequences/") + subFolder.str() + string("/");
////    createFolder(folder);

////    //#define out cout
////    stringstream out;
////    out << "Testing the properties ot the affinity function for the following parameters : \n";
////    out << "   ->     L= " << L << "\t(Size of sequences)" << endl;
////    out << "   ->     R= " << R << "\t(specificity parameter)" << endl;
////    out << "   -> maxCl= " << maxClusters << "\t(cluster size scale)" << endl;
////    switch (typeAffinityFunction) {
////    case seqAff: {out << "   -> Using standard affinity (Saham's)\n"; break;}
////    case seqAffNorm: {out << "   -> Using standard affinity normalized by maxCl^r\n"; break;}
////    case seqAffWindow: {out << "   -> Using the maximum affinity of a sliding window\n"; break;}
////    }
////    out << "==== Part 1 : enumerates all (if possible), or a lot of sequences and sort them by affinity to get the best ones : ====" << endl;

////#define resolutiondistrib 100
////#define maxSequencesToEnumeate 500000

////    // The reference antigen is 00000...
////    binarySequence * ref = new binarySequence(L);

////    // will store a large list of random sequences, with their affinity to ref
////    vector<pair<double, binarySequence*> > store;


////    out << "1 ----------- Distribution of affinities ----------------------" << endl;
////    vector<double> distribution;
////    vector<double> logDistrib;
////    distribution.resize(resolutiondistrib + 1);

////    // in case L is big, will only sample maxSequencesToEnumeate sequences.
////    int total = 0;
////    int maxim = pow(2, L);
////    if (L > 26) {
////        maxim = maxSequencesToEnumeate + 1;           // to avoid it to become negative ...
////    }

////    // if possible to enumerate all, do it one by one.
////    bool enumerateAll = (maxim < maxSequencesToEnumeate);
////    if (enumerateAll) {
////        for (int i = 0; i < maxim; ++i) {
////            binarySequence * a = new binarySequence(L, (long) i);
////            double affi = binarySequence::affinity(a, ref, R, maxClusters, typeAffinityFunction);
////            distribution[(int) (((double) resolutiondistrib) * affi)] += 1.0; // put into the histogram
////            store.push_back(pair<double, binarySequence*> (affi,a));
////            logDistrib.push_back(log10(affi+1e-6));
////            total++;
////        }
////    } else {
////        // if not, than sample randomly
////        for (int i = 0; i < maxSequencesToEnumeate; ++i) {
////            if(((i%100000) == 0) || (i == maxSequencesToEnumeate-1)) cout << i << "/" << maxSequencesToEnumeate << endl;
////            binarySequence * a = new binarySequence(L);
////            a->randomize();
////            double affi = binarySequence::affinity(a, ref, R, maxClusters, typeAffinityFunction);
////            distribution[(int) (((double) resolutiondistrib) * affi)] += 1.0;
////            store.push_back(pair<double, binarySequence*> (affi,a));
////            logDistrib.push_back(log10(affi+1e-6));
////            total++;
////        }
////    }

////    out << "Distribution of affinities\n";
////    for (int i = 0; i < (int) distribution.size(); ++i) {
////        distribution[i] /= (double) total;
////        out << i << "\t"  << double (i) * (1.0 / (double) resolutiondistrib)
////            << "\t" << double (i + 1) * (1.0 / (double) resolutiondistrib)
////            << "\t" << distribution[i] <<  endl;
////    }
////    out << "Distribution of affinities in LOG scale\n";
////    histogramFromDistrib v(logDistrib, 100);
////    out << v.print() << endl;

////    out << "\nSequences and affinity, "
////        << ((enumerateAll) ? " in the order of ID\n" : " randomly generated\n");
////    for (int i = 0; i < 200; ++i) {
////        out << i << "\t" << store[i].second->print() << "\t" << store[i].first << "\n";
////    }

////    out << "\nSequences, sorted from the best, out of the "
////        << min(maxim,(int) maxSequencesToEnumeate) << " evaluated sequences\n";
////    std::sort(store.begin(), store.end(), compSequences);
////    for (int i = 0; i < 200; ++i) {
////        out << i << "\t" << store[i].second->print() << "\t" << store[i].first << "\n";
////    }

////    if (enumerateAll) {
////        out << "\nAffinity of sequences taken randomly\n";
////        for (int i = 0; i < 100; ++i) {
////            binarySequence * seqtmp = new binarySequence(L);
////            seqtmp->randomize();
////            out << i << "\t" << binarySequence::affinity(seqtmp,ref,R,maxClusters,typeAffinityFunction) << "\t" << seqtmp->print() << "\n";
////        }
////    }

////    out << "2 ------------------------- Mutation histograms -----------------------" << endl;

////    // Cuts the distribution of affinities in blocs of 5% of the sequences
////#define sizeClasses 0.05
////    out << "\nEqual classes of affinities (best percentiles)" << endl;
////    vector<double> classes;
////    vector<double> valuesClasses;

////    std::reverse(store.begin(), store.end()); // now to put it increasing

////    // every 5% of sequences browsed (fraction), store the value of affinity
////    double fraction = 0.00;
////    for(unsigned int i = 0; i < store.size(); ++i){
////        //out << store[i].second->print() << " Aff " << store[i].first << endl;
////        if((((double) i / (double) (store.size())) >= fraction - 1e-6) || (i == store.size() - 1)){
////            //out << "]" << fraction - sizeClasses << "," << fraction << "]\t" << store[i].first << endl;
////            valuesClasses.push_back(store[i].first);
////            classes.push_back(fraction);
////            fraction += sizeClasses;
////        }
////    }
////    out << "Affinity classes: cell at percent XX has the affinity value of: (last line=log, 0->-10)" << endl;
////    out << printVec(classes) << endl;
////    out << printVec(valuesClasses) << endl;
////    for(unsigned int i = 0; i < valuesClasses.size(); ++i){
////        out << "\t" << log10(valuesClasses[i] + 1e-10);
////    } out << endl;

////    out << " Now, for each class of sequences, makes the distribution of affinities." << endl;

////    vector<double> newAffinitiesByMutation;      // collecting new affinities only inside one class
////    vector<double> foldAffinitiesByMutation;
////    vector< vector<double>> tableAbsMutations;   // to store the histograms inside each class.
////    vector< vector<double>> tableFoldMutations;
////    vector<double> TOTALnewAffinitiesByMutation; // collecting new affinities for all sequences.
////    vector<double> TOTALfoldAffinitiesByMutation;
////    // For making histograms in term of fold induction, the following groups/classes will be used to make histograms
////    vector<double> classesFoldInd = {0,0.1,0.2,0.4,0.6,0.8,0.9,0.99,1.01,1.1,1.2,1.4,2,2.5,3.25,5,7.5,10};
////    // For making histograms in term of affinity, the classes inside valueclasses will be used.

////    classes.push_back(1.0); // to avoid seg fault

////    // Now will browse again, class by class. Values of sequences inside valuesClasses[currentClass-1] and valuesClasses[currentClass]
////    fraction = sizeClasses;
////    int currentClass = 1;

////    for(unsigned int i = 0; i < store.size(); ++i){
////        binarySequence* thisSeq = store[i].second;
////        double oldAff = binarySequence::affinity(thisSeq, ref, R, maxClusters, typeAffinityFunction);
////        if((oldAff > valuesClasses[currentClass] + 1e-6) || (oldAff < valuesClasses[currentClass-1] - 1e-6) ) cerr << "Class problems, sequence outside its class" << endl;

////        // within one class, add the new possible affinities
////        for(int j = 0; j < thisSeq->size * 2; ++j){ // might get sequence of different lengths !
////            binarySequence seq(thisSeq);
////            seq.mutateOnePosition();
////            double newAff = binarySequence::affinity(&seq, ref, R, maxClusters, typeAffinityFunction);
////            newAffinitiesByMutation.push_back(newAff);
////            if(oldAff > 0) foldAffinitiesByMutation.push_back(newAff / oldAff);
////            TOTALnewAffinitiesByMutation.push_back(newAff);
////            if(oldAff > 0) TOTALfoldAffinitiesByMutation.push_back(newAff / oldAff);
////            //out << thisSeq->print() << " Aff " << oldAff << "->" << seq.print() << " Newaff " << newAff << " Fold " << newAff / oldAff << endl;
////        }
////        if((((double) i / (double) (store.size())) >= fraction - 1e-6) || (i == store.size() - 1)){
////            //out << "Histograms for mutation for the class [" << valuesClasses[currentClass-1] << "," << valuesClasses[currentClass] << "]" << endl;
////            histogramFromDistrib resAbsForThisClass(newAffinitiesByMutation,valuesClasses);
////            tableAbsMutations.push_back(resAbsForThisClass.densities);
////            //out << "New affinities" << endl;
////            //out << resAbsForThisClass.print(true);
////            histogramFromDistrib resFoldForThisClass(foldAffinitiesByMutation,classesFoldInd);
////            tableFoldMutations.push_back(resFoldForThisClass.densities);
////            //out << "Fold increase in the affinity" << endl;
////            //out << resFoldForThisClass.print(true);
////            //out << endl << endl << endl;
////            newAffinitiesByMutation.clear();
////            foldAffinitiesByMutation.clear();
////            currentClass++;
////            fraction += sizeClasses;
////        }

////    }
////    out << "Outputing the affinity changes from all sequences" << endl;
////    histogramFromDistrib TOTALresAbsForThisClass(TOTALnewAffinitiesByMutation,valuesClasses);
////    out << "New affinities" << endl;
////    out << TOTALresAbsForThisClass.print(true);
////    histogramFromDistrib TOTALresFoldForThisClass(TOTALfoldAffinitiesByMutation,classesFoldInd);
////    out << "Fold increase in the affinity" << endl;
////    out << TOTALresFoldForThisClass.print(true);

////    out << "Outputing the affinity changes (new affinity value), from each class of affinity" << endl;
////    out << "\t";
////    for(unsigned int j = 0; j < tableAbsMutations.size(); ++j){
////        out << "[" << valuesClasses[j] << "," << valuesClasses[j+1] << "]\t";
////    }
////    out << "\n";
////    for(unsigned int i = 0; i < TOTALresAbsForThisClass.densities.size(); ++i){
////        out << "[" << TOTALresAbsForThisClass.lowBoundsXs[i] << "," << TOTALresAbsForThisClass.highBoundsXs[i] << "]";
////        for(unsigned int j = 0; j < tableAbsMutations.size(); ++j){
////            out << "\t" << tableAbsMutations[j][i];
////        }
////        out << endl;
////    }

////    out << "Outputing the affinity changes (new affinity value), from each class of affinity" << endl;
////    out << "\t";
////    for(unsigned int j = 0; j < tableFoldMutations.size(); ++j){
////        out << "[" << valuesClasses[j] << "," << valuesClasses[j+1] << "]\t";
////    }
////    out << "\n";
////    for(unsigned int i = 0; i < TOTALresFoldForThisClass.densities.size(); ++i){
////        out << "[" << TOTALresFoldForThisClass.lowBoundsXs[i] << "," << TOTALresFoldForThisClass.highBoundsXs[i] << "]";
////        for(unsigned int j = 0; j < tableFoldMutations.size(); ++j){
////            out << "\t" << tableFoldMutations[j][i];
////        }
////        out << endl;
////    }




////    out << "==== Part 2 : Evaluating cross-reactivity in the system : ====" << endl;

////    for(int k = 0; k < 10; ++k){

////        binarySequence ref2(ref);
////        ref2.mutateOnePosition();
////        ref2.mutateOnePosition();

////        stringstream fname;
////        fname << folder << "DotPlot2CloseAG-L" << L << "_r" << R << nameTypeAff(typeAffinityFunction) << "_C" << maxClusters << "rep" << k << "_" << ref->print() << "to" << ref2.print() << ".txt";
////        fstream of(fname.str().c_str(), ios::out);

////        vector<double> affRef;
////        vector<double> affRef2;

////        for(int i = 0; i < 20000; ++i){
////            binarySequence test(ref2);
////            test.randomize();
////            double aff1 = binarySequence::affinity(&test, ref, R, maxClusters, typeAffinityFunction);
////            double aff2 = binarySequence::affinity(&test, &ref2, R, maxClusters, typeAffinityFunction);
////            affRef.push_back(aff1);
////            affRef2.push_back(aff2);
////            of << aff1 << "\t" << aff2 << endl;
////        }
////        of.close();
////    }

////    for(int k = 0; k < 10; ++k){

////        binarySequence ref2(ref);
////        ref2.mutateOnePosition();
////        ref2.mutateOnePosition();
////        ref2.mutateOnePosition();
////        ref2.mutateOnePosition();
////        ref2.mutateOnePosition();
////        ref2.mutateOnePosition();
////        ref2.mutateOnePosition();

////        stringstream fname;
////        fname << folder << "DotPlot6CloseAG-L" << L << "_r" << R << nameTypeAff(typeAffinityFunction) << "_C" << maxClusters << "rep" << k << "_" << ref->print() << "to" << ref2.print() << ".txt";
////        fstream of(fname.str().c_str(), ios::out);

////        vector<double> affRef;
////        vector<double> affRef2;

////        for(int i = 0; i < 20000; ++i){
////            binarySequence test(ref2);
////            test.randomize();
////            double aff1 = binarySequence::affinity(&test, ref, R, maxClusters, typeAffinityFunction);
////            double aff2 = binarySequence::affinity(&test, &ref2, R, maxClusters, typeAffinityFunction);
////            affRef.push_back(aff1);
////            affRef2.push_back(aff2);
////            of << aff1 << "\t" << aff2 << endl;
////        }
////        of.close();
////    }

////    for(int k = 0; k < 10; ++k){
////        binarySequence ref2(ref);
////        ref2.randomize();

////        stringstream fname;
////        fname << folder << "DotPlotRandCloseAG-L" << L << "_r" << R << nameTypeAff(typeAffinityFunction) << "_C" << maxClusters << "rep" << k << "_" << ref->print() << "to" << ref2.print() << ".txt";
////        fstream of(fname.str().c_str(), ios::out);


////        vector<double> affRef;
////        vector<double> affRef2;

////        for(int i = 0; i < 20000; ++i){
////            binarySequence test(ref2);
////            test.randomize();
////            double aff1 = binarySequence::affinity(&test, ref, R, maxClusters, typeAffinityFunction);
////            double aff2 = binarySequence::affinity(&test, &ref2, R, maxClusters, typeAffinityFunction);
////            affRef.push_back(aff1);
////            affRef2.push_back(aff2);
////            of << aff1 << "\t" << aff2 << endl;
////        }
////        of.close();
////    }

////    ofstream ffin(folder + string("Output.txt"));
////    ffin << out.str();
////    ffin.close();
////    return out.str();


////    int nbAntigens = 10;
////    out << "Generating randomly " << nbAntigens << " antigens " << endl;
////    vector<binarySequence*> ags;
////    for (int i = 0; i < nbAntigens; ++i) {
////        binarySequence * seq = new binarySequence(L);
////        seq->randomize();
////        ags.push_back(seq);
////        out << "\tAg nr " << i << "\t" << seq->print() << endl;
////    }
////    out << "\nNumber of antigens recognized by randomly generated sequences, based on threshold\n";

////    out << "  -> (for the first 100 sequences : ) In the case of random sequences" << endl;
////    total = 0;
////#define thresholdRecoAg 0.1
////    int nbDiscardedSeq = 0;  // sequences that don't recognize anything
////    int countprint = 0;
////    for (int k = 0; k < min(maxim, (int) maxSequencesToEnumeate); ++k) {
////        if (k == 100) {
////            out
////                    <<
////                       "  -> (for the remaining sequences) for sequences recognizing at least an antigen with affinity 0.1"
////                    << endl;
////        }
////        total++;

////        // for each sequence,
////        bool recoAtLeastOne = false;
////        vector<double> nbRecoDepThresh(10, 0.0);
////        vector<double> affinityEach(nbAntigens, 0.0);
////        binarySequence * seqtmp = new binarySequence(L);
////        seqtmp->randomize();
////        for (int j = 0; j < nbAntigens; ++j) {
////            double thisAff = binarySequence::affinity(seqtmp,
////                                                      ags[j],
////                                                      R,
////                                                      maxClusters,
////                                                      typeAffinityFunction);
////            if ((thisAff > thresholdRecoAg) || (k < 100)) {
////                recoAtLeastOne = true;
////            } else { nbDiscardedSeq++; }
////            affinityEach[j] = thisAff;
////            for (int i = 0; i <= (int) (9.99 * thisAff); ++i) {
////                if (i < 10) { nbRecoDepThresh[i]++; }
////            }
////        }
////        if (recoAtLeastOne && (countprint < 5000)) {
////            countprint++;
////            out << "RandSeq " << k << ", " << seqtmp->print() << " ";
////            out << "nbAgPerThreshold:";
////            for (int i = 0; i < 10; ++i) {
////                out << "\t" << nbRecoDepThresh[i];
////            }
////            out << "\taffPerAg:";
////            for (int i = 0; i < nbAntigens; ++i) {
////                out << "\t" << affinityEach[i];
////            }
////            out << endl;
////        }
////        delete seqtmp;
////    }
////    out << "   ... Nb of sequences analyzed: " << total << endl;
////    out << "   ... Nb of sequences discarded: " << nbDiscardedSeq
////        << "(except the 100 first ones, i.e. among the :" << total - 100 << " remaining)" << endl;

////    out << "==== Part 3 : Evaluating the effect of mutations : ====" << endl;

////    binarySequence * start = new binarySequence(L);    // starting by '0000' : the best sequence
////    out << "NbMut\tsequence\taffinity\n";
////    for (int i = 0; i < 2 * L; ++i) {
////        out << i << "\t" << start->print() << "\t" << binarySequence::affinity(start,
////                                                                               ref,
////                                                                               R,
////                                                                               maxClusters,
////                                                                               typeAffinityFunction)
////            << endl;
////        start->mutateOnePosition();
////    }

////    out << "\tReaching a good affinity\t" << endl;
////    start->randomize();
////    double prevaff = binarySequence::affinity(start, ref, R, maxClusters, typeAffinityFunction);

////    bool stop = false;
////    for (int i = 0; (i < L) && (!stop); ++i) {
////        out << "sequence : " << start->print() << "\tAff:\t" << prevaff << "\t";
////        out << "PossibleMut:";
////        vector<int> posGoodMutations;
////        for (int i = 0; i < L; ++i) {
////            binarySequence stmp = binarySequence(start);
////            stmp.content[i] = !stmp.content[i];
////            double newaff
////                    = binarySequence::affinity(&stmp, ref, R, maxClusters, typeAffinityFunction);
////            out << "\t" << newaff;
////            if (newaff > prevaff) { posGoodMutations.push_back(i); }
////        }
////        out << endl;
////        if (posGoodMutations.size() > 0) {
////            int nextmut = random::uniformInteger(0,posGoodMutations.size() - 1);
////            start->content[posGoodMutations[nextmut]] = !start->content[posGoodMutations[nextmut]];
////            prevaff = binarySequence::affinity(start, ref, R, maxClusters, typeAffinityFunction);
////        } else {
////            stop = true;
////        }
////    }


////    for (int i = 0; i < (int) store.size(); ++i) {
////        delete store[i].second;
////    }


////    out << "\tReaching a good affinity\t" << endl;
////    start->randomize();
////    /*double prevaff = binarySequence::seq_affinity(start, ref, R, maxClusters, typeAffinityFunction);

////       bool stop = false;
////       for (int i = 0; (i < L) && (!stop); ++i) {
////          out << "sequence : " << start->print() << "\tAff:\t" << prevaff << "\t";
////          out << "PossibleMut:";
////          vector<int> posGoodMutations;
////          for (int i = 0; i < L; ++i) {
////             binarySequence stmp = binarySequence(start);
////             stmp.content[i] = !stmp.content[i];
////             double newaff
////                = binarySequence::seq_affinity(&stmp, ref, R, maxClusters, typeAffinityFunction);
////             out << "\t" << newaff;
////             if (newaff > prevaff) { posGoodMutations.push_back(i); }
////          }
////          out << endl;
////          if (posGoodMutations.size() > 0) {
////             int nextmut = irandom(posGoodMutations.size() - 1);
////             start->content[posGoodMutations[nextmut]] = !start->content[posGoodMutations[nextmut]];
////             prevaff = binarySequence::seq_affinity(start, ref, R, maxClusters, typeAffinityFunction);
////          } else {
////             stop = true;
////          }
////       }*/



////    out << "2D plot of affinities to 2 antigens" << endl;
////    out << "Case 1 : two close antigens, one mutation away" << endl;

////    return out.str();
////}





functions from quality.cpp
 * // mask=1 = mutable
string mutateNPositions(string content, int nMut, string mask) {
    size_t size = content.size();
    if(mask.size() != size) {
        cerr << "ERR:mutateOnePosition, wrong mask size " << endl;
        return string("Wrong Mask");
    }
   int cpt = 0;
   int success = 0;
   while(cpt < 25){
       char test = randomAA();
       size_t position = random::uniformInteger(0,size-1);
       if ((position >= size) || (position < 0)) {
          cerr << "Random does shit, sequencespace.cpp::mutateOnePosition";
          return content; // fail
       }
       if((test != content[position]) && (mask[position] == '1')){
            content[position] = test;
            success++;
            mask[position] = '0'; // do not mutate again there
            if(success == nMut){
                return content; // success
            }
       }
       cpt++;
   }
   cerr << "ERR: foldedFree::mutateOnePosition, got trouble to mutate an AA" << endl;
   return string("MutateFailed"); // silent
}


// starts from an AA sequence, test the quality, then looks as N-points mutants and looks at the quality as well,
// Mask = 0100011 with 1 are positions that can be mutated.
// note: the random sequence is NOT following the mask, it's completely independent of the original one
// use NbMuts < 0 tto completely randomize each time
void evaluateAntigen(string antigenID, int nrAGs, int nBCRseqPerAG, bool startFromRandomAAseq, int NbMuts, string maskMut, string initialSeq){
    cout << "Now real mask" << endl;

    precise_stopwatch stopwatch; // careful, doesnt work for more than 2 hours

    std::pair<superProtein*, vector<int> > AG = getAntigen(antigenID);

    string AAseq = AG.first->getAAseq();

    // if gave initial AG sequence, takes it
    if(initialSeq.size() > 0) AAseq = initialSeq;
    size_t nAAs = AAseq.size();

    // checks if mask has good size. If not put  everywhere
    if(maskMut.size() != nAAs){
        cerr << "ERR: evaluateAntigen(" << antigenID << ", mask= " << maskMut << " is not same size as AA sequence from antigen " << AAseq << endl;
        maskMut = string(nAAs, '1');
    }

    // if want random sequence, randomizes only mutable positions
    if(startFromRandomAAseq){
        string newProt = randomProt(nAAs);
        for(size_t i = 0; i < newProt.size(); ++i){
            if(maskMut[i] == '1'){
                AAseq[i] = newProt[i];
            }
        }
    }

    int sizeReceptors = 10;
    int minInteract = 11;

    // store the original AA sequence
    string originalSeq = AAseq;


    for(int KL = 0; KL < nrAGs; ++KL){
        // the first time, use the original sequence
        if(KL == 0){
            AAseq = originalSeq;
        } else {
            // each time want a new sequence, will do mutations around
            if(NbMuts > 0)
                AAseq = mutateNPositions(originalSeq, NbMuts, maskMut);
            else {
                // randomizes all mutable positions
                string newProt = randomProt(nAAs);
                for(size_t i = 0; i < newProt.size(); ++i){
                    if(maskMut[i] == '1'){
                        AAseq[i] = newProt[i];
                    }
                }
            }
        }


        AG.first->setAAs(AAseq);
        affinityOneLigand T1 = affinityOneLigand(AG.first, sizeReceptors, minInteract, -1, 1, AG.second);

        unsigned int actual_wait_time = stopwatch.elapsed_time<unsigned int, std::chrono::microseconds>();
        cerr << actual_wait_time / 1000. << "ms Elapsed to compute structures" << endl;

        cout << "Computing 1000 affinities for ligand " << antigenID << " with sequence " << AAseq << " as reference, receptorSize " << sizeReceptors << " minI=" << minInteract << endl;

        vector<double> listAff;
        double vmin = 0;
        double avg = 0;
        int cpt = 0;
        for(int i = 0; i < nBCRseqPerAG; ++i){
            string Px = randomProt(sizeReceptors+1);
            std::pair<double, double> res = T1.affinity(Px);
            vmin = min(vmin, res.first);
            avg += res.first;
            listAff.push_back(res.first);
            //cout << "affinity(" << Px << ") = \t" << res.first << "\t" << res.second << endl;
        }
        cout << "AG:" << AAseq << ", avg binding " << avg / (double) cpt << endl;

        unsigned int actual_wait_time2 = stopwatch.elapsed_time<unsigned int, std::chrono::microseconds>();
        cerr << actual_wait_time2 / 1000.<< " ms Elapsed to compute 1000 affinities" << endl;

//        vector<double> intBoundaries;
//        for(double d = -120; d < -10; ++d){
//            intBoundaries.push_back(d-0.5);
//        }
//        cout << "   -> Distributions of best energies" << endl;
//        histogramFromDistrib h1(listAff, intBoundaries);
//        cout << h1.print(true) << endl;
    }

    //    cout << "Best Aff(type Best) so far:" << vmin << endl;
    //    for(int i = 0; i < 1; ++i){
    //        string Px = randomProt(sizeReceptors+1);
//        cout << "Details of the structures and affinities" << endl; // for " << simpleAccessible << " (Ag=" << AAsimple << ", BCR=" << Px << "), receptors " << receptorSize << " minI=4" << endl;
//        std::pair<double, double> res = T1->affinity(Px, true);
//        cout << "affinity(" << Px << ") = \t" << res.first << "\t" << res.second << endl;
//    }

}

#endif
