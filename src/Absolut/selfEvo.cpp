#include "selfEvo.h"
#include "../Ymir/ymir.h"
#include "../Tools/nucleotides.h"
#include "../Tools/zaprandom.h"
#include "motifFeatures.h"
#include "antigenLib.h"
#include "../Ymir/proteins.h"
#include "../Ymir/plot3d.h"

// checks if a position is surrounded by least 4 positions around it in the same plane (donut)
bool isDonut(int position, set<int>& occupiedPositions){
    vector<int> pos0 = lattice::positionFromID(position);
    int pos1 = lattice::idFromPosisition({pos0[0] + 1, pos0[1]    , pos0[2]});
    int pos2 = lattice::idFromPosisition({pos0[0] - 1, pos0[1]    , pos0[2]});
    int pos3 = lattice::idFromPosisition({pos0[0]    , pos0[1] + 1, pos0[2]});
    int pos4 = lattice::idFromPosisition({pos0[0]    , pos0[1] - 1, pos0[2]});
    int pos5 = lattice::idFromPosisition({pos0[0]    , pos0[1]    , pos0[2] + 1});
    int pos6 = lattice::idFromPosisition({pos0[0]    , pos0[1]    , pos0[2] - 1});
    bool t1 = (occupiedPositions.find(pos1) != occupiedPositions.end());
    bool t2 = (occupiedPositions.find(pos2) != occupiedPositions.end());
    bool t3 = (occupiedPositions.find(pos3) != occupiedPositions.end());
    bool t4 = (occupiedPositions.find(pos4) != occupiedPositions.end());
    bool t5 = (occupiedPositions.find(pos5) != occupiedPositions.end());
    bool t6 = (occupiedPositions.find(pos6) != occupiedPositions.end());
    if(t1 && t2 && t3 && t4) return true;
    if(t1 && t2 && t5 && t6) return true;
    if(t3 && t4 && t5 && t6) return true;
    return false;
}
//  technique to fill the embarrassing points by forbidden positions.
//  alreadyBlocked starts empty, and each recursive call will fill the embarassing points (i.e. dead ends = points with 5 contacts)
//  then the function is called again as long as it finds more points to fill. the return will be the list of embarrassing points.
set<int> listEmbarrasingPoints(superProtein* ligand, set<int> alreadyBlocked, bool silent){
    if(!silent) cout << "listEmbarrasingPoints function called with " << alreadyBlocked.size() << " elements" << endl;
    if(ligand->structure == nullptr) cerr << "ERR: Choucroute garnie, listEmbarrassing points called with a superProtein of NULL structure" << endl;
    set<int> enlargedPoints = alreadyBlocked;
    set<int> largerProtein = union_sets(ligand->structure->occupiedPositions, alreadyBlocked);
    set<int> neigh = neighborPositions(largerProtein);

    for(set<int>::iterator it = neigh.begin(); it != neigh.end(); ++it){
        int nContacts = nbTouchPoints(largerProtein, *it);
        bool goingThrough = isDonut(*it, ligand->structure->occupiedPositions);
        if((nContacts >= 5) || ((nContacts >= 4) && goingThrough)){
            if(!silent) cout << "Position " << *it << " -> " << printVector(lattice::positionFromID(*it)) << " has " << nContacts << " neighbors " << endl;
            enlargedPoints.insert(*it);
        }
    }
    if(enlargedPoints.size() != alreadyBlocked.size()){
        return listEmbarrasingPoints(ligand, enlargedPoints, silent);
    }
    if(!silent) {
        cout << "Finished layer with " << enlargedPoints.size() << " elements " << endl;
        cout << "blockV = {";
        int i = 0;
        for(set<int>::iterator it = enlargedPoints.begin(); it != enlargedPoints.end(); ++it){
            if(i > 0) cout << ", ";
            cout << *it;
            i++;
        }
        cout << "};" << endl;
    }
    return enlargedPoints;
}

void testEmbarrassing(){

//    superProtein ligand;
//    ligand = superProtein();

//    return;

    superProtein* P1 = new superProtein("SUSRRDDRRUUDDRLULDRLRDRRUURRLLSUUDLDUUDUDDRUDURDUUDRRULRRDLDRLRLRLULRDRLLDSDLUDSSLURRLDULDSSL", 133152);
    superProtein* P2 = new superProtein(insert(P1, "BDSUURRLLRRDLDDRURSLDRRDDLRRDDURUSDDLDLSSLDSRLSUDSSLLDUDDRRUDLRSSSLLURUSUUDSUUDDSURSLDURDLLRL", 132903, 1000));
    P2->setAAs("VVKFMDVYQRSYCHPIETLVDIFQEYPDEIEYIFKPSCVPLMRCGGCCNDEGLECVPTEESNITMQIMRIKPHQGQHIGEMSFLQHNKCECRPKVVKFMDVYQRSYCHPIETLVDIFQEYPDEIEYIFKPSCVPLMRCGGCCNDEGLECVPTEESNITMQIMRIKPHQGQHIGEMSFLQHNKCECRPK");
    affinityOneLigand T1 = affinityOneLigand(P2, 4, 8, -1, 1);

    T1.setUltraFast(false);
    vector<string> bestStructures;
    std::pair<double, double> res = T1.affinity("ACDFH", false, &bestStructures);
    cout << "Got:" << printVector(bestStructures);

    set<int> test = listEmbarrasingPoints(P2);
    cout << print(test) << endl;

#ifdef ALLOW_GRAPHICS
    glDisplay();
    addToDisplay(P2);
    addToDisplay(new set<int>(test)); // note: do not give test, it will be destroyed when the function ends
    glutMainLoop();
    return;
#endif
}

// try batch possible mutations, keeps the best provided it's an increase.
// stops when no mutation increases (try 5x batch)
// gets the receptor size from affinityOneligand
string oneStochasticNuclMutation(affinityOneLigand* T1, string refSequence = "", int batch = 20){
    // will also make the graph of mutations (possible each time?)
    int size = (T1->sizeReceptors+1)*3;
    string currentSeq = ((refSequence.size() > 0) ? refSequence : clearStops(randomDNA(size)));
    if(currentSeq.size() != static_cast<size_t>(size)) cerr << "ERR: oneStochasticNuclMutation, the sequence does not match the size of receptors inside T1" << endl;

    double currentAff = T1->affinity(convertToProtein(currentSeq)).first;
    bool stop = false;
    while(!stop){
        double bestAffBatch = 1e6;
        string bestBatch = currentSeq;
        for(int i = 0; i < batch; ++i){
            string test = mutateDNA(currentSeq, true, true);
            string asProt = convertToProtein(test);
            //cout << "Prot seq: " << asProt << endl;
            double newAff = T1->affinity(asProt).first;
            if(newAff < bestAffBatch){
                bestBatch = test;
                bestAffBatch = newAff;
            }
        }
        if(bestAffBatch >= currentAff) stop = true; // ie. if it is equal or worst
        else { // do not taje bestAffBatch if it's worse ...
            //cout << currentSeq << "\t" << convertToProtein(currentSeq) << "\t" << currentAff << endl;
            currentAff = bestAffBatch;
            currentSeq = bestBatch;
        }
    }
    return currentSeq;
}

// test all possible mutations around a sequence, and takes it.
// this is working as AA sequences
string oneOptimalMutation(affinityOneLigand* T1, string refSequence = ""){
    // will also "L3paper" the graph of mutations (possible each time?)
    int size = (T1->sizeReceptors+1);
    string currentSeq = ((refSequence.size() > 0) ? refSequence : convertToProtein(clearStops(randomDNA(3*size))));
    if(currentSeq.size() != static_cast<size_t>(size)) cerr << "ERR: oneStochasticNuclMutation, the sequence does not match the size of receptors inside T1" << endl;
    double currentAff = T1->affinity(currentSeq).first;
    bool stop = false;
    while(!stop){
        // one round
        double bestAffBatch = 1e6;
        string bestBatch = currentSeq;
        vector<char> AAs = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};
        vector<int> pos(size, 0);
        std::iota(pos.begin(), pos.end(), 0.);
        //cout << printVector(pos);
        random::shuffle(pos);
        for(size_t k = 0; k < size; ++k){
            for(size_t j = 0; j < 20; ++j){
                size_t i = pos[k];
                string test = currentSeq;
                test[i] = AAs[j];
                double newAff = T1->affinity(test).first;
                if(newAff < bestAffBatch){
                    bestBatch = test;
                    bestAffBatch = newAff;
                }
            }
        }
        if(bestAffBatch >= currentAff) stop = true; // ie. if it is equal or worst
        else {
            cout << currentSeq << "\t" << currentAff << endl;
            currentAff = bestAffBatch;
            currentSeq = bestBatch;
        }
    }
    cout << currentSeq << "\t" << T1->affinity(currentSeq).first << endl;
    return currentSeq;
}

void testIncreasingMutations(){

    //std::pair<superProtein*, vector<int> > res = getAntigen("L3paper");
    string antigenName = "L3paper";
    int sizeRec = 7;
    affinityOneLigand* T1 = getAffinityAntigen(antigenName, sizeRec, sizeRec*2);

    string AAseq;
    cout << "Test randomly increasing DNA mutations around ligand " << antigenName << " (function oneStochasticNuclMutation)" << endl;
    cout << "Startaff " << T1->affinity(string(sizeRec+1, 'L')).first << endl;

    string AllLeu;
    for(int k = 0; k < sizeRec+1;++k){
        AllLeu = AllLeu + string("CTC");
    }
    AAseq = convertToProtein(oneStochasticNuclMutation(T1, AllLeu, 50)); // try with only leucine
    for(int i = 0; i < 25; ++i){
        AAseq = convertToProtein(oneStochasticNuclMutation(T1, "", 50));
        cout << "Result=" << AAseq << "\t" << T1->affinity(AAseq).first << endl;
    }
    cout << "Test all optimal mutations one by one to increasing affinity around ligand " << antigenName << " (function oneOptimalMutation)" << endl;
    for(int i = 0; i < 25; ++i){
        AAseq = oneOptimalMutation(T1);
        cout << "Result=" << AAseq << "\t" << T1->affinity(AAseq).first << endl;
    }

//   set<int> block = struct3D("SSSSSUUSSSS", UnDefined, lattice::idFromPosisition(31,34,33)).occupiedPositions;
//   for(set<int>::iterator it = block.begin(); it != block.end(); ++it){
//       cout << *it << "\t";
//   }
//   set<int> block2 = struct3D("RSSSUUSSS", UnDefined, lattice::idFromPosisition(36,33,33)).occupiedPositions;
//   for(set<int>::iterator it = block.begin(); it != block.end(); ++it){
//       cout << *it << "\t";
//   }
}

// from features.h
// double affinityCodeTot(string receptor, string codeInters){
// string interactions(struct3D& ligand, struct3D & s2, string AAsequenceLigand);
// codeSelfInteractions()
// string codeSelfInteractions(struct3D& s);
// need a function for affinity between protein and structure.
//double affinity(superProtein* ligand, string recept){
//    string code = interactions(ligand->structure, struct3D(recept), string AAsequenceLigand);

//}


double selfTotAffinity(string currentSeq, string AAseqToFold){
    struct3D currentStruct = struct3D(currentSeq);
    if(!currentStruct.properlyFolded) return +1e6;
    string codeSelf = codeSelfInteractions(currentStruct);
    double currentAff = affinityCodeTot(AAseqToFold, codeSelf);
    return currentAff;
}

string stochasticSelfFolding(int sizeReceptors, string AAseqToFold, string baseStructure = ""){
    size_t size = static_cast<size_t>(sizeReceptors);
    string currentSeq = ((baseStructure.size() > 0) ? baseStructure : string(size, 'S'));
    if(currentSeq.size() != static_cast<size_t>(size)) cerr << "ERR: oneStochasticNuclMutation, the sequence does not match the size of receptors inside T1" << endl;

    double currentAff = selfTotAffinity(currentSeq, AAseqToFold);

    bool stop = false;
    while(!stop){
        // one round
        double bestAffBatch = 1e6;
        string bestBatch = currentSeq;
        vector<char> Muts = {'S','U','D','R','L'};
        for(size_t i = 0; i < size; ++i){
            for(size_t j = 0; j < 20; ++j){
                string test = currentSeq;
                test[i] = Muts[j];
                double newAff = selfTotAffinity(test, AAseqToFold);
                if(newAff < bestAffBatch){
                    bestBatch = test;
                    bestAffBatch = newAff;
                }
            }
        }
        if(bestAffBatch >= currentAff) stop = true; // ie. if it is equal or worst
        cout << currentSeq << "\t" << currentAff << endl;
        currentAff = bestAffBatch;
        currentSeq = bestBatch;
    }
    return currentSeq;
}

//
string selfFolding(){

    return string();
}
