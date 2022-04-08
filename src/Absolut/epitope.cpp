#include "epitope.h"
#include "fileformats.h"
#include <map>
#include <algorithm>
#include "../Ymir/compact.h"
#include "../Ymir/plot3d.h"
#include "../Ymir/receptorligand.h"

#include "motifFeatures.h"
#include "poolstructs.h"
#include "antigenLib.h"
#include <sstream>



static vector< std::map< set<int>, set< set<int> > *  >  > dynamicMemory;

void clearMemory(){
    for(size_t i = 0; i < dynamicMemory.size(); ++i){
        for(std::map< set<int>, set< set<int> > *  >::iterator it = dynamicMemory[i].begin(); it != dynamicMemory[i].end(); ++it){
            delete (it->second);
        }
    }
    dynamicMemory.clear();
}

// Aim: input = positions bound by a binding, output: all subsets of these positions of size k
// note, if the input is smaller than k, returns only one set with the positions, and fills remaining with -1, -2, -3...,
// so at the end every binding has at least one set (later they can be discarded).
set<set<int>> generateSubsets(set<int> & input, int k){


    //cout << " Called k=" << k << " set=" << print(input) << endl;
    set<set<int>> res;
    int S = static_cast<int>(input.size());
    if(k < 1) return res;
    if(k == 1){
        for(set<int>::iterator it = input.begin(); it != input.end(); ++it){
            set<int> Uni = {*it};
            res.insert(Uni);
        }
        return res;
    }
    if(S == k) {
        res.insert(input);
        return res;
    }
    if(S < k){
        set<int> toAugment = input;
        for(int i = 0; i < (k-S); ++i){
            toAugment.insert(-1-i);
        }
        res.insert(toAugment);
        return res;
    }
//    if(S - k > 10){
//        cerr << "WRN: generating subsets of size " << k << " from a set of size " << S << " Is likely to take forever ... " << endl;
//    }


    // Make sure to not recompute twice the same subsets (would become exponential in the stack)
    // for each k separately
    // defined outside static vector< std::map< set<int>, set< set<int> > *  >  > dynamicMemory;
    dynamicMemory.resize(1000);
    std::map<set<int>, set<set<int>>* >::iterator found = dynamicMemory.at(static_cast<size_t>(k)).find(input);
    if(found != dynamicMemory.at(static_cast<size_t>(k)).end()){
        //cout << "Recall" << endl;
        return *(found->second);
    }



    for(set<int>::iterator it = input.begin(); it != input.end(); ++it){
        set<int> minusOneElement = input;
        minusOneElement.erase(*it);

        set< set<int> > recursiveSmallerSets = generateSubsets(minusOneElement, k-1);
        for(set< set<int> >::iterator j = recursiveSmallerSets.begin(); j != recursiveSmallerSets.end(); ++j){
            set<int> got = *j;
            got.insert(*it);
            //cout << "Before " << res.size() << "\t";
//            if(res.size() == 1){
//                set<int> bef = *(res.begin());
//                cout << print(bef);
//            }
            res.insert(got);
            //cout << " insert " << print(got) << endl;
            //cout << "After " << res.size() << "\n";

        }
    }
    //cout << " Level " << S << " return " << res.size() << " elements" << endl;

    // storing in memory - needs a limitation because memory explodes easily
    if(res.size() < 10000){
        set<set<int>>* res4memory = new set<set<int>>(res);
        dynamicMemory.at(static_cast<size_t>(k))[input] = res4memory;
    }
    if(res.size() > 50000) cerr << "¤";
    else if(res.size() > 10000) cerr << "!";
    else if(res.size() > 100) cerr << ".";
    return res;
}

inline int Factorial(int x) {
  return (x == 1 ? x : (x < 1 ? 1 : x * Factorial(x - 1)));
}

void testGenerateSubsets(){
    set<int> S1 = {1,3,5,8,9,10, 13, 14, 15}; //, 200, 3000, 4000};
    int S = static_cast<int>(S1.size());
    for(int k = 1; k < S+1; ++k){
        cout << "Generating subsets of size " << k << " of set " << print(S1) << endl;
        set< set<int> > res = generateSubsets(S1, k);
        for(set< set<int> > ::iterator it = res.begin(); it != res.end(); ++it){
            set<int> got = *it;
            cout << "Got: " << print(got) << "\n";
        }
        int expected = Factorial(S)/(Factorial(k)*Factorial(S-k));
        if(expected != res.size()) cerr << "ERR: generateSubsets failed to get the good number of subsets " << endl;
        cout << "Expected nb: " << expected << " and got " << res.size() << " subsets " << endl;
        cout << endl;
    }
}




// comparison function for sorting
bool compSetsToSeqs(std::pair< set<int>, set<string>> a, std::pair< set<int>, set<string>> b){
    return a.second.size() < b.second.size();
}

vector< std::pair<string, string> > setCoveringStructures( dataset< analyzedBinding >& annotatedDataset, int sizeSet){
    // aim: implement set cover that cover all binding structures with sets of k positions on the antigen.

    // First, for each binding, compute all positions subsets of size k of its binding positions,
    // then make two maps: subset => (remaining) Structures that cover it
    //                and: structure => bool assigned to set (or not yet)

    // For each set, the list of AA sequences that fit into this subset;
    std::map< set<int>, set<string> > seqsPerSubset;

    //std::map< string, bool > sequenceIsAssociatedToSubset;

    size_t S = annotatedDataset.nLines();
    for(size_t i = 0; i < S; ++i){
        // try to avoid to copy the structure binding, so directly get the fields.

        analyzedBinding* elem = annotatedDataset.getElement(i);
        string AAseq = elem->AAseq;
        set<int> bindingPos = elem->bindingAGresiduesID; //bindingAGPositionsInSpace;

//        // 1- start a list of sequences saying if they have been treated
//        std::map< string, bool >::iterator foundSeq = sequenceIsAssociatedToSubset.find(AAseq);
//        if(foundSeq != sequenceIsAssociatedToSubset.end()){
//            cout << "Note: sequence " << AAseq << " is given multiple times" << endl;
//        } else {
//            sequenceIsAssociatedToSubset[AAseq] = false;
//        }

        // 2- Generate all possible subsets of the bound positions, and associate the AA sequence to all these subsets
        // cerr << "Generating K=" << sizeSet << "-subsets from " <<  bindingPos.size() << " elements" << endl;
        cerr << "Treating " << elem->AAseq << "->" << print(bindingPos) << endl;
        set< set<int> > possibleSubsets = generateSubsets(bindingPos, sizeSet);
        for(set< set<int> >::iterator it = possibleSubsets.begin(); it != possibleSubsets.end(); ++it){
            std::map< set<int>, set<string> >::iterator found = seqsPerSubset.find(*it);
            if(found != seqsPerSubset.end()){
                found->second.insert(AAseq);
                // note: if the sequence is already there we don't care, it will not be inserted twice (set)
            } else {
                set<string> Uni = {AAseq};
                seqsPerSubset[*it] = Uni;
            }
        }
        clearMemory();
    }

    // 3- Now, we are in the minilal set cover problem: find a minimal number of sets (from the list) that cover all sequences.
    // Here, we take a greedy algorithm to start: take the best set, assign all its structures to him, then remove all those structures from the sets,
    // then sort again, to cover the remaining structures.


    // Output: the list of subsets and their sequences, that would cover the most sequences
    vector< std::pair< set<int>, set<string>> > winners;


    // Let's do it iteratively, later will be a while loop
    int maxCoveringSets = 25;
    for(int ci = 0; ci < maxCoveringSets; ++ci){

        cout << " =========== Round number " << ci << " ===============" << endl;
        // Stopping condition: We will discard subsets, and at the end there will be no more subsets in the map.
        if(seqsPerSubset.size() == 0) break;

        // 3a- sort the subset by number of sequences they cover
        // trick: transforms into a vector in order to use the std::sort function => subsets sorted by the number of sequences they explain
        vector< std::pair< set<int>, set<string>> > vecSorted;

        // copy key-value pairs from the map to the vector
        std::copy(seqsPerSubset.begin(), seqsPerSubset.end(), std::back_inserter<std::vector<std::pair< set<int>, set<string>>>>(vecSorted));

        // Sort according to the number of structures that the subset explains... (see compSetsToSeqs function)
        std::sort(vecSorted.begin(), vecSorted.end(),compSetsToSeqs);

        // Print, for debugging
        //for (auto const &pair: vecSorted) {
        //    set<int> subset = pair.first;
        //    std::cout << '{' << print(subset) << "," << pair.second.size() << '}' << '\n';
       // }

        // 3b- Take the best subset,
        winners.push_back(vecSorted.back());
        set<int> bestSet = vecSorted.back().first;
        set<string> coveredSequences = vecSorted.back().second;
        seqsPerSubset.erase(bestSet);

        // 3c- Remove the covered sequences from all other subsets.
        // trick: browse the list of sets inside the vector, but when empty, removes the sets from the dictionnary,
        //        it wouldnt work to remove elements from a structure (the map or the vector) while looping of it.
        size_t nS = vecSorted.size();
        for(size_t i = 0; i < nS; ++i){
            set<int> subset = vecSorted[i].first;
            if(seqsPerSubset.size() > 0){
                std::map< set<int>, set<string> >::iterator it = seqsPerSubset.find(subset);
                if(it != seqsPerSubset.end()){
                    for(set<string>::iterator it2 = coveredSequences.begin(); it2 != coveredSequences.end(); ++it2){
                        if(it->second.find(*it2) != it->second.end()){
                            it->second.erase(*it2);
                        }
                    }
                    if(it->second.empty()){
                        seqsPerSubset.erase(subset);
                    }
                }
            }
        }
        cout << " =========== Round Finished ===============" << endl;
    }


    cout << "Greedy covering of " << sizeSet << " antigen positions " << endl;
    vector< std::pair<string, string> > res;
    for(size_t i = 0; i < winners.size(); ++i){
        set<int> subset = winners[i].first;
        set<string> sequences = winners[i].second;
        stringstream pulledSeqs;
        pulledSeqs << sequences.size();
        for(set<string>::iterator it = sequences.begin(); it != sequences.end(); ++it){
            pulledSeqs << " " << *it;
        }
        cout << print(subset) << " => " << sequences.size() << endl; // << pulledSeqs.str() << endl;
        res.push_back(std::pair<string, string>(print(subset), pulledSeqs.str()));
    }
    cerr << "Covering finished" << endl;
    return res;
}



std::map<int, int> posAntigens(superProtein *P){
    std::map<int, int> res;
    if(!P) return res;
    int S = P->size();
    for(int i = 0; i < S; ++i){
        int ID = (*P)[i].IDresidue;
        int pos = (*P)[i].IDposition;
        //cout << ID << "=>" << pos << endl;
        if(res.find(ID) != res.end()){
            cerr << "ERR: posAntigens, a protein residue ID is defined twice: " << ID << endl;
        } else {
            res.insert(std::pair<int,int>(ID, pos));
        }
    }
    return res;
}




dataset<analyzedBinding> exampleAnalyzedDatasetFor1FBI(){
    cerr << "Start generating dataset " << endl;
    // Here is a curated example for antigen FB1, of a list of different high affinity sequences binding on it,
    // I removed many sequences with the same binding pattern, to be shorter and with equilibrated number of sequences for each pattern. Just for testing,
    // in real life, some patterns will be much more frequent.
    string exData = ""
                    "	CARLLSSITTVRDYW	CARLLSSITTV	-82.5	133154	RLUSDLSDUL	"
                    "	CARSLIGYYYAMDYW	ARSLIGYYYAM	-81.56	133154	RLUSDLSDUL	"
                    "	CARSLIGYYYAMDYW	ARSLIGYYYAM	-81.56	133154	SRUSLUSLRU	"
                    "	CARIFITTVVAKYFDVW	IFITTVVAKYF	-80.59	137187	UDDLLRRLDU	"
                    "	CARLLLRRYYFDYW	LLLRRYYFDYW	-82.29	137187	UDDLLRRLDU	"
                    "	CTRLLLRLYYFDYW	LLLRLYYFDYW	-84.39	137187	UDDLLRRLDU	"
                    "	CTRWLLRSHWYFDVW	LLRSHWYFDVW	-82.32	137187	UDDLLRRLDU	"
                    "	CARHEDRGIYYDYEWAMDYW	IYYDYEWAMDY	-83.34	137187	UDRUSUDRDU	"
                    "	CARHEESLITTVVATDYFDYW	LITTVVATDYF	-88.74	137187	UDRUSUDRDU	"
                    "	CARHEGGYYYGKLAMDYW	GYYYGKLAMDY	-81.33	137187	UDRUSUDRDU	"
                    "	CARDKGIYYGNYVGYAMDYW	IYYGNYVGYAM	-84.15	137187	UDRUSUDRDU	"
                    "	CARDLALLGQGMDYW	DLALLGQGMDY	-83.47	137187	UDRUSUDRDU	"
                    "	CARDLYGSRGDWYFDVW	DLYGSRGDWYF	-82.5	137187	UDRUSUDRDU	"
                    "	CARDLYYGSNDWYFDVW	DLYYGSNDWYF	-83.71	137187	UDRUSUDRDU	"
                    "	CARDLYYGSNDWYFDVW	YYGSNDWYFDV	-80.88	137187	UDRUSUDRDU	"
                    "	CARDPFYGSSYAYWYFDVW	PFYGSSYAYWY	-80.85	137187	UDRUSUDRDU	"
                    "	CARDPITTVVARIGYAMDYW	PITTVVARIGY	-81.04	137187	UDRUSUDRDU	"
                    "	CARDPITTVVARIGYAMDYW	VVARIGYAMDY	-83.35	137187	UDRUSUDRDU	"
                    "	CARDPIVYDYDYAMDYW	IVYDYDYAMDY	-83.3	137187	UDRUSUDRDU	"
                    "	CTTSLLWLRRRRYFDYW	SLLWLRRRRYF	-81.36	137188	BUURSUDSRR	"
                    "	CARSLLYYYGSGYFDYW	SLLYYYGSGYF	-80.6	137188	BUURSUDSRR	"
                    "	CARRIGGLLVYYFDYW	GGLLVYYFDYW	-81.43	137252	RRUURSUDRL	"
                    "	CARRRWLLPFGYFDVW	RRWLLPFGYFD	-80.69	137252	RRUURSUDRL	"
                    "	CARSASLLWFFFDYW	ASLLWFFFDYW	-86.13	137252	RRUURSUDRL	"
                    "	CARYGLLDVYWYFDVW	YGLLDVYWYFD	-83.54	137252	RRUURSUDRL	"
                    "	CASRALLRLDWYFDVW	RALLRLDWYFD	-80.51	137252	RRUURSUDRL	"
                    "	CVGWLLLWFAYW	VGWLLLWFAYW	-84.57	137252	RRUURSUDRL	"
                    "	CVRDGLLDWYFDVW	DGLLDWYFDVW	-81.33	137252	RRUURSUDRL	"
                    "	CVRDTLLVVHRYFDVW	DTLLVVHRYFD	-80.99	137252	RRUURSUDRL	"
                    "	CVRVFMVLYAMDYW	VRVFMVLYAMD	-82.01	137252	RRUURSUDRL	"
                    "	CARSLNLLLRPGLFAYW	LNLLLRPGLFA	-84.08	137379	SRUDDURUUL	"
                    "	CAIRSIYDGYLWFAYW	IRSIYDGYLWF	-82.18	137379	SRUDDURUUL	"
                    "	CARSIRLRLLPYFDYW	IRLRLLPYFDY	-81.31	137379	SRUDDURUUL	"
                    "	CARQLRLRFYWYFDVW	LRLRFYWYFDV	-81.84	137379	SRUDDURUUL	"
                    "	CAIGLLRPPWFAYW	IGLLRPPWFAY	-80.62	137379	SRUDDURUUL	"
                    "	CARNSLIYYYGSSFFAYW	IYYYGSSFFAY	-81.11	137379	SRUDDURUUL	"
                    "	CARPFHYYYGPWFAYW	FHYYYGPWFAY	-81.38	137379	SRUDDURUUL	"
                    "	CIYYYGSSLFAYW	IYYYGSSLFAY	-81.15	137379	SRUDDURUUL	"
                    "	CTSLLLAWFAYW	CTSLLLAWFAY	-81.99	137379	SRUDDURUUL	"
                    "	CTALLYPWFAYW	CTALLYPWFAY	-81.8	137379	SRUDDURUUL	"
                    "	CTPLYYDFFAYW	CTPLYYDFFAY	-80.97	137379	SRUDDURUUL	"
                    "	CTRFRDSLLLFDYW	CTRFRDSLLLF	-81.33	137379	SRUDDURUUL	"
                    "	CTTMGWLLRDWFAYW	MGWLLRDWFAY	-83.46	137379	SRUDDURUUL	"
                    "	CTWLLRDWFAYW	CTWLLRDWFAY	-82.03	137379	SRUDDURUUL	"
                    "	CVRILLHGLMAMDYW	VRILLHGLMAM	-84.94	137379	SRUDDURUUL	"
                    "	CARDLDGYYGVFFAYW	LDGYYGVFFAY	-81.57	137379	SRUDDURUUL	"
                    "	CTRRILLLDVW	CTRRILLLDVW	-81.15	137379	SRUDSRDDLD	"
                    "	CARSLRSYYHFFDYW	LRSYYHFFDYW	-83.97	137379	SRUDSRDDLD	"
                    "	CTRFRDSLLLFDYW	FRDSLLLFDYW	-84.52	137379	SRUDSRDDLD	"
                    "	CARLEGLTVVAPYFDVW	LEGLTVVAPYF	-80.68	137379	SRUDSSDDRD	"
                    "	CARSGLQELWVVDWYFDVW	LQELWVVDWYF	-82.56	137379	SRUDSSDDRD	"
                    "	CVRRLIYYQAFFDVW	VRRLIYYQAFF	-84.47	137379	SRUDSSDDRD	"
                    "	CARMDLGQLRLMGYFDYW	DLGQLRLMGYF	-80.86	137443	RLRSSRUURS	"
                    "	CTPLWLRRGPWFAYW	WLRRGPWFAYW	-81.03	137443	RLRSSRUURS	"
                    "	CTRYLLWVGGWFAYW	LLWVGGWFAYW	-85.19	137443	RLRSSRUURS	"
                    "	CTTIWLRRGGFFDYW	WLRRGGFFDYW	-81.89	137443	RLRSSRUURS	"
                    "	CARSLSLLLLYYFDYW	CARSLSLLLLY	-81.24	137443	RLRSSRUURS	"
                    "	CARSRGQLRLQAWFAYW	QLRLQAWFAYW	-81.29	137443	RLRSSRUURS	"
                    "	CASYFYGSAWFAYW	YFYGSAWFAYW	-80.74	137443	RLRSSRUURS	"
                    "	CARALYYFGYFDVW	ALYYFGYFDVW	-80.79	137443	RLRSSRUURS	"
                    "	CARDRLLRWPYYFDYW	LLRWPYYFDYW	-81.91	137443	RLRSSRUURS	"
                    "	CARFLTAQATLAWFAYW	FLTAQATLAWF	-81.55	137443	RLRSSRUURS	"
                    "	CARHEERGWLLRDGAWFAYW	LLRDGAWFAYW	-81.83	137443	RLRSSRUURS	"
                    "	CASVWLRQRAWFAYW	WLRQRAWFAYW	-81.05	137443	RLRSSRUURS	"
                    "	CTINWDALFAYW	TINWDALFAYW	-80.74	137443	RLRSSRUURS	"
                    "	CACARGEMGLRRLFAYW	EMGLRRLFAYW	-81.06	137443	RLRSSRUURS	"
                    "	CTVRQLRLRAWFAYW	QLRLRAWFAYW	-81.22	137443	RLRSSRUURS	"
                    "	CARDFSHLTWFAYW	DFSHLTWFAYW	-80.62	137443	RLRSSRUURS	"
                    "	CARDLGIRPFFAYW	DLGIRPFFAYW	-83.38	137443	RLRSSRUURS	"
                    "	CTTMGWLLRDWFAYW	TMGWLLRDWFA	-81.15	137443	RLRUDDURUU	"
                    "	CVRILLHGLMAMDYW	CVRILLHGLMA	-83.56	137443	RLRUDDURUU	"
                    "	CARSLNLLLRPGLFAYW	SLNLLLRPGLF	-86.38	137443	RLRUDDUURR	"
                    "	CARSLSLLLLYYFDYW	SLSLLLLYYFD	-83.15	137443	RLRUDDUURR	"
                    "	CARNSLIYYYGSSFFAYW	SLIYYYGSSFF	-81.22	137443	RLRUDDUURR	"
                    "	CARSASLLWFFFDYW	SASLLWFFFDY	-81.07	137443	RLRUDSRDDR	"
                    "	CARPLLLFLDYAMDYW	CARPLLLFLDY	-82.8	137443	RLRUDSRDDR	"
                    "	CAYLLYYYGILYAMDYW	LLYYYGILYAM	-84.39	137443	RLRUDSRDDR	"
                    "	CSLFITTVVAFYYYAMDYW	FITTVVAFYYY	-80.97	137443	RLRUDSRDDR	"
                    "	CVGWLLLWFAYW	CVGWLLLWFAY	-84.88	137443	RLRUDSRDDR	"
                    "	CVRFRLWFYAMDYW	CVRFRLWFYAM	-80.63	137443	RLRUDSRDDR	"
                    "	CVRGLLRFYAMDYW	CVRGLLRFYAM	-80.6	137443	RLRUDSRDDR	"
                    "	CVRVFMVLYAMDYW	CVRVFMVLYAM	-83.29	137443	RLRUDSRDDR	"
                    "	CATGWLLPYWYFDVW	ATGWLLPYWYF	-80.66	141220	DLDLLUSLRU	"
                    "	CAEGWLLWYFDVW	AEGWLLWYFDV	-80.51	141220	DLDLLUSLRU	"
                    "	CTRGEGWLLRLAWFAYW	GEGWLLRLAWF	-80.87	141220	DLDLLUSLRU	"
                    "	CVRDSLLAWFAYW	RDSLLAWFAYW	-80.8	141220	DLDLLUSLRU	"
                    "	CARRRWLLPFGYFDVW	RWLLPFGYFDV	-80.6	141220	LLDDRDUUDR	"
                    "	CARSLLRSFYWYFDVW	RSLLRSFYWYF	-83.98	141220	LLDDRDUUDR	"
                    "	CARSLNLLLRPGLFAYW	LLLRPGLFAYW	-82.76	141220	LLDDRDUUDR	"
                    "	CASGLLLRSYFDYW	GLLLRSYFDYW	-81.22	141220	LLDDRDUUDR	"
                    "	CATGWLLPYWYFDVW	GWLLPYWYFDV	-82.49	141220	LLDDRDUUDR	"
                    "	CATWWLLRDAWFAYW	WWLLRDAWFAY	-80.52	141220	LLDDRDUUDR	"
                    "	CTRDPLLLRPYWYFDVW	LLLRPYWYFDV	-81.15	141220	LLDDRDUUDR	"
                    "	CTRLLLRLYYFDYW	RLLLRLYYFDY	-81.57	141220	LLDDRDUUDR	"
                    "	CTRTWWLLQAWFAYW	WWLLQAWFAYW	-84.84	141220	LLDDRDUUDR	"
                    "	CTTLLWLRIWYMDYW	LLWLRIWYMDY	-81.82	141220	LLDDRDUUDR	"
                    "	CTTMGWLLRDWFAYW	GWLLRDWFAYW	-82.02	141220	LLDDRDUUDR	"
                    "	CTVEALLLRYPWTYFDYW	ALLLRYPWTYF	-82.25	141220	LLDDRDUUDR	"
                    "	CTWLLRDWFAYW	TWLLRDWFAYW	-82.3	141220	LLDDRDUUDR	"
                    "	CTRLLATVVAYYAMDYW	LLATVVAYYAM	-81.22	141283	DURUUSSUDR	"
                    "	CTGPPFITTVVPFAYW	FITTVVPFAYW	-80.56	141283	DURUUSSUDR	"
                    "	CARLLPSYYSISWFAYW	LLPSYYSISWF	-80.8	141283	DUULSLRSUU	"
                    "	CARLEYGNFPWYFDVW	LEYGNFPWYFD	-86.03	141283	SRUSUDRDUL	"
                    "	CARSELGLGLDYW	CARSELGLGLD	-81.09	141283	SRUSUDRDUL	"
                    "	CARIFITTVVAKYFDVW	IFITTVVAKYF	-80.59	141284	BDDRDUUDRL	"
                    "	CARLLLRRYYFDYW	LLLRRYYFDYW	-82.29	141284	BDDRDUUDRL	"
                    "	CARRRWLLPFGYFDVW	WLLPFGYFDVW	-83.42	141284	BDDRDUUDRL	"
                    "	CARSLLRSFYWYFDVW	SLLRSFYWYFD	-84.67	141284	BDDRDUUDRL	"
                    "	CASPLLRRWYFDVW	PLLRRWYFDVW	-82.13	141284	BDDRDUUDRL	"
                    "	CATGWLLPYWYFDVW	WLLPYWYFDVW	-85.11	141284	BDDRDUUDRL	"
                    "	CATLFSLYYFDYW	TLFSLYYFDYW	-81.33	141284	BDDRDUUDRL	"
                    "	CAIDDLLLRWVYFDYW	LLLRWVYFDYW	-84.67	141284	BDDRDUUDRL	"
                    "	CAILYSNYLFAYW	ILYSNYLFAYW	-81.31	141284	BDDRDUUDRL	"
                    "	CARSVLLWELGVYFDYW	VLLWELGVYFD	-83.01	141284	BDDRDUUDRL	"
                    "	CARTGWLLRDWYFDVW	WLLRDWYFDVW	-84.03	141284	BDDRDUUDRL	"
                    "	CARTLLRSPLYWYFDVW	LLRSPLYWYFD	-81.41	141284	BDDRDUUDRL	"
                    "	CARWLLLGDWYFDVW	LLLGDWYFDVW	-85.33	141284	BDDRDUUDRL	"
                    "	CARRIFPYWYFDVW	RIFPYWYFDVW	-81.32	141284	BDDRDUUDRL	"
                    "	CASRLLPYWYFDVW	RLLPYWYFDVW	-83	141284	BDDRDUUDRL	"
                    "	CATPFYYGYLYWYFDVW	FYYGYLYWYFD	-81.05	141284	BDDRDUUDRL	"
                    "	CTTQGLLLRSYYFDYW	LLLRSYYFDYW	-82.23	141284	BDDRDUUDRL	"
                    "	CARVLYDGYLAWFAYW	VLYDGYLAWFA	-81.43	141284	BDDRDUUDRL	"
                    "	CARWVFGYLYFDYW	WVFGYLYFDYW	-81.24	141284	BDDRDUUDRL	"
                    "	CTGLIYYYVYWYFDVW	LIYYYVYWYFD	-82.61	141284	BDDRDUUDRL	"
                    "	CTGWLLRSWYFDVW	WLLRSWYFDVW	-83.93	141284	BDDRDUUDRL	"
                    "	CTRLLLRLYYFDYW	LLLRLYYFDYW	-84.39	141284	BDDRDUUDRL	"
                    "	CVRRLIYYQAFFDVW	LIYYQAFFDVW	-80.6	141284	BDDRDUUDRL	"
                    "	CVRSLWSFWYFDVW	SLWSFWYFDVW	-81.51	141284	BDDRDUUDRL	"
                    "	CASRALLRLDWYFDVW	LLRLDWYFDVW	-80.57	141412	RSRDDRSSRL	"
                    "	CATLFSLYYFDYW	CATLFSLYYFD	-81.08	141412	RSRDDRSSRL	"
                    "	CTGGVLLLRYFDVW	GVLLLRYFDVW	-82.5	141412	RSRDDRSSRL	"
                    "	CAILLLPFYAMDYW	AILLLPFYAMD	-81.86	141412	RSRDDRSSRL	"
                    "	CARSVDFALFITTADFDYW	FALFITTADFD	-82.44	141412	RSRDDRSSRL	"
                    "	CARWLLLGDWYFDVW	RWLLLGDWYFD	-83.18	141412	RSRDDRSSRL	"
                    "	CALLWLTTWFAMDYW	LLWLTTWFAMD	-82.11	141412	RSRDDRSSRL	"
                    "	CALSLLRAYYFDYW	ALSLLRAYYFD	-80.8	141412	RSRDDRSSRL	"
                    "	CTAITTVVAFYWYFDVW	VVAFYWYFDVW	-80.55	141412	RSRDDRSSRL	"
                    "	CTPYLLLRYFDVW	PYLLLRYFDVW	-81.15	141412	RSRDDRSSRL	"
                    "	CTSLDLLWLRRGAYFDYW	LLWLRRGAYFD	-81.18	141412	RSRDDRSSRL	"
                    "	CTTLLWLRIWYMDYW	TLLWLRIWYMD	-80.55	141412	RSRDDRSSRL	"
                    "	CVRDTLLVVHRYFDVW	LLVVHRYFDVW	-80.92	141412	RSRDDRSSRL	"
                    "	CARRGLRYLFDYW	CARRGLRYLFD	-81.24	145633	SUUSRDUSUL	"
                    "	CARSDFAVLLRPHYAMDYW	CARSDFAVLLR	-81.25	145633	SUUSRDUSUL	"
                    "	CARSGLNLFFDYW	CARSGLNLFFD	-83.09	145633	SUUSRDUSUL	"
                    "	CAINTLRGYFDYW	AINTLRGYFDY	-81.93	145634	BUDLUDSDRU	"
                    "	CASLRDYDWYFDVW	ASLRDYDWYFD	-80.52	145635	BSUDLUDSDR	"
                    "	CAVRSVNWYFDVW	CAVRSVNWYFD	-82	145635	BSUDLUDSDR	"
                    "	CAINTLRGYFDYW	CAINTLRGYFD	-83.22	145635	BSUDLUDSDR	"
                    "	CAIRTLNWAFAYW	CAIRTLNWAFA	-85.72	145635	BSUDLUDSDR	"
                    "	CAISNLPYYFDYW	CAISNLPYYFD	-84.39	145635	BSUDLUDSDR	"
                    "	CASSGGLYYDYDWFAYW	GGLYYDYDWFA	-80.84	149471	LRULSLRDLR	"
                    "	CATPTLYYDYSWFAYW	PTLYYDYSWFA	-80.84	149471	LRULSLRDLR	"
                    "	CPTLYYDFPWFAYW	PTLYYDFPWFA	-83.78	149471	LRULSLRDLR	"
                    "	CATPTLYYDYSWFAYW	TLYYDYSWFAY	-80.95	149535	SULSLRDLRR	"
                    "	CARSDRIYYDYDGFAYW	RIYYDYDGFAY	-80.54	149535	SULSLRDLRS	"
                    "	CARAFYYDYEGFAYW	AFYYDYEGFAY	-82.41	149535	SULSLRDLRS	"
                    "	CARKILLDYEGFAYW	KILLDYEGFAY	-85.63	149535	SULSLRDLRS	"
                    "	CARSGFYYGYDEFAYW	GFYYGYDEFAY	-81.47	149535	SULSLRDLRS	"
                    "	CSLPIYYDYGGFAYW	PIYYDYGGFAY	-81.3	149535	SULSLRDLRS	"
                    "	CAKTHGNLLWDYAMDYW	NLLWDYAMDYW	-81.68	149535	SULSLRDSLL	"
                    "	CGNLLWDYAMDYW	NLLWDYAMDYW	-81.68	149535	SULSLRDSLL	"
                    "	CARPQIYDGYYPFAYW	IYDGYYPFAYW	-80.67	149536	USRDSURDRD	"
                    "	CARWSFYDGYYDFDYW	FYDGYYDFDYW	-81.17	149536	USRDSURDRD	"
                    "	CVREGYDFLFAYW	CVREGYDFLFA	-81.85	149536	USRDSURDRD	"
                    "	CARHEDEDDYDFPWFAYW	DDYDFPWFAYW	-83.38	153506	SLRDRLSLDR	"
                    "	CARHEDREGDYYGSFYYFDYW	YYGSFYYFDYW	-86.82	153506	SLRDRLSLDR	"
                    "	CARHEDSTMVTTTFSPFAYW    VTTTFSPFAYW	-84.54	153506	SLRDRLSLDR	";
    stringstream toread(exData);

    std::pair<superProtein*, vector<int> > AG = getAntigen("1FBI_X");

    // Will transform the data into a 'binding dataset'
    dataset< analyzedBinding > annotatedDataset;
    string longseq, CDR3;
    int ID = 0;
    cout << "List sequences / structures with their bound positions" << endl;
    while(toread >> longseq){
        toread >> CDR3;
        ID++;

        double energy;
        int pos;
        string structure;
        toread >> energy >> pos >> structure;
        cerr << longseq << " " << CDR3 << " " << energy << " " << pos << " " << structure << " ";
        struct3D s = struct3D(structure, UnDefined, pos);
        superProtein s2(s);
        string structureID = getUniqueIDStructure(s);

        //enum {interCodeWithIDpos, listAAPairs, AAcompoAGEpitope, AAcompoABParatope, seqAGEpitope, seqABParatope, motifAGEpitope, motifABParatope, motifsSizeGapsLigand, motifsSizeGapsRec, motifsChemicalLig, motifsChemicalRec, agregatesAGEpitope, agregatesABParatope, chemicalAGEpitope, chemicalABParatope, segmentedABParatope, segmentedAGEpitope, interMaskABParatope, interMaskAGEpitope, positionsBound, NB_features}; //nbneighbors, distChem, selfFolding,
        int degree = 1;
        vector<string> analyzedFeatures = structuralFeatures(*(AG.first), s2, degree);

        set<int> IDresonAntigen = stringToSet(analyzedFeatures[positionsBound]);

        analyzedBinding* bd = new analyzedBinding(CDR3, energy, structureID, "noInterCodeYet", IDresonAntigen, "UnknownYet", vector<string>());
        cerr << analyzedFeatures[positionsBound] << endl;

        stringstream genID; genID << "ID" << ID;
        annotatedDataset.addLine(genID.str(), longseq, bd, true);

    }
    cerr << "Finished generating dataset " << endl;
    return annotatedDataset;
}

// a bit symmetrical function, shows two proteins and the interface of binding.
void showParatopeEpitope(superProtein* prot1, superProtein* prot2, vector<int> forbiddenPos){
    set<int> PosProt1 = getOccupiedPositions(prot1);
    set<int> PosProt2 = getOccupiedPositions(prot2);
    set<int> vicinity1 = neighborPositions(PosProt1);
    set<int> vicinity2 = neighborPositions(PosProt2);

    // epitope
    set<int> interface1 = intersection_sets<int>(PosProt1, vicinity2);
    // paratope
    set<int> interface2 = intersection_sets<int>(PosProt2, vicinity1);

    cout << "Interface 1 is: " << print(interface1) << endl;
    cout << "Interface 2 is: " << print(interface2) << endl;
    #ifdef ALLOW_GRAPHICS
    setHeatmapColorZero({0.2, 0.2, 0.2});
    setHeatmapColorOne({1.0, 0.1, 0.1});
    #endif

    vector<std::pair<int, double> > heatmap1;
    for(set<int>::iterator it = interface1.begin(); it != interface1.end(); ++it){
        heatmap1.push_back(std::pair<int, double> (*it, 1.0));
    }

    vector<std::pair<int, double> > heatmap2;
    for(set<int>::iterator it = interface2.begin(); it != interface2.end(); ++it){
        heatmap2.push_back(std::pair<int, double> (*it, 1.0));
    }

    #ifdef ALLOW_GRAPHICS
    //Note: can not delete these pointers, they are needed for plotting after the function is finished
    if(forbiddenPos.size() > 0){
        addToDisplay(new set<int>(forbiddenPos.begin(), forbiddenPos.end()));
    }
    addToDisplay(new vector<std::pair<int, double> >(heatmap1));
    addToDisplay(new vector<std::pair<int, double> >(heatmap2));
    #endif
}

void showBindingHotspots(dataset<analyzedBinding>& annotatedDataset, string antigenID, int sizeSet){

    vector< std::pair<string, string> > res = setCoveringStructures(annotatedDataset, sizeSet);



    std::pair<superProtein*, vector<int> > AG = getAntigen(antigenID);

    std::map<int, int> posMap = posAntigens(AG.first);

    #ifdef ALLOW_GRAPHICS
    setHeatmapColorZero({0.2, 0.2, 0.2});
    setHeatmapColorOne({1.0, 0.1, 0.1});
    //setHeatmapColorOne({149./256., 196./256., 234./256.});  1FBI

    #endif

    // Shows first cluster:
    if(res.size() == 0) cerr << "NO Cluster" << endl;
    vector< set<int> > fullyBoundPerHotspot;
    cout << "Found " << res.size() << " different clusters " << endl;
    for(size_t i = 0; i < res.size(); ++i){
        cout << "Cluster is " << res[i].first << endl;
        set<int> S = stringToSet(res[i].first);

        // Option 1: show only the common positions
        vector<std::pair<int, double> > heatmap;
        for(set<int>::iterator it = S.begin(); it != S.end(); ++it){
            //cout << "Heatmap pos " << *it << endl;
            if(posMap.find(*it) != posMap.end()){
                heatmap.push_back(std::pair<int, double> (posMap[*it], 1.0));
            } else {
                cerr << "ERR: residue ID " << *it << " in antigen is no known" << endl;
            }
        }
        //addToDisplay(heatmap);

        size_t DS = annotatedDataset.nLines();
        // Option 2: Show a density of binding for the sequences in this 'epitope'
        int cptIncl = 0;
        for(size_t j = 0; j < DS; ++j){
            set<int> bindingPos = annotatedDataset.getElement(j)->bindingAGresiduesID; //bindingAGPositionsInSpace;
            if(isIncluded(S, bindingPos)){
                cptIncl++;
            }
        }
        // list of sequences included:

        std::map<int, double> heatmapDensity;
        for(size_t j = 0; j < DS; ++j){
            // try to avoid to copy the structure binding, so directly get the fields.
            string AAseq = annotatedDataset.getElement(j)->AAseq;
            set<int> bindingPos = annotatedDataset.getElement(j)->bindingAGresiduesID; //bindingAGPositionsInSpace;
            if(isIncluded(S, bindingPos)){
                for(set<int>::iterator itk = bindingPos.begin(); itk != bindingPos.end(); ++itk){
                    if(posMap.find(*itk) != posMap.end()){
                        if(heatmapDensity.find(*itk) != heatmapDensity.end()){
                            heatmapDensity[posMap[*itk] ] = 1. / static_cast<double>(cptIncl);
                        } else {
                             heatmapDensity[posMap[*itk] ] += 1. / static_cast<double>(cptIncl);
                        }
                    } else {
                        cerr << "ERR: residue ID " << *itk << " in antigen is no known" << endl;
                    }
                }
                cout << AAseq << "\t" << i+1 << "/" << res.size() << "\t" << res[i].first << "\tPositions\t" << print(bindingPos) << endl; //"\tStructures\t" << res[i].second << endl;
            }
        }

        // Show the 100% occupied positions

        set<int> fullyBound;
        for(std::map<int, double>::iterator itk = heatmapDensity.begin(); itk != heatmapDensity.end(); ++itk){
            if(itk->second >= 0.999){
                // reverts position in lattice => position on antigen
                bool found = false;
                for(std::map<int,int>::iterator it2 = posMap.begin(); it2 != posMap.end(); ++it2){
                    if(it2->second == itk->first){
                        fullyBound.insert(it2->first);
                        found = true;
                    }
                }
                if(!found) {
                    cerr << "ERR: lattice position " << itk->first << " not found on antigen" << endl;
                }
            }
        }
        fullyBoundPerHotspot.push_back(fullyBound);
        cout << "100% occupied positions for cluster " << i+1 << "\t" << print(fullyBound) << endl;



        // transforms the map into a vector,
        vector< std::pair<int, double> > HeatmapDensityForPlotting;
        std::copy(heatmapDensity.begin(), heatmapDensity.end(), std::back_inserter<std::vector<std::pair< int, double > > > (HeatmapDensityForPlotting));

        #ifdef ALLOW_GRAPHICS
        addToDisplay(new vector< std::pair<int, double> >(HeatmapDensityForPlotting));
        #endif
    }


    // ========= Now outputs the sum up information: ===============
    ofstream hotspots_report(string("Hotspots") + antigenID + string(".txt"));
    hotspots_report << "#Antigen " << antigenID << endl;
    hotspots_report << "#hotspotsCore = { ";
    for(size_t i = 0; i < res.size(); ++i){
        if(i > 0) hotspots_report << ",";
        hotspots_report << "{" << setToString(stringToSet(res[i].first, ' '),',') << "}  ";
    }
    hotspots_report << "};" << "\t";
    hotspots_report << "hotspotsLarge = { ";
    for(size_t i = 0; i < res.size(); ++i){
        if(i > 0) hotspots_report << ",";
        hotspots_report << "{" << setToString(fullyBoundPerHotspot[i], ',') << "}  ";
    }
    hotspots_report << "};" << endl;

    hotspots_report << "Structure\tHotspotID\tHospotCore\tHospot100pcts\tnBound100pcts\tBoundPos\n";
    for(size_t i = 0; i < res.size(); ++i){
        set<int> S = stringToSet(res[i].first);
        size_t DS = annotatedDataset.nLines();
        for(size_t j = 0; j < DS; ++j){
            string AAseq = annotatedDataset.getElement(j)->AAseq;
            set<int> bindingPos = annotatedDataset.getElement(j)->bindingAGresiduesID; //bindingAGPositionsInSpace;
            if(isIncluded(S, bindingPos)){
                hotspots_report << AAseq << "\t" << antigenID << "_H" << i+1 << "\t" << setToString(stringToSet(res[i].first, ' '),'-') << "\t" << setToString(fullyBoundPerHotspot[i], '-') <<  "\t" << fullyBoundPerHotspot[i].size() << "\t" << setToString(bindingPos, '-') << endl;
            }
        }
    }
    hotspots_report.close();



    #ifdef ALLOW_GRAPHICS
    // Note: Control will be given to openGL. Don't put code afterwards.
    // we let the user decide to show the ligand, because displayLigand destroys much and does glutmainloop
    //displayLigand(AG.first, AG.second, true);
    //addToDisplay(AG.first, true);
    //addToDisplay(AG.second);
    #endif
}

void testSetCoveringStructures(){

}

//epitope::epitope()
//{
//    // the data structure will be: AA sequence - structureID - set of occupied positions - affinity



//}
