#include "motifFeatures.h"
#include "../Ymir/ymir.h"
#include "../Ymir/plot3d.h"
#include "../Ymir/proteins.h"
#include "fileformats.h"
#include "poolstructs.h"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <regex>
#include "antigenLib.h" // for testing only
#include "../Tools/dirent.h"
using namespace std;

vector<string> split(const string& str, const string& delim)
{
    vector<string> tokens;
    size_t prev = 0, pos = 0;
    do
    {
        pos = str.find(delim, prev);
        if (pos == string::npos) pos = str.length();
        string token = str.substr(prev, pos-prev);
        if (!token.empty()) tokens.push_back(token);
        prev = pos + delim.length();
    }
    while (pos < str.length() && prev < str.length());
    return tokens;
}

void testFeatures(){
    {
        superProtein v0 = superProtein();
        struct3D S01 = struct3D("SURRDRDDSDRLLDDRRUUL", UnDefined, 133152);
        superProtein v1 = insert(v0, S01, 0);
        struct3D S02 = struct3D("LUSDDURDRDDLLSURRUULDDUDSRLRR", UnDefined, 128867);
        superProtein v2 = insert(v1, S02, 1000);
        v2.setAAs(string("GIVEQCCTSICSLYQLENYCNFVNQHLCGSHLVEALYLVCGERGFFYTPKA"));

        {
// This is old code. renew
            //            features f(&v2, 1);
//            ofstream fwr("C:/Users/pprobert/Desktop/PulledDataset/Analyzed.txt");
//            fwr << f.getPropertiesDeprecated("C:/Users/pprobert/Desktop/PulledDataset/BigDataset.txt");
//            fwr.close();
        }
    }
}
//4	Subchains
//133152	SUSSDUSSDSSDLLUSSLSDDSLSSSDRLSUDDLURDSLDLUURRLUSRLSSSUDLSRDRLSRDDLUSLSRRDLLRSULRUSLRLRLRDURDUURDSURDDSLSLDSLRSSLRSSDSUUDRRSSLURRDDLSRLLUURLUURDSDSSSSRDRSRLRSURRLDUSSUUDUSLRLSDDSUSLSRDRRLRDLRLRLSDSSUURSRSSSLSUDLRULRLDUSDSURRSRULLRRDDLLSULLRSDUSSLSUUDSUULLRSDSRDLLULRLURDSSLURDDSURLRSRDSSUSSRSSUSDUUDLLSRDLLURRLLDDUDSDDUUDLSURR?
//161709	RDSURDSSUSSLDDSRUDDUURUSSDUUSDRUDSUURRSSRSLUURLDUDUDLUDSRSSDRSRDDLRSDLLSURDLSRDSRLSRSDLSRSUDRDSLUURDSLRLLUDDUSDURLDSRLURURLRDLLUURRSLRDDLSLLRSLSSDSSLDUDLDSLRDLSRSDLRSRSRDULRSLLSDUUDSDSRUDDSSSLRSLDSURDLRLDRUDDRLRLSSSRLSSRUSSRSRLLRRDLLRLULSSSSSSRSRUSRLDDRRURLLRRDLLSSLSURDUDULSLSRLSDLLUDSURRSSUDSRLLSDDUUDDLLUURDLULRLDLLUURSLUR?
//189672	LLSSRLRLRDLDURDLSUUDDLDSULLDUSRLRSSLLUUDUUDLUSLSDUSDSSSLSRRSLDRDLRSLLURRUUDLSRDSRSLRSDSLRULLDLURDDLSURDURLUSRLSRLRLSDSURRSURRDSLRDUUDDLLUULRRLSUSULDURLUULSRLLUSSDULSDLLDRSLDRSRDSSSRLSLLRLSSSUDUDUUDSDSULURDUDSDRLRLUDSRLSRUDURSURRDDUUDDLSDSSSRLSSUULUDSURSDUURRDDUULLSUSSDUDSUDUURLSSUDLUUDRLLSSLRLRRDSLLUUDDLURDLLUUSLSDLLUURRDDL?
//152796	SDURLSSSSUSURRDUDSRURSUDRLRUDSRLRDDLULDRDLULLSLSURLDUSSSSLDSSLRULLRURRDDURLLRDSLDSSSDURLULRDLLSURRDLRRLSRSRSLRUDLRLUDUDSDDLRLURDDLLUDSUURRLRSUSRLUSSDUULUSSSLLRSLSURLDLDULDUDSURUDRSLRULSLLURLSDUSLSSRRLLSSDURDLDDRSLRUDLSRSSUDDRDLURSDSLULUSDUSDSLUSULLSSRUSDURRDDLURDRDUDSUSSSSLDLSSRSDSLSDUDSLSDSSSURDRUURSDLUUDLDDRRDUSUURRDLLRRD




// function to get affinities
//double val = AAaffinity(AA_ID(receptorAASeq[i]), AA_ID(receptorAASeq[j]));

// Recalculates the affinity from an interaction code (aAbGaH...)
// note: this is the total affinity, but we want to maximize the binding...
double affinityCodeTot(string receptor, string codeInters){
    double E = 0;
    if((codeInters.size() % 2) != 0) cerr << "ERR: affinityCode, codeInters should have even size " << codeInters << endl;
    size_t receptorSizeAA = receptor.size();
    for(size_t j = 0; j < codeInters.size() / 2; ++j){
        char A1 = codeInters[2*j];
        char A2 = codeInters[2*j+1];
        if((A1 >= 'a') && (A1 <= 'z')){
            size_t AApos1 = static_cast<size_t>(A1 - 'a');
            if(AApos1 >= receptorSizeAA) cerr << "ERR: affinityCode, the position " << A1 << " is out of bounds of receptor sequence " << receptor << endl;
            char AAchar1 = receptor[AApos1];
            if((A2 >= 'A') && (A2 <= 'Z')){
                E += AAaffinity(AA_ID(AAchar1), AA_ID(A2));
            } else {
                if((A2 >= 'a') && (A2 <= 'z')){
                    size_t AApos2 = static_cast<size_t>(A2 - 'a');
                    if(AApos2 >= receptorSizeAA) cerr << "ERR: affinityCode, the position " << A2 << " is out of bounds of receptor sequence " << receptor << endl;
                    char AAchar2 = receptor[AApos2];
                    E += AAaffinity(AA_ID(AAchar1), AA_ID(AAchar2));
                } else {
                    cerr << "ERR: affinityCode, code " << A1 << A2 << " incorrect" << endl;
                }
            }
        } else {
            cerr << "ERR: affinityCode, code " << A1 << A2 << " incorrect" << endl;
        }
    }
    return E;
}

// ----------------- pieces of codes to try to find the best antibody ---------------------------
// Produces all possible masks with N positions
vector<string> recursiveMask(vector<bool> positionsToTryAll, int positionToDo){
    vector<string> result;
    if(static_cast<size_t>(positionToDo) == positionsToTryAll.size()) return result;
    if(positionToDo < 0) {cerr << "ERR: recursiveMask, the positionToDo is getting negative. Should not happen." << endl; return result;}
    vector<string> previous = recursiveMask(positionsToTryAll, positionToDo+1);
    bool tryAll = positionsToTryAll[static_cast<size_t>(positionToDo)];
    // then need to enumerate all AAs
    size_t nEnum = previous.size();
    for(size_t i = 0; i < nEnum; ++i){
        string smaller = previous[i];
        if(tryAll){
            for(int k = static_cast<int>(AA::Cys); k < static_cast<int>(AA::NB_AAs); ++k){
                result.push_back(smaller + string(1, AAname(static_cast<AA>(k))));
            }
        } else {
            result.push_back(smaller + string(1, '-'));
        }
    }
    return result;
}



// input: masked AAs = either - = inknown or a capital letter: fixed.
// Note: this is optimizing both the binding and total energy, because the self folding AAs are masked (fixed)
string bestEnergyWithReceptorMask(string maskedAAs, string codeInters){
    string res;
    if((codeInters.size() % 2) != 0) cerr << "ERR: bestEnergyWithReceptorMask, codeInters should have even size " << codeInters << endl;
    vector<string> listOptCharsPerPosition = vector<string>(maskedAAs.size(), string(""));
    // each position is coded as a lower case letter.
    for(size_t i = 0; i < maskedAAs.size(); ++i){
        if(maskedAAs[i] == '-'){
            char thisPosition = 'a' + static_cast<char>(i);
            // now tries all AAs at this position and keeps the best.
            string bestAAsHere = "";
            double bestEnerg = +1e6;
            for(int k = static_cast<int>(AA::Cys); k < static_cast<int>(AA::NB_AAs); ++k){
                double thisAAEnerg = 0;
                // calculates the contribution of this AA in the interaction
                for(size_t j = 0; j < codeInters.size() / 2; ++j){
                    char A1 = codeInters[2*j];
                    char A2 = codeInters[2*j+1];
                    if(A1 == thisPosition){
                        if((A2 >= 'A') && (A2 <= 'Z')){
                           thisAAEnerg += AAaffinity(static_cast<AA>(k), AA_ID(A2));
                        } else {
                            cerr << "ERR: bestEnergyWithReceptorMask, all the AAs that make self-folding should be fixed in the mask. Here, got code " << A1 << A2 << ", so can not find the optimal AAs for this position" << endl;
                        }
                    }
                }
                // and only keeps if improving
                if(thisAAEnerg < bestEnerg){
                    bestAAsHere = string(1, AAname(static_cast<AA>(k)));
                }
                if(fabs(thisAAEnerg - bestEnerg) < 1e-6){
                    bestAAsHere = bestAAsHere + string(1, AAname(static_cast<AA>(k)));
                }
                bestEnerg = min(bestEnerg, thisAAEnerg);
            }
            listOptCharsPerPosition[i] = bestAAsHere;
        } else {
            listOptCharsPerPosition[i] = string(1, maskedAAs[i]);
        }
        res[i] = listOptCharsPerPosition[i].at(0);
    }
    cout << "Mask " << maskedAAs << "Optimal = " << res << ", possibles = ";
    for(size_t i = 0; i < listOptCharsPerPosition.size(); ++i){
        cout << listOptCharsPerPosition[i] << ",";
    }
    cout << endl;
    return res;
}

// Pieces of code to get the best possible
std::pair<string, double> bestSequenceForStructure(string codeInters, int sizeReceptors){
    // 1- find positions that are self-binding
    if(sizeReceptors < 0) cerr << "ERR: bestSequenceForStructure, screw Santa, you gave a negative sizeReceptors!" << endl;
    vector<bool> selfInterPos = vector<bool>(static_cast<size_t>(sizeReceptors), false);
    for(size_t j = 0; j < codeInters.size() / 2; ++j){
        char A1 = codeInters[2*j];
        char A2 = codeInters[2*j+1];
        if((A1 < 'a') || (A1 > 'z')) cerr << "ERR: bestSequenceForStructure, the codes of interactions should start by position (lower case) then position or AA. Code " << A1 << A2 << " is incorrect." << endl;
        size_t AApos1 = static_cast<size_t>(A1 - 'a');
        if((A2 >= 'a') | (A2 <= 'z')){
            int insidePos = static_cast<int>(A2 - 'a');
            if((insidePos < 0) || (insidePos >= sizeReceptors)){
                cerr << "bestSequenceForStructure(" << codeInters << ", sizeRec=" << sizeReceptors << "), a code position " << A2 << " is out of bounds " << endl;
                return std::pair<string,double> (string("Fail"), NAN);
            }
            selfInterPos[static_cast<size_t>(insidePos)] = true;
            selfInterPos[static_cast<size_t>(AApos1)] = true;
        }
    }
    vector<string> toDos = recursiveMask(selfInterPos);
    vector<string> bestMasks; // to store all optimal masks, so it would be possible to recompute all possible sequences.
    double bestE = 1e6;
    string globalBest = "None";
    for(size_t i = 0; i < toDos.size(); ++i){
        string maskedAAs = toDos[i];
        string localBest = bestEnergyWithReceptorMask(maskedAAs, codeInters);
        double localE = affinityCodeTot(localBest, codeInters);
        if(localE < bestE){
            bestMasks.clear();
        }
        if( localE <= bestE){
            bestMasks.push_back(maskedAAs);
            globalBest = localBest;
        }
        bestE = min(localE, bestE);
    }
    return std::pair<string,double>(globalBest, bestE);
}

void testBestSequencesPerStructure(){
    // testing double affinityCode(string receptor, string codeInters){
    // example
    // NTAHCSQEPKM	-71.82	iGjIkCeLgLdYkQhLjLkYjCbFdFkVbLcYaVkYjTiA	1	120675	UUDUSLSSLS
    cout << "Expect -71.82 : " << affinityCodeTot("NTAHCSQEPKM", "iGjIkCeLgLdYkQhLjLkYjCbFdFkVbLcYaVkYjTiA");
    cout << "Expect -3.39 : " << affinityCodeTot("MGH", "ab");

    // testing string bestEnergyWithReceptorMask(string maskedAAs, string codeInters){

    // testing vector<string> recursiveMask(vector<bool> positionsToTryAll, int positionToDo = 0){

    // testing std::pair<string, double> bestSequenceForStructure(string codeInters, int sizeReceptors){

}





// Chemical properties
char chemCode(char AAletter){
    switch(AAletter){
    case 'D': { return 'c';}
    case 'E': { return 'c';}
    case 'H': { return 'c';}
    case 'R': { return 'c';}
    case 'K': { return 'c';}
    case 'S': { return 'p';}
    case 'T': { return 'p';}
    case 'N': { return 'p';}
    case 'Q': { return 'p';}
    case 'G': { return 'n';}
    case 'A': { return 'n';}
    case 'V': { return 'n';}
    case 'L': { return 'n';}
    case 'I': { return 'n';}
    case 'P': { return 'n';}
    case 'M': { return 'n';}
    case 'C': { return 'n';}
    case 'F': { return 'r';}
    case 'W': { return 'r';}
    case 'Y': { return 'r';}
    case 'B': { return '!';}
    case 'J': { return '!';}
    case 'O': { return '!';}
    case 'U': { return '!';}
    case 'X': { return '!';}
    case '-': { return '-';}
    case '0': { return '0';}
    case '1': { return '1';}
    case '2': { return '2';}
    case '3': { return '3';}
    case '4': { return '4';}
    case '5': { return '5';}
    case '6': { return '6';}
    case '7': { return '7';}
    case '8': { return '8';}
    case '9': { return '9';}
    default: return '!';
    }
}



// this function is too complex!
// input: vector size 26, for each letter => returns the number of each AA:   A01B02C10D01
string codeAAcompo(vector<int> numberPerLetter){
    stringstream res;
    for(int i = static_cast<int>(Cys); i < static_cast<int>(NB_AAs); ++i){
        int pos = static_cast<int>(AAname(static_cast<AA>(i)) - 'A');
        //cout << AAname(static_cast<AA>(i)) << "->" << pos << endl;
        //cout << numberPerLetter.size() << endl;
        if((pos >= 0) && (pos < 26)) {
            int nb = numberPerLetter.at(static_cast<size_t>(pos));
            res << AAname((AA) i);
            if(nb < 9) res << '0';
            if(nb < 99) res << '0';
            res << nb;
        }
    }
    return res.str();
}

// input: vector size 26, for each letter => returns the number of each AA:   A01C10D01 (only with good letters)
string shortCodeAAcompo(vector<int> numberPerLetter){
    stringstream res;
    for(size_t i = 0; i < numberPerLetter.size(); ++i){
        if(chemCode('A'+i) != '!'){
            int nb = numberPerLetter[i];
            if(nb < 10) res << static_cast<char>('A' + i) << '0' << numberPerLetter[i];
            else res << ('A' + i) << numberPerLetter[i];
        }
    }
    return res.str();
}

string fullCodeAAcompo(vector<int> numberPerLetter, char sep){
    double sum = 0;
    for(size_t i = 0; i < numberPerLetter.size(); ++i){
        sum += static_cast<double>(numberPerLetter[i]);
    }
    stringstream res;
    for(size_t i = 0; i < numberPerLetter.size(); ++i){
        if(chemCode('A'+i) != '!'){
            if(i > 0) res << sep;
            if(sum < 0.5){
                res << 0;
            } else {
                res << static_cast<double>(numberPerLetter[i]) / sum;
            }
        }
    }
    return res.str();
}

string convertAAtoAAusage(string AAseq, char sep){
   vector<int> AAcompLig(26,0);
   for(size_t i = 0; i < AAseq.size(); ++i){
       char AAlig = AAseq[i];
       if(chemCode(AAlig) != '!'){
           AAcompLig.at(static_cast<size_t>(AAlig - 'A'))++;
       }
   }
   return fullCodeAAcompo(AAcompLig, sep);
}

string convertAAtoAAcode(string AAseq){
   vector<int> AAcompLig(26,0);
   for(size_t i = 0; i < AAseq.size(); ++i){
       char AAlig = AAseq[i];
       if(chemCode(AAlig) != '!'){
           AAcompLig.at(static_cast<size_t>(AAlig - 'A'))++;
       }
   }
   return shortCodeAAcompo(AAcompLig);
}

void testsAAcompo(){
    string s1= "CAIDGYSLYWYFDVW";
    cout << convertAAtoAAusage(s1) << endl;
    cout << convertAAtoAAcode(s1) << endl;
}


// transforms an integer into string with total of nChars (zeros before)
string fixedInt(int n, int nChars){
    stringstream res;
    res << n;
    if(res.str().size() > nChars) cerr << "ERR: fixedInt(" << n << ", nChars=" << nChars << "), number too big for this number of digits!" << endl;
    return string(nChars - res.str().size(), '0') + res.str();
}


// converts a string of AAs into a string of chemical properties
string convertChemical(string whatever){
    string res = string(whatever.size(), '?');
    for(size_t i = 0; i < whatever.size(); ++i){
        res[i] = chemCode(whatever[i]);
    }
    return res;
}

// A minimal function if needs to be performed veryfast, but essentially redundant with the next function
std::pair<string , string> minimalFeaturesDegreeOne(superProtein& ligand, superProtein & s2){
    std::pair<string , string> res("","");

    stringstream interCodeID;
    if((!ligand.structure) || (!s2.structure) || (ligand.points.size() == 0) || (s2.points.size() == 0)) {return res;}
    if(!ligand.structure->properlyFolded){cerr << "! unfolded - "; return res;}

    size_t SL = ligand.points.size();
    size_t S2L = s2.points.size();
    vector<bool> interactingPositionsLig(SL, false);
    vector<bool> interactingPositionsRec(S2L, false);
    vector<int> degreesPositionsLig(SL, 0);
    vector<int> degreesPositionsRec(S2L, 0);
    for(size_t i = 0; i < SL; ++i){ // lig
        for(size_t j = 0; j < S2L; ++j){ // rec
            if(lattice::areNeighbors(ligand.points[i].IDposition, s2.points[j].IDposition)){
                // this is degree 1 binding
                interactingPositionsLig[i] = true;
                interactingPositionsRec[j] = true;
                degreesPositionsLig[i]++;
                degreesPositionsRec[j]++;
                interCodeID << static_cast<char>('a' + j) << fixedInt(ligand.points[i].IDresidue,4);
            }
        }
    }
    stringstream posBound;
    for(size_t i = 0; i < SL; ++i){ // lig
        int IDposInLigand = ligand[static_cast<int>(i)].IDresidue;
        if(interactingPositionsLig[i]){
            posBound << IDposInLigand << " ";
        }
    }
    res.first = interCodeID.str();
    res.second = posBound.str();
    return res;
}



// actually we only need the features of the best folding, so it's OK.
// This time, both
//enum {interCodeWithIDpos, listAAPairs, AAcompoAGEpitope, AAcompoABParatope, seqAGEpitope, seqABParatope, motifAGEpitope, motifABParatope, motifsSizeGapsLigand, motifsSizeGapsRec, motifsChemicalLig, motifsChemicalRec, agregatesAGEpitope, agregatesABParatope, chemicalAGEpitope, chemicalABParatope, positionsBound, segmentedABParatope, segmentedAGEpitope, interMaskABParatope, interMaskAGEpitope, NB_features}; //nbneighbors, distChem, selfFolding,
vector<string> structuralFeatures(superProtein& ligand, superProtein & s2, int minDegreeInteract, bool includeDegree){
    //minDegreeInteract = 1;
    //includeDegree = true;
    vector<string> res = vector<string>(NB_features+10, string(""));
    size_t SL = ligand.points.size();
    size_t S2L = s2.points.size();

    string listPairs;
    vector<int> AAcompLig(26,0);
    vector<int> AAcompRec(26,0);
    vector<bool> interactingPositionsLig(SL, false);
    vector<bool> interactingPositionsRec(S2L, false);
    string interMaskLig = string(SL, '0');
    string interMaskRec = string(S2L, '0');
    vector<int> degreesPositionsLig(SL, 0);
    vector<int> degreesPositionsRec(S2L, 0);
    stringstream interCodeID;
    stringstream interCodeIDinternal;

    stringstream codeInters;
    if((!ligand.structure) || (!s2.structure) || (ligand.points.size() == 0) || (s2.points.size() == 0)) {
        return res;
    }
    if(!ligand.structure->properlyFolded){
        cerr << "! unfolded - ";
        return res;
    }

    vector<double> contribEpi = vector<double>(SL, 0.0);
    vector<double> contribPara = vector<double>(S2L, 0.0);
    vector<double> maskcontribEpi = vector<double>(SL, 0.0);
    vector<double> maskcontribPara = vector<double>(S2L, 0.0);

    // 1- transforms the interaction code into two binary vectors (the interaction masks)
    for(size_t i = 0; i < SL; ++i){ // lig
        for(size_t j = 0; j < S2L; ++j){ // rec
            if(lattice::areNeighbors(ligand.points[i].IDposition, s2.points[j].IDposition)){
                double interEner = AAaffinity(ligand[static_cast<int>(i)].TypeResidue, s2[static_cast<int>(j)].TypeResidue);
                contribEpi[i] += interEner;
                contribPara[j] += interEner;
                maskcontribEpi[i] += 1;
                maskcontribPara[j] += 1;

                char AAlig = AAname(ligand[static_cast<int>(i)].TypeResidue);
                char AArec = AAname(s2[static_cast<int>(j)].TypeResidue);

                listPairs = listPairs + string(1, AAlig) + string(1, AArec);

                // these 2 lines are useless because overwritten by degree just after the loops
                interactingPositionsLig[i] = true;
                interactingPositionsRec[j] = true;
                interMaskLig[i] = (char) interMaskLig[i] + 1; // chars
                interMaskRec[j] = (char) interMaskLig[j] + 1; // chars
                degreesPositionsLig[i]++;
                degreesPositionsRec[j]++;

                // should not count an AA twice!
                //res.push_back(std::pair<int,int>(j,-NB_AAs-1+(int) ligand.points[i].TypeResidue));
                // it means residue i in ligand is interacting with residue j in

                //codeInters << (char) ('a' + (char) j) << AAname(ligand.points[i].TypeResidue);
                interCodeID << (char) ('a' + (char) j) << fixedInt(ligand.points[i].IDresidue,4);

            }
        }
    }

    // there will be at least one interaction.
    stringstream scontrib1;
    for(size_t i = 0; i < SL; ++i){ // lig
        if(maskcontribEpi[i] > 0){
            scontrib1 << AAname(ligand[static_cast<int>(i)].TypeResidue) << i << ":" << maskcontribEpi[i] << "_" << contribEpi[i] << ",";
        }
    }
    res[contribPerAAepi] = scontrib1.str();

    stringstream scontrib2;
    for(size_t j = 0; j < S2L; ++j){ // rec
        if(maskcontribPara[j] > 0){
            scontrib2 << AAname(s2[static_cast<int>(j)].TypeResidue) << static_cast<char>('a' + j) << ":" << maskcontribPara[j] << "_" << contribPara[j] << ",";
        }
    }
    res[contribPerAAparaBind] = scontrib2.str();


    vector<double> contribParaFold = vector<double>(S2L, 0.0);
    vector<double> maskcontribParaFold = vector<double>(S2L, 0.0);
    for(size_t j1 = 0; j1 < S2L; ++j1){ // rec
        for(size_t j2 = j1+2; j2 < S2L; ++j2){ // rec
            if(lattice::areNeighbors(s2.points[j1].IDposition, s2.points[j2].IDposition)){
                interCodeIDinternal << (char) ('a' + (char) j1) << (char) ('a' + (char) j2) << "-";
                double interEner = AAaffinity(s2[static_cast<int>(j1)].TypeResidue, s2[static_cast<int>(j2)].TypeResidue);
                contribParaFold[j1] += interEner;
                contribParaFold[j2] += interEner;
                maskcontribParaFold[j1] += 1;
                maskcontribParaFold[j2] += 1;
            }
        }
    }

    stringstream scontrib3;
    for(size_t j1 = 0; j1 < S2L; ++j1){ // rec
        if(maskcontribParaFold[j1] > 0){
            scontrib3 << AAname(s2[static_cast<int>(j1)].TypeResidue) << static_cast<char>('a' + j1) << ":" << maskcontribParaFold[j1] << "_" << contribParaFold[j1] << ",";
        }
    }
    res[contribPerAAparaFold] = scontrib3.str();


    for(size_t i = 0; i < SL; ++i){
        interactingPositionsLig[i] = (degreesPositionsLig[i] >= minDegreeInteract);
    }
    for(size_t j = 0; j < S2L; ++j){
        interactingPositionsRec[j] = (degreesPositionsRec[j] >= minDegreeInteract);
    }
    string seqLig;
    string motifLig;
    string motifGapLig;
    string motifGapSizLig;
    string aggrLig;
    string segmentedLig;
    bool alreadyInGap = false;
    int gapsize = 0;
    for(size_t i = 0; i < SL; ++i){ // lig
        char AAlig = AAname(ligand[static_cast<int>(i)].TypeResidue);
        string degreeSuffix = ((includeDegree) ? string(1, '1' + degreesPositionsLig[i] - 1) : string(""));
        if(interactingPositionsLig[i]){
            if(i == 0 || !interactingPositionsLig[i-1]) segmentedLig = segmentedLig + string("["); //end of gap,
            alreadyInGap = false;   //end of gap, start of binding region
            if(gapsize > 0) motifGapSizLig = motifGapSizLig + fixedInt(gapsize, 3);
            gapsize = 0;


            seqLig = seqLig + string(1, AAlig) + degreeSuffix;
            segmentedLig = segmentedLig + string(1, AAlig) + degreeSuffix;
            motifLig = motifLig + string(1, 'X') + degreeSuffix;
            aggrLig = aggrLig + string(1, AAlig) + degreeSuffix;
            motifGapSizLig = motifGapSizLig + string(1, 'X'); //string(1, AAlig);
            if(AAlig != '?'){
                //cerr << AAlig << "+";
                AAcompLig.at(static_cast<size_t>(AAlig - 'A'))++;
            }
        } else {
            gapsize++;  //extending gap
            if((!alreadyInGap) && (seqLig.size() > 0)){ // start of a gap, ending of binding
                motifLig = motifLig + string(1, '-');
                aggrLig = aggrLig + string(1, '-');
                if(includeDegree) { // makes it multiple of 2 if including degree
                    motifLig = motifLig + string(1, '-');
                    aggrLig = aggrLig + string(1, '-');
                }
                segmentedLig = segmentedLig + string("]");
            }
            segmentedLig = segmentedLig + string(1, AAlig) + degreeSuffix;
            if(seqLig.size() > 0) alreadyInGap = true; // start or continuation of gap
        }
    }
    // sometimes there is one or two '-'s at the end, removes them
    if(!alreadyInGap) segmentedLig = segmentedLig + string("]");
    if(motifLig.back() == '-') motifLig = motifLig.substr(0, motifLig.size()-1);
    if(motifLig.back() == '-') motifLig = motifLig.substr(0, motifLig.size()-1);
    if(aggrLig.back() == '-') aggrLig = aggrLig.substr(0, aggrLig.size()-1);
    if(aggrLig.back() == '-') aggrLig = aggrLig.substr(0, aggrLig.size()-1);

    string seqRec;
    string motifRec;
    string motifRecLig;
    string motifGapSizRec;
    string aggrRec;
    string segmentedRec;
    alreadyInGap = false;
    gapsize = 0;
    for(size_t j = 0; j < S2L; ++j){ // rec
        char AArec = AAname(s2[static_cast<int>(j)].TypeResidue);
        string degreeSuffix = ((includeDegree) ? string(1, '1' + degreesPositionsRec[j] - 1) : string(""));
        if(interactingPositionsRec[j]){
            if(j == 0 || !interactingPositionsRec[j-1]) segmentedRec = segmentedRec + string("["); //end of gap,
            alreadyInGap = false;
            if(gapsize > 0) motifGapSizRec = motifGapSizRec + fixedInt(gapsize, 3);
            gapsize = 0;

            seqRec = seqRec + string(1, AArec) + degreeSuffix;
            segmentedRec = segmentedRec + string(1, AArec) + degreeSuffix;
            motifRec = motifRec + string(1, 'X') + degreeSuffix;
            aggrRec = aggrRec + string(1, AArec) + degreeSuffix;
            motifGapSizRec = motifGapSizRec + string(1, 'X'); //+ string(1, AArec);
            //cerr << AArec << "-";
            if(AArec != '?'){
                AAcompRec.at(static_cast<size_t>(AArec - 'A'))++;
            }

        } else {
            gapsize++;
            if((!alreadyInGap) && (seqRec.size() > 0)){
                motifRec = motifRec + string(1, '-');
                aggrRec = aggrRec + string(1, '-');
                if(includeDegree) { // makes it multiple of 2 if including degree
                    motifRec = motifRec + string(1, '-');
                    aggrRec = aggrRec + string(1, '-');
                }
                segmentedRec = segmentedRec + string("]");
            }
            segmentedRec = segmentedRec + string(1, AArec) + degreeSuffix;
            if(seqRec.size() > 0) alreadyInGap = true;
        }
    }
    if(!alreadyInGap) segmentedRec = segmentedRec + string("]");
    if(motifRec.back() == '-') motifRec = motifRec.substr(0, motifRec.size()-1);
    if(motifRec.back() == '-') motifRec = motifRec.substr(0, motifRec.size()-1);
    if(aggrRec.back() == '-') aggrRec = aggrRec.substr(0, aggrRec.size()-1);
    if(aggrRec.back() == '-') aggrRec = aggrRec.substr(0, aggrRec.size()-1);

    stringstream posBound;
    for(size_t i = 0; i < SL; ++i){ // lig
        int IDposInLigand = ligand[static_cast<int>(i)].IDresidue;
        if(interactingPositionsLig[i]){
            posBound << IDposInLigand << " ";
        }
    }

//    int selffolding = 0;
//    for(size_t i = 0; i < S2L; ++i){
//        for(size_t j = i+2; j < S2L; ++j){ // j=i1 are neighbors by definition !
//            if(lattice::areNeighbors(s2.points[i].IDposition, s2.points[j].IDposition)){
//                res.push_back(std::pair<int,int>(i,j));
//                codeInters << (char) ('a' + (char) i) << (char) ('a' + (char) j);
//            }
//        }
//    }

    res[interCodeWithIDpos] = interCodeID.str();
    res[listAAPairs] = listPairs;
    res[seqAGEpitope] = seqLig;
    res[seqABParatope] = seqRec;
    res[motifAGEpitope] = motifLig;
    res[motifABParatope] = motifRec;
    res[motifsSizeGapsLigand] = motifGapSizLig;
    res[motifsSizeGapsRec] = motifGapSizRec;
    res[motifsChemicalLig] = convertChemical(motifLig);
    res[motifsChemicalRec] = convertChemical(motifRec);
    res[agregatesAGEpitope] = aggrLig;
    res[agregatesABParatope] = aggrRec;
    res[chemicalAGEpitope] = convertChemical(aggrLig);
    res[chemicalABParatope] = convertChemical(aggrRec);
    res[positionsBound] = posBound.str();
    res[AAcompoAGEpitope] = codeAAcompo(AAcompLig);
    res[AAcompoABParatope] = codeAAcompo(AAcompRec);
    res[AAcompoFullSlice] = convertAAtoAAusage(s2.getAAseq(), '_');
    res[AAcompoFullCDR] = "UnknownYet-Compute outside structural features please";
    res[interMaskAGEpitope] = interMaskLig;
    res[interMaskABParatope] = interMaskRec;
    res[segmentedAGEpitope] = segmentedLig;
    res[segmentedABParatope] = segmentedRec;
    res[interCodeInternal] = interCodeIDinternal.str();

    return res;
}

//features::features(superProtein* p, int _minDegree, ) : rep(){
//    minDegree = _minDegree;
//    if(p == nullptr) {cerr << "ERR: creating a features() class from NULL superprotein pointer." << endl; return;}
//    currentLigand = p;
//    hotspotsLarge = _hotspotsLarge;

//}

// this is preparing the infos on this antigen before generating features of sequences one by one
features::features(string _antigenID, int _minDegree, int _includeDegreeInMotifs) : rep(){
    minDegree = _minDegree;
    antigenID = _antigenID;
    includeDegreeInMotifs = _includeDegreeInMotifs;

    antigenInfo AG = getAntigenInfos(_antigenID);
    currentLigand = AG.first;
    hotspotsLarge.clear(); //resize(AG.hotspotsLarge.size());
    for(size_t i = 0; i < AG.hotspotsLarge.size(); ++i){
        set<int> toPush = set<int>(AG.hotspotsLarge[i].begin(), AG.hotspotsLarge[i].end());
        hotspotsLarge.push_back(toPush);
    }
    cout << "   ... Preparing features analysis for " << antigenID << " (" << hotspotsLarge.size() << " hotspots known),  with minDegree=" << minDegree << ", " << ((includeDegreeInMotifs)? "": "NOT-") << "including degree in motifs" << endl;
}

string features::titleFeatures(){
     return string("seqAGEpitope\tseqABParatope\tmotifAGEpitope\tmotifABParatope\tagregatesAGEpitope\tagregatesABParatope\tchemicalAGEpitope\tchemicalABParatope\tmotifsSizeGapsLigand\tmotifsSizeGapsRec\tinterMaskAGEpitope\tinterMaskABParatope\tsegmentedAGEpitope\tsegmentedABParatope\tAAcompoFullSlice\tAAcompoFullCDR\tsizeCDR3\tInterCodeInternal\tcontribPerAAepi\tcontribPerAAparaBind\tcontribPerAAparaFold"); // \thotspot_ID");
}

string features::printLine(vector<string> res){
    //enum {interCodeWithIDpos, listAAPairs, AAcompoAGEpitope, AAcompo , positionsBound, NB_features}; //nbneighbors, distChem, selfFolding,
    stringstream res2;

    res2 << res[seqAGEpitope] << "\t" << res[seqABParatope] << "\t" << res[motifAGEpitope] << "\t" << res[motifABParatope] << "\t" << res[agregatesAGEpitope] << "\t" << res[agregatesABParatope] << "\t" << res[chemicalAGEpitope] << "\t" << res[chemicalABParatope] << "\t" << res[motifsSizeGapsLigand] << "\t" << res[motifsSizeGapsRec] << "\t" << res[interMaskAGEpitope] << "\t" << res[interMaskABParatope] << "\t" << res[segmentedAGEpitope] << "\t" << res[segmentedABParatope] << "\t" << res[AAcompoFullSlice] << "\t" << res[AAcompoFullCDR] << "\t" << res[sizeCDR3] << "\t" << res[interCodeInternal] << "\t" << res[contribPerAAepi] << "\t" << res[contribPerAAparaBind] << "\t" << res[contribPerAAparaFold]; // << "\t" << res[hotspot_ID];
    return res2.str();
}

vector<string>features::getProperties(binding &b, string fullSequence){
    std::pair<int,string> str = retrieveStructureFromPosAndStructure(b.structureID, '-');
    //cout << b.structureID << "->" << str.first << "," << str.second << endl;
    return getProperties(b.AAseq,  str.first, str.second, fullSequence);
}

// the file will read, for each AA sequence, the structure and AA acid content of the receptor.
// this function also doesn't give the AA composition of the CDR3 sequence, need to be done outside
vector<string> features::getProperties(string AAligand, int posStart, string structLigand, string fullSequence){
    //cout << "Now analyzing " << AAligand << " and structure " << posStart << "-" << structLigand << endl;
    if(AAligand.size() != structLigand.size() + 1){
        cerr << "ERR: features::getProperties(" << AAligand << ", pos=" << posStart << ", struct=" << structLigand << "),\n     non amtching sizes between structure size (nb of bonds), given=" << structLigand << ", size=" << structLigand.size() << ", and peptide size (nb of AAs, that should bd nb bonds+1) given=" << AAligand.size() << endl;
    }
    struct3D S = struct3D(structLigand, UnDefined, posStart);
    superProtein Rec = superProtein(S);
    Rec.setAAs(AAligand);
    //cout << "The antigen is " << endl;
    //cout << print(*currentLigand);
    //cout << "The receptor is " << endl;
    //cout << print(Rec);
    vector<string> res = structuralFeatures(*currentLigand, Rec, minDegree, includeDegreeInMotifs);
    res[AAcompoFullCDR] = convertAAtoAAusage(fullSequence, '_');
    stringstream toStr; toStr << fullSequence.size();
    res[sizeCDR3] = toStr.str();

    // Hotspots are defined for degree 1... [sorry, bit of waisting time...]
    vector<string> resForHotspots = ((minDegree == 1) ? res : structuralFeatures(*currentLigand, Rec, 1));

    res[hotspot_ID] = "Unknown";
    if(hotspotsLarge.size() > 0){
        bool foundHotspot = false;
        for(size_t i = 0; i < hotspotsLarge.size(); ++i){
            set<int> bound = stringToSet(resForHotspots[positionsBound]);
            set<int> SetHot = hotspotsLarge[i];
            if(isIncluded(SetHot, bound)){
                //cout << print(SetHot) << " included in " << print(bound) << endl;
                stringstream hotname; hotname << antigenID << "_H" << i+1;
                if(foundHotspot){
                    res[hotspot_ID] += hotname.str();
                } else {
                    res[hotspot_ID] = hotname.str();
                }
                foundHotspot = true;
            }
        }
    }
    return res;
}

string concat(vector<string> v, string sep){
    stringstream res;
    for(size_t i = 0; i < v.size(); ++i){
        if(i > 0) res << sep;
        res << v.at(i);
    }
    return res.str();
}



//With C++11 compiler, for non-negative integers I would use something like this (note the :: instead of std::):

bool is_number(const std::string &s) {
    if(s.size() == 0) return false;
    bool has_only_digits = (s.find_first_not_of( "0123456789" ) == string::npos);
    return has_only_digits;
    //return !s.empty() && std::all_of(s.begin(), s.end(), ::isdigit);
}

// this is an outdated function, no hotspot
string features::getPropertiesDeprecated(string fileName, bool display, int nToDisplay){
    cout << "WRN: the function features::getProperties(string fileName) is deprecated. please read the file from the fileformats functions and use getProperties(line) - see in delimain.cpp option getFeatures, for how to do." << endl;
    exampleSeq = "";
    exampleEnergy = 0;
    exampleFirstInterCode = "";
    exampleFirstStruct = "";

    // file format to be read: AAsequence, Energy, InteractionProfile (for checking), nbStructures, StartPos-Structure1, StartPos-Structure2 ...
    stringstream res;
    string sep = "\t";
    ifstream fr;
    fr.open(fileName);
    if(!fr) cerr << "ERR: " << fileName << ", file not found." << endl;
    int lineNr = 0;

    while(fr.good()){
        lineNr++;
        string AAseq;

        double E;
        string interCode;
        int nS;
        int startPos;
        string structure;
        if(fr >> AAseq){

            rep.addSequence(AAseq, -1, true, true);

            if(is_number(AAseq)){
                //cerr << AAseq << " is number " << endl;
                int testFormat = stoi(AAseq);
                if((lineNr == 1) && (testFormat > -1) && (testFormat < 1e10)){
                    cout << "   ... Recognized alternate format 2 for binding datasets\n";
                    return string();
                } else {
                    cout << "   ... ERR in the formatting. got a weird integer as AA sequence in the first line instead of an AA sequence (format 1)" << endl;
                    return string();
                }
            }

            fr >> E >> interCode >> nS;
            //cout << AAseq << "\t" << E << "\t" << interCode << "\t" << nS << endl;
            if((nS < 1) || (nS > 1000)){cerr << "ERR: features::getProperties, getting a too high number of structures (" << nS <<"), around line " << lineNr << ", from reading file " << fileName << endl; return res.str();}
            for(int i = 0; i < nS; ++i){
                fr >> startPos >> structure;

                // for last line, stupid iostream reads empty stuff
                if(AAseq.size() > 0){
#ifdef ALLOW_GRAPHICS
                if(display && (lineNr <= nToDisplay) && (lineNr < 10000)){
                    addToDisplay(new struct3D(structure, UnDefined, startPos));
                }
#endif
                // normally, will just take the last structure, because the motifs will be the same.
                //res << AAseq << sep << E << sep << concat(getProperties(AAseq, startPos, structure), sep) << endl;
                // 4 Valentin:
                vector<string> vres = getProperties(AAseq, startPos, structure); // now we can add full sequence here

                if((lineNr == 1) && (i == 0)){
                    exampleSeq = AAseq;
                    exampleEnergy = E;
                    exampleFirstInterCode = interCode;
                    exampleFirstStruct = structure;
                }


                if(i == 0) res << AAseq << sep << E << sep << interCode << sep << vres[interCodeWithIDpos] << sep << startPos << sep << structure << endl;
                }
                //cout << vres[interCodeWithIDpos] << endl;
                //cout << AAseq << sep << E << sep << concat(getProperties(AAseq, startPos, structure), sep) << endl;
            }
        }
    }
    return res.str();
}

string features::getPropertiesFormat2Deprecated(string fileName, bool display, int nToDisplay){
    exampleSeq = "";
    exampleEnergy = 0;
    exampleFirstInterCode = "";
    exampleFirstStruct = "";

    // file format to be read: AAsequence, Energy, InteractionProfile (for checking), nbStructures, StartPos-Structure1, StartPos-Structure2 ...
    stringstream res;
    string sep = "\t";
    ifstream fr(fileName);
    if(!fr) cerr << "ERR: " << fileName << ", file not found." << endl;
    int lineNr = 0;
    while(fr.good()){
        lineNr++;
        string AAseq;
        double E;
        string interCode;
        int nS;
        int startPos;
        string structure;
        int ID;
        fr >> ID;
        string largeStr;
        fr >> largeStr;
        //cout << "ID=" << ID << endl;
        int lineNr = 0;
        if(fr >> AAseq){


            rep.addSequence(AAseq, ID, true, true);
            lineNr++;
            fr >> E >> nS;

            if(is_number(AAseq)) {
                cerr << "ERR comes from:" << nS << " structures, ID=" << ID << ", AAseq " << AAseq << endl;
            }

            //cout << AAseq << "\t" << E << "\t" << interCode << "\t" << nS << endl;
            if((nS < 1) || (nS > 1000)){cerr << "ERR: features::getProperties, getting a too high number of structures (" << nS <<"), around line " << lineNr << ", from reading file " << fileName << endl; return res.str();}
            for(int i = 0; i < nS; ++i){
                string interCode;
                fr >> startPos >> structure >> interCode;
                //cerr << AAseq << " " << startPos << " " << structure << " " << nS << endl;
#ifdef ALLOW_GRAPHICS
                if(display && (lineNr <= nToDisplay) && (lineNr < 10000)){
                    addToDisplay(new struct3D(structure, UnDefined, startPos));
                }
#endif

                if((lineNr == 1) && (i == 0)){
                    exampleSeq = AAseq;
                    exampleEnergy = E;
                    exampleFirstInterCode = interCode;
                    exampleFirstStruct = structure;
                }


                if(nS >= 2){
                    fr >> interCode; // bug in the dataset
                    if(is_number(interCode)) {
                        cout << "Wrong formatting for sequence " << AAseq << ", aborting alternative structures " << endl;
                        i = nS;
                        char buf[10001];
                        fr.getline(buf, 10000);
                    }
                }
                if(nS >= 3){
                    fr >> interCode; // bug in the dataset
                    if(is_number(interCode)) {
                        cout << "Wrong formatting for sequence " << AAseq << ", aborting alternative structures " << endl;
                        i = nS;
                        char buf[10001];
                        fr.getline(buf, 10000);
                    }
                }
                if(nS >= 4){
                    fr >> interCode; // bug in the dataset
                    if(is_number(interCode)) {
                        cout << "Wrong formatting for sequence " << AAseq << ", aborting alternative structures " << endl;
                        i = nS;
                        char buf[10001];
                        fr.getline(buf, 10000);
                    }
                }
                if(nS >= 5){
                    fr >> interCode; // bug in the dataset
                    if(is_number(interCode)) {
                        cout << "Wrong formatting for sequence " << AAseq << ", aborting alternative structures " << endl;
                        i = nS;
                        char buf[10001];
                        fr.getline(buf, 10000);
                    }
                }
                if(nS >= 6){
                    fr >> interCode; // bug in the dataset
                    if(is_number(interCode)) {
                        cout << "Wrong formatting for sequence " << AAseq << ", aborting alternative structures " << endl;
                        i = nS;
                        char buf[10001];
                        fr.getline(buf, 10000);
                    }
                }
                // normally, will just take the last structure, because the motifs will be the same.
                //res << AAseq << sep << E << sep << concat(getProperties(AAseq, startPos, structure), sep) << endl;
                // 4 Valentin:

                vector<string> vres = getProperties(AAseq, startPos, structure);
                //if(i == 0)
                res << ID << sep << largeStr << sep << AAseq << sep << E << sep << interCode << sep << vres[interCodeWithIDpos] << sep << startPos << sep << structure; //<< endl;
                res << "\t" << vres[seqAGEpitope] << "\t" << vres[seqABParatope] << "\t" << vres[motifAGEpitope] << "\t" << vres[motifABParatope] << "\t" << vres[motifsSizeGapsLigand] << "\t" << vres[motifsSizeGapsRec] << "\t" << vres[motifsChemicalLig] << "\t" << vres[motifsChemicalRec] << "\t" << vres[agregatesAGEpitope] << "\t" << vres[agregatesABParatope] << "\t" << vres[chemicalAGEpitope] << "\t" << vres[chemicalABParatope] << "\t" << vres[positionsBound] << endl;
                //cout << vres[interCodeWithIDpos] << endl;
                //cout << AAseq << sep << E << sep << concat(getProperties(AAseq, startPos, structure), sep) << endl;
            }
        }
    }
    return res.str();
}


/*
 * #include <iostream>
#include <regex>
using namespace std;

int main() {
    std::string pattern("A|D");         // Regex expression
    std::regex rx(pattern);             // Getting the regex object

    std::string s("ABCDEABCABD");       // Defining the string input
    std::ptrdiff_t number_of_matches = std::distance(  // Count the number of matches inside the iterator
        std::sregex_iterator(s.begin(), s.end(), rx),
        std::sregex_iterator());

    std::cout << number_of_matches << std::endl;  // Displaying results
    return 0;
}

#include <iostream>
#include <dirent.h>
#include <sys/types.h>

using namespace std;
void list_dir(const char *path) {
   struct dirent *entry;
   DIR *dir = opendir(path);

   if (dir == nullptr) {
      return;
   }
   while ((entry = readdir(dir)) != nullptr) {
   cout << entry->d_name << endl;
   }
   closedir(dir);
}
int main() {
   list_dir("/home/username/Documents");
}
*/




