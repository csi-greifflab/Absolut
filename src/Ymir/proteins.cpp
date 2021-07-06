

#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cmath>
#include <set>
#include <map>
#include <cstring> // strcpy
#include "../Tools/zaprandom.h"
using namespace std;

#include "proteins.h"
#include "lattice.h"

#include "plot3d.h"

double AAaffinity(AA a1, AA a2){
        // Reference: Miyazawa 1996: Residue-Residue potentials with a favorable contact pair term and unfavorable high packing density term
        // in this matrix, the upper half (right) are the energies in RT units. The lower part (eij) is the difference energy from eii+eji to eij only.
        // We only use the upper half part of the matrix.
//        Cys   Met  Phe Ile  Leu  Val  Trp  Tyr  Ala  Gly  Thr  Ser  Asn  Gln  Asp  Glu  His  Arg  Lys  Pro
//    Cys −5.44 −4.99 −5.80 −5.50 −5.83 −4.96 −4.95 −4.16 −3.57 −3.16 −3.11 −2.86 −2.59 −2.85 −2.41 −2.27 −3.60 −2.57 −1.95 −3.07 Cys
//    Met 0.46 −5.46 −6.56 −6.02 −6.41 −5.32 −5.55 −4.91 −3.94 −3.39 −3.51 −3.03 −2.95 −3.30 −2.57 −2.89 −3.98 −3.12 −2.48 −3.45 Me
//    Phe 0.54 −0.20 −7.26 −6.84 −7.28 −6.29 −6.16 −5.66 −4.81 −4.13 −4.28 −4.02 −3.75 −4.10 −3.48 −3.56 −4.77 −3.98 −3.36 −4.25 Phe
//    Ile 0.49 −0.01 0.06 −6.54 −7.04 −6.05 −5.78 −5.25 −4.58 −3.78 −4.03 −3.52 −3.24 −3.67 −3.17 −3.27 −4.14 −3.63 −3.01 −3.76 Ile
//    Leu 0.57 0.01 0.03 −0.08 −7.37 −6.48 −6.14 −5.67 −4.91 −4.16 −4.34 −3.92 −3.74 −4.04 −3.40 −3.59 −4.54 −4.03 −3.37 −4.20 Leu
//    Val 0.52 0.18 0.10 −0.01 −0.04 −5.52 −5.18 −4.62 −4.04 −3.38 −3.46 −3.05 −2.83 −3.07 −2.48 −2.67 −3.58 −3.07 −2.49 −3.32 Val
//    Trp 0.30 −0.29 0.00 0.02 0.08 0.11 −5.06 −4.66 −3.82 −3.42 −3.22 −2.99 −3.07 −3.11 −2.84 −2.99 −3.98 −3.41 −2.69 −3.73 Trp
//    Tyr 0.64 −0.10 0.05 0.11 0.10 0.23 −0.04 −4.17 −3.36 −3.01 −3.01 −2.78 −2.76 −2.97 −2.76 −2.79 −3.52 −3.16 −2.60 −3.19 Tyr
//    Ala 0.51 0.15 0.17 0.05 0.13 0.08 0.07 0.09 −2.72 −2.31 −2.32 −2.01 −1.84 −1.89 −1.70 −1.51 −2.41 −1.83 −1.31 −2.03 Ala
//    Gly 0.68 0.46 0.62 0.62 0.65 0.51 0.24 0.20 0.18 −2.24 −2.08 −1.82 −1.74 −1.66 −1.59 −1.22 −2.15 −1.72 −1.15 −1.87 Gly
//    Thr 0.67 0.28 0.41 0.30 0.40 0.36 0.37 0.13 0.10 0.10 −2.12 −1.96 −1.88 −1.90 −1.80 −1.74 −2.42 −1.90 −1.31 −1.90 Thr
//    Ser 0.69 0.53 0.44 0.59 0.60 0.55 0.38 0.14 0.18 0.14 −0.06 −1.67 −1.58 −1.49 −1.63 −1.48 −2.11 −1.62 −1.05 −1.57 Ser
//    Asn 0.97 0.62 0.72 0.87 0.79 0.77 0.30 0.17 0.36 0.22 0.02 0.10 −1.68 −1.71 −1.68 −1.51 −2.08 −1.64 −1.21 −1.53 Asn
//    Gln 0.64 0.20 0.30 0.37 0.42 0.46 0.19 −0.12 0.24 0.24 −0.08 0.11 −0.10 −1.54 −1.46 −1.42 −1.98 −1.80 −1.29 −1.73 Gln
//    Asp 0.91 0.77 0.75 0.71 0.89 0.89 0.30 −0.07 0.26 0.13 −0.14 −0.19 −0.24 −0.09 −1.21 −1.02 −2.32 −2.29 −1.68 −1.33 Asp
//    Glu 0.91 0.30 0.52 0.46 0.55 0.55 0.00 −0.25 0.30 0.36 −0.22 −0.19 −0.21 −0.19 0.05 −0.91 −2.15 −2.27 −1.80 −1.26 Glu
//    His 0.65 0.28 0.39 0.66 0.67 0.70 0.08 0.09 0.47 0.50 0.16 0.26 0.29 0.31 −0.19 −0.16 −3.05 −2.16 −1.35 −2.25 His
//    Arg 0.93 0.38 0.42 0.41 0.43 0.47 −0.11 −0.30 0.30 0.18 −0.07 −0.01 −0.02 −0.26 −0.91 −1.04 0.14 −1.55 −0.59 −1.70 Arg
//    Lys 0.83 0.31 0.33 0.32 0.37 0.33 −0.10 −0.46 0.11 0.03 −0.19 −0.15 −0.30 −0.46 −1.01 −1.28 0.23 0.24 −0.12 −0.97 Lys
//    Pro 0.53 0.16 0.25 0.39 0.35 0.31 −0.33 −0.23 0.20 0.13 0.04 0.14 0.18 −0.08 0.14 0.07 0.15 −0.05 −0.04 −1.75 Pro

    // The ENUM in proteins.h has been defined in the same order, very important !!
    //                             Cys    Met    Phe    Ile    Leu    Val    Trp    Tyr    Ala    Gly    Thr    Ser    Asn    Gln    Asp    Glu    His    Arg    Lys    Pro
    static vector<double> cCys = {-5.44, -4.99, -5.80, -5.50, -5.83, -4.96, -4.95, -4.16, -3.57, -3.16, -3.11, -2.86, -2.59, -2.85, -2.41, -2.27, -3.60, -2.57, -1.95, -3.07}; // Cys
    static vector<double> cMet = {-4.99, -5.46, -6.56, -6.02, -6.41, -5.32, -5.55, -4.91, -3.94, -3.39, -3.51, -3.03, -2.95, -3.30, -2.57, -2.89, -3.98, -3.12, -2.48, -3.45}; // Me
    static vector<double> cPhe = {-5.80, -6.56, -7.26, -6.84, -7.28, -6.29, -6.16, -5.66, -4.81, -4.13, -4.28, -4.02, -3.75, -4.10, -3.48, -3.56, -4.77, -3.98, -3.36, -4.25}; // Phe
    static vector<double> cIle = {-5.50, -6.02, -6.84, -6.54, -7.04, -6.05, -5.78, -5.25, -4.58, -3.78, -4.03, -3.52, -3.24, -3.67, -3.17, -3.27, -4.14, -3.63, -3.01, -3.76}; // Ile
    static vector<double> cLeu = {-5.83, -6.41, -7.28, -7.04, -7.37, -6.48, -6.14, -5.67, -4.91, -4.16, -4.34, -3.92, -3.74, -4.04, -3.40, -3.59, -4.54, -4.03, -3.37, -4.20}; // Leu
    static vector<double> cVal = {-4.96, -5.32, -6.29, -6.05, -6.48, -5.52, -5.18, -4.62, -4.04, -3.38, -3.46, -3.05, -2.83, -3.07, -2.48, -2.67, -3.58, -3.07, -2.49, -3.32}; // Val
    static vector<double> cTrp = {-4.95, -5.55, -6.16, -5.78, -6.14, -5.18, -5.06, -4.66, -3.82, -3.42, -3.22, -2.99, -3.07, -3.11, -2.84, -2.99, -3.98, -3.41, -2.69, -3.73}; // Trp
    static vector<double> cTyr = {-4.16, -4.91, -5.66, -5.25, -5.67, -4.62, -4.66, -4.17, -3.36, -3.01, -3.01, -2.78, -2.76, -2.97, -2.76, -2.79, -3.52, -3.16, -2.60, -3.19}; // Tyr
    static vector<double> cAla = {-3.57, -3.94, -4.81, -4.58, -4.91, -4.04, -3.82, -3.36, -2.72, -2.31, -2.32, -2.01, -1.84, -1.89, -1.70, -1.51, -2.41, -1.83, -1.31, -2.03}; // Ala
    static vector<double> cGly = {-3.16, -3.39, -4.13, -3.78, -4.16, -3.38, -3.42, -3.01, -2.31, -2.24, -2.08, -1.82, -1.74, -1.66, -1.59, -1.22, -2.15, -1.72, -1.15, -1.87}; // Gly
    static vector<double> cThr = {-3.11, -3.51, -4.28, -4.03, -4.34, -3.46, -3.22, -3.01, -2.32, -2.08, -2.12, -1.96, -1.88, -1.90, -1.80, -1.74, -2.42, -1.90, -1.31, -1.90}; // Thr
    static vector<double> cSer = {-2.86, -3.03, -4.02, -3.52, -3.92, -3.05, -2.99, -2.78, -2.01, -1.82, -1.96, -1.67, -1.58, -1.49, -1.63, -1.48, -2.11, -1.62, -1.05, -1.57}; // Ser
    static vector<double> cAsn = {-2.59, -2.95, -3.75, -3.24, -3.74, -2.83, -3.07, -2.76, -1.84, -1.74, -1.88, -1.58, -1.68, -1.71, -1.68, -1.51, -2.08, -1.64, -1.21, -1.53}; // Asn
    static vector<double> cGln = {-2.85, -3.30, -4.10, -3.67, -4.04, -3.07, -3.11, -2.97, -1.89, -1.66, -1.90, -1.49, -1.71, -1.54, -1.46, -1.42, -1.98, -1.80, -1.29, -1.73}; // Gln
    static vector<double> cAsp = {-2.41, -2.57, -3.48, -3.17, -3.40, -2.48, -2.84, -2.76, -1.70, -1.59, -1.80, -1.63, -1.68, -1.46, -1.21, -1.02, -2.32, -2.29, -1.68, -1.33}; // Asp
    static vector<double> cGlu = {-2.27, -2.89, -3.56, -3.27, -3.59, -2.67, -2.99, -2.79, -1.51, -1.22, -1.74, -1.48, -1.51, -1.42, -1.02, -0.91, -2.15, -2.27, -1.80, -1.26}; // Glu
    static vector<double> cHis = {-3.60, -3.98, -4.77, -4.14, -4.54, -3.58, -3.98, -3.52, -2.41, -2.15, -2.42, -2.11, -2.08, -1.98, -2.32, -2.15, -3.05, -2.16, -1.35, -2.25}; // His
    static vector<double> cArg = {-2.57, -3.12, -3.98, -3.63, -4.03, -3.07, -3.41, -3.16, -1.83, -1.72, -1.90, -1.62, -1.64, -1.80, -2.29, -2.27, -2.16, -1.55, -0.59, -1.70}; // Arg
    static vector<double> cLys = {-1.95, -2.48, -3.36, -3.01, -3.37, -2.49, -2.69, -2.60, -1.31, -1.15, -1.31, -1.05, -1.21, -1.29, -1.68, -1.80, -1.35, -0.59, -0.12, -0.97}; // Lys
    static vector<double> cPro = {-3.07, -3.45, -4.25, -3.76, -4.20, -3.32, -3.73, -3.19, -2.03, -1.87, -1.90, -1.57, -1.53, -1.73, -1.33, -1.26, -2.25, -1.70, -0.97, -1.75}; // Pro
    static vector<vector<double>> matrix = {cCys, cMet, cPhe, cIle, cLeu, cVal, cTrp, cTyr, cAla, cGly, cThr, cSer, cAsn, cGln, cAsp, cGlu, cHis, cArg, cLys, cPro};
    if((a1 == UndefinedYet) || (a1 == NB_AAs) || (a2 == UndefinedYet) || (a2 == NB_AAs)) return NAN;
    return matrix[min(a1, a2)][max(a1,a2)]; // the upper right triangle
}

void testAAaffinities(){
    cout << "MetIle - -6.02 " << AAaffinity(Met, Ile) << " - " << AAaffinity(Ile, Met) << endl;
    cout << "LeuLeu - -7.37 " << AAaffinity(Leu, Leu) << endl;
    cout << "ArgPhe - -3.98 " << AAaffinity(Arg, Phe) << " - " << AAaffinity(Arg, Phe) << endl;
}

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
    default: {
        //cerr << "Char " << AA_name << " is not an AA" << endl;
        return NB_AAs;
        break;
    }
    }
    return NB_AAs;
}

char AAname(AA aa){
    switch(aa){
    case Cys: return 'C';
    case Met: return 'M';
    case Phe: return 'F';
    case Ile: return 'I';
    case Leu: return 'L';
    case Val: return 'V';
    case Trp: return 'W';
    case Tyr: return 'Y';
    case Ala: return 'A';
    case Gly: return 'G';
    case Thr: return 'T';
    case Ser: return 'S';
    case Asn: return 'N';
    case Gln: return 'Q';
    case Asp: return 'D';
    case Glu: return 'E';
    case His: return 'H';
    case Arg: return 'R';
    case Lys: return 'K';
    case Pro: return 'P';
    case UndefinedYet: return '?';
    case NB_AAs: return '!';
    }
    return ' ';;
}


bool correctAA(char AA_name){
    AA res = AA_ID(AA_name);
    return(!((res == UndefinedYet) || (res == NB_AAs)));
}

bool checkAAseq(string s){
    for(unsigned int i = 0; i < s.size(); ++i){
        if((AA_ID(s[i]) ==  UndefinedYet) || (AA_ID(s[i]) == NB_AAs)){
            cerr << "Incorrect AA sequence: " << s << ", because contains " << s[i] << endl;
            return false;
        }
    }
    return true;
}

string randomProt(int size){
    string allAAs = "CMFILVWYAGTSNQDEHRKP";
    string res = string(size, '?');
    for(int i = 0; i < size; ++i){
        int newAAid = random::uniformInteger(0,19);
        res[i] = allAAs[newAAid];
    }
    return res;
}

char randomAA(){
    string allAAs = "CMFILVWYAGTSNQDEHRKP";
    int newAAid = random::uniformInteger(0,19);
    return allAAs[newAAid];
}






int getNewID(int channel){
    static vector<int> cpt = vector<int>(1, -1);
    if(channel < 0) return -1;
    if(channel >= (int) cpt.size()){
        cpt.resize(channel+1, -1); // hope it fills -1 to only the new positions
    }
    return ++cpt[channel];
}

string print(residue &a){
    stringstream res;
    vector<int> position = lattice::positionFromID(a.IDposition);
    res << a.IDresidue << "\t" << AAname(a.TypeResidue) << " (" << a.TypeResidue << ")\t" << a.IDposition << "\t" <<  position[0] << "\t" << position[1] << "\t" << position[2];
    return res.str();
}

#define startingPos1 0
#define startingPos2 1

// the important storage inside the structure is occupiedPositions.
superProtein::superProtein() {
    structure = new struct3D(); // note: never ever copy an external pointer into structure because the destructor would delete the external pointer then.
    ID = getNewID();
    points.clear();
}

superProtein::~superProtein(){
    //cout << "Deleting " << structure << " Prot of size " << this->size() << endl; //print(*this) << endl;

    // I still don't understand why this creates a segfault when calling:
    if(structure) delete structure; // for instance if say superProtein p1, then p1 = p2, p1 will be destroyed. might have mo structure?
}

superProtein::superProtein(const struct3D& source){
    //cout << "Copy called" << endl;
    structure = new struct3D(source);
    ID = getNewID();
    createPoints();
}

superProtein::superProtein(string absoluteSequence, int initialPos) {
    structure = new struct3D(absoluteSequence, UnDefined, (initialPos == -1) ? lattice::centralPosition() : initialPos);
    ID = getNewID();
    createPoints();
    // reconstitutes the structure point by point. This is redundant with the struct3D constructor (sorry!)
    // I wanted to keep the struct3D constructor as simple as possible for computations of millions of structures
    // where we don't really care about the order of positions, just the starting and ending ones.
}

// creates points from the structure to the list of residues points[]
void superProtein::createPoints(){
    if(structure == nullptr) {cerr << "superProtein::createPoints, structure = NULL" << endl; return;}
    int cptIDresidue = 0;
    points.clear();
    points.push_back(residue(structure->startingPosition, cptIDresidue));
    vector<int> currentPos = lattice::positionFromID(structure->startingPosition);
    if(structure->sequence.size() == 0) return;     // this is fine, no need to raise an error.

    cptIDresidue++;
    moveDirection currentDir = (moveDirection) charMoveToInt(structure->sequence[0]);
    moveDirection currentYaxis = initialYaxis((moveDirection) charMoveToInt(structure->sequence[0]));
    vector<int> actualMove = moveVector(currentDir);
    currentPos[0] += actualMove[0];
    currentPos[1] += actualMove[1];
    currentPos[2] += actualMove[2];
    int IDcurrentPos = lattice::idFromPosisition(currentPos);
    points.push_back(residue(IDcurrentPos, cptIDresidue));
    cptIDresidue++;

    for(int i = 1; i < (int) structure->sequence.size(); ++i){
        std::pair<moveDirection, moveDirection> next = nextAbsoluteMove(currentDir, currentYaxis, (moveDirection) charMoveToInt(structure->sequence[i]));
        currentDir = next.first;
        currentYaxis = next.second;
        vector<int> actualMove = moveVector(currentDir);
        currentPos[0] += actualMove[0];
        currentPos[1] += actualMove[1];
        currentPos[2] += actualMove[2];
        IDcurrentPos = lattice::idFromPosisition(currentPos);
        points.push_back(residue(IDcurrentPos, cptIDresidue));
        cptIDresidue++;
    }
}


int superProtein::size(){
    return points.size();
}

bool superProtein::isFree(int idPosition){
    if((!structure) && (points.size() > 0)) cerr << "superProtein::isFree(), protein with empty structure but has points." << endl;
    if(structure) return((structure->occupiedPositions.find(idPosition)) == structure->occupiedPositions.end());
    return true;
}

residue superProtein::operator[](int _IDresidue){
    if((_IDresidue < 0) || (_IDresidue >= (int) points.size())) {cerr << "superProtein::getResidue(" << _IDresidue << "), out of bounds - maybe empty protein (size " << size() << ")" << endl; return residue(0);}
    return points[_IDresidue];
}
residue* superProtein::end(){
    if(points.size() == 0) return nullptr;
    return &(points[points.size()-1]);
}
residue* superProtein::begin(){
    if(points.size() == 0) return nullptr;
    return &(points[0]);
}

void superProtein::push_back(AA newAA, int IDnewPos){
    if((end() != nullptr) && (!lattice::areNeighbors(end()->IDposition, IDnewPos))) {cerr << "ERR: superProtein::push_back (AA,pos), Inserting a new AA at a non-neighboring position : " << AAname(newAA) << " at position ID " << IDnewPos << endl; return;}
    points.push_back(residue(IDnewPos, points.size(), newAA));
    //points[points.size()-1].IDresidue = points.size()-1;
    if(!structure) cerr << "superProtein::push_back(), protein has structure=NULL" << endl;
    structure->occupiedPositions.insert(IDnewPos);
}

void superProtein::push_back(residue &toAddToTail){ // ID residue automatically computed.
    if((end() != nullptr) && (!lattice::areNeighbors(end()->IDposition, toAddToTail.IDposition))) {cerr << "ERR: superProtein::push_back(residue) Inserting a new AA at a non-neighboring position : " << print(toAddToTail) << " while ending position was " << end()->IDposition << ":" << printVector(lattice::positionFromID(end()->IDposition)) << endl; return;}
    points.push_back(toAddToTail);
    points.back().IDresidue = points.size()-1;
    if(!structure) cerr << "superProtein::push_back(), protein has structure=NULL" << endl;
    structure->occupiedPositions.insert(toAddToTail.IDposition);
}

bool superProtein::contiguous(){
    // cout << sequence << endl;
    if(structure->sequence.find('+') != std::string::npos) return false; // i.e found a '+'
    // cout << "IS CONTIGUOUS" << endl;
    int NP = points.size();
    int lastID = -1e6;
    for(int i = 0; i < NP; ++i){
        int newID = points[i].IDresidue;
        if((i > 0) && (newID != lastID + 1)){
            return false;
        }
        lastID = newID;
    }
    return true;
}

// this adds + when residue IDs are not continuous
string superProtein::getAAseq(){
    stringstream res;
    int NR = points.size();
    int lastID = -1e6;
    for(int i = 0; i < NR; ++i){
        int newID = points[i].IDresidue;
        //if((i > 0) && (newID != lastID + 1)){
            // res << "+";, no, I prefer making sure that the number of residue is the same as the number of AAs in the AAsequence. More consistent
        //}
        res << AAname(points[i].TypeResidue);
        lastID = newID;
    }
    return res.str();
}

vector<int> startingPositions(superProtein* P){
    vector<int> res = vector<int>();
    if((P == nullptr) || (P->structure == nullptr)) return res;
    res.push_back(P->structure->startingPosition);
    if(P->contiguous()) return res;

    size_t NP = P->points.size();
    int lastID = -1e6;
    for(size_t i = 0; i < NP; ++i){
        int newID = P->points[i].IDresidue;
        if((i > 0) && (newID != lastID + 1)){
            res.push_back(P->points[i].IDposition);
        }
        lastID = newID;
    }
    return res;
}

set<int> getOccupiedPositions(superProtein *P){
    set<int> res;
    if(P == nullptr) return res;
    size_t NS = P->points.size();
    for(size_t i = 0; i < NS; ++i){
        res.insert(P->points[i].IDposition);
    }
    return res;
}


// note: this function is not efficient, this is more for outputing
vector<std::pair<int, string> > getSubChains(superProtein* P){

    vector<std::pair<int, string>> res;
    if((P == nullptr) || (P->structure == nullptr)) return res;
    //
    if(P->contiguous()) {
        res.push_back(std::pair<int, string>(P->structure->startingPosition, P->structure->sequence));
        return res;
    }

    size_t NP = P->points.size();

    string fullSeq = P->structure->sequence;
//    if(P->structure->sequence.size() <= NP) {
//        cerr << "ERR: getSubChains(...), the stucture " << P->structure->sequence << ", is too short compared to the nb of AAs " << P->getAAseq() << endl;
//        return res;
//    }
    int pos = 0;
    string buffer = "";
    int lastStartPos = P->structure->startingPosition;

    int lastID = -1e6;
    for(size_t i = 0; i < NP; ++i){


        int newID = P->points[i].IDresidue;
        if((i > 0) && (newID != lastID + 1)){
            res.push_back(std::pair<int,string>( lastStartPos, buffer));
            buffer = "";
            //pos++; // for the + => There are always one more AA.
            lastStartPos = P->points[i].IDposition;
        }
        lastID = newID;
        if(pos < P->structure->sequence.size()) buffer += string(1, P->structure->sequence[pos]);
        pos++;
    }
    res.push_back(std::pair<int,string>( lastStartPos, buffer));
    return res;
}

/*void superProtein::push_front(residue &toAddToFront){ // computationally expensive.
    if(lattice::isNeighbor(toAddToFront.IDposition, begin()->IDposition)) {cerr << "ERR: Inserting a new AA at a non-neighboring position : " << AAname(begin()->TypeResidue) << " at position ID " << begin()->IDposition << endl; return;}
    points.push_front(toAddToFront);
    for(int i = 0; i < points.size(); ++i){
        points[i].IDresidue = i;
    }
    occupiedPositions.insert(toAddToFront.IDposition);
}*/

//void superProtein::push_back(AA newAA, moveDirection newMove){
    // need to get the current plane ...

//        occupiedPositions.insert(IDnewPos);
//}



superProtein::superProtein (const superProtein& toCopy) {
    //: struct3D(toCopy.sequence, (toCopy.listYAxis.size() > 0) ? (moveDirection) charMoveToInt(toCopy.listYAxis[0]) : UnDefined, toCopy.startingPosition)
    structure = new struct3D(*(toCopy.structure));
    ID = getNewID();
    points = toCopy.points;
    // should be done occupiedPositions = toCopy.occupiedPositions;
}

// stupid thing, if I call the function print(superProtein) from inside the function print() of another class,
// it gets confused, so I had to add this function, very stupid.
string printProtein(superProtein &P){
    return print(P);
}

string print(superProtein &P){
    stringstream res;
    res << print(*P.structure);
    res << "Position and ID residues one by one: \n";
    int NP = P.size();
    for(int i = 0; i < NP; ++i){
        res << P[i].IDresidue << "\t" << AAname(P[i].TypeResidue) << "\t" << P[i].IDposition << "(" << printVector(lattice::positionFromID(P[i].IDposition)) << ")\n";
    }
    return res.str();
}


//// Here, use lattice functions, and remove -Width from here
bool excludeCriterion(superProtein* toCheck){
    if(toCheck == nullptr) return true;

    for(int i = 0; i < toCheck->size(); ++i){
        int IDcoord = toCheck->operator [](i).IDposition;
        vector<int> v = lattice::positionFromID(IDcoord);
        int z = v[2]; //IDcoord / ZWidth;
        int y = v[1]; //(IDcoord - ZWidth*z)/YWidth;
        int x = v[0]; //(IDcoord - ZWidth*YWidth*z - YWidth*y);
        //if((abs((XWidth / 2) -x) > 1) || (abs((YWidth / 2)-y) > 1) || (abs((ZWidth / 2)-z) > 1)) return true;
        if((abs(50 -x) > 1) || (abs(50-y) > 1) || (abs(50-z) > 1)) return true;
    }
    return false;
}

// This function doesn't work ???
/*bool excludeCriterion(int IDcoord){
    vector<int> v = lattice::positionFromID(IDcoord);
    int z = v[2]; //IDcoord / ZWidth;
    int y = v[1]; //(IDcoord - ZWidth*z)/YWidth;
    int x = v[0]; //(IDcoord - ZWidth*YWidth*z - YWidth*y);
    //if((abs((XWidth / 2) -x) > 1) || (abs((YWidth / 2)-y) > 1) || (abs((ZWidth / 2)-z) > 1)) return true;
    if((abs(50 -x) > 1) || (abs(50-y) > 1) || (abs(50-z) > 1)) return true;
}*/

ensProts getStartingStructures(int currentSize, vector<AA> sequence){
    ensProts resGroup = ensProts();

    // for a cube, either start at 50 50 50 or at 49 49 49
    {
        superProtein* start = new superProtein(string("S"), lattice::idFromPosisition(50,50,50));
        //nbProts++;
        if(currentSize > 2){
            for(int i = 0; i < currentSize-1; ++i){
                residue sr = residue(lattice::idFromPosisition(50+i,50,50), i, sequence[i]);
                start->push_back(sr);
            }
            residue sr = residue(lattice::idFromPosisition(50+currentSize-2,51,50), currentSize-1, sequence[currentSize-1]);
            start->push_back(sr);
        } else {
            for(int i = 0; i < currentSize; ++i){
                residue sr = residue(lattice::idFromPosisition(50+i,50,50), i, sequence[i]);
                start->push_back(sr);
            }
        }
        resGroup.add(start);
    }

    {
        superProtein* start = new superProtein(string("S"), lattice::idFromPosisition(49,49,49));
        //nbProts++;
        if(currentSize > 2){
            for(int i = 0; i < currentSize-1; ++i){
                residue sr = residue(lattice::idFromPosisition(49+i,49,49), i, sequence[i]);
                start->push_back(sr);
            }
            residue sr = residue(lattice::idFromPosisition(49+currentSize-2,50,49), currentSize-1, sequence[currentSize-1]);
            start->push_back(sr);
        } else {
            for(int i = 0; i < currentSize; ++i){
                residue sr = residue(lattice::idFromPosisition(49+i,49,49), i, sequence[i]);
                start->push_back(sr);
            }
        }
        resGroup.add(start);
    }

    //cout << "Start for Size = " << currentSize << "\t" << resGroup.size() << endl;
    //cout << print(resGroup) << endl;
    return resGroup;
}

ensProts getAllStructures(int _size, vector<AA> sequence){
    ensProts resGroup = ensProts();
    if((int) sequence.size() < _size) {cerr << "making proteins with size " << _size << "with incomplete sequence " << endl; return resGroup;}
    if(_size < 1) return resGroup;

    if(_size > 3){
        ensProts smallerGroup = getAllStructures(_size-1, sequence);
        int smallerSize = smallerGroup.size();
        for(int i = 0; i < smallerSize; ++i){
            superProtein* currentBase = smallerGroup[i];
            if(!currentBase->structure) cerr << "getAllstructures, got proteins with empty structure inside" << endl;
            vector<int> freeNeighb = lattice::idFreeNeighbors(currentBase->end()->IDposition, currentBase->structure->occupiedPositions);
            int nbFree = freeNeighb.size();
            for(int j = 0; j < nbFree; ++j){
                    superProtein* newProt = new superProtein(*currentBase);
                    //nbProts++;
                    newProt->push_back(sequence[_size-1], freeNeighb[j]);
                    if(!excludeCriterion(newProt)) resGroup.add(newProt);
                    else {
                        delete newProt;
                    }

            }
        }
    }

    // In any case the starting structures can be enumerated. Will depend on the stopping constraint
    ensProts newStructToAdd = getStartingStructures(_size, sequence);
    int ns = newStructToAdd.size();
    for(int i = 0; i < ns; ++i){
        superProtein* examinedProt =newStructToAdd[i];
        if(!excludeCriterion(examinedProt)) resGroup.add(examinedProt);
        else {
            //cout << "Executed" << print(*examinedProt) <<  endl;
            delete examinedProt;
            //nbProts--;
        }
    }
    cout << "Size(AAs) = " << _size << "\t" << resGroup.size() << endl;
    return resGroup;
}





/*vector<struct3D*> getStartingStructures(int currentSize){


    // for a cube, either start at 50 50 50 or at 49 49 49
    {
        superProtein* start = new superProtein(string("S"), lattice::idFromPosisition(50,50,50));
        //nbProts++;
        if(currentSize > 2){
            for(int i = 0; i < currentSize-1; ++i){
                residue sr = residue(lattice::idFromPosisition(50+i,50,50), i, sequence[i]);
                start->push_back(sr);
            }
            residue sr = residue(lattice::idFromPosisition(50+currentSize-2,51,50), currentSize-1, sequence[currentSize-1]);
            start->push_back(sr);
        } else {
            for(int i = 0; i < currentSize; ++i){
                residue sr = residue(lattice::idFromPosisition(50+i,50,50), i, sequence[i]);
                start->push_back(sr);
            }
        }
        resGroup.add(start);
    }

    {
        superProtein* start = new superProtein(string("S"), lattice::idFromPosisition(49,49,49));
        //nbProts++;
        if(currentSize > 2){
            for(int i = 0; i < currentSize-1; ++i){
                residue sr = residue(lattice::idFromPosisition(49+i,49,49), i, sequence[i]);
                start->push_back(sr);
            }
            residue sr = residue(lattice::idFromPosisition(49+currentSize-2,50,49), currentSize-1, sequence[currentSize-1]);
            start->push_back(sr);
        } else {
            for(int i = 0; i < currentSize; ++i){
                residue sr = residue(lattice::idFromPosisition(49+i,49,49), i, sequence[i]);
                start->push_back(sr);
            }
        }
        resGroup.add(start);
    }

    //cout << "Start for Size = " << currentSize << "\t" << resGroup.size() << endl;
    //cout << print(resGroup) << endl;
    return resGroup;
}

ensProts getAllStructures(int _size, vector<AA> sequence){
    ensProts resGroup = ensProts();
    if((int) sequence.size() < _size) {cerr << "making proteins with size " << _size << "with incomplete sequence " << endl; return resGroup;}
    if(_size < 1) return resGroup;

    if(_size > 3){
        ensProts smallerGroup = getAllStructures(_size-1, sequence);
        int smallerSize = smallerGroup.size();
        for(int i = 0; i < smallerSize; ++i){
            superProtein* currentBase = smallerGroup[i];
            vector<int> freeNeighb = lattice::idFreeNeighbors(currentBase->end()->IDposition, *currentBase);
            int nbFree = freeNeighb.size();
            for(int j = 0; j < nbFree; ++j){
                    superProtein* newProt = new superProtein(*currentBase);
                    //nbProts++;
                    newProt->push_back(sequence[_size-1], freeNeighb[j]);
                    if(!excludeCriterion(newProt)) resGroup.add(newProt);
                    else {
                        delete newProt;
                    }

            }
        }
    }

    // In any case the starting structures can be enumerated. Will depend on the stopping constraint
    ensProts newStructToAdd = getStartingStructures(_size, sequence);
    int ns = newStructToAdd.size();
    for(int i = 0; i < ns; ++i){
        superProtein* examinedProt =newStructToAdd[i];
        if(!excludeCriterion(examinedProt)) resGroup.add(examinedProt);
        else {
            //cout << "Executed" << print(*examinedProt) <<  endl;
            delete examinedProt;
            //nbProts--;
        }
    }
    cout << "Size(AAs) = " << _size << "\t" << resGroup.size() << endl;
    return resGroup;
}*/







//string getAAs(struct3D* s){
//    if(!s) return string("");
//    int L = s->points.size();
//    string res = string(L, ' ');
//    for(int i = 0; i < L; ++i){
//        res[i] = s->points[i].TypeResidue;
//    }
//    return res;
//}

void superProtein::setAAs(string seqAAs){
    int L = points.size();
    if((int) seqAAs.size() != L){
        cerr << "ERR: setAAs, not same length ! Points are " << points.size() << " while you give " << seqAAs.size() << " AAs \n  Details: " << seqAAs << " on structure " << print(*this) << endl;
        return;
    }
    for(int i = 0; i < L; ++i){
        if(correctAA(seqAAs[i])){
            points[i].TypeResidue = AA_ID(seqAAs[i]);
        } else {
            cerr << "ERR: setAAs(" << seqAAs << "), incorrect AA: " << seqAAs[i] << endl;
        }
    }
}





string print(ensProts ep){
    stringstream res;
    for(int i = 0; i < ep.size(); ++i){
        res << print(*(ep.list[i])) << "\n\n";
    }
    return res.str();
}



bool collide(superProtein& s1, struct3D& s2){
    std::vector<int> common_points;
    if(!s1.structure){
        cerr << "ERR: collide(s1,s2), one protein has NULL structure" << endl;
        return false;
    }
    set_intersection(s1.structure->occupiedPositions.begin(),s1.structure->occupiedPositions.end(),s2.occupiedPositions.begin(),s2.occupiedPositions.end(), std::back_inserter(common_points));
    return(common_points.size() > 0);
}

bool collide(superProtein& s1, superProtein& s2){
    std::vector<int> common_points;
    if((!s1.structure) || (!s2.structure)){
        cerr << "ERR: collide(s1,s2), one protein has NULL structure" << endl;
        return false;
    }
    set_intersection(s1.structure->occupiedPositions.begin(),s1.structure->occupiedPositions.end(),s2.structure->occupiedPositions.begin(),s2.structure->occupiedPositions.end(), std::back_inserter(common_points));
    return(common_points.size() > 0);
}


void testProteins(){

    #ifdef ALLOW_GRAPHICS
    if(false){
        cout << "Test copy and pushBackAbsoluteMove" << endl;
        char *c[] = {(char*)"Hello",nullptr};
        glDisplay(0,c);

        struct3D* a = new struct3D(string("DDLUSRSLDURLSLL"), Left, lattice::centralPosition());
        cerr << print(*a);
        addToDisplay(a, false); // need to be a new XXX to be still existing after the function ends.
        for(int i = 0; i < Nb_Moves_Absolute; ++i){
            struct3D* b = new struct3D(*a);
            if(i == 0) {
                cerr << "Copy ..." << endl;
                cerr << print(*b) << endl;
            }
            if(b->pushBackAbsoluteMove((moveDirection) i)){
                if(b->properlyFolded){
                    cerr << "Successfully added absolute " << intToMoveChar(i) << endl;
                    cerr << print(*b) << endl;
                    addToDisplay(b, false); // need to be a new XXX to be still existing after the function ends.
                }
            }
        }
        glutMainLoop();
    }
    #endif

    cout << "Test revert " << endl;
    int pos = lattice::idFromPosisition(30,30,30);
    struct3D* a = new struct3D(string("DDLUSRSLDURLS"), UnDefined, pos);
    string test = revert(a->sequence);
    int endPos = a->endingPosition;
    struct3D* b = new struct3D(test, UnDefined, endPos);
    cout << "Reverting " << a->sequence << ", with ABS " << a->separatedSingleMoves << " from pos " << printVector(lattice::positionFromID(pos));
    cout << "\n     into " << b->sequence << ", with ABS " << b->separatedSingleMoves << " from pos " << printVector(lattice::positionFromID(endPos)) << endl;
    #ifdef ALLOW_GRAPHICS
    if(true){
        char *c[] = {(char*)"Hello",nullptr};
        glDisplay(0,c);
        addToDisplay(a, false); // need to be a new XXX to be still existing after the function ends.
        addToDisplay(b, false); // need to be a new XXX to be still existing after the function ends.
        glutMainLoop();
    }
    #endif

    return;
    // only on freeglut: glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,GLUT_ACTION_GLUTMAINLOOP_RETURNS);
//    cout << "Test struct3D" << endl;
//    struct3D* s1 = new struct3D(string("BRRURD"));
//    struct3D* s1b = new struct3D(string("BRRURDULLSRLSDUSSSUD"), Nb_Moves, lattice::idFromPosisition(54,54,54));
//    struct3D* s1c = new struct3D(string("BRRRRRR"));
//    cout << print(*s1) << endl;
//    glDisplay(0,c);
//    addToDisplay(s1); // need to be a new XXX to be still existing after the function ends.
//    addToDisplay(s1b); // need to be a new XXX to be still existing after the function ends.
//    glutMainLoop();

    struct3D s = struct3D(string("BRRURDULLSRLSDUSSSUD"));
    cout << "Sequence s:" << print(s) << endl;
    cout << "Testing contains() for structure s1" << endl;
    vector<int> positions = {11201, lattice::centralPosition(), lattice::idFromPosisition(31,32,33), lattice::idFromPosisition(31,32,32)};
    for(unsigned int i = 0; i < positions.size(); ++i){
        cout << ((contains(s.occupiedPositions, positions[i])) ? "Contains         " : "Does not contain ") << positions[i] << " -> " << printVector(lattice::positionFromID(positions[i])) << endl;
    }
    cout << endl;

    // Examples for collide, with s1
    // SSURLDS from 31 32 31 => collide, 1 pt
    // SSRRUDLSDS from 31 32 31=> do not collide
    // SUULUDLSDS from 31 32 32 => collides 2 times
    struct3D s2 = struct3D(string("SUULUDLSDS"), UnDefined, lattice::idFromPosisition(31,32,31));
    cout << "Sequence s2:" << print(s2) << endl;
    cout << "These sequences " << ((collide(s,s2)) ? "Collide " : "Do not collide ") << endl;
    vector<int> res = intersection(s, s2);
    cout << "Points of intersection :" << res.size() << endl;
    for(unsigned int i = 0; i < res.size(); ++i){
        cout << "   " << res[i] << " -> " << printVector(lattice::positionFromID(res[i])) << endl;
    }

    cout << "Test of the touch function (for s2)" << endl;
    vector<int> positions2 = {11201, lattice::centralPosition(), lattice::idFromPosisition(31,32,33), lattice::idFromPosisition(31,32,32), lattice::idFromPosisition(32,31,30)};
    for(unsigned int i = 0; i < positions2.size(); ++i){
        cout << ((touch(s2, positions2[i])) ? "Touches         " : "Does not touch  ") << positions2[i] << " -> " << printVector(lattice::positionFromID(positions2[i])) << endl;
    }

    cout << "Test neighbor function from structure s2" << endl;
    struct3D* s3 = new struct3D(string("SUULUDLSDS"), UnDefined, lattice::idFromPosisition(31,32,31));
    set<int>* res3 = new set<int>(neighborPositions(*s3));
    cout << print(*res3) << endl;
//    glDisplay(0,c);
//    addToDisplay(s3); // need to be a new XXX to be still existing after the function ends.
//    addToDisplay(res3); // need to be a new XXX to be still existing after the function ends.
//    glutMainLoop();


    cout << "test superProtein" << endl;
    superProtein p = superProtein(string("SUULUDLSDS"));
    cout << print(p) << endl;


    cout << "Test fuse proteins" << endl;
    struct3D* s4 = new struct3D(string("RULSDDL"), UnDefined, lattice::idFromPosisition(21,22,21));
    struct3D* s5 = new struct3D(string("DDLSSURLSSSS"), UnDefined, lattice::idFromPosisition(31,32,31));
    struct3D* s6 = new struct3D(fuse(string("RULSDDL"), string("DDLSSURLSSSS")), UnDefined, lattice::idFromPosisition(41,42,41));
    cout << print(*s4) << endl;
    cout << print(*s5) << endl;
    cout << print(*s6) << endl;

    #ifdef ALLOW_GRAPHICS
    if(true){
        char *c[] = {(char*)"Hello",nullptr};
        glDisplay(0,c);
        addToDisplay(s4); // need to be a new XXX to be still existing after the function ends.
        addToDisplay(s5); // need to be a new XXX to be still existing after the function ends.
        addToDisplay(s6); // need to be a new XXX to be still existing after the function ends.
        glutMainLoop();
    }
    #endif
}


int nbProts;
void testEnsProts() {
    nbProts = 0;

    vector<AA> sequence;
    sequence.push_back(Met);
    sequence.push_back(Cys);
    sequence.push_back(Arg);
    sequence.push_back(Thr);
    sequence.push_back(Gln);
    sequence.push_back(Asp);
    sequence.push_back(His);
    sequence.push_back(Arg);
    sequence.push_back(Trp);
    sequence.push_back(Lys);
    sequence.push_back(Gly);
    sequence.push_back(Gly);
    sequence.push_back(Gly);
    sequence.push_back(Gly);
    sequence.push_back(Gly);
    sequence.push_back(Gly);
    sequence.push_back(Gly);
    sequence.push_back(Gly);
    sequence.push_back(Gly);
    sequence.push_back(Gly);
    sequence.push_back(Gly);
    sequence.push_back(Gly);
    sequence.push_back(Gly);
    sequence.push_back(Gly);
    sequence.push_back(Gly);
    sequence.push_back(Gly);
    sequence.push_back(Gly);
    sequence.push_back(Gly);

    //cout << 5*5*5*5*5*5*5*5 << endl;

    // Numbers : without contraint
    //  1 -> 1
    //  2 -> 6  (backwards)
    //  3 -> 30
    //  4 -> 150
    //  5 -> 726
    //  6 -> 3534
    //  7 -> 16926
    //  8 -> 81390
    //  9 -> 387966
    //  10-> too big
    //  11->
    //  12->
    //  13->
    //  14->
    //  15->
    /*    Staying at x,y,z +/- 3 each direction,
    Size = 1	1
    Size = 2	6
    Size = 3	30
    Size = 4	150
    Size = 5	722
    Size = 6	3450
    Size = 7	15930
    Size = 8	72522
    Size = 9	321794
    Size = 10	1415450
    */
    /*    Staying at x,y,z +/- 2 each direction,
    Size = 1	1
    Size = 2	6
    Size = 3	30
    Size = 4	146
    Size = 5	658
    Size = 6	2858
    Size = 7	11802
    Size = 8	48282
    Size = 9	193714
    Size = 10	781114
    */
    /*    Staying at x,y,z +/- 1 each direction,
    Size = 1	1
    Size = 2	6
    Size = 3	26
    Size = 4	98
    Size = 5	330
    Size = 6	1130
    Size = 7	3746
    Size = 8	12802
    Size = 9	42498
    Size = 10	143610
    Size = 11	472242

    and on top, starting from angle,
    Start for Size = 1	1
    Size = 1	1	nP1	nIDs1
    Start for Size = 2	1
    Size = 2	1	nP2	nIDs3
    Start for Size = 3	1
    Size = 3	1	nP3	nIDs5
    Start for Size = 4	1
    Size = 4	3	nP6	nIDs12
    Start for Size = 5	1
    Size = 5	9	nP15	nIDs28
    Start for Size = 6	1
    Size = 6	31	nP46	nIDs73
    Start for Size = 7	1
    Size = 7	105	nP151	nIDs217
    Start for Size = 8	1
    Size = 8	373	nP524	nIDs718
    Start for Size = 9	1
    Size = 9	1277	nP1801	nIDs2442
    Start for Size = 10	1
    Size = 10	4380	nP6181	nIDs8438
    Start for Size = 11	1
    Size = 11	14518	nP20699	nIDs28468
    Start for Size = 12	1
    Size = 12	48444	nP69143	nIDs95516
    Start for Size = 13	1
    Size = 13	158234	nP227377	nIDs314946
    Start for Size = 14	1
    Size = 14	520420	nP747797	nIDs1037906
    */

    ensProts P = getAllStructures(27, sequence);
    //ensProts P = getStartingStructures(10, sequence);
    cerr << P.size();

    #ifdef ALLOW_GRAPHICS
    if(false){
        char *c[] = {(char*)"Hello",nullptr};
        glDisplay(0,c);
        for(int i = 0; i < min(10000, P.size()); ++i){
            addToDisplay(P[i], false); // need to be a new XXX to be still existing after the function ends.
        }
        glutMainLoop();
    }
    #endif


    return;
}








// reasoning: now, I want
//  -1/ a function to show proteins with AAs ?
//  0/ with an example, check that the affinity is the good one.
//  1/ distribution of folding energies for one sequence, along all structures (make it for 100 sequences)
//  2/ distribution of affinities for one fixed prot (ligand - take 10) to random proteins (take 10000 receptors)
//  3/ can we build manually an optimal structure ?
//  4/ profiles of mutations
//  5/ how do we get the dt ?
//  6/ What are always these 2 or 3 folded optimal ?
//          are they symmetries ?
//          and check that the optimal are not because of errors of double.
// idea for paralellization: will need the affinity far later than the mutation happens -> put a wish list, computes in parallel, and then
// when the cell enters selection it knows the affinity.


superProtein insert(superProtein* existingProt, struct3D* toAdd, int IDfirstResidue){
    return insert(*existingProt, *toAdd, IDfirstResidue);
}

superProtein insert(superProtein* existingProt, string sequence, int startPos, int IDfirstResidue){
    struct3D newChain = struct3D(sequence, UnDefined, startPos);
    return insert(*existingProt, newChain, IDfirstResidue);
}

superProtein insert(superProtein& existingProt, struct3D& toAdd, int IDfirstResidue){
    superProtein P = superProtein(); // make a copy! Would be bad to return a modifidied input...
    if(!existingProt.structure) {cerr << "superProtein insert(), empty structure" << endl; return P;}
    P.points = existingProt.points;

    P.structure->occupiedPositions = existingProt.structure->occupiedPositions;

    if(collide(P, toAdd)) {
        cerr << "ERR: insert(" << print(existingProt) << " with structure  << print(toAdd) <<  and IDresidue=" << IDfirstResidue << ", they collide! => insertion denied" << endl;
        return P;
    }
    P.structure->sequence = existingProt.structure->sequence + string("+") + toAdd.sequence;
    P.structure->separatedSingleMoves = existingProt.structure->separatedSingleMoves + string("+") + toAdd.separatedSingleMoves;
    P.structure->listYAxis = existingProt.structure->listYAxis + string("+") + toAdd.listYAxis;
    P.structure->occupiedPositions = union_sets(P.structure->occupiedPositions, toAdd.occupiedPositions);

    // Now creates the list of residues in space from this structure
    superProtein ProtToAdd = superProtein(toAdd);
    if(existingProt.size() == 0) {
        P.structure->startingPosition = ProtToAdd.structure->startingPosition;
    } else {
        P.structure->startingPosition = existingProt.structure->startingPosition;
    }
    int nR = ProtToAdd.size();          // careful, a struct3D.size() and a protein.size() are different
    for(int i = 0; i < nR; ++i){
        residue toInsert = residue(ProtToAdd.points[i]);
        toInsert.IDresidue = i + IDfirstResidue;
        P.points.push_back(toInsert);
    }
    P.structure->endingPosition = ProtToAdd.structure->endingPosition;
    return P;
}

void testMultichain(){
    struct3D a = struct3D("SULLD");
    superProtein b = superProtein(a);
    cout << print(b) << endl;
    cout << b.getAAseq() << endl;
    //b.setAAs("ABCDE");
    b.setAAs("FGHIJK");
    cout << print(b) << endl;
    cout << b.getAAseq() << endl;

//  << " Hash=" << hashCode(P);
}
