#include "dlab.h"
// 0049 176 6871 9497
dlab::dlab()
{

}

#include "../Ymir/ymir.h"
#include "../Tools/md5.h"
#include "motifFeatures.h"
#include "quality.h"
#include "selfEvo.h"
#include "../Tools/zaprandom.h"
#include "../Tools/dirent.h"
#include "antigenLib.h"
#include "importrepertoire.h"
#include "poolstructs.h"
#include "epitope.h"
#include "html.h"
#include "fileformats.h"
#include "topology.h"
#include <iostream>
#include <fstream>

#include <cmath>

// we will have a 3D tensor here...

// Options:

void testDlab(std::string antigenID, std::string AAsequence);


/* ID_slide_Variant	CDR3	Best	Slide	Energy	Structure
1_00a	CAGPSTTVPYYFDYW	false	CAGPSTTVPYY	-70.09	132966-BDUDDRDURR
1_01a	CAGPSTTVPYYFDYW	false	AGPSTTVPYYF	-77.71	132966-BDUDDRDURR
1_02a	CAGPSTTVPYYFDYW	false	GPSTTVPYYFD	-74.03	132965-DDUULULDDS
1_03a	CAGPSTTVPYYFDYW	true	PSTTVPYYFDY	-78.32	132966-BDSRRLUUSL
1_04a	CAGPSTTVPYYFDYW	false	STTVPYYFDYW	-73.95	132966-BDUDDRDURR
2_00a	CARAYYSNDYW	true	CARAYYSNDYW	-67.54	128929-BLLUUSDLLD
3_00a	CARWDDYDDWFAYW	true	CARWDDYDDWF	-81.87	128996-RLLRLUDUUS
3_01a	CARWDDYDDWFAYW	false	ARWDDYDDWFA	-71.97	124837-RRUURLRRSR
3_02a	CARWDDYDDWFAYW	false	RWDDYDDWFAY	-75.78	132959-ULLDLUURRU
3_03a	CARWDDYDDWFAYW	false	WDDYDDWFAYW	-77.12	124837-RUULRRSRSU
4_00a	CARESSGYGYW	true	CARESSGYGYW	-76.38	137190-LSDDUDLRLL
5_00a	CARYNYGPMDYW	false	CARYNYGPMDY	-79.54	128996-RLLRLUSURD
5_01a	CARYNYGPMDYW	true	ARYNYGPMDYW	-81.4	141287-DDLUSURLRR
7_00a	CARVPNWDVNWSFDVW	false	CARVPNWDVNW	-80.33	128996-RLLRLUDUUS
7_01a	CARVPNWDVNWSFDVW	false	ARVPNWDVNWS	-79.75	157796-SRLUUSRDUL
7_02a	CARVPNWDVNWSFDVW	false	RVPNWDVNWSF	-69.38	132896-LDLLULUDLD
7_03a	CARVPNWDVNWSFDVW	false	VPNWDVNWSFD	-75.46	124837-RRURSUUSLD
7_04a	CARVPNWDVNWSFDVW	false	PNWDVNWSFDV	-87.47	141281-BDUDLDRLUL
7_05a	CARVPNWDVNWSFDVW	true	NWDVNWSFDVW	-91.71	137248-RRDLDUDRDR */

vector <vector <int> > possibleRoundings(vector<double> c);

vector< std::pair<string, string> > latticesComplexes(superProtein* antigen, int startPos, string structure, string AAseq, size_t latticeSize){
    vector< std::pair<string, string> > res;

    //voxelGrid g1 = voxelGrid("RRDLDUDRDR", 137248, "NWDVNWSFDVW").reshapeAroundCenter(20, 20, 20);
    //std::pair<superProtein*, vector<int>> AG = getAntigen("1ADQ_A");
    voxelGrid g1 = voxelGrid(structure, startPos, AAseq).reshapeAroundCenter(20, 20, 20);
    voxelGrid g2 = voxelGrid(antigen).reshapeAroundCenter(20, 20, 20);

    voxelGrid filtered1 = g1.interfaceWith(g2);
    voxelGrid filtered2 = g2.interfaceWith(g1);
    voxelGrid complexOnlyInterfaces = filtered1 + filtered2;

    vector<double> c1 = filtered1.getCentre();
    vector<double> c2 = filtered2.getCentre();
    vector<double> c3 = {0.5*(c1[0] + c2[0]), 0.5*(c1[1] + c2[1]), 0.5*(c1[2] + c2[2])};

    vector<voxelGrid*> centered1 = filtered1.PossibleCenterings(c3, latticeSize, latticeSize, latticeSize);
    vector<voxelGrid*> centered2 = filtered2.PossibleCenterings(c3, latticeSize, latticeSize, latticeSize);

    for(size_t i = 0; i < centered1.size(); ++i){
        res.push_back(std::pair<string,string>(centered1[i]->asText(false), centered2[i]->asText(false)));
        delete centered1[i];
        delete centered2[i];
    }

    return res;
}


void testCentering(){

    //vector <vector <int> > test1 = possibleRoundings({1,2,3});
    //vector <vector <int> > test2 = possibleRoundings({1.5,2.5,3.5});
    //vector <vector <int> > test3 = possibleRoundings({1,2,3.5});


    cout << "Empty constructor" << endl;
    voxelGrid g = voxelGrid(8, 8, 8);
    g.content[3][3][3] = 1;
    cout << g.asText() << "\n\n";



    vector<double> c3 = {3, 3, 3};
    vector<voxelGrid*> centered1 = g.PossibleCenterings(c3, 6, 6, 6);

    for(size_t i = 0; i < centered1.size(); ++i){
        cout << centered1[i]->asText() << "\n\n";
        delete centered1[i];
    }
}

void testCentering2(){
    voxelGrid g = voxelGrid(20, 20, 20);
    g.content[15][15][15] = 1;
    g.content[14][15][15] = 18;
    g.content[16][15][15] = 1;
    g.content[15][14][15] = 1;
    g.content[15][16][15] = 1;

    // old center = where the old grid was centered
    vector<double> c3 = g.getCentre();
    cout << printVector(c3) << endl;

    //
    vector<voxelGrid*> centered1 = g.PossibleCenterings(c3, 6, 6, 6);

    for(size_t i = 0; i < centered1.size(); ++i){
        cout << centered1[i]->asText() << "\n\n";
        delete centered1[i];
    }
}

void testVoxelGrid(){

    cout << "Empty constructor" << endl;
    cout << voxelGrid(8, 8, 8).asText() << endl;

    cout << "From structure: " << endl;
    voxelGrid g1 = voxelGrid("RRDLDUDRDR", 137248, "NWDVNWSFDVW").reshapeAroundCenter(8, 8, 8);
    cout << printVector(g1.getCentre()) << endl;
    cout << g1.asText() << endl;
    cout << "Antigen 1ADQ_A" << endl;
    std::pair<superProtein*, vector<int>> AG = getAntigen("1ADQ_A");
    voxelGrid g2 = voxelGrid(AG.first).reshapeAroundCenter(8, 8, 8);
    cout << g2.asText() << endl;
    cout << printVector(g2.getCentre()) << endl;
    cout << "Creating the complex" << endl;
    voxelGrid complex = g1 + g2;
    cout << complex.asText() << endl;
    cout << printVector(complex.getCentre()) << endl;


    // Note: we are supposed to keep only interacting residues, damned!
    voxelGrid filtered1 = g1.interfaceWith(g2);
    voxelGrid filtered2 = g2.interfaceWith(g1);
    voxelGrid complexOnlyInterfaces = filtered1 + filtered2;

    cout << "Now showing only residues involved in the interface" << endl;
    cout << complexOnlyInterfaces.asText() << endl;

    cout << "Recentering the complex with the average of the two INTERFACES centres" << endl;

    cout << "Getting the centres" << endl;
    vector<double> c1 = filtered1.getCentre();
    vector<double> c2 = filtered2.getCentre();
    vector<double> c3 = {0.5*(c1[0] + c2[0]), 0.5*(c1[1] + c2[1]), 0.5*(c1[2] + c2[2])};
    cout << "C1: " << printVector(c1) << " C2: " << printVector(c2) << " C3; " << printVector(c3) << endl;

    {
        vector<voxelGrid*> centered1 = filtered1.PossibleCenterings(c3);
        vector<voxelGrid*> centered2 = filtered2.PossibleCenterings(c3);

        cout << "Pairs of centered complexes:" << endl;
        for(size_t i = 0; i < centered1.size(); ++i){
            cout << centered1[i]->asText() << "\n" << centered2[i]->asText() << "\n\n";
            delete centered1[i];
            delete centered2[i];
        }
    }

    cout << "Now creating a smaller centered lattice" << endl;
    {
        vector<voxelGrid*> centered1 = filtered1.PossibleCenterings(c3, 6, 6, 6);
        vector<voxelGrid*> centered2 = filtered2.PossibleCenterings(c3, 6, 6, 6);

        cout << "Pairs of centered complexes:" << endl;
        for(size_t i = 0; i < centered1.size(); ++i){
            cout << centered1[i]->asText() << "\n" << centered2[i]->asText() << "\n\n";
            delete centered1[i];
            delete centered2[i];
        }
    }


}

void testDlab(std::string antigenID, std::string AAsequence){
    // this became option11 in Delimain

    // Questions: Data enhancement? or take the best?
    // How many pairs
    // fnat - do we pick a random one?

    // Step : pick binders pairs

    // Create non-binders from them

    // Do affinity calculation and get binding poses

    // get the fnat score of each binding pose (if they are binder)

    // generate the lattice version  (+/- data augmentation)

    // Data augmentation: Pick non-binder pairs (and we KNOW they are non-binders)





    // what happens if we take suboptimal binders actually?




    // show we can do polyreactivity? hehe

}


set<string> split5(string s){
    size_t L = s.size();
    if(L % 5 != 0) cerr << "ERR: split5, Code size not multiple of 5" << s << endl;
    set<string> res;
    for(size_t i = 0; i < L/5; ++i){
        string pattern = s.substr(i*5, 5);
        if(res.find(pattern) != res.end()){
            cerr << "ERR: duplicate pattern inside sring " << s << ", pattern " << pattern << endl;
        }
        res.insert(pattern);
    }
    return res;
}

/* tests:
 * cout << fnat("abc", "h13564") << endl;
    cout << fnat("", "") << endl;
    cout << fnat("i0042i0042", "j0039") << endl;
    cout << fnat("j0039j0041i0042", "j0039j0041i0042g0044") << endl;
    cout << fnat("j0039i0042", "j0039j0041i0042g0044") << endl;
    cout << fnat("j0039j0043i0042", "j0039j0041i0042g0044") << endl;
    cout << fnat("j0039j0041i0042g0044i0044j0045g0046k0046a0048k0050k0070h0071i0074j0081", "j0039j0041i0042g0044i0044j0045k0046a0048k0050k0070h0073i0074j0081") << endl;

    gives
        ERR: split5, Code size not multiple of 5abc
        ERR: split5, Code size not multiple of 5h13564
        0
        ERR: duplicate pattern inside sring i0042i0042, pattern i0042
        0
        0
        3
        2
        2
        12
    */


int fnat(string code1, string code2){
    // splits by 5 characters each
    set<string> c1 = split5(code1);
    set<string> c2 = split5(code2);
    set<string> inter = intersection_sets<string>(c1, c2);
    return inter.size();
}



voxelGrid::voxelGrid(size_t x, size_t y, size_t z) : sizeX(x), sizeY(y), sizeZ(z) {
    // https://slaystudy.com/initialize-3d-vector-in-c/
    content = vector< vector< vector<int> > > (x , vector< vector<int> > (y, vector<int> (z , 0) ) );
}

voxelGrid::voxelGrid(const voxelGrid& toCopy) : sizeX(toCopy.sizeX), sizeY(toCopy.sizeY), sizeZ(toCopy.sizeZ) {
    content = vector< vector< vector<int> > > (sizeX , vector< vector<int> > (sizeY, vector<int> (sizeZ , 0) ) );
    for(size_t i = 0; i < sizeX; ++i){
        for(size_t j = 0; j < sizeY; ++j){
            for(size_t k = 0; k < sizeZ; ++k){
                content[i][j][k] = toCopy.content[i][j][k];
            }
        }
    }
}

voxelGrid::voxelGrid(superProtein* P){
    sizeX = XWidth;
    sizeY = YWidth;
    sizeZ = ZWidth;
    content = vector< vector< vector<int> > > (sizeX , vector< vector<int> > (sizeY, vector<int> (sizeZ , 0) ) );
    for(int i = 0; i < P->size(); ++i){
        //int IDresidue = P[i].IDresidue;
        int IDposition = (*P)[i].IDposition;
        AA thisAA = (*P)[i].TypeResidue;
        vector<int> pos = lattice::positionFromID(IDposition);
        content[pos[0]][pos[1]][pos[2]] = int(thisAA);
        //cout << pos[0] << "," << pos[1] << "," << pos[2] << " (" << int(thisAA) << " ) | ";
    }
}


voxelGrid::voxelGrid(string structure, int position, string AAs){
    sizeX = XWidth;
    sizeY = YWidth;
    sizeZ = ZWidth;
    content = vector< vector< vector<int> > > (sizeX , vector< vector<int> > (sizeY, vector<int> (sizeZ , 0) ) );
    superProtein P(structure, position);
    P.setAAs(AAs);
    for(int i = 0; i < P.size(); ++i){
        //int IDresidue = P[i].IDresidue;
        int IDposition = P[i].IDposition;
        AA thisAA = P[i].TypeResidue;
        vector<int> pos = lattice::positionFromID(IDposition);
        content[pos[0]][pos[1]][pos[2]] = thisAA;
    }
}

voxelGrid voxelGrid::operator + (const voxelGrid &toAdd){

     if(abs(int(sizeX - toAdd.sizeX)) + abs(int(sizeY - toAdd.sizeY)) + abs(int(sizeZ - toAdd.sizeZ)) > 0){
         cerr << "ERR: voxelGrid::operator +, the two lattices have different dimensions" << endl;
         return voxelGrid(sizeX, sizeY, sizeZ);
     }
     voxelGrid res(*this);
     for(size_t i = 0; i < sizeX; ++i){
         for(size_t j = 0; j < sizeY; ++j){
             for(size_t k = 0; k < sizeZ; ++k){
                 if(content[i][j][k] * toAdd.content[i][j][k] != 0){
                     cerr << "ERR: operator voxelGrid '+', combining two overlapping grids at position " << i << "," << j << "," << k << endl;
                 }
                 res.content[i][j][k] = content[i][j][k] + toAdd.content[i][j][k];
             }
         }
     }
     return res;
}

vector<double> voxelGrid::getCentre(){
    double sumX = 0;
    double sumY = 0;
    double sumZ = 0;
    double totPoints = 0;
    for(size_t i = 0; i < sizeX; ++i){
        for(size_t j = 0; j < sizeY; ++j){
            for(size_t k = 0; k < sizeZ; ++k){
                double used = ((content[i][j][k] != 0) ? 1 : 0);
                sumX += static_cast<double>(i) * used;
                sumY += static_cast<double>(j) * used;
                sumZ += static_cast<double>(k) * used;
                totPoints += used;
            }
        }
    }
    return {sumX / totPoints, sumY / totPoints, sumZ / totPoints};
}

double voxelGrid::access(int x, int y, int z){
    if((x < 0) || (y < 0) || (z < 0)){return 0;}
    if((x >= sizeX) || (y >= sizeY) || (z >= sizeZ)){return 0;}
    return content[x][y][z];
}

string filterPrint(int value){
    if(value == 0) return "_";
    stringstream res;
    res << AAname(AA(value));
    return res.str();
}
string voxelGrid::asText(bool dim2){
    stringstream res;
    if(dim2) res << "    Y=0" << string(max(0, (int) sizeZ-3), ' ') << " Y=1" << string(max(0, (int)sizeZ-3), ' ') << " Y=2" << string(max(0, (int)sizeZ-3), ' ') << " ...\n"; //_____ ________ ________ ________ ________ ________ ________
    if(!dim2) res << "\"";
    for(size_t i = 0; i < sizeX; ++i){
        if(dim2) res << "X=" << i << " ";
        for(size_t j = 0; j < sizeY; ++j){
            for(size_t k = 0; k < sizeZ; ++k){
                //if(i + j + k != 0) res << ", ";
                res << filterPrint(content[i][j][k]);
            }
            res << " ";
        }
        if(dim2) res << "\n";
        else res << ",";
    }
    if(!dim2) res << "\"";
    if(dim2) res << printVector(getCentre()) << endl;
    return res.str();
}

vector <vector <int> > possibleRoundings(vector<double> c){
    if(c.size() != 3) {cerr << "ERR: PossibleRoundings, wrong size" << endl; return vector< vector<int> > ();}
    vector<int> cl = {(int) std::floor(c[0]), (int) std::floor(c[1]), (int) std::floor(c[2])};
    vector<int> ch = {(int) std::ceil(c[0]), (int) std::ceil(c[1]), (int) std::ceil(c[2])};

    vector< vector<int> > res = {{cl[0], cl[1], cl[2]}};
    if(cl[0] != ch[0]){
        size_t existing = res.size();
        for(size_t t = 0; t < existing; ++t){
            vector<int> newV = res[t];
            newV[0] = ch[0];
            res.push_back(newV);
        }
    }
    if(cl[1] != ch[1]){
        size_t existing = res.size();
        for(size_t t = 0; t < existing; ++t){
            vector<int> newV = res[t];
            newV[1] = ch[1];
            res.push_back(newV);
        }
    }
    if(cl[2] != ch[2]){
        size_t existing = res.size();
        for(size_t t = 0; t < existing; ++t){
            vector<int> newV = res[t];
            newV[2] = ch[2];
            res.push_back(newV);
        }
    }
    //cout << "--- initial= " << printVector(c) << endl;
    //for(size_t t = 0; t < res.size(); ++t){
    //    cout << printVector(res[t]) << endl;
    //}
    return res;
    // old way: vector< vector<int> > centres = { {cl[0], cl[1], cl[2]}, {cl[0], cl[1], ch[2]}, {cl[0], ch[1], cl[2]}, {cl[0], ch[1], ch[2]}, {ch[0], cl[1], cl[2]}, {ch[0], cl[1], ch[2]}, {ch[0], ch[1], cl[2]}, {ch[0], ch[1], ch[2]}};
}

// according to a desired (floating centre)
// The trick is that the new grid can be smaller!
vector<voxelGrid*> voxelGrid::PossibleCenterings(vector<double> c, size_t newSizeX, size_t newSizeY, size_t newSizeZ){
    if(newSizeX == 0) newSizeX = this->sizeX;
    if(newSizeY == 0) newSizeY = this->sizeY;
    if(newSizeZ == 0) newSizeZ = this->sizeZ;

    vector<double> newCenter = {static_cast<double>(newSizeX-1) / 2., static_cast<double>(newSizeY-1) / 2., static_cast<double>(newSizeZ-1) / 2.};
    vector<double> tr = {newCenter[0] - c[0], newCenter[1] - c[1], newCenter[2] - c[2]};
    vector< vector<int> > possibleTrs = possibleRoundings(tr);

    //cout << printVector(c) << endl;
    //cout << printVector(newCenter) << endl;

    vector<voxelGrid*> res;
    for(size_t s = 0; s < possibleTrs.size(); ++s){
        voxelGrid* g = new voxelGrid(newSizeX, newSizeY, newSizeZ);
        for(int i = 0; i < static_cast<int>(newSizeX); ++i){
            for(int j = 0; j < static_cast<int>(newSizeY); ++j){
                for(int k = 0; k < static_cast<int>(newSizeZ); ++k){
                    g->content[i][j][k] = this->access(i - possibleTrs[s][0], j - possibleTrs[s][1], k - possibleTrs[s][2]);
                }
            }
        }
        res.push_back(g);
    }
    return res;
}


voxelGrid voxelGrid::reshapeAroundCenter(size_t newSizeX, size_t newSizeY, size_t newSizeZ) {
    vector<size_t> oldCentre = {sizeX/2, sizeY/2, sizeZ/2};

    voxelGrid g = voxelGrid(newSizeX, newSizeY, newSizeZ);

    vector<int> newCenter = {static_cast<int>(newSizeX-1) / 2, static_cast<int>(newSizeY-1) / 2, static_cast<int>(newSizeZ-1) / 2};
    int tr_x = newCenter[0] - static_cast<int>(oldCentre[0]); // can be negative so not a size_t type
    int tr_y = newCenter[1] - static_cast<int>(oldCentre[1]);
    int tr_z = newCenter[2] - static_cast<int>(oldCentre[2]);

    for(int i = 0; i < static_cast<int>(newSizeX); ++i){
        for(int j = 0; j < static_cast<int>(newSizeY); ++j){
            for(int k = 0; k < static_cast<int>(newSizeZ); ++k){
                g.content[i][j][k] = this->access(i - tr_x, j - tr_y, k - tr_z);
            }
        }
    }
    return g;
}




voxelGrid voxelGrid::interfaceWith(voxelGrid &other){
    voxelGrid res(*this);
    if(abs(int(sizeX - other.sizeX)) + abs(int(sizeY - other.sizeY)) + abs(int(sizeZ - other.sizeZ)) > 0){
        cerr << "ERR: voxelGrid::interface, the two lattices have different dimensions" << endl;
        return voxelGrid(sizeX, sizeY, sizeZ);
    }
    for(size_t i = 0; i < sizeX; ++i){
        for(size_t j = 0; j < sizeY; ++j){
            for(size_t k = 0; k < sizeZ; ++k){
                int sum = 0;
                sum += abs(other.access(i-1, j, k));
                sum += abs(other.access(i+1, j, k));
                sum += abs(other.access(i, j-1, k));
                sum += abs(other.access(i, j+1, k));
                sum += abs(other.access(i, j, k-1));
                sum += abs(other.access(i, j, k+1));
                if(sum == 0) res.content[i][j][k] = 0;
            }
        }
    }
    return res;
}


