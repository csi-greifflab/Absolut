//#include <vector>
//#include <string>
//#include <sstream>
//#include <fstream>
//#include <iostream>
//#include <cmath>
//#include <set>
//#include <map>
//#include <cstring> // strcpy
//#include <getopt.h>     // to parse main arguments

#include "Ymir.h"
#include <string>
using namespace std;

int main(void){

    string simpleStructure = string("DUUSURUSLDDUUDDUDDSUUSRDDRUDDUUDDUDDSUULS");
    string AAsimple = "YFHGCARRATLNTTISWEYVSVDMEKIRVGGNEWFNHTMYVT";     //or: randomProt(simpleAccessible.size()+1);
    int receptorSize = 8;
    int minNbContacs = 4;
    double kT = 1.0;

    vector<int> forbOptions = {POSITIONS_SQUARE_BLOCKED, lattice::idFromPosisition(30,22,22)};

    /* Possibility to display:
    glDisplay(argc, argv);
    addToDisplay(ligand, true);
    addToDisplay(forb);
    glutMainLoop();
    */

    affinityOneLigand T1 = affinityOneLigand(simpleStructure, AAsimple, lattice::idFromPosisition(31,30,34), receptorSize, minNbContacs, -1, kT, forbOptions);

    cout << "Details of the structures and affinities for " << simpleStructure << " (" << AAsimple << "), receptors " << receptorSize << " minI=" << minNbContacs << endl;
    for(int i = 0; i < 1; ++i){
        string Px = randomProt(receptorSize+1);
        std::pair<double, double> res = T1.affinity(Px, true);
        cout << "affinity(" << Px << ") = \t" << res.first << "\t" << res.second << endl;
    }

}

