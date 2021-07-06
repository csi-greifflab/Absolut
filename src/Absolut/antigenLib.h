#ifndef ANTIGENLIB_H
#define ANTIGENLIB_H

#include <string>
#include <vector>
#include "../Ymir/ymir.h"
#include "../Tools/distribution.h"
#include "../Tools/zaprandom.h"

void showListIDs();
std::vector<std::string> listIDs();

// Antigens are defined as 1/ a protein (i.e. multiple structures + position of each AA + AA sequence)
//                     and 2/ a list of forbidden positions (possibly glycans);
std::pair<superProtein*, vector<int> >  getAntigen(std::string ID, bool display = false);

// I initially programmed an antigen to be a pair <protein, vector> like the getAntigen() funcion returns.
// Anthough only these two properties are sufficient for simulations, with time it became important to get 
// more info on an antigen
// => I created this structure as a future replacement for the type "std::pair<superProtein*, vector<int> >"
// => In order for the code to be compatible, I named first and second the two equivalent fields, 
struct antigenInfo {
    antigenInfo();
    string ID;
    superProtein* first;
    vector<int> second;    
    vector<int> glycans; // position of glycans
    
    // note: the following fields were calculated from the 7M murine CDR3
    vector< vector<int> > hotspotsCore; // list of positions 100% shared in a hotspot
    vector< vector<int> > hotspotsLarge; // list of positions 100% shared in a hotspot
    vector<double> thresholds; // list of affinity thresholds from the murine dataset.
    string antibodyChains;

    string print();
};

antigenInfo getAntigenInfos(std::string ID);

affinityOneLigand* getAffinityAntigen(std::string ID, int sizeReceptors, int minInteract);

string IDshortcut(int ID);

string IDshortcut(string justPDB);

#endif // ANTIGENLIB_H
