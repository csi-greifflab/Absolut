#include "poolstructs.h"
#include "ymir.h"
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>  // also for string.replace
using namespace std;


bool areStructEqual(string struc1, string struc2, set<int>& pos1, set<int>& pos2){
    if(struc1.size() != struc2.size()) return false;
    if(pos1.size() != pos2.size()) return false;

    // rotates both structures to start by S*U and compare
    string rel1 = normalizeAbsolute(struc1);
    string rel2 = normalizeAbsolute(struc2);
    // also rotates if takes the other structure from the other end
    string rel3 = normalizeAbsolute(revert(struc2)); // from the other side
    // if both of them are different, then different
    if(rel1.compare(rel2) && (rel1.compare(rel3))) return false;

    // If the structures are compatible, check that they are in the same space points (like not rotated)
    // does it by checking they occupy the exact same positions in space.
    // (technically, we could also have just compared 3 positions but would be more complex)
    // trick: if the union of two sets is the same amount of elements, then the two sets were identical. Maybe there are better functions...
    set<int> pooledPos = union_sets(pos1, pos2);
    if(pooledPos.size() != pos1.size()) return false;
    return true;
}

bool areStructEqual(int pos1, string struc1, int pos2, string struc2){
    struct3D s1 = struct3D(struc1, UnDefined, pos1);
    struct3D s2 = struct3D(struc2, UnDefined, pos2);
    return areStructEqual(s1.sequence, s2.sequence, s1.occupiedPositions, s2.occupiedPositions);
}

bool areStructEqual(struct3D& s1, struct3D& s2){
    return areStructEqual(s1.sequence, s2.sequence, s1.occupiedPositions, s2.occupiedPositions);
}



// Idea: there are always two ways to describe the same structure (from both ends)
// so we would like to always use the same ID. Easy, from the two starting positions, take the lowest.

// The getinverted returns that the structure was inverted
std::pair<string, int> uniqueStructure(struct3D& s1){ //, bool& gotInverted){
    if(s1.endingPosition < s1.startingPosition){
        string rev = revert(s1.sequence);

        // test, can be ignored later
        struct3D fromOtherEnd(rev, UnDefined, s1.endingPosition);
        if(!areStructEqual(s1, fromOtherEnd)){
            cerr << "ERR: getUniqueIDStructure, for whatever reason, the reversing of structure failed.\n"
                    "Origin: " << s1.sequence << "-" << s1.startingPosition << " and raised " <<  fromOtherEnd.sequence << "-" << fromOtherEnd.startingPosition << endl;
            return std::pair<string, int>("",0);
        }
        //gotInverted = true;
        return std::pair<string, int>(rev,s1.endingPosition);
    }
    //gotInverted = false;
    return std::pair<string, int>(s1.sequence, s1.startingPosition);
}


// The getinverted returns that the structure was inverted
std::pair<string, int> oppositeEqualStructure(struct3D& s1){ //, bool& gotInverted){
    if(s1.endingPosition >= s1.startingPosition){
        string rev = revert(s1.sequence);

        // test, can be ignored later
        struct3D fromOtherEnd(rev, UnDefined, s1.endingPosition);
        if(!areStructEqual(s1, fromOtherEnd)){
            cerr << "ERR: getUniqueIDStructure, for whatever reason, the reversing of structure failed.\n"
                    "Origin: " << s1.sequence << "-" << s1.startingPosition << " and raised " <<  fromOtherEnd.sequence << "-" << fromOtherEnd.startingPosition << endl;
            return std::pair<string, int>("",0);
        }
        //gotInverted = true;
        return std::pair<string, int>(rev,s1.endingPosition);
    }
    //gotInverted = false;
    return std::pair<string, int>(s1.sequence, s1.startingPosition);
}

// Format: SSULLUD 315865
string getUniqueIDStructure(struct3D& s1){
    std::pair<string, int> res = uniqueStructure(s1);
    stringstream pool;
    pool << res.first << "\t" << res.second;
    return pool.str();
}

// Format: 315865-SSULLUD (no space in written files)
string getWrittenUniqueIDStructure(struct3D& s1){
    std::pair<string, int> res = uniqueStructure(s1);
    stringstream pool;
    pool << res.second << "-" << res.first;
    return pool.str();
}


std::pair<string, int> retrieveStructureFromID(string ID, char sep){
    if(sep != '\t'){
        std::replace( ID.begin(), ID.end(), sep, '\t');
    }
    stringstream read(ID);
    string structure;
    read >> structure;
    int startpos;
    read >> startpos;
    return std::pair<string,int>(structure, startpos);
}
// opposite order
std::pair<int, string> retrieveStructureFromPosAndStructure(string ID, char sep){
    if(sep != '\t'){
        std::replace( ID.begin(), ID.end(), sep, '\t');
    }
    stringstream read(ID);
    int startpos;
    read >> startpos;
    string structure;
    read >> structure;
    return std::pair<int,string>(startpos, structure);
}


bool fastAreStructEqual(struct3D& s1, struct3D& s2){
    return (!getUniqueIDStructure(s1).compare(getUniqueIDStructure(s2)));
}

string testEqual(struct3D& s1, struct3D& s2){
    stringstream res;
    res << "Structures " << s1.startingPosition << "-" << s1.sequence << " and ";
    res << s2.startingPosition << "-" << s2.sequence;

    if(areStructEqual(s1, s2)) res << " are EQUAL";
    else res << " are DIFFERENT";
    if(fastAreStructEqual(s1, s2)) res << " and EQUAL by fast ";
    else res << " and DIFFERENT by fast ";
    return res.str();
}

void testAreStructEqual(){
    // example of structures that are equal:
    // 116835 BRSRUUDUDR
    // 116835 BRSRUUDUDR
    {
        struct3D s1("BRSRUUDUDR", UnDefined, 116835);
        struct3D s2("BRSRUUDUDR", UnDefined,116835 );
        cout << testEqual(s1,s2) << " should be EQUAL" << endl;
    }

    // 116835 BRSRUUDUDR
    // 129185 LDRLRLLDSD
    {
        struct3D s1("BRSRUUDUDR", UnDefined, 116835);
        struct3D s2( "LDRLRLLDSD", UnDefined,129185);
        cout << testEqual(s1,s2) << " should be EQUAL" << endl;
    }
    // exaple of structures that are not equal:
    // 116835 BRSRUUDUDR
    // 111111 LDRLRLLDSD
    {
        struct3D s1("BRSRUUDUDR", UnDefined, 116835);
        struct3D s2("LDRLRLLDSD", UnDefined, 111111);
        cout << testEqual(s1,s2) << " should be DIFFERENT" << endl;
    }

    // 116835 BRSRUUDUDR
    // 116836 BRSRUUDUDR
    {
        struct3D s1("BRSRUUDUDR", UnDefined, 116835);
        struct3D s2("BRSRUUDUDR", UnDefined, 116836);
        cout << testEqual(s1,s2) << " should be DIFFERENT" << endl;
    }

    // 116835 BSRSRUUDUDR
    // 116835 SBRSRUUDUDR
    {
        struct3D s1("BSRSRUUDUDR", UnDefined, 116835);
        struct3D s2("SSRSRUUDUDR", UnDefined, 116835);
        cout << testEqual(s1,s2) << " should be DIFFERENT" << endl;
    }

    cout << "Test IDs" << endl;
    // 116835 BRSRUUDUDR
    // 129185 LDRLRLLDSD
    // 129185 LDRLRLLDSD
    {
        struct3D s1("BRSRUUDUDR", UnDefined, 116835);
        struct3D s2("LDRLRLLDSD", UnDefined, 129185);
        struct3D s3("LDRLRLLDSU", UnDefined, 129186);
        cout << "BRSRUUDUDR-116835 (end" << s1.endingPosition << ") => ID: " << getUniqueIDStructure(s1) << " returns to " << retrieveStructureFromID(getUniqueIDStructure(s1)).first << "-" << retrieveStructureFromID(getUniqueIDStructure(s1)).second << endl;
        cout << "LDRLRLLDSD-129185 (end" << s2.endingPosition << ")=> ID: " << getUniqueIDStructure(s2) << " returns to " << retrieveStructureFromID(getUniqueIDStructure(s2)).first << "-" << retrieveStructureFromID(getUniqueIDStructure(s2)).second << endl;
        cout << "LDRLRLLDSU-129186 (end" << s3.endingPosition << ")=> ID: " << getUniqueIDStructure(s3) << " returns to " << retrieveStructureFromID(getUniqueIDStructure(s3)).first << "-" << retrieveStructureFromID(getUniqueIDStructure(s3)).second << endl;
    }
}

// returns a structure and number of equal structures in the list
std::map<string, int> groupStructuresInClasses(vector<struct3D*> toParse){
    std::map<string, int> res;
    size_t NS = toParse.size();
    for(size_t i = 0; i < NS; ++i){
        struct3D* sNew = toParse[i];
        string ID = getUniqueIDStructure(*sNew);

        // look for this ID in the dictionnary
        std::map<string, int>::iterator it = res.find(ID);
        if(it != res.end()){

            it->second++;
        } else {
            res.insert(std::pair<string, int>(ID, 1));
        }
    }
    return res;
}


// the data comes as: CDR-sequence-structure-affinity
// how do we sort? pre-decide per CDR or per sliding window => pre processing before.
// we want to:
//  Affinity -> position
//  - know where high affinity sequences bind (on the antigen, and where are they in space)
//      can pre-sort high affinity sequences, and make a heatmap of their number => only need a list of structures
//
//  Position -> affinity
//  - for an epitope, know which sequences bind there and what are their affinity
//  Correlations: which group of positions on the antigen makes epitopes?


// For each antigen position, list of
std::map<int, int> distributionStructuresAroundAntigen(vector<struct3D*> toParse){


}


