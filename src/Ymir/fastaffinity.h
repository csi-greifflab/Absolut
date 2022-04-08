#ifndef FASTAFFINITY_H
#define FASTAFFINITY_H

#include <map>
#include <set>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
using namespace std;

// to include the full library
#include "receptorligand.h"

#include "pthread.h"
extern pthread_mutex_t lockAccessPrecompAffinities;

// All info of a ligand pose, if one wants to output the top binding poses
struct pose {
    pose(int _startPos, double _affInteract, double _totalAff, string _structure, string _interactionCode) :
        startPos(_startPos), affInteract(_affInteract), totalAff(_totalAff), structure(_structure), interactionCode(_interactionCode) {}
    int startPos;
    double affInteract;
    double totalAff;
    string structure;
    string interactionCode;
    string print();
};

class affinityOneLigand
{
public:
    //#define showBestStructures

    affinityOneLigand(string _ligandStructureSeq, string _ligandAAseq, int startPosition, int _sizeReceptors, int _minimalNInteract, int _minSelfFoldings, double _KT, vector<int> listForbiddenPositions = vector<int>());
    affinityOneLigand(superProtein* _ligand, int _sizeReceptors, int _minimalNInteract, int _minSelfFoldings, double _KT, vector<int> listForbiddenPositions = vector<int>());

    // 1 - list of options requested
    string ligandSeq;
    string ligandAAseq;
    superProtein* ligand;
    int sizeReceptors;
    int minimalNInteract;
    int nInterCodes;
    int minNrSelfInteractions;
    int nFoldingCodes;
    double KT;
    string fileStructures;
    string fileSelfFoldings;

    bool modeUltraFast; // will only return best affinity, not the statistical
    void setUltraFast(bool v){modeUltraFast = v;}

    // 2 - Pre-calculated list of structures (actually not the structures themselves, but their interaction profiles)
    // and stored / loaded into files
    vector<string>* interactions;
    vector<int> nbRepeats;
    vector<string>* selfInteractions;
    vector<int> nbSelfRepeats;

    // 3 - to be called for getting an affinity
    // returns both best and statistical energies.
    std::pair<double, double> affinity(string receptorAASeq, bool showStructures = false, vector<string> *returnBestStructures = nullptr, size_t nPoses = 0, vector<pose> *returnTopPoses = nullptr);

    // just simpler way, just returns the bestEnergy
    double bestAffinity(string receptorAASeq, bool showStructures = false, vector<string> *returnBestStructures = nullptr){
        return affinity(receptorAASeq, showStructures, returnBestStructures).first;}
    // Dictionary of previously calculated affinities
    //map<int, double> affSingleInteractions; // for a receptor, the dictionnary of possible single interactions (position in receptor - AA in ligand)
                                            // will be regenerated for each new receptorAASeq, and will be applied to all interaction profiles

    // 4 - memory of the previously called affinities.
    map<string, double> knownBestAffinities;
    map<string, double> knownStatisticalAffinities;
    map<string, string> knownBestInteractions;

    void printInfos();

    std::map<string, vector<std::pair<int,string> > > profileToStructure;

    affinityOneLigand(const affinityOneLigand& copy){cerr << "NO copy of affinityOneLigand accepted" << endl;}
    // common function called by both constructors. I tried by one constructor to tall another, but it creates a duplicate struct instead and returns an uninstantiated one.
    void initialize(superProtein* _ligand, int _sizeReceptors, int _minimalNInteract, int _minSelfFoldings, double _KT, vector<int> listForbiddenPositions = vector<int>());

    // name of files corresponding to the parameters given in the constructor.
    string fStruct; string fAll; string fCompact;
};




struct receptorLibrary {
    receptorLibrary(int LseqBCRs,int minContacts,string agStruct,string agSeq,double threshold,vector<int> forbiddenPos);
    int LseqBCRs; // in AAs
    int minContacts;
    string agStruct;
    string agSeq;
    vector<int> forbiddenPos;
    double threshold;
    vector<std::pair<double, string>> content; // list of energy - BCRs
    string getRandBCR(double minEnergyRequested);
    void generateLibFile();
    void readLibFile();
    bool libFileExists();
};


set<int>* generateForbidden(vector<int> listForbiddenPositions);
void testFastAffinity();

// Question : does the interaction changes by cutting the sequences with less than 5 affinities ...

#endif // FASTAFFINITY_H
