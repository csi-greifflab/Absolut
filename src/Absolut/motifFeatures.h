#ifndef Motif_FEATURES_H
#define Motif_FEATURES_H

#include "../Ymir/ymir.h"
#include "importrepertoire.h"

vector<string> split(const string& str, const string& delim);


double affinityCodeTot(string receptor, string codeInters);
vector<string> recursiveMask(vector<bool> positionsToTryAll, int positionToDo = 0);
string bestEnergyWithReceptorMask(string maskedAAs, string codeInters);
std::pair<string, double> bestSequenceForStructure(string codeInters, int sizeReceptors);
void testBestSequencesPerStructure();
void testFeatures();

// This is the main function that calculates the features.
enum {interCodeWithIDpos, listAAPairs, seqAGEpitope, seqABParatope, motifAGEpitope, motifABParatope, motifsSizeGapsLigand, motifsSizeGapsRec, motifsChemicalLig, motifsChemicalRec, agregatesAGEpitope, agregatesABParatope, chemicalAGEpitope, chemicalABParatope, positionsBound, hotspot_ID, segmentedABParatope, segmentedAGEpitope, interMaskABParatope, interMaskAGEpitope, AAcompoAGEpitope, AAcompoABParatope, AAcompoFullSlice, AAcompoFullCDR, sizeCDR3, interCodeInternal, contribPerAAepi, contribPerAAparaBind, contribPerAAparaFold, NB_features}; //nbneighbors, distChem, selfFolding,
vector<string> structuralFeatures(superProtein& ligand, superProtein & s2, int minDegreeInteract = 1, bool includeDegree = false);


std::pair<string , string> minimalFeaturesDegreeOne(superProtein& ligand, superProtein & s2);

struct binding;

class features
{
    superProtein* currentLigand;
    string antigenID;
    vector< set<int> > hotspotsLarge;

    int minDegree;
    bool includeDegreeInMotifs;

public:
    // tracks list sequences already read for this ligand.
    importRepertoire rep;


    //features(superProtein* p, int _minDegree = 1, string _antigenID);
    features(string _antigenID, int _minDegree = 1, int _includeDegreeInMotifs = false);

    // the file will read, for each AA sequence, the structure and AA acid content of the receptor.
    vector<string> getProperties(string AAligand, int posStart, string structLigand, string fullSequence = "");
    vector<string> getProperties(binding &b, string fullSequence = "");
    static string titleFeatures();
    static string printLine(vector<string> res);

    // file format to be read: AAsequence, Energy, InteractionProfile (for checking), nbStructures, StartPos-Structure1, StartPos-Structure2 ...
    // use display = true only when using openGL libraries, then it shows
    string getPropertiesDeprecated(string fileName, bool display = false, int nToDisplay = 0);
    string getPropertiesFormat2Deprecated(string fileName, bool display = false, int nToDisplay = 0);
    //void poolBindingsFiles(string directory, string filePatterns, bool daughterFolders = false);

    // Every time getProperties (or format2) are called, they remember the first element, as an example, for testing
    string exampleSeq;
    double exampleEnergy;
    string exampleFirstInterCode;
    string exampleFirstStruct;
};

string shortCodeAAcompo(vector<int> numberPerLetter);
string fullCodeAAcompo(vector<int> numberPerLetter, char sep = '\t');
string convertAAtoAAusage(string AAseq, char sep = '\t');
string convertAAtoAAcode(string AAseq);
void testsAAcompo();

void testBestSequencesPerStructure();

#endif // FEATURES_H
