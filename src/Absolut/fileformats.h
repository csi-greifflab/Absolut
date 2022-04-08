#ifndef FILEFORMATS_H
#define FILEFORMATS_H

/// @brief This file defines the main files formats handling sequences. Briefly,
/// (file format)
/// - a repertoire is a 2-column text table without header, one line is ID_sequence Sequence  (like a dictionary, IDs should be unique)
///        I/O: repertoire a("myFile.txt");  a.write("newFile.txt")
///
/// (line format)
/// - a binding is the raw structural annotation of how one sequence binds an antigen. It includes the 11-AA Slice, energy and structure
/// - an analyzedBinding is the structural annotation + binding features on how one sequence binds an antigen.
///        I/O: (analyzed)binding a(string lineToRead); string theLine = a.printLine();
///
/// (file formats)
/// - a dataset<binding> is like a repertoire, but one line is ID_sequence, sequence, isThisBindingTheBest?, and the infos of the binding
/// - a dataset<analyzedBinding> is like a repertoire, but one line is ID_sequence, sequence, isThisBindingTheBest?, and the infos of the analyzedBinding
///        I/O:  I/O: repertoire a("myFile.txt");  a.write("newFile.txt", (option) sorted=true);
///
/// note: a dataset file has headers, can be accessed with
///     printHeaders<repertoire>()
///     printHeaders<binding>()
///     printHeaders<analyzedBinding>()
///
/// note: since dataset<> is a template, all its code is inside the .h file, sorry for less visibility

#include <vector>
#include <string>
#include <map>
#include <set>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <type_traits>      // to compare types for printHeaders

#include "ymir.h"           // for printVector
#include "motifFeatures.h"
using namespace std;

/// A - Tool functions

/// @brief Tool function that reads a dataset<binding> and extracts only structures between two affinity thresholds.
/// takes only the best binding for each CDR3.
vector<struct3D*> filterStructuresFromAbsolutRawData(string bindingDatasetFile, double affMin, double affMax);

/// @brief Tool function Parses a file and gives its type. possible answers: "notFound", "repertoire", "rawBindings", "analyzedBindings", "notRecognized"
string getFormat(string fileName);

/// B - Classes to read or write files and handle file formats (all text files, tab separated)

/// 1 - ‘repertoire’: A list of CDR3 sequences with an ID and additional information as 'tag' (no headers)
/// Ex:
///      CDR_3073    CARRDDYDGFDYW
///      CDR_3078    CARRGGSSHWYFDVW
///  ...

// Note: for sorting repertoires, a max size for the IDs in number of characters is decided. See compFirst function for details
#define maxSizeIDs 250

struct repertoire{
    repertoire(){}
    repertoire(string fileName, bool skipHeader = false);
    string originalFileName;
    vector<string> listIDs;
    vector<string> sequences;
    virtual bool write(string fileName);

    // tool functions
    virtual size_t nLines();
    vector< pair<string, string > > getLines(int beginLine, int endingLine, set<string> IDstoIgnore = set<string>());
    pair<string, string > getLine(int line);
    vector< pair<string, string > > getRandomLines(int nWantedLines);
    pair<string, string > getRandomLine();

    bool containsSequence(string sequence);
    bool containsID(string ID);
    string findSeqFromID(string ID);
    vector<string> findIDfromSeq(string sequence);

protected:
    // Internal functions to make sure an ID is not used 2 times
    // Maps where are each sequence and ID seen so far (to make search much faster)
    map<string, size_t > listPositionsIDs;
    map<string, vector<size_t> > listPositionsSequences;
public:
    bool addSequence(string ID, string sequence);

public:
    virtual string showMappings();
    virtual ~repertoire(){}
    vector<size_t> sortedPositionsPerID();
};

/// @brief test function / example with a repertoire
void testRepertoire();

/// Function that the returns the header of a file type, defined at the end
template<typename T>
string printHeaders();




/// 2- Binding dataset: One CDR3 is cut into ‘slides’, and one line is ID of this slide, CDR3, and a 'binding' containing info on this slide
/// Two cases:     raw binding dataset: Just minimal information per line (class binding)
/// Ex:
///        ID_slide_Variant        CDR3                Best?        Slide            Energy    Structure
///        3073_00a        CARRDDYDGFDYW        true        CARRDDYDGFD        -67.81        145697-DDRSUDUURD
///        3073_01a        CARRDDYDGFDYW        false        ARRDDYDGFDY        -66.83        141666-BLRLSSDDUL
///
///        analyzed binding dataset: All features of binding (class analyzedBinding)

/// 2a-  Slides: Cutting a CDR3 into slides
/// do not confuse with function slice(seq, delim) in Ymir/fastaffinity that cuts a string by a delimiter
vector<string> slides(string largeSeq, int size);

/// 2b-  ID encoding: Each line will get a different ID, that will be coded as:
///      IDCDR3_SlideNr[xx]VariantStructure[a-z]
///      ex: CDR3Number4_01a
string combinedID(string ID, int nrSlide, int nrStructure);

/// 2c-  One RAW binding (information how one slide binds to the antigen):
///     Slide         Energy     StructureOfBinding (position-structure)
///     CARRDDYDGFD    -67.81    145697-DDRSUDUURD
struct binding {
    binding(string lineToRead);
    binding(string _AAseq, double _bindEnergy, string _structureID); //,  string _interactionCode, set<int>& _bindingAGresiduesID);

    // binding(const binding & ){cerr << "ERR: Are you sure the copy constructor of Binding should be called? " << endl;}

    string AAseq;
    double bindEnergy;
    string structureID;

    virtual string printLine();

    virtual ~binding(){}
};

/// 2d-  One Analyzed binding (information how one slide binds to the antigen):
struct analyzedBinding : public binding {

    analyzedBinding(string lineToRead) : binding (lineToRead) {
        // Note: the three first elements are already read by binding() constructor...
        string _AAseq;
        double _bindEnergy;
        string _structureID;
        string _interactionCode;
        string _bindingAGresiduesID;
        string _hotspotID;

        stringstream parseRestLine(lineToRead);

        // tries to parse the line. At least one integer should be read for binding residues. If still buffer = -1, it means the line is incomplete (didnt reach this part)
        parseRestLine >> _AAseq >> _bindEnergy >> _structureID >> _interactionCode >> _bindingAGresiduesID >> _hotspotID;

        if(_bindingAGresiduesID.size() > 0){
            AAseq = _AAseq;
            bindEnergy = _bindEnergy;
            structureID = _structureID;
            interactionCode = _interactionCode;
            bindingAGresiduesID = stringToSet(_bindingAGresiduesID, '-');
            hotspotID = _hotspotID;
            char buf[10001];
            parseRestLine.getline(buf, 10000);
            string _remainingFeatures = string(buf);
            remainingFeatures =  split(_remainingFeatures, "-");
        } else {
            cerr << "ERR: analuzedBinding(read from line), incomplete line: " << lineToRead << endl;
        }
    }
    analyzedBinding(string _AAseq, double _bindEnergy, string _structureID,  string _interactionCode, set<int> _bindingAGresiduesID, string _hotspotID, vector<string> _remainingFeatures)
        : binding(_AAseq, _bindEnergy, _structureID) {
        interactionCode = _interactionCode;
        bindingAGresiduesID = _bindingAGresiduesID;
        hotspotID = _hotspotID;
        remainingFeatures = _remainingFeatures;
    }

    string interactionCode;
    set<int> bindingAGresiduesID;
    string hotspotID;
    vector<string> remainingFeatures;


    //set<int> getOccupiedPositions();
    //set<int> getBindingAGPositionsInSpace();

    string printLine() {
        stringstream res;
        res << binding::printLine() << "\t" << interactionCode << "\t" << setToString(bindingAGresiduesID, '-') << "\t" << hotspotID << "\t" << features::printLine(remainingFeatures);
        return res.str();
    }

    ~analyzedBinding(){}
};

template<typename T>
struct dataset : public repertoire {

    string headers(){
        return string("ID_slide_Variant\tCDR3\tBest\t") + printHeaders<T>();
    }
    dataset(string fileName) : repertoire() {

        ifstream f;
        f.open(fileName.c_str());
        if(!f){cerr << "ERR::rawBindingsRepertoire::rawBindingsRepertoire, file not found: " << fileName << endl; return;}

        cout << "   ... Opening repertoire file: " << fileName << endl;
        originalFileName = fileName;



        // Then should be the headers:
        // ID_slide_Variant\tCDR3\tBest\t + header of the type T, for instance, binding is: Slide\tBest\tEnergy\tStructure
        char bufHeader[10001];
        f.getline(bufHeader, 10000); // finishes current line
        //cout << "Line 1" << bufHeader << endl;

        // should start by '#Antigen ' then the antigen
        stringstream line1(bufHeader);
        string start;
        line1 >> start >> nameAntigen;
        if(start.compare("#Antigen")) {
            cerr << "ERR: reading dataset " << fileName << ", this file should start by #Antigen nameAntigen" << endl;
            return;
        }
        if(nameAntigen.size() == 0){
            cerr << "ERR: reading dataset " << fileName << ", the file starts witn #Antigen but without antigen name" << endl;
            return;
        }

        f.getline(bufHeader, 10000); // Now gets the next line
        //cout << "Line 2" << bufHeader << endl;
        // Philippe: I have inconsistencies with older formats with 'Best?' and new format 'Best'... so not checking now
        //if(string(bufHeader).substr(0, headers().size()).compare(headers())){
        //    cerr << "WRN: reading dataset " << fileName << ", the second line should be exactly:\n" << headers() << "\n while was " << bufHeader << endl;
        //    return;
        //}

        int nrLinesERR = 0;
        // this part is similar as the repertoire read function
        while(f.good()){
            bool badLine = false;
            string ID;
            string sequence;
            f >> ID >> sequence;
            if(ID.size() > 0){      // last line / EOF always creates problems
                //listIDs.push_back(ID);
                //sequences.push_back(sequence);
                addSequence(ID, sequence);

                bool bestPerCDR3 = false;
                string readNext; // to check this is '0' or '1', or 'true' or 'false', other choices means problem of parsing (like missing 3rd column)
                f >> readNext;
                if(!(readNext.compare("0") && readNext.compare("false") && readNext.compare("FALSE"))){
                    bestPerCDR3 = false;
                } else if(!(readNext.compare("1") && readNext.compare("true") && readNext.compare("TRUE"))){
                    bestPerCDR3 = true;
                } else {
                    cerr << "ERR: rawBindingsRepertoire, the file should have the following columns: ID\tCDR3Sequence\tBestForThisCDR?('true'/'false')\t" << printHeaders<binding>() << endl;
                    badLine = true;
                }

                // I decide to parse line by line, in case there are arguments, to not start reading the next line,
                // and the binding class takes care of it
                char buf[10001];
                f.getline(buf, 10000);
                T* b = new T(string(buf));
                if(b->AAseq.size() == 0) badLine = true;
                bestForThisCDR3.push_back(bestPerCDR3);
                storage.push_back(b);

                // updates the mapping where each Slided Sequence is (multiple times is ok)
                // this assumes that there was no ID conflict, because the mother function addSequence would have detected it anyways.

                string slidedSeq = b->AAseq;
                std::map<string, vector<size_t> >::iterator it = listPositionsSlidedSequences.find(slidedSeq);
                if(it == listPositionsSlidedSequences.end()){
                    listPositionsSlidedSequences.insert(std::pair<string, vector<size_t> >(slidedSeq, vector<size_t>(1, sequences.size() - 1)));
                } else {
                    listPositionsSlidedSequences[slidedSeq].push_back(sequences.size() - 1);
                }

                // just compting how many lines had problem, for info
                if(badLine == true) nrLinesERR++;
            }
        }
        f.close();

        cout << "   ... Got : " << listIDs.size() << " lines, among which " << nrLinesERR << " had a parsing problem " << endl;
    }
    dataset(){}

    bool write(string fileName, bool sorted = true){

        // creates a vector<string> with all IDs, then will print them in this order.
        vector<size_t> orderPositionsToPrint;

        if(sorted){
            orderPositionsToPrint = sortedPositionsPerID();
            if(orderPositionsToPrint.size() != nLines()) cerr << "ERR (Internal) dataset<binding>:write, Wrong size of sortedPositionsPerID" << endl;
        } else {
            orderPositionsToPrint.resize(nLines(), 0);
            std::iota(orderPositionsToPrint.begin(),orderPositionsToPrint.end(),0); // 0,1,2,3,4,5,6
        }



    // thorough documentaion of ifstream/ofstream: http://umich.edu/~eecs381/handouts/filestreams.pdf
//    template<typename T>
//    bool dataset<T>::write(string fileName){
        ofstream f;
        f.open(fileName);
        if(!f) {cerr << "ERR: repertoire::write, could not find file: " << fileName << endl; return false;}
        if(nameAntigen.size() == 0) {
            nameAntigen = "UNKNOWN?";
            cerr << "WRN: Writing a dataset<> with no known nameAntigen => using UNKNOWN? instead";
        }
        f << "#Antigen " << nameAntigen << "\n";
        f << headers() << "\n";
        size_t nL = nLines();
        for(size_t j = 0; j < nL; ++j){
            size_t i = orderPositionsToPrint[j];
            if(storage[i] == nullptr) cerr << "ERR: empty binding inside rawBindingsRepertoire, line " << i << " starting at 0" << endl;
            f << listIDs[i] << "\t" << sequences[i] << "\t" << ((bestForThisCDR3[i]) ? "true" : "false") << "\t" << storage[i]->printLine() << "\n";
        }
        f.close();
        return true;
    }


    string nameAntigen;

    void setNameAntigen(string _nameAntigen){nameAntigen = _nameAntigen;}
    bool addLine(string ID, string longSequence, T* b, bool bestBindingForThisSequence = false){
//    template<typename T>
//    bool dataset<T>::addLine(string ID, string longSequence, T* b, bool bestBindingForThisSequence){

        if(b == nullptr) {cerr << "ERR: adding empty binding to rawBindingsRepertoire" << endl;
            return false;}
        if(b->AAseq.size() == 0) {cerr << "ERR: adding empty binding to rawBindingsRepertoire" << endl;
            return false;}

        // if addSequence says false, it means the ID existed already => did not add the line, so we don't add anything
        if(!addSequence(ID, longSequence)) return false;

        bestForThisCDR3.push_back(bestBindingForThisSequence);
        storage.push_back(b);

        // updates the mapping where each Slided Sequence is (multiple times is ok)
        // this assumes that there was no ID conflict, because the mother function addSequence would have detected it anyways.
        string slidedSeq = b->AAseq;
        std::map<string, vector<size_t> >::iterator it = listPositionsSlidedSequences.find(slidedSeq);
        if(it == listPositionsSlidedSequences.end()){
            listPositionsSlidedSequences.insert(std::pair<string, vector<size_t> >(slidedSeq, vector<size_t>(1, sequences.size() - 1)));
        } else {
            listPositionsSlidedSequences[slidedSeq].push_back(sequences.size() - 1);
        }
        return true;
    }


    // inherited fields from the mother class repertoire:
    //    string originalFileName;
    //    vector<string> listIDs;
    //    vector<string> sequences; // this will be the CDR3 sequence + the ID
    // Note: the binding includes the slide sub sequence of the CDR3

    vector<bool> bestForThisCDR3;       // just to make it easy to extract the best binding per CDR3
    vector<T*> storage;

    T* getElement(size_t i){
        if((i >= 0) && (i < nLines())) return storage[i];
        cerr << "ERR: dataset<>, getLine(" << i << " out of bounds " << endl;
        return nullptr;
    }
    size_t nLines(){
        size_t nIDs = listIDs.size();
        if(sequences.size() != nIDs) cerr << "ERR: rawBindingsRepertoire, nLines, it occurs that there are " << nIDs << " IDs while " << sequences.size() << " sequences... this is not good!" << endl;
        if(bestForThisCDR3.size() != nIDs) cerr << "ERR: rawBindingsRepertoire, nLines, it occurs that there are " << nIDs << " IDs while " << bestForThisCDR3.size() << " bestForThisCDR3... this is not good!" << endl;
        if(storage.size() != nIDs) cerr << "ERR: rawBindingsRepertoire, nLines, it occurs that there are " << nIDs << " IDs while " << storage.size() << " bindings... this is not good!" << endl;
        return nIDs;
    }


    // the functions containsSequence, findSeq from ID, and IDfromSeq all work with the CDR3. Here you can also look for slided sequences,
    bool containsSlidedSequence(string sequence);
    vector<string> findIDfromSlidedSeq(string slidedSequence);
    string findSlidedSeqFromID(string ID);

    string showMappings(){
        stringstream res;
        res << repertoire::showMappings();

        for(std::map<string, vector<size_t> >::iterator it2 = listPositionsSlidedSequences.begin(); it2 != listPositionsSlidedSequences.end(); ++it2){
            res << "slidedSequence " << it2->first << " is at position " << printVector(it2->second) << "\n";
        }
        res << "\n";
        return res.str();
    }
protected:
     map<string, vector<size_t> > listPositionsSlidedSequences;
public:
    ~dataset(){
         listPositionsSlidedSequences.clear();
        size_t nD = storage.size();
        //cerr << "Delete content" << endl;
        for(size_t i = 0; i < nD; ++i){
            if(storage[i]){
                delete storage[i];
                storage[i] = nullptr;
            }
        }
        //cerr << "Finished deleting content" << endl;
     }

};

void testRawBindingsRepertoire();

struct analyzedBinding;

template<typename T>
string printHeaders(){
    if (std::is_same<T, binding>::value){
        return string("Slide\tEnergy\tStructure");
    } else if (std::is_same<T, analyzedBinding>::value){
        return string("Slide\tEnergy\tStructure\tinteractionCode\tAGboundPositions\thotspot_ID\t") + features::titleFeatures();
    } else {
        return string("");
    }
}

void testSlide();


#endif // FILEFORMATS_H
