#include "fileformats.h"
#include "antigenLib.h"
#include "poolstructs.h"
#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;

vector<struct3D*> filterStructuresFromAbsolutRawData(string bindingDatasetFile, double affMin, double affMax);
vector<struct3D*> filterStructuresFromRawDataOldFormat(string bindingDatasetFile, double affMin, double affMax, int maxStructures);



// possible answers: "notFound", "repertoire", "rawBindings", "analyzedBindings", "notRecognized"
string getFormat(string fileName){
    {
        ifstream test(fileName);
        if(!test) return string("notFound");
    }
    {
        dataset<analyzedBinding> test(fileName);
        if(test.nLines() > 0) return string("analyzedBindings");
    }
    {
        dataset<binding> test(fileName);
        if(test.nLines() > 0) return string("rawBindings");
    }
    {
        repertoire test(fileName);
        if(test.nLines() > 0) return string("repertoire");
    }
    return string("notRecognized");
}


// 1- A list of CDR3 sequences with an ID is called a 'repertoire'"
// Format requirements:
//      if the first character is a '>' then will recognize as fasta format, the line with '>' is the ID,
//      expect two columns per line, one for the ID, then sequence, then everything else in that line will be discarded
//          (so the file can actually have more sequences, they will just be ignored)
repertoire::repertoire(string fileName, bool skipHeader){
    ifstream f;
    f.open(fileName.c_str());
    if(!f){cerr << "ERR::repertoire::repertoire, file not found: " << fileName << endl; return;}

    //cout << "   ... Opening repertoire file: " << fileName << endl;
    originalFileName = fileName;

    char bufHeader[10001];
    if(skipHeader) {f.getline(bufHeader, 10000);}

    while(f.good()){
        string ID;
        string sequence;
        f >> ID >> sequence;
        if(ID.size() > 0){      // last line / EOF always creates problems
            //listIDs.push_back(ID);
            //sequences.push_back(sequence);
            addSequence(ID, sequence);
        }
        char buf[10001];
        f.getline(buf, 10000);
    }
    f.close();

    //cout << "   ... Got : " << listIDs.size() << " lines " << endl;
}

size_t repertoire::nLines(){
    size_t nIDs = listIDs.size();
    if(sequences.size() != nIDs) cerr << "ERR: repertoire, nLines, it occurs that there are " << nIDs << " IDs while " << sequences.size() << " sequences... this is not good!" << endl;
    return nIDs;
}

// beginLine and endingLine are both included and start at 0
// returns <ID, sequence>
vector< pair<string, string > > repertoire::getLines(int beginLine, int endingLine, set<string> IDstoIgnore){
    size_t nL = nLines();
    vector< pair<string, string > > res;
    if(beginLine > endingLine) return res;
    if(beginLine < 0) cerr << "ERR: repertoire::getLines(" << beginLine << ", " << endingLine << "), negative line value" << endl;
    if(endingLine >= static_cast<int>(nL)) cerr << "ERR: repertoire::getLines(" << beginLine << ", " << endingLine << "), the repertoire only has " << nL << " lines " << endl;
    endingLine = min(endingLine, 1000000000);

    // THERE WAS A BUG HEEEEERE!!! Last line was omitted (changed bountary up there to -1 for ending line
    for(size_t i = static_cast<size_t>(beginLine); i <= static_cast<size_t>(endingLine); ++i){
        // only adds this sequence if it is not in the list to ignore
        if(IDstoIgnore.find(listIDs[i]) == IDstoIgnore.end()){
            res.push_back(std::pair<string, string>(listIDs[i], sequences[i]));
        }
    }
    return res;
}

pair<string, string > repertoire::getLine(int line){
    if((line < 0) || (line >= static_cast<int>(nLines()))){
        cerr << "ERR: repertoire::getLine(" << line << ", out of bounds, " << nLines() << " lines." << endl;
        return pair<string, string >();
    }
    return std::pair<string, string>(listIDs[line], sequences[line]);
}

pair<string, string > repertoire::getRandomLine(){
    int line = random::uniformInteger(0, nLines()-1);
    if((line < 0) || (line >= nLines())){cerr << "ERR:repertoire::getRandomLine(), choucroute in random generator " << endl; }
    return std::pair<string, string>(listIDs[line], sequences[line]);
}

vector< pair<string, string > > repertoire::getRandomLines(int nWantedLines){
    if(nWantedLines >= nLines()){
        return getLines(0, nLines()-1);
    }
    vector<size_t> remainingLinesToPick(nLines(), 0);
    std::iota(remainingLinesToPick.begin(), remainingLinesToPick.end(), 0);
    printVector(remainingLinesToPick);
    random::shuffle(remainingLinesToPick);

    vector< pair<string, string > > res;
    if(remainingLinesToPick.size() < nWantedLines){cerr << "ERR: Shuffled vector lost weight " << endl;}
    for(size_t i = 0; i < nWantedLines; ++i){
        res.push_back(getLine(remainingLinesToPick[i]));
    }
    return res;
}



bool compFirst(std::pair<string, size_t> a , std::pair<string, size_t> b){
    if(a.first.length() >= maxSizeIDs){
        cerr << "ERR: the IDs are too big to be sorted. Maybe change the size in fileFormats::compFirst function" << endl;
    }
    string s1 = string(maxSizeIDs - a.first.length(), ' ') + a.first;
    string s2 = string(maxSizeIDs - b.first.length(), ' ') + b.first;
    return (s1 < s2); // this is the lexicographic order from string
}


vector<size_t> repertoire::sortedPositionsPerID(){

    vector< std::pair< string, size_t> > vecSorted = vector< std::pair< string, size_t> >();

    std::copy(listPositionsIDs.begin(), listPositionsIDs.end(), std::back_inserter<std::vector<std::pair< string, size_t > > > (vecSorted));

    //cout << vecSorted.size() << endl;
    std::sort(vecSorted.begin(), vecSorted.end(), compFirst);
    //cout << vecSorted.size() << endl;

    size_t N = vecSorted.size();
    vector<size_t> res(N, 0);
    for(size_t i = 0; i < N; ++i){
        res[i] = vecSorted[i].second;
    }
    return res;

}


// thorough documentaion of ifstream/ofstream: http://umich.edu/~eecs381/handouts/filestreams.pdf
bool repertoire::write(string fileName){
    ofstream f;
    f.open(fileName);
    if(!f) {cerr << "ERR: repertoire::write, could not find file: " << fileName << endl; return false;}
    size_t nL = nLines();
    for(size_t i = 0; i < nL; ++i){
        f << listIDs[i] << "\t" << sequences[i] << "\n";
    }
    f.close();
}

bool repertoire::containsSequence(string sequence){
    //std::vector<string>::iterator it = find (sequences.begin(), sequences.end(), sequence);
    std::map<string, vector<size_t> >::iterator it = listPositionsSequences.find(sequence);
    return (it != listPositionsSequences.end());
}

bool repertoire::containsID(string ID){
    //std::vector<string>::iterator it = find (listIDs.begin(), listIDs.end(), sequence);
    std::map<string, size_t >::iterator it = listPositionsIDs.find(ID);
    return (it != listPositionsIDs.end());
}

string repertoire::findSeqFromID(string ID){
    std::map<string, size_t >::iterator it = listPositionsIDs.find(ID);
    if(it != listPositionsIDs.end()){
        size_t pos = it->second;
        if(pos >= nLines()) {
            cerr << "ERR: findSeqFromID: ID " << ID << ", the sequence is supposed to be at position " << pos << " which is not compatible with the number of lines " << nLines() << "."
            "This is an internal error. Maybe the list of IDs or positions was changed without changing the maps of where are the sequences / IDs" << endl;
            return string("");
        }
        return sequences[pos];
    }
    // in case of not found.
    return string("");
}

vector<string> repertoire::findIDfromSeq(string sequence){
    vector<string> res;
    std::map<string, vector< size_t > >::iterator it = listPositionsSequences.find(sequence);
    if(it != listPositionsSequences.end()){
        vector<size_t> listPos = it->second;
        for(size_t i = 0; i < listPos.size(); ++i){
            size_t pos = listPos[i];
            if(pos >= nLines()) {
                cerr << "ERR: findIDfromSeq: sequence " << sequence << ", the sequence is supposed to be (in particular) at position " << pos << " which is not compatible with the number of lines " << nLines() << "."
                "This is an internal error. Maybe the list of IDs or positions was changed without changing the maps of where are the sequences / IDs" << endl;
                res.clear();
                return res;
            }
            res.push_back(listIDs[pos]);
        }
    }
    // in case of not found.
    return res;
}

bool repertoire::addSequence(string ID, string sequence){
    if(containsID(ID)){
        cerr << "ERR: repertoire::addSequence(" << sequence << ", " << ID << "), this ID is already taken. The sequence/ID is NOT added" << endl;
        if(sequence.compare(findSeqFromID(ID))){
            cerr << "   ... even worse, the original sequence with the same ID was different!:" << findSeqFromID(ID) << endl;
        }
        return false;
    }

    listIDs.push_back(ID);
    sequences.push_back(sequence);
    listPositionsIDs.insert(std::pair<string,size_t>(ID, listIDs.size() - 1));
    std::map<string, vector<size_t> >::iterator it = listPositionsSequences.find(sequence);
    if(it == listPositionsSequences.end()){
        listPositionsSequences.insert(std::pair<string, vector<size_t> >(sequence, vector<size_t>(1, sequences.size() - 1)));
    } else {
        listPositionsSequences[sequence].push_back(sequences.size() - 1);
    }
    return true;
}

string repertoire::showMappings(){
    stringstream res;
    size_t nL = nLines();
    for(size_t i = 0; i < nL; ++i){
        res << "position " << i << "\tID " << listIDs[i] << "\tSeq " << sequences[i] << "\n";
    }
    for(std::map<string,size_t>::iterator it1 = listPositionsIDs.begin(); it1 != listPositionsIDs.end(); ++it1){
        res << "ID " << it1->first << " is at position " << it1->second << "\n";
    }
    for(std::map<string, vector<size_t> >::iterator it2 = listPositionsSequences.begin(); it2 != listPositionsSequences.end(); ++it2){
        res << "sequence " << it2->first << " is at position " << printVector(it2->second) << "\n";
    }
    res << "\n";
    return res.str();
}




void testRepertoire(){
    ofstream f("TestRep.txt");
    f << "SequenceOne\tHKISGTW" << endl;
    f << "SequenceTwo\tHKISGTW" << endl;
    f << "SequenceTwo\tHKISGTW" << endl;
    f << "SequenceThree\tXKALWNGT" << endl;
    f << "SequenceTthree\tXKALWNGTK" << endl;
    f.close();

    repertoire T = repertoire("TestRep.txt");
    cout << "Lines:" << T.nLines() << endl;
    cout << "ID from sequence HKISGTW, found " << T.findIDfromSeq("HKISGTW").size() << "IDs " << endl;
    cout << "sequence from ID SequenceThree, found " << T.findSeqFromID("SequenceThree") <<  endl;
    T.write("ReconstitutedRep.txt");
}


// do not mix with function slice(seq, delim)
vector<string> slides(string largeSeq, int size){
    vector<string> res;
    size_t tot = largeSeq.size();
    if(size <= 0) return res;
    if(static_cast<size_t>(size) > tot) {
        //cout << "WRN: slide, sequence " << largeSeq << " discarded because too small (slices expected of size " << size << ")" << endl;
        return res;
    }
    for(size_t i = 0; i <= tot - size; ++i){
        res.push_back(largeSeq.substr(i, size));
    }
    return res;
}

void testSlide(){
    cout << printVector(slides("ABCDEFGHI",3));
    cout << printVector(slides("AB",3));
    cout << printVector(slides("ABCDEFGHI",1));
    cout << printVector(slides("ABCDEFGHI",0));
    cout << printVector(slides("ABCDEFGHI",5));
}





string combinedID(string ID, int nrSlide, int nrStructure){
    stringstream res;
    res << ID << "_";
    if(nrSlide < 0) cerr << "ERR: combinedID, negative nrSlide" << endl;
    if(nrSlide < 10) res << "0";
    if(nrSlide > 99) cerr << "ERR: combinedID, too high nrSlide" << nrSlide << ", max 99" <<endl;
    res << nrSlide;
    if(nrStructure < 0) cerr << "ERR: combinedID, negative nrStructure" << endl;
    if(nrStructure > 25) cerr << "ERR: combinedID, only max 26 authorized nrStructures (0..25), here using " << nrStructure << endl;
    res << static_cast<char>('a' + nrStructure);
    return res.str();
}








binding::binding(string _AAseq, double _bindEnergy, string _structureID){ //, string _interactionCode, set<int>& _bindingAGresiduesID){
    AAseq = _AAseq;
    bindEnergy = _bindEnergy;
    structureID = _structureID;
    //interactionCode = _interactionCode;
    //bindingAGresiduesID = _bindingAGresiduesID;
}



string binding::printLine(){
    stringstream res;
    res << AAseq << "\t" << bindEnergy << "\t" << structureID;
             // << "\t" << bindingAGresiduesID.size() << "\bindingAGresiduesID\t" << print(bindingAGresiduesID) << endl;
    return res.str();
}


// reads one line
binding::binding(string lineToRead){

    string _AAseq;
    double _bindEnergy;
    string _structureID;
    // string _interactionCode;
    // set<int> _bindingAGresiduesID;

    stringstream parseRestLine(lineToRead);

    // tries to parse the line. At least one integer should be read for binding residues. If still buffer = -1, it means the line is incomplete (didnt reach this part)
    parseRestLine >> _AAseq >> _bindEnergy >> _structureID;

//    >> _interactionCode;
//    int residueID = -1;
//    do {
//        residueID = -1;
//        parseRestLine >> residueID;
//        if(_bindingAGresiduesID.find(residueID) != _bindingAGresiduesID.end()) cerr << "ERR: binding(reading file line) for slide " << AAseq << ", a bound residue ID is repeated at least twice(" << residueID << ")" << ", from line " << lineToRead << endl;
//        if(residueID >= 0) _bindingAGresiduesID.insert(residueID);
//    }
//    while((!parseRestLine.eof()) && (residueID != -1));
//    if(_bindingAGresiduesID.size() == 0){
//        cerr << "ERR: binding(read from line), incomplete line, missing ID bound residues. Will return empty binding. Line was: " << lineToRead << endl;
//    } else {

    if(_structureID.size() > 0){
        AAseq = _AAseq;
        bindEnergy = _bindEnergy;
        structureID = _structureID;
    } else {
        cerr << "ERR: binding(read from line), incomplete line: " << lineToRead << endl;
    }

//        interactionCode = _interactionCode;
//        bindingAGresiduesID = _bindingAGresiduesID;

}

//template<typename T>
//dataset<T>::dataset(string fileName)
/* : repertoire() {

    ifstream f;
    f.open(fileName.c_str());
    if(!f){cerr << "ERR::rawBindingsRepertoire::rawBindingsRepertoire, file not found: " << fileName << endl; return;}

    cout << "   ... Opening repertoire file: " << fileName << endl;
    originalFileName = fileName;

    f >> nameAntigen;

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
            if(!(readNext.compare("0") && readNext.compare("false"))){
                bestPerCDR3 = true;
            } else if(!(readNext.compare("1") && readNext.compare("true"))){
                bestPerCDR3 = false;
            } else {
                cerr << "ERR: rawBindingsRepertoire, the file should have the following columns: ID\tCDR3Sequence\tBestForThisCDR?('true'/'false')\t" << binding::printHeaders() << endl;
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
}*/





template<typename T>
bool dataset<T>::containsSlidedSequence(string sequence){
    //std::vector<string>::iterator it = find (sequences.begin(), sequences.end(), sequence);
    std::map<string, vector<size_t> >::iterator it = listPositionsSlidedSequences.find(sequence);
    return (it != listPositionsSlidedSequences.end());
}

template<typename T>
vector<string> dataset<T>::findIDfromSlidedSeq(string slidedSequence){
    //std::vector<string>::iterator it = find (sequences.begin(), sequences.end(), sequence);
    vector<string> res;
    std::map<string, vector< size_t > >::iterator it = listPositionsSlidedSequences.find(slidedSequence);
    if(it != listPositionsSlidedSequences.end()){
        vector<size_t> listPos = it->second;
        for(size_t i = 0; i < listPos.size(); ++i){
            size_t pos = listPos[i];
            if(pos >= nLines()) {
                cerr << "ERR: findIDfromSeq: sequence " << slidedSequence << ", the sequence is supposed to be (in particular) at position " << pos << " which is not compatible with the number of lines " << nLines() << "."
                "This is an internal error. Maybe the list of IDs or positions was changed without changing the maps of where are the sequences / IDs" << endl;
                res.clear();
                return res;
            }
            res.push_back(listIDs[pos]);
        }
    }
    // in case of not found.
    return res;
}

template<typename T>
string dataset<T>::findSlidedSeqFromID(string ID){
    std::map<string, size_t >::iterator it = listPositionsIDs.find(ID);
    if(it != listPositionsIDs.end()){
        size_t pos = it->second;
        if(pos >= nLines()) {
            cerr << "ERR: findSeqFromID: ID " << ID << ", the sequence is supposed to be at position " << pos << " which is not compatible with the number of lines " << nLines() << "."
            "This is an internal error. Maybe the list of IDs or positions was changed without changing the maps of where are the sequences / IDs" << endl;
            return string("");
        }
        return storage[pos]->AAseq;
    }
    // in case of not found.
    return string("");
}




void testRawBindingsRepertoire(){
    ofstream f("TestRawBindRep.txt");
    f << "#Antigen 1EGJ\n"
    "ID_slide_Variant\tCDR3\tBest\t" << printHeaders<binding>() << "\n"
    "Seq1000_01a\tCARGGDHYYGSYWYF\tfalse\tCARGGDHYYGS\t-79.89\t145380-UDDSLLRDUU\n"
    "Seq1000_02a\tCARGGDHYYGSYWYF\tfalse\tARGGDHYYGSS\t-73.91\t145380-SDSLLSRDUU\n"
    "Seq1000_03a\tCARGGDHYYGSYWYF\tfalse\tRGGDHYYGSSY\t-77\t145380-SDSLLSRSDD\n"
    "Seq1000_04a\tCARGGDHYYGSYWYF\tfalse\tGGDHYYGSSYW\t-81.37\t145380-SDSLLSRSDD\n"
    "Seq1000_05a\tCARGGDHYYGSYWYF\ttrue\tGDHYYGSSYWY\t-84.4\t145380-SDSLLSRDUU\n"
    "Seq1000_05b\tCARGGDHYYGSYWYF\ttrue\tGDHYYGSSYWY\t-84.4\t145380-SDSLLSRSDD\n";
    f.close();

    dataset<binding> T = dataset<binding>("TestRawBindRep.txt");
    cout << "Lines:" << T.nLines() << endl;
    cout << "ID from sequence CARGGDHYYGSYWYF, found " << T.findIDfromSeq("CARGGDHYYGSYWYF").size() << "IDs " << endl;
    cout << "sequence from ID Seq1000_02a, found " << T.findSeqFromID("Seq1000_02a") <<  endl;

    cout << "ID from slidedSequence GDHYYGSSYWY, found " << T.findIDfromSlidedSeq("GDHYYGSSYWY").size() << "IDs " << endl;
    cout << "slidedSequence from ID Seq1000_04a, found " << T.findSlidedSeqFromID("Seq1000_04a")  << endl;

    T.write("ReconstitutedRawRep.txt");    

    cout << "Internal position mappings for info" << endl;
    cout << T.showMappings() << endl;
}







vector<struct3D*> filterStructuresFromAbsolutRawData(string bindingDatasetFile, double affMin, double affMax){
    dataset<binding> data = dataset<binding>(bindingDatasetFile);
    string AntigenName = data.nameAntigen;
    cout << "   ... loading antigen " << AntigenName << " from the library" << endl;
    std::pair<superProtein*, vector<int> > AG = getAntigen(AntigenName);
    vector<struct3D*> res;
    for(size_t i = 0; i < data.nLines(); ++i){
        bool best = data.bestForThisCDR3[i];
        double E = data.storage[i]->bindEnergy;
        if(best && ( E <= affMax) && (E >= affMin)){
            string structureID = data.storage[i]->structureID;
            std::pair<int,string> str = retrieveStructureFromPosAndStructure(structureID, '-');
            struct3D* S = new struct3D(str.second, UnDefined, str.first);
            //cout << data.sequences[i] << "->" << data.storage[i]->AAseq << ",E=" << E << "\t" << str.first << " " << str.second << endl;
            res.push_back(S);
        }
    }
    return res;
}



vector<struct3D*> filterStructuresFromRawDataOldFormat(string bindingDatasetFile, double affMin, double affMax, int maxStructures){ //, bool keepMultiBind = false){
    ifstream rdfile(bindingDatasetFile);
    if(!rdfile){
        cerr << "ERR: filterStructuresFromRawDataOldFormat(), couldn't open " << bindingDatasetFile << ", sorry" << endl;
        return vector<struct3D*>();
    }
    vector<double> energiesSoFar;
    int ID;
    string largeStr;
    string slide;
    double E;
    int posStart;
    string structure;
    vector<struct3D*> res;
    int cpt = 0;
    while(rdfile.good()){
        ID = -1;
        rdfile >> ID >> largeStr >> slide >> E;
        //cout << ID << " " << E << endl;
        energiesSoFar.push_back(E);
        if(ID == -1) {break;} //rdfile.close(); return res;} // end of file
        if(( E <= affMax) && (E >= affMin)){
            cout << "- took:E=" << E << "- ";
            cpt++;
            string trashCodes;
            rdfile >> trashCodes >> trashCodes;
            //if(nS >= 1){
                rdfile >> posStart >> structure;
                cout << largeStr << "->" << slide << ",E=" << E << "\t" << posStart << " " << structure << endl;
                struct3D* toAdd = new struct3D(structure, UnDefined, posStart);
                res.push_back(toAdd);
            //}
            if(cpt >= maxStructures) {break;}//rdfile.close(); return res;}
        }
        char trash[10001];
        rdfile.getline(trash, 10000);
    }
    cout << "Found total of " << cpt << "sequences" << endl;
    histogramFromDistrib view = histogramFromDistrib(energiesSoFar, 100);
    cout << view.print() << endl;
    rdfile.close();
    return res;
}
