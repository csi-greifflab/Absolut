#include "importrepertoire.h"
#include "antigenLib.h"
#include "../Ymir/ymir.h"
#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;

//// do not mix with function slice(seq, delim)
//vector<string> slides(string largeSeq, int size){
//    vector<string> res;
//    size_t tot = largeSeq.size();
//    if(size <= 0) return res;
//    if(static_cast<size_t>(size) > tot) {
//        cout << "WRN: slide, sequence " << largeSeq << " discarded because too small (slices expected of size " << size << ")" << endl;
//        return res;
//    }
//    for(size_t i = 0; i <= tot - size; ++i){
//        res.push_back(largeSeq.substr(i, size));
//    }
//    return res;
//}

//void testSlide(){
//    cout << printVector(slides("ABCDEFGHI",3));
//    cout << printVector(slides("AB",3));
//    cout << printVector(slides("ABCDEFGHI",1));
//    cout << printVector(slides("ABCDEFGHI",0));
//    cout << printVector(slides("ABCDEFGHI",5));
//}


//void analyzeRepertoire(string AntigenName, int receptorSize, int minInteract, int nParallel, bool restartFromScratch, int IDprocess){
//    std::pair<superProtein*, vector<int> > AG = getAntigen(AntigenName);
//    affinityOneLigand T3 = affinityOneLigand(AG.first, receptorSize, minInteract, -1, 1, AG.second);
//    T3.setUltraFast(true);

//    string progressFile = "OngoingRepertoireFor" + AntigenName + ".txt";
//    string resultFile = "OngoingResultsFor" + AntigenName + ".txt";

//    // Test if the 'progressed repertoire' exists
//    ifstream f(progressFile.c_str());
//    if(restartFromScratch || (!f)){ // if not, just starts it with treated=false for each sequence
//        importRepertoire rep = importRepertoire("C:/Users/pprobert/Desktop/Main/B-CurrentZapotect/UniqueCDR3s.txt", false);
//        rep.writeUpdatedRepertoire(progressFile);
//        cout << "   ... Generated a local copy of the repertoire with tagging of already computed sequences" << endl;
//    }
//    if(f) f.close();

//    if(restartFromScratch) {
//        if(remove(resultFile.c_str()) != 0){
//            // erll, couldn't remove
//        }
////        ofstream erasing(resultFile, ios::out);
////        erasing << "-- FreshStart --";
////        erasing.close();
//    }
//    // Then gets sequences...
//    //cout << "   ... Now opening repertoire ";
//    importRepertoire rep2 = importRepertoire(progressFile, true);
//    int n = rep2.nLines;

//    // will slice it into slices of size nParallel
//    int nBlock = (n / nParallel) + 1; // integer division
//    rep2.currentIndex = nBlock * IDprocess;
//    nBlock = min(nBlock, 1000000); // Mas 1e6 sequences per run.
//    rep2.writeUpdatedRepertoire("", nBlock); // this will block the sequences in the repertoire file eith a '2' even if not treated yet
//    cout << "   ... Have read " << n << " sequences. will process a block of " << nBlock << " sequences " << endl;

//    cout << "   ... Will now save results in " << resultFile << endl;
//    ofstream fout(resultFile, ios::app);
//    if(!fout) {cerr << "ERR: Couldn't open output file: " << resultFile << endl; return;}

//    // Will write output per blocks, and from time-to-time append on the file, so it avoids mixed lines from different parallel simulations
//    stringstream blocksOut;

//    string nextSeq;
//    int cpt = 0;
//    do {
//        std::pair<int, string> next = rep2.nextSequenceAndID(true);
//        nextSeq = next.second; //rep2.nextSequence(true);
//        int ID = next.first;
//        cout << "treat" << ID << endl;
//        if(nextSeq.size() > 0){
//            //nextSeq = "CAGPSTTVPYYFDYW";
//            //cout << "Next=" << nextSeq << endl;
//            vector<string> cutSlides = slides(nextSeq, receptorSize+1);
//            size_t nSl = cutSlides.size();
//            for(size_t i = 0; i < nSl; ++i){
//                string subSeq = cutSlides[i];
//                vector<string> optimalStruct;
//                std::pair<double, double> affs = T3.affinity(subSeq, false, &optimalStruct);
//                size_t nS = optimalStruct.size();
//                if(cpt < 100) cout << ID << "\t" << nextSeq << "\t" << subSeq << "\t" << affs.first << "\t" << nS;
//                blocksOut << ID << "\t" << nextSeq << "\t" << subSeq << "\t" << affs.first << "\t" << nS;
//                for(size_t j = 0; j < nS; ++j){
//                    if(cpt < 100) cout << "\t" << optimalStruct[j];
//                    blocksOut << "\t" << optimalStruct[j];
//                }
//                if(cpt < 100) cout << endl;
//                blocksOut << endl;
//            }
//            if(cpt == 99) cout << "   ... now continues in silent mode ... " << endl;
//            cpt++;
//            if((cpt % 100) == 0){
//                cout << "-- write " << cpt << "/" << nBlock << endl;
//                fout << blocksOut.str();
//                blocksOut.str("");
//                blocksOut.clear(); // obviously these two functions fail to clear the stirngstream (wtf). Clear seems to remove the errs. str() doesnt do anything.
//                blocksOut = stringstream(); // so, use brute force rplacement
//            }
//        }
//    } while((nextSeq.size() > 0) && (cpt < 1e7));

//    cout << "-- write " << cpt << "/" << nBlock << " -- FINISHED!" << endl;
//    fout << blocksOut.str();
//    blocksOut.str("");
//    blocksOut.clear();

//    fout.close();
//    // now says that the sequences are treated
//    rep2.writeUpdatedRepertoire();
//}


// Adding a sequence, like for merging repertoires together, it will check if the same ID or same sequence exist.
// Adds a sequence. If a sequence is given whose ID already exists, tests the sequence. If sequence the same, do not add. If sequence is
// different, then raises an error. If ID=-1, will find a new ID for it (caution, might collide with IDs of added sequences later that would have
// this ID. If mergeIDsWhenSequenceIsEqual, a sequence already existing will not be added even if it has a non-used ID.
void importRepertoire::addSequence(std::string sequence, int ID, bool treated, bool mergeIDsWhenSequenceIsEqual){
//    static int cptErrors = 0;

    std::vector<int>::iterator it = find (listIDs.begin(), listIDs.end(), ID);
    // if the ID is found
    if (it != listIDs.end()){
         size_t position = std::distance(listIDs.begin(), it);
         if(listIDs.at(position) != ID) {
             cerr << "ERR: importRepertoire::addSequence(), ERR when finding ID position " << position << " in sequence vector" << endl;
            return;
         }
         string foundSeq = sequences.at(position);
//         if(foundSeq != sequence){
//             cptErrors++;
//             if(cptErrors < 20) cerr << "ERR: importRepertoire::addSequence(), adding an entry with same ID but difference sequence. Will shut up at 20x this error" << endl;
//            return;
//         }
     }
     if(mergeIDsWhenSequenceIsEqual){
        std::vector<string>::iterator it = find (sequences.begin(), sequences.end(), sequence);
        // sequence already found
        // adds the ID, but not sequence and says treated
        if (it != sequences.end()){
             listIDs.push_back(ID);
             sequences.push_back("");
             alreadyTreated.push_back(true);
             nLines++;
             return;
        }
     }
    listIDs.push_back(ID);
    sequences.push_back(sequence);
    alreadyTreated.push_back(true);
    nLines++;
}



importRepertoire::importRepertoire(string fileName, bool thirdColumn) : currentIndex(0)
{
    nLines = 0;
    currentIndex = -1;
    firstBlockedLine = -1;
    lastBlockedLine = -1;

    ifstream f;
    f.open(fileName.c_str());
    if(!f){cerr << "ERR::importRepertoire::importRepertoire, file not found: " << fileName << endl; return;}
    cout << "   ... Opening repertoire file: " << fileName << endl;
    originalFileName = fileName;
    while(f.good()){
        int ID = -1;
        string sequence;
        bool treated = false;
        f >> ID >> sequence;
        if(thirdColumn) f >> treated;
        if(ID < 0) return; // stupid iostream continues reading one line even if file is finished.
        listIDs.push_back(ID);
        sequences.push_back(sequence);
        alreadyTreated.push_back(treated);
        nLines++;
    }
    f.close();
    if(sequences.size() != alreadyTreated.size()){
        cerr << "ERR: Pêche melba, importRepertoire::importRepertoire has read a different number of sequences and treated values" << endl;
    }
}
std::pair<int, string> importRepertoire::nextSequenceAndID(bool setTreated, int lastLine){
    string res = nextSequence(setTreated, lastLine);
    int ID = -1;
    if(res.size() > 0){
        if(currentIndex <= lastLine){
            if(currentIndex < listIDs.size()){
                ID = this->listIDs[currentIndex];
            }
        }
    }
    return std::pair<int, string>(ID, res);
}

string importRepertoire::nextSequence(bool setTreated, int lastLine){
    currentIndex++;
    while((currentIndex <= lastLine) && (currentIndex < sequences.size()) && (alreadyTreated[currentIndex])){
        currentIndex++;
    }
    if((currentIndex == sequences.size()) || (currentIndex == lastLine+1)){
        return string("");
    }
    if((lastBlockedLine >= 0) && (currentIndex > lastBlockedLine)){
        return string("");
    }
    if(setTreated) alreadyTreated.at(currentIndex) = true;
    return sequences.at(currentIndex);
}


bool importRepertoire::writeUpdatedRepertoire(string fileName){
    ofstream f;
    //if(fileName.size() == 0) fileName = originalFileName;

    if((sequences.size() != alreadyTreated.size()) || (sequences.size() != listIDs.size())){
        cerr << "ERR: 'Île flottante', importRepertoire::writeUpdatedRepertoire has read a different number of sequences and treated values" << endl;
        return false;
    }
    stringstream res;
//    if(blockLines > 0){
//        string useless = nextSequence(false); // just to get the current Index right
//        firstBlockedLine = static_cast<int>(currentIndex);
//        //int remainToBlock = blockLines;
//        for(size_t i = 0; i < sequences.size(); ++i){
//            if((i >= static_cast<size_t>(firstBlockedLine)) && (remainToBlock > 0)){
//                lastBlockedLine = static_cast<int>(i);
//                remainToBlock--;
//                res << listIDs[i] << "\t" << sequences[i] << "\t" << ((alreadyTreated[i]) ? 1 : 2) << "\n";
//            } else {
//                res << listIDs[i] << "\t" << sequences[i] << "\t" << alreadyTreated[i] << "\n";
//            }
//        }
//    } else {
        for(size_t i = 0; i < sequences.size(); ++i){
            res << listIDs[i] << "\t" << sequences[i] << "\t" << alreadyTreated[i] << "\n";
        }
//    }


    // I prefer first putting things into a sstream and then writing to avoid that two programs would write lines at the same time.
    f.open(fileName);
    if(!f){cerr << "ERR::importRepertoire::writeUpdatedRepertoire, file not found: " << fileName << endl; return false;}
    f << res.str();
    f.close();

    return true;
}

// test:
void testImportRepertoire(){
    ofstream f("TestRep.txt");
    f << "1\tSEQUENCE_UNE" << endl;
    f << "3\tSEQUENCE_TROIS" << endl;
    f << "2\tSEQUENCE_DEUX" << endl;
    f << "4\tSEQUENCE_QUATRE" << endl;
    f << "5\tSEQUENCE_CINQ" << endl;
    f << "6\tSEQUENCE_SIX" << endl;
    f << "7\tSEQUENCE_SEPT" << endl;
    f << "8\tSEQUENCE_HUIT" << endl;
    f.close();

    importRepertoire T = importRepertoire("TestRep.txt");
    cout << "Picked " << T.nextSequence(true) << endl;
    T.writeUpdatedRepertoire("TestRep2.txt");
    importRepertoire T2 = importRepertoire("TestRep2.txt");
    cout << "Picked " << T.nextSequence(true) << endl;
    cout << "Picked " << T.nextSequence(true) << endl;
    cout << "Picked " << T.nextSequence(true) << endl;
    cout << "Picked " << T.nextSequence(true) << endl;
    cout << "Picked " << T.nextSequence(true) << endl;
    cout << "Picked " << T.nextSequence(true) << endl;
}




