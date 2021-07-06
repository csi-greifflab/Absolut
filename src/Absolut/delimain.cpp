// This is Absolut (main file), C++ code to discretize PDB antigens and generate bindings of short peptides around it.

// Three ways to compile it:
//      Using Absolut.pro, it includes latFit library, Qt user interface, openGL visualisation of antigen-receptor structures
//          this .pro file automatically includes "#define ALLOW_GRAPHICS", that will include the openGL code (see plot3d2h).
//      Using AbsolutNoLib.pro, it includes only the code to generate bindings. No GUI, no openGL. Can be compiled on any computer without any pre-installed library.
//          this .pro file automally includes "#define NO_QT" that will exclude any Qt dependent code, so pdb.h is not included and neither latFit.
//          ALLOW_GRAPHICS is node defined, so the openGL code is also excluded
//      Using AbsolutBoLibMPI.pro, it includes only the code to generate bindings, but allows the use of MPI.
//          additionnaly, it includes "#define ALLOW_GRAPHICS", and the MPI library is included/linked, and some scripts are made parallel.
//
// Note, even without MPI, Absolut can be used multi-thread up to 50 threads. The MPI version is more intended for High Performance Computing architectures.
// See common.h for more options on including / excluding libraries

#include "../Ymir/ymir.h"
#include "../Tools/md5.h"
#include "motifFeatures.h"
#include "quality.h"
#include "selfEvo.h"
#include "../Tools/zaprandom.h"
#include "../Tools/dirent.h"
#include <regex>
#include "antigenLib.h"
#include "importrepertoire.h"
#include "poolstructs.h"
#include "epitope.h"
#include "html.h"
#include "fileformats.h"
#include "topology.h"
#include <iostream>
#include <fstream>

// for absolut this is 10 and 11
#define DefaultReceptorSizeBonds 10
#define DefaultContactPoints 11

// This option includes the crystal structure from the PDB antibody that was binding the antigen
#define defaultIncludePDBantibody false

// for sleep
#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif

#include <cstdlib>

#ifndef NO_LIBS
#include "discretize.h"
#endif


using namespace std;

#ifndef NOQT
#include <QApplication>
#include "pdb.h"
#endif

#ifdef USE_MPI
#include "mpi.h"
// to run from qtreator,
// add a new executable in the run panel => mpiexec.exe
// command: -n 2 release/AbsolutNoLib.exe repertoire 1ADQ ShortSequences.txt
// and add the working folder as the main compiling one
//https://stackoverflow.com/questions/31783667/set-number-of-processes-mpi-in-cmake-project-at-qtcreator
#endif

#include "pthread.h"

static pthread_mutex_t lockSaveCommonDataset =  PTHREAD_MUTEX_INITIALIZER;

// The main.cpp is below. First some separate uses of the code are defined / shown when absolut is called without argument.

string getHelp(int argc, char** argv){
    if(argc < 1) return string("Executable with no name?");
    stringstream res;
    string programName = string(argv[0]).substr(string(argv[0]).find_last_of("/\\") + 1);
    res << " -- How to use Absolut -- \n";
    res << " \n";
    #ifndef NOQT
    res << " ====== !!! Note: Option 1 and 7 are not available because they need libraries (Qt, latfit) - use Absolut.pro for that ====== !!!" << endl;
    #endif
    #ifdef USE_MPI
    res << " ====== Note: you are using the parallelized version, using MPI (interesting for option 2 only). For the other options, only one CPU will be used ===== " << endl;
    #endif
    res << " Option 1: Open graphical interface to discretize an antigen using latfit. By giving more arguments, precomputations are done automatically.\n";
    res << " " << programName << " discretize  [PDB_ID Chains_of_interest]  [Resolution (5.25 default)]  [type_Discretization (CA / CoM / FuC default)]\n";
    res << "           ex: " << programName << " discretize\n";
    res << "           ex: " << programName << " discretize 1CZ8 VW\n";
    res << "           ex: " << programName << " discretize 1CZ8 VW 5.25\n";
    res << "           ex: " << programName << " discretize 1CZ8 VW 5.25 FuC\n";
    res << "           - please do not run this option in debug, the call to python pdb-tools scripts using system() usually fails (no idea why)\n";
    res << "           - please keep this order of options -\n";
    res << " \n";
    res << " Option 2: Generate a dataset of bindings of a repertoire of CDR3 sequences around an antigen (from the library)\n";
    res << " " << programName << " repertoire ID_antigen RepertoireFile [nThreads] [startIDSeq endIDseq]\n";
    res << "           ex: " << programName << " repertoire 1CZ8_VW sequences.txt\n";
    res << "           ex: " << programName << " repertoire 1CZ8_VW sequences.txt 10\n";
    res << "           ex: " << programName << " repertoire 1CZ8_VW sequences.txt 10 IDinFileNames\n";
    res << "           ex: " << programName << " repertoire 1CZ8_VW sequences.txt 10 15000 50000\n";
    res << "           ex: " << programName << " repertoire 1CZ8_VW sequences.txt 10 IDinFileNames 15000 50000\n";
    res << "           - you are allowed to use the short PDB ID of an antigen, provided there is only one 'version' of it in the system -\n";
    res << "           ex: " << programName << " repertoire 1CZ8 sequences.txt\n";
    res << "           - you are allowed to use the the integer ID of the antigen instead of its name. use listAntigens to get the IDs. -\n";
    res << "           ex: " << programName << " repertoire 5    sequences.txt\n";
    res << "           - Note: we have decided to use as default, sliding windows of 11 Amino Acids and at least 11 interactions for the structures. -\n";
    res << " \n";
    res << " Option 3: See the list of antigens in the library\n";
    res << " " << programName << " listAntigens\n";
    res << " \n";
    res << " Option 4: Get the binding information of one CDR3 around one antigen\n";
    res << " " << programName << " singleBinding ID_antigen CDR3Sequence\n";
    res << "           ex: " << programName << " singleBinding 1CZ8_VW CAGPSTTVPYYFDY\n";
    res << "           - you are allowed to use the shorter PDB ID antigen name or its ID in the library \n";
    res << "           - Note: we have decided to use as default: sliding windows of 11 Amino Acids and at least 11 interactions for the structures. -\n";
    res << " \n";
    res << " Option 5: Analyse the binding features of a binding dataset\n";
    res << " " << programName << " getFeatures ID_antigen bindingDatasetFile outputFeaturesFile [degree=1] [includeDegree=false] \n";
    //res << "           - note: the dataset file should have the following format:\n";
    res << "           ex: " << programName << " getFeatures 1FBI_X 1FBI_X_Mascotte.txt AnalyzedFile.txt\n";
    res << "           ex: " << programName << " getFeatures 1FBI_X 1FBI_X_Mascotte.txt AnalyzedFile.txt 2\n";
    res << "           ex: " << programName << " getFeatures 1FBI_X 1FBI_X_Mascotte.txt AnalyzedFile.txt 1 true\n";
    //res << " Option 5b: POOL the binding features of binding datasets (from a folder and the file pattern, like BindingMyAntigen*.txt or so)\n";
    //res << " " << programName << " poolFeatures ID_antigen folderFiles patternFileNames outputFeaturesFile outputListTreatedSequencesFile\n";
    //res << "          - note:  Directory should finish by /, very important\n";
    //res << " \n";
    //res << " Option 6: Check how many sequences have been successfully analysed inside a dataset / remaining sequences.\n";
    //res << " " << programName << "checkBindings ID_antigen bindingDatasetFile\n";
    res << " \n";
    res << " Option 7: Visualization tools: antigen; list of structures around it and binding hotspots.\n";
    res << " " << programName << "visualize ID_antigen\n";
    res << " " << programName << "visualize ID_antigen bindingDatasetFile [maxEnergy=1e6]\n";
    res << " " << programName << "hotspots ID_antigen bindingDatasetFile [maxEnergy=1e6] [sizeKsets=4] [degree=1]\n";
    res << " " << programName << "interface ID_antigen CDRH3seq\n";
    res << " \n";
    res << " Option 8: Generate Batch scripts to run repertoire analyses on a slurm cluster for all library antigens.\n";
    res << "          - Generates per antigen: BatchXXX.sh (72 / 120 hrs max), SlurmTestXXX.sh using short file (5 mins max)\n";
    res << "          - Together with GetXXXX_X_Structures.sh to download the structures from the database\n";
    res << " " << programName << "batch FRAM/SAGA/FLORIDA FileCDR3s [nbCores=4] [ShortTestFileCDR3s]\n";
    res << " \n";
    res << " Option 9: Generate html report, for one antigen. If bindings are available, give folder.\n";
    res << " " << programName << "html ID_antigen [outputFolfer] [folderWithBindings]\n";
    res << " \n";
    res << " Option 10: Develop affinities for random sequences if they would bind with a specific structure.\n";
    res << "          - either by giving the antigen and the structure explicitely:\n";
    res << " " << programName << "develop ID_antigen startingPosStructure structure nrSequencesToGenerate\n";
    res << "            ex: develop 1FBI_X 141220 LLDDRDUUDR 20\n";
    res << "          - or by directly giving the interaction code of a structure (doesn't need to give the antigen / structure):\n";
    res << " " << programName << "develop interactionCode sizeRecetors nrSequencesToGenerate\n";
    res << "            ex: develop iFkFhEcSgSdNfNdNdTaNcNaNcGcWbVgAhKhDjDkVkIbebgfi 11 20\n";
    res << " \n";
    res << " Others: Simple tools.\n";
    res << " " << programName << "info_position   X Y Z  => Transform into lattice position\n";
    res << " " << programName << "info_position   N      => Transform lattice position N into X,Y,Z \n";
    res << " " << programName << "info_structure  pos codeStructure [AAs] \n";
    res << " " << programName << "info_structure  pos1 codeStructure1 pos2 codeStructure2 ... [AAs] \n";
    res << " " << programName << "info_antigens \n";
    res << " " << programName << "info_antigen    AntigenID \n";
    res << " " << programName << "info_fileNames  AntigenID \n";
    res << " " << programName << "info_IDstructure  pos codeStructure \n";
    return res.str();
}

// Only functions headers here, the functions' content is defined at the bottom of the file.

// discretize
int option1(string PDB_ID = "", string chains = "", double resolution = 5.25, string typeDiscrete = "FuC");

// repertoire
void option2(string ID_antigen, string repertoireFile, int nThreads = 1, string prefix = string(""), int startingLine = 0, int endingLine = 1000000000);

// listAntigens
void option3(){showListIDs();}

// singleBinding
void option4(string ID_antigen, string CDR3);

// getFeatures
void option5(string ID_antigen, string bindingDatasetFile, string outputFeaturesFile, int minDegree = 1, bool includeDegreeInMotifs = false);

// checkBindings - this option is obsolete
void option6(string ID_antigen, string bindingDatasetFile){cout <<"Sorry, option 6 is not available" << endl;}

// visualize and visualize_hotspot
void option7(string ID_antigen, bool generateHotspots = false, string bindingDatasetFile = "", double maxE = 1e+6, int sizeKsets = 4, int degree = 1);

// visualize only paratope and epitope of one binding
void option7b(string ID_antigen, string sequence);

// batch
void option8(string fileCDR3s, int nbNodes, string shortTestFile, string architecture);

// html
void option9(string ID_antigen, string outputFolder = "", string folderWithBindings = "", string folderWithStructures = "");

// info_position
void info_position(int X, int Y, int Z){cout << lattice::idFromPosisition(X, Y, Z) << endl;}
void info_position(int N){cout << printVector(lattice::positionFromID(N)) << endl;}
void info_fileNames(string ID_antigen);

// info_antigens
void infosAllAntigens(){
    vector<string> toParse = listIDs();
    stringstream res;
    res << "no\tID\tAAseq\tnAAs\tsurfaceAAs\tcore\tlarge\n";
    for(size_t i = 0; i < toParse.size(); ++i){
        //std::pair<superProtein*, vector<int> > AGi = getAntigen(toParse[i]);
        antigenInfo AGi = getAntigenInfos(toParse[i]);

        vector<vector<int> > hotspotsCore = AGi.hotspotsCore;
        set<int> accuCore;
        for(size_t k = 0; k < hotspotsCore.size(); ++k){
            set<int> thisone = set<int>(hotspotsCore[k].begin(), hotspotsCore[k].end());
            accuCore.insert(thisone.begin(), thisone.end());
        }
        size_t residuesCore = accuCore.size();

        vector<vector<int> > hotspotsLarge = AGi.hotspotsLarge;
        set<int> accuLarge;
        for(size_t k = 0; k < hotspotsLarge.size(); ++k){
            set<int> thisone = set<int>(hotspotsLarge[k].begin(), hotspotsLarge[k].end());
            accuLarge.insert(thisone.begin(), thisone.end());
        }
        size_t residuesLarge = accuLarge.size();

        if(true){
            int receptorSize = DefaultReceptorSizeBonds; //10;
            int minInteract = DefaultContactPoints; //11;

            std::pair<superProtein*, vector<int> > AG = getAntigen(toParse[i]);
            res << fnameStructures(AG.first, receptorSize, minInteract, AG.second) << "\t";
        }

        res << i << "\t" << toParse[i] << "\t" <<  AGi.first->getAAseq() << "\t" << AGi.first->getAAseq().size() << "\t" << getSurfaceAAs(AGi.first) << "\t" << residuesCore << "\t" << residuesLarge << "\n";
    }
    cout << res.str();
}

// These functions will generate many affinities but only for one structure.
void option10a(string interactionCode, int sizeReceptors, int nrSequences);
void option10b(string antigenID, string structureID, int nrSequences); // the sizeReceptors will be taken from the structureID



void option10a(string interactionCode, int sizeReceptors, int nrSequences){
    if(nrSequences > 50000000) {cerr << "ERR: option10a, max 50 000 000 sequences, requested " << nrSequences <<"\n"; return;}
    string res;
    for(int i = 0; i < nrSequences; ++i){
        string receptor = randomProt(sizeReceptors);
        double aff = affinityCodeTot(receptor, interactionCode);
        stringstream buf; buf << receptor << "\t" << aff << "\n";
        res.append(buf.str());

        if(i % 1000000 == 0){
            cout << res;
            res.clear();
        }
    }
    cout << res;
}


void option10b(string ID_antigen, int startingPos, string structure, int nrSequences){

    string AntigenName = "";
    if(ID_antigen.size() < 4){
        AntigenName = IDshortcut(std::stoi(ID_antigen));
    } else {
        AntigenName = IDshortcut(ID_antigen);
    }

    cout << "   ... loading antigen " << AntigenName << " from the library" << endl;
    std::pair<superProtein*, vector<int> > AG = getAntigen(AntigenName);

    int nrAAs = static_cast<int>(structure.size() + 1);
    superProtein* receptor = new superProtein(structure, startingPos);
    string interactionCode = interactions(*(AG.first), *receptor);
    cout << "Interaction Code is " << interactionCode << endl;
    option10a(interactionCode, nrAAs, nrSequences);
}


void infoOneAntigen(string ID_antigen){
    string AntigenName = "";
    if(ID_antigen.size() < 4){
        AntigenName = IDshortcut(std::stoi(ID_antigen));
    } else {
        AntigenName = IDshortcut(ID_antigen);
    }

    antigenInfo AGinfos = getAntigenInfos(AntigenName);
    cout << AGinfos.print() << endl;
}






































int main(int argc, char** argv){

    // If using MPI, only option 2 can be run multi-processes, the other ones will only be started by the main process, the other ones will do nothing
    int nJobs = 1;
    int rankProcess = 0;
    #ifdef USE_MPI
//  If MPI is used (amd compiled with), it will just start independent Jobs, with a certain ID (to split sequences to treat)
//  MPI_Init(nullptr, nullptr);
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    if (provided != MPI_THREAD_FUNNELED) {
        cerr << "ERR: The MPI installed in this computer seems to be incompatible with threads" << endl;
        exit(-1);
    }
    MPI_Comm_size(MPI_COMM_WORLD, &nJobs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankProcess);

    #endif

    if(argc <= 1){
        cout << getHelp(argc, argv);
        return 0;
    }

    bool success = false;
    string command = string(argv[1]);

    if(!command.compare("repertoire")){
        success = true;
        switch(argc){
        //" repertoire 1CZ8_VW sequences.txt\n";
        //" repertoire 1CZ8_VW sequences.txt 10\n";
        //" repertoire 1CZ8_VW sequences.txt 10 IDinFileNames\n";
        //" repertoire 1CZ8_VW sequences.txt 10 15000 50000\n";
        //" repertoire 1CZ8_VW sequences.txt 10 IDinFileNames 15000 50000\n";
        // void option2(string ID_antigen, string repertoireFile, int nThreads = 1, string prefix = string(""), int startingLine = 0, int endingLine = 1000000000)
        case 4: option2(argv[2], argv[3]); break;
        case 5: option2(argv[2], argv[3], atoi(argv[4])); break;
        case 6: option2(argv[2], argv[3], atoi(argv[4]), argv[5]); break; // default 1 thread
        case 7: option2(argv[2], argv[3], 1, "", atoi(argv[4]), atoi(argv[5])); break; // default 1 thread
        case 8: option2(argv[2], argv[3], atoi(argv[4]), argv[5], atoi(argv[6]), atoi(argv[7])); break;
        default: success = false; cerr << "ERR: couldn't get the correct number of arguments. Expected 'repertoire + 4 arguments. got in total " << argc-1 << endl;
        }
    }

    // these options are not made for multi process, so only one process will do it
    if(rankProcess == 0){
        if(!command.compare("discretize")){
            success = true;
            switch(argc){
            case 2: option1(); break;
            case 3: option1(argv[2]); break;
            case 4: option1(argv[2], argv[3]); break;
            case 5: option1(argv[2], argv[3], atof(argv[4])); break;
            case 6: option1(argv[2], argv[3], atof(argv[4]), argv[5]); break;
            default: success = false; cerr << "ERR: couldn't get the correct number of arguments. Expected 'discretize + 1,2,3 or 4 arguments. got in total " << argc-1 << endl;
            }
        }

        if(!command.compare("listAntigens")){
            success = true;
            option3();
        }

        if(!command.compare("singleBinding")){
            success = true;
            if(argc == 4) option4(argv[2], argv[3]);
            else cerr << "ERR: couldn't get the correct number of arguments. Expected 'listAntigens + 2 arguments. got in total " << argc-1 << endl;
        }

        if(!command.compare("getFeatures")){
            success = true;
            switch(argc){
            //case 4: option5(argv[2], argv[3]); break;
            case 5: option5(argv[2], argv[3], argv[4]); break;
            case 6: option5(argv[2], argv[3], argv[4], atoi(argv[5])); break;
            case 7: {
                string incl = string(argv[6]);
                if(incl.compare("true") && incl.compare("True") && incl.compare("TRUE") && incl.compare("false") && incl.compare("False") && incl.compare("FALSE")){
                    cerr << "ERR: getFeatures, last argument " << argv[6] << " should be true/True/TRUE/false/False/FALSE. Please retry" << endl;
                    return -1;
                }
                option5(argv[2], argv[3], argv[4], atoi(argv[5]), ((!incl.compare("true")) || (!incl.compare("True")) || (!incl.compare("TRUE")))); break;
            }
            default: success = false; cerr << "ERR: couldn't get the correct number of arguments. Expected 'getFeatures + 3, 4, or 5 arguments. got in total " << argc-1 << endl;
            }
        }

        if(!command.compare("poolFeatures")){
            success = true;
            switch(argc){
            //case 7: option5b(argv[2], argv[3], argv[4], argv[5], argv[6]); break;
            default: success = false; cerr << "ERR: couldn't get the correct number of arguments. Expected 'poolFeatures + 5 arguments. got in total " << argc-1 << endl;
            }
        }

//        if(!command.compare("checkBindings")){
//            success = true;
//            if(argc == 4) option6(argv[2], argv[3]);
//            else cerr << "ERR: couldn't get the correct number of arguments. Expected 'checkBindings + 2 arguments. got in total " << argc-1 << endl;
//        }

        if(!command.compare("visualize")){
            success = true;
            switch(argc){
            case 3: option7(argv[2], false); break;
            case 4: option7(argv[2], false, argv[3]); break;
            case 5: option7(argv[2], false, argv[3], atof(argv[4])); break;
            default: success = false; cerr << "ERR: couldn't get the correct number of arguments. visualize + 1,2 or 3 arguments. got in total " << argc-1 << endl;
            }
        }

        if(!command.compare("interface")){
            success = true;
            switch(argc){
            case 4: option7b(argv[2], argv[3]); break;
            default: success = false; cerr << "ERR: couldn't get the correct number of arguments. visualize + 2 argument. got in total " << argc-1 << endl;
            }
        }

        //./Absolut visualize ID_antigen\n";
        //./Absolut visualize ID_antigen bindingDatasetFile [maxEnergy=1e6]\n";
        //./Absolut hotspots ID_antigen bindingDatasetFile [maxEnergy=1e6] [sizeKsets=4] [degree=1]\n";

        if(!command.compare("hotspots")){
            success = true;
            switch(argc){
            case 4: option7(argv[2], true, argv[3]); break;
            case 5: option7(argv[2], true, argv[3], atof(argv[4])); break;
            case 6: option7(argv[2], true, argv[3], atof(argv[4]), atoi(argv[5])); break;
            case 7: option7(argv[2], true, argv[3], atof(argv[4]), atoi(argv[5]), atoi(argv[6])); break;
            default: success = false; cerr << "ERR: couldn't get the correct number of arguments. visualize_hotspots + 1,2 or 3 arguments. got in total " << argc-1 << endl;
            }
        }

        if(!command.compare("info_position")){
            success = true;
            switch(argc){
            case 3: info_position(atoi(argv[2])); break;
            case 5: info_position(atoi(argv[2]), atoi(argv[3]), atoi(argv[4])); break;
            default: success = false; cerr << "ERR: couldn't get the correct number of arguments. info_positions + 1 or 3 arguments. got in total " << argc-1 << endl;
            }
        }

        //batch FRAM/SAGA/FLORIDA FileCDR3s [nbCores=4] [ShortTestFileCDR3s]
        if(!command.compare("batch")){
            success = true;
            switch(argc){
            case 4: option8(argv[3], 4, string(""), argv[2]); break;
            case 5: option8(argv[3], atoi(argv[4]), string(""), argv[2]); break;
            case 6: option8(argv[3], atoi(argv[4]), argv[5], argv[2]); break;
            default: success = false; cerr << "ERR: couldn't get the correct number of arguments. batch + 1 or 2 arguments. got in total " << argc-1 << endl;
            }
        }
//        if(!command.compare("info_IDstructure")){
//            success = true;
//            switch(argc){
//                case 3: option8(argv[2], 4, argv[2]); break;
//                default: success = false; cerr << "ERR: couldn't get the correct number of arguments. info_IDstructure + 1 argument. got in total " << argc-1 << endl;
//            }
//        }
        if(!command.compare("info_antigen")){
            success = true;
            switch(argc){
                case 3: infoOneAntigen(argv[2]); break;
                default: success = false; cerr << "ERR: couldn't get the correct number of arguments. info_antigen + 1 arguments. got in total " << argc-1 << endl;
            }
        }
        if(!command.compare("info_antigens")){
            success = true;
            switch(argc){
                case 2: infosAllAntigens(); break;
                default: success = false; cerr << "ERR: couldn't get the correct number of arguments. info_antigens + 0 arguments. got in total " << argc-1 << endl;
            }
        }
        if(!command.compare("html")){
            success = true;
            switch(argc){
                case 3: option9(argv[2]); break;
                case 4: option9(argv[2], argv[3]); break;
                case 5: option9(argv[2], argv[3], argv[4]); break;
                case 6: option9(argv[2], argv[3], argv[4], argv[5]); break;
                default: success = false; cerr << "ERR: couldn't get the correct number of arguments. html + 1 to 4 arguments. got in total " << argc-1 << endl;
            }
        }

        if(!command.compare("develop")){
            success = true;
            switch(argc){
                // option10a string interactionCode, int sizeReceptors, int nrSequences
                case 5: option10a(argv[2], atoi(argv[3]), atoi(argv[4])); break;
                // option10b string ID_antigen, int startingPos, string structure, int nrSequences
                case 6: option10b(argv[2], atoi(argv[3]), argv[4], atoi(argv[5])); break;
                default: success = false; cerr << "ERR: couldn't get the correct number of arguments. html + 1 to 4 arguments. got in total " << argc-1 << endl;
            }
        }

        if(!command.compare("info_fileNames")){
            success = true;
            switch(argc){
                case 3: info_fileNames(argv[2]); break;
                default: success = false; cerr << "ERR: couldn't get the correct number of arguments. html + 1 to 4 arguments. got in total " << argc-1 << endl;
            }
        }

        if(!success) cerr << "Sorry, command " << command << " could not be found/understood. Run without argument for help" << endl;
    }

#ifdef USE_MPI
    cerr << "Job" << rankProcess + 1 << " Ends up MPI " << endl;
MPI_Finalize();
#endif

    return 0;
}






































// ================= Option 1: Launching the User Interface to discretize antigens =====================.
//int option1(string PDB_ID = "", string chains = "", double resolution = 5.25, string typeDiscrete = "FuC"){
int option1(string PDB_ID, string chains, double resolution, string typeDiscrete){

    #ifndef NOQT
    char* argv[] = {(char*) "Nothing!"};
    int argc = 1;
    QApplication appl(argc, argv);
    bool pipeline = false; // this option would trigger an automatic comparison/iteration of lots of lattice resolutions
    PDB* a2 = new PDB(PDB_ID, chains, pipeline, resolution, typeDiscrete);
    // use example: PDB* a2 = new PDB("1CZ8", "VW");
    a2->show();
    appl.exec();
    #else
    cerr << "You are intending to run the graphical interface of Absolut. Please make sure to use Zapotec.pro and that NOQT / NO_LIBS are not defined." << endl;
    #endif
    return 0;
}





// ================= Option 2: Multithread repertoire binding calculation =====================.

// In order to use threads, the computation function needs to receive only one argument as void*
// therefore, we pool all the arguments into one structure.
// http://www.cse.cuhk.edu.hk/~ericlo/teaching/os/lab/9-PThread/Pass.html
struct argsThread {
    argsThread(){}
    affinityOneLigand* T3;
    vector< std::pair<string, string> >* listToProcess;
    string resultFile;
    string identificationThread;
    string antigenName;
    int receptorSize;
    dataset<binding>* savingLocation;
};



// Main function for the repertoire option, that will be run on each single thread separately.
void *oneThreadJob(void *input){
    // 1 - Reconstituting the arguments from the pointer to the class argsThreads
    string resultFile = ((struct argsThread*)input)->resultFile;
    int receptorSize = ((struct argsThread*)input)->receptorSize;
    affinityOneLigand* T3 =  ((struct argsThread*)input)->T3;
    vector< std::pair<string, string> >* listToProcess = ((struct argsThread*)input)->listToProcess;
    string antigenName =  ((struct argsThread*)input)->antigenName;
    string identificationThread = ((struct argsThread*)input)->identificationThread;
    dataset<binding>* savingLocation = ((struct argsThread*)input)->savingLocation;

    // 2- Preparing output by two ways:
    // a/ writing in a file. It was problems by keeping the ofstream opened, so I will close it each time and do append
    // first, erases the file if existed
    {
        ofstream fout(resultFile.c_str()); //, ios::app);
        if(!fout) {
            cerr << "ERR: Couldn't open output file: " << resultFile << endl;
            return nullptr;
            //pthread_exit(nullptr); => This is bad because if only one thread, it would destroy the main thread...
        }
        fout << antigenName << "\n";
        fout.close();
    }
    string blockToWrite;    // stores output by blocs to write 100 lines by 100 lines into the result files
                            // Note: it ended up impossible to clear a stringstream, that's why I store as a string,
                            // but i build by blocks using stringstream (see blocksOut)

    // Treats each CDR3 sequence (from the list to process by this thread)
    size_t Nseq = listToProcess->size();
    cout << "Thread " << identificationThread << " got " << Nseq << " sequences to proceed " << endl;

    for(size_t i = 0; i < Nseq; ++i){
        std::pair<string, string> next = listToProcess->at(i);
        string CDR3seq = next.second;
        string ID = next.first;

        // cout << "Process " << next.first << "=> " << next.second << endl;
        if(CDR3seq.size() > 0){
            stringstream blocksOut;
            //cout << "Thread " << identificationThread << " treats ID " << ID << endl;
            vector<string> cutSlides = slides(CDR3seq, receptorSize+1);
            size_t nSl = cutSlides.size();

            // 1st step: calculate the affinities of each sliding window and finds the best(or equally best) slides
            // Note: it's not very optimal to duplicate affinities, Luckily, the affinity() function remembers previously calculated ones
            double bestEnergy = +1000000;
            for(size_t si = 0; si < nSl; ++si){
                string subSeq = cutSlides[si];
                std::pair<double, double> affs = T3->affinity(subSeq, false);
                bestEnergy = min(affs.first, bestEnergy);
            }

            // 2: Now, knowing who is the best, re-requests affinity of each sliding windows and says 'best' for the equally best ones
            for(size_t si = 0; si < nSl; ++si){

                string subSeq = cutSlides[si];
                vector<string> optimalStructures;

                // Gets affinity and all optimal bindings (sometimes more than one)
                std::pair<double, double> affs = T3->affinity(subSeq, false, &optimalStructures);

                // Now we know what is the best binding energy for this CDR3
                bool isBest = (fabs(affs.first - bestEnergy) < 1e-6);

                #define commonSaving true
                if(commonSaving) pthread_mutex_lock(&lockSaveCommonDataset);

                size_t nS = optimalStructures.size();
                // For each optimal binding structure
                for(size_t j = 0; j < nS; ++j){

                    string IDthisSlide = combinedID(ID, static_cast<int>(si), static_cast<int>(j));

                    // only says what is happening for the top 10 sequences.
                    if(i < 10) cout << "Thread " << identificationThread << " " << IDthisSlide << "\t" << CDR3seq << "\t" << ((isBest) ? "true" : "false") << "\t" << subSeq << "\t" << affs.first << "\t" << optimalStructures[j] << "\n";

                    blocksOut << IDthisSlide << "\t" << CDR3seq << "\t" << ((isBest) ? "true" : "false") << "\t" << subSeq << "\t" << affs.first << "\t" << optimalStructures[j] << "\n";

                    if(commonSaving) if(savingLocation){
                        binding* b = new binding(subSeq, affs.first, optimalStructures[j]);
                        savingLocation->addLine(IDthisSlide, CDR3seq, b, isBest);
                    }
                }
                if(nS == 0) {
                     blocksOut << combinedID(ID, static_cast<int>(si), 0) << "\t" << CDR3seq << "\t" << ((isBest) ? "true" : "false") << "\t" << subSeq << "\t" << affs.first << "ERR:No_optimal_structure_???" << "\n";
                }
                if(commonSaving) pthread_mutex_unlock(&lockSaveCommonDataset);

            } // end for each sliding window

            blockToWrite += blocksOut.str();
            if(i == 10) cout << "Thread " << identificationThread << "   ... now continues in silent mode ... " << endl;


            // Writing in output file every 100 CDR3 sequences
            if((i % 100) == 0){
                cout << "Thread " << identificationThread << "-- write " << i << "/" << Nseq << endl;
                ofstream fout(resultFile.c_str(), ios::app);
                if(fout){
                    fout << blockToWrite;
                    fout.close();
                    blockToWrite.clear();
                } else {
                    // decided to put errors in cout in threads, because cerr cuts the couts into pieces
                    cout << "Thread " << identificationThread << ", had problems to write into file, keeping into memory till next try" << endl;
                    cout << "Thread " << identificationThread << ", file was:" << resultFile << endl;
                }
            }
        }
    }

    cout << "Thread " << identificationThread << "-- write FINAL " << Nseq << endl;

    {
        ofstream fout(resultFile.c_str(), ios::app);
        if(fout){
            fout << blockToWrite;
            fout.close();
        } else {
            // decided to put errors in cout in threads, because cerr cuts the couts into pieces
            cout << "Thread " << identificationThread << ", had problems to write into file, WILL NOW OUTPUT RESULT " << endl;
            cout << "Thread " << identificationThread << ", file was:" << resultFile << endl;
            cout << blockToWrite;
        }
    }

    //pthread_exit(nullptr);

    // normally this is done by pthreads, but compiler complains reaching end without return
    //cerr << "Will exit thread" << endl;
    return nullptr;
}

// the IDs start at 0
#define MAX_ALLOWED_THREADS 50
//void option2(string ID_antigen, string repertoireFile, int nThreads = 1, string prefix = string(""), int startingLine = 0, int endingLine = 1000000000){
void option2(string ID_antigen, string repertoireFile, int nThreads, string prefix, int startingLine, int endingLine){

    // Default options for our foldings:
    int receptorSize = DefaultReceptorSizeBonds; //10; // defined in number of bounds, add +1 to get the number of AAs
    int minInteract = DefaultContactPoints; //11;

    // nJobs will be the ID of this process (if MPI is used, each job will get a different rank), if not, there is only one job with rank 0
    int nJobs = 1;
    int rankProcess = 0;
    #ifdef USE_MPI
        // If MPI is used (amd compiled with), it will just start independent Jobs, with a certain ID (to split sequences to treat)
        //MPI_Init(nullptr, nullptr);
        MPI_Comm_size(MPI_COMM_WORLD, &nJobs);
        MPI_Comm_rank(MPI_COMM_WORLD, &rankProcess);
        cout << "MPI started!" << endl;
    #endif

    // This mutex control the access to shared memory inside affinityOneLigand::affinity()
    if (pthread_mutex_init(&lockAccessPrecompAffinities, nullptr) != 0){
        cerr << "\nERR: option repertoire, mutex lockAccessPrecompAffinities init failed, problem with pthreads?" << endl;
    }
    // This mutex controls writing the results on the common memory
    if (pthread_mutex_init(&lockSaveCommonDataset, nullptr) != 0){
        cerr << "\nERR: option repertoire, mutex lockSaveCommonDataset init failed, problem with pthreads?" << endl;
    }

    // each process gets an ID that will be used for every communication, because they will output at any time/order
    stringstream IDjob; IDjob << "Job" << rankProcess+1 << "/" << nJobs;
    string myID = IDjob.str();


    // => Checking inputs one by one

    // 0 - Basics
    if(nThreads < 0) {cerr << "ERR: option repertoire, wrong number of threads, will take 1" << nThreads << endl; nThreads = 1;}
    if(nThreads > MAX_ALLOWED_THREADS) {cerr << "ERR: option repertoire, Does not allow more than " << MAX_ALLOWED_THREADS << " threads (requested: " << nThreads << "), will take 50" << endl; nThreads = MAX_ALLOWED_THREADS;}
    if(startingLine > endingLine){cerr << "ERR: option repertoire, the ending line " << endingLine << " is before the starting line " << startingLine << " -> ending" << endl; return;}

    // 1 - Reading the list of sequences to process:
    cout << myID << "   ... loading repertoire " << repertoireFile << "\n";
    //<< "\n       expected to have 2 or 3 columns: ID , CDR3 sequence , [optional tag]" << endl;
    repertoire rep = repertoire(repertoireFile);
    cout << myID << "   ... Found " << rep.nLines() << " lines/sequences " << endl;
    if(rep.nLines() < 1) return;

    // 2a - Loading the antigen from the library. Either call with a number, or with a substring of its name.
    string AntigenName = (ID_antigen.size() < 4) ? IDshortcut(std::stoi(ID_antigen)) : IDshortcut(ID_antigen);
    cout << myID << "   ... loading antigen " << AntigenName << " from the library " << endl;
    std::pair<superProtein*, vector<int> > AG = getAntigen(AntigenName);

    // 2b - Preparing or loading the structures for this antigen
    // Note: we make each process check the file exists, in case problems of copying/file access between different machines
    // first, check if the structures are available (or the prepared compact interaction codes for this AA sequence):
    string fStruct = fnameStructures(AG.first, receptorSize, minInteract, AG.second);
    ifstream f(fStruct.c_str());
    if(f.good()){f.close();}
    else {
        string fCompact =   fileNameCompactForAASeqLigand(AG.first, receptorSize, minInteract, AG.second);
        ifstream f(fStruct.c_str());
        if(f.good()){f.close();}
        else {
            if(rankProcess == 0){
                cout << "\nERR: the list of binding structures for this antigen could not been found in this folder..." << endl;
                cout << "     its calculation can take typically 12 to 50 hours, so we do not recompute structures inside " << endl;
                cout << "     the 'repertoire' option, which is made to treat lots of sequences in multithreads, and the " << endl;
                cout << "     calculation of structures is not multithreaded, so it would waste resources." << endl;
                cout << "     => Please either find the structures file on the Absolut server, " << endl;
                cout << "     or run this program with the option 'singleBinding' and one CDR3 AA sequence, it will compute " << endl;
                cout << "     the structures and save them in the current folder." << endl;
                cout << "\n";
                cout << "     For information, the lacking file is:" << endl;
                cout << "     " << fStruct << endl;
                cout << "\n";
                cout << "     Or, alternately, the file with precomputed interaction codes:" << endl;
                cout << "     " << fStruct << endl;
                cout << "Bye!" << endl;
            }
            // each process will close
            return;
        }
    }



    // One process only makes an affinity calculation to regenerate the compact file if necessary
    #ifdef USE_MPI
    if(rankProcess == 0){
        cout << "The main process will check everything is ready to calculate affinities, and regenerate compact files if necessary" << endl;
        affinityOneLigand Ttest = affinityOneLigand(AG.first, receptorSize, minInteract, -1, 1, AG.second);
        Ttest.affinity(string(receptorSize+1, 'A'));
        cout << "   -> Everything ready!" << endl;
    }
    int res2 = MPI_Barrier(MPI_COMM_WORLD);
    if (res2 != MPI_SUCCESS) cout << myID << ", problem with MPI_Barrier inside option2()" << endl;
    #endif

    affinityOneLigand T3 = affinityOneLigand(AG.first, receptorSize, minInteract, -1, 1, AG.second);
    // This shortcuts a lots of documenting/debugging calculations, and only calculates best energy
    T3.setUltraFast(true);

    // 3 - Separating the job:
    //      - between the processes (1 ... nJobs) => startingLineThisJob ... endingLineThisJob (both included)
    //      - then between the threads insied each process => Will be done later
    int minLine = max(0, startingLine);
    int maxLine = min(static_cast<int>(rep.nLines())-1, endingLine);
    int nToProcess = maxLine - minLine + 1;

    int nPerJob = nToProcess / nJobs;   // number per MPI process

    int startingLineThisJob = minLine + rankProcess * nPerJob; // + 10000;
    int endingLineThisJob = minLine + (rankProcess + 1) * nPerJob - 1;
    if((endingLineThisJob - startingLineThisJob) > 1e6){
        cout << "WRN: Be prepared, very high number of sequences for this process. Might take forever...!\n" << endl;
    }
    if(rankProcess == nJobs - 1) endingLineThisJob = maxLine; // just for rounding it up
    cout << myID << "    ... will process lines " << startingLineThisJob << " to " << endingLineThisJob << " and then split into " << nThreads << " threads " << endl;


    dataset<binding>* commonSavingLocation = new dataset<binding>();
    commonSavingLocation->setNameAntigen(AntigenName);

    pthread_t tid[MAX_ALLOWED_THREADS]; // we allow max 50 thread
    for(size_t i = 0; i < static_cast<size_t>(nThreads); ++i){

        // cut into equal blocks between each thread
        int nToProcessThreads = endingLineThisJob - startingLineThisJob + 1;
        int nPerThread = nToProcessThreads / nThreads;
        int startingLineThisThread = startingLineThisJob + static_cast<int>(i) * nPerThread;
        int endingLineThisThread = startingLineThisJob + (static_cast<int>(i) + 1) * nPerThread - 1;
        if(static_cast<int>(i) == nThreads - 1) endingLineThisThread = endingLineThisJob;
        //cout << myID << "t" << i << "    ... will process lines " << startingLineThisThread << " to " << endingLineThisThread << "" << endl;

        stringstream resultFile;  resultFile << prefix << "TempBindingsFor" << AntigenName << "_t" << i << "_Part" << rankProcess+1 <<"_of_" << nJobs << ".txt";
        vector< std::pair<string, string> >* listToProcess = new vector< std::pair<string, string> >(rep.getLines(startingLineThisThread, endingLineThisThread));


        stringstream threadName; threadName << "Proc" << rankProcess << "_t" << i << "/" << nThreads;
        cout << threadName.str() << ", will treat " << listToProcess->size() << ", i.e. lines " << startingLineThisThread << " to " << endingLineThisThread << " included,  and save result in " << endl;
        cout << threadName.str() << ", " << resultFile.str() << endl;

        //struct argsThread *arguments = (struct argsThread *) malloc(sizeof(struct argsThread)); // this did segfault, C way of doing, bad bad bad
        argsThread *arguments = new argsThread();
        arguments->resultFile = resultFile.str();
        arguments->receptorSize = receptorSize;
        arguments->T3 = &T3;
        arguments->listToProcess = listToProcess;
        arguments->identificationThread = threadName.str();
        arguments->antigenName = ID_antigen;
        arguments->savingLocation = commonSavingLocation;

        if(nThreads > 1){
            int err = pthread_create(&(tid[i]), nullptr, &oneThreadJob, (void *) arguments);
            if (err != 0) cerr << "\ncan't create thread :" << strerror(err) << endl;
        } else {
            // in this case we don't need to start a thread, just call the function
            oneThreadJob((void*) arguments);
        }
    }

    if(nThreads > 1){
        cerr << myID << ", waiting for all threads to complete" << endl;
        for(size_t i = 0; i < static_cast<size_t>(nThreads); ++i){
            pthread_join(tid[i], nullptr);
            cout << myID << ", thread " << i << "/" << nThreads << " has finished " << endl;
        }
    }


    stringstream fName; fName << prefix << AntigenName << "FinalBindings_Process_" << rankProcess+1 << "_Of_" << nJobs << ".txt";
    cout << myID << "Saving pooled results into " << fName.str() << endl;
    commonSavingLocation->write(fName.str());



    // This mutex control the access to shared memory, can be cleared now
    pthread_mutex_destroy(&lockAccessPrecompAffinities);
    pthread_mutex_destroy(&lockSaveCommonDataset);

    cout << myID << ", =============== process finished! ==============" << endl;


    #ifdef USE_MPI
// This was a code to force MPI processes to wait for each-other using MPI barrier. Apparently, it is not necessary.
//    double centroid[3];/*ignore this array*/
//    if (rankProcess != 0) {
//        MPI_Recv(&centroid, 3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//        // sleep(1); /*This represent many calculations that will happen here later, instead of sleep*/
//        cout << myID << "=> Will now wait for other process to complete " << endl;
//        int res = MPI_Barrier(MPI_COMM_WORLD);
//        if (res != MPI_SUCCESS) cout << myID << ", problem with MPI_Barrier" << endl;
//    } else {
//        for (int i=0; i<nJobs-1; i++)        {
//            MPI_Send(&centroid, 3, MPI_DOUBLE, i+1, 0, MPI_COMM_WORLD);
//        }
//        int res = MPI_Barrier(MPI_COMM_WORLD);
//        if (res != MPI_SUCCESS) cout << myID << ", problem with MPI_Barrier" << endl;
//        cout << myID << " => All MPI processes have been completed! Main prosses will take over " << endl;
//    }

#ifdef WIN32
    Sleep(2);
    if(commonSavingLocation) delete commonSavingLocation;
    Sleep(2);
#else
    sleep(2);
    if(commonSavingLocation) delete commonSavingLocation;
    sleep(2);
#endif
    #endif
}


// ================================= Option 4: Single Binding (does all precalculations) ====================================

void option4(string ID_antigen, string CDR3){
    int receptorSize = DefaultReceptorSizeBonds; //10;
    int minInteract = DefaultContactPoints; //11;

    // Loading antigen
    string AntigenName = "";
    if(ID_antigen.size() < 4){
        AntigenName = IDshortcut(std::stoi(ID_antigen));
    } else {
        AntigenName = IDshortcut(ID_antigen);
    }

    cout << "   ... loading antigen " << AntigenName << " from the library" << endl;
    std::pair<superProtein*, vector<int> > AG = getAntigen(AntigenName);

    // Loading or calculating receptor structures around it
    affinityOneLigand T3 = affinityOneLigand(AG.first, receptorSize, minInteract, -1, 1, AG.second);
    // This shortcuts a lots of documenting/debugging calculations, and only calculates best energy
    T3.setUltraFast(true);

    // Getting sliding windows of it
    vector<string> cutSlides = slides(CDR3, receptorSize+1);
    size_t nSl = cutSlides.size();

    // For each sliding windows
    for(size_t i = 0; i < nSl; ++i){

        string subSeq = cutSlides[i];
        vector<string> optimalStruct;
        std::pair<double, double> affs = T3.affinity(subSeq, false, &optimalStruct);
        size_t nS = optimalStruct.size();

        // only says what is happening for the top 100 sequences.
        cout << CDR3 << "\t" << subSeq << "\t" << affs.first << "\t" << nS;

        // For each optimal binding structure
        for(size_t j = 0; j < nS; ++j){
            cout << "\t" << optimalStruct[j];
        }

        cout << endl;
    }
}

// very similar than singleBinding, except we will take the best binding and show paratope and epitope
void option7b(string ID_antigen, string CDR3){
    int receptorSize = DefaultReceptorSizeBonds; //10;
    int minInteract = DefaultContactPoints; //11;

    // Loading antigen
    string AntigenName = "";
    if(ID_antigen.size() < 4){
        AntigenName = IDshortcut(std::stoi(ID_antigen));
    } else {
        AntigenName = IDshortcut(ID_antigen);
    }

    cout << "   ... loading antigen " << AntigenName << " from the library" << endl;
    std::pair<superProtein*, vector<int> > AG = getAntigen(AntigenName);

    // Loading or calculating receptor structures around it
    affinityOneLigand T3 = affinityOneLigand(AG.first, receptorSize, minInteract, -1, 1, AG.second);
    // This shortcuts a lots of documenting/debugging calculations, and only calculates best energy
    T3.setUltraFast(true);

    // Getting sliding windows of it
    vector<string> cutSlides = slides(CDR3, receptorSize+1);
    size_t nSl = cutSlides.size();

    #ifdef ALLOW_GRAPHICS
    glDisplay();
    addToDisplay(AG.first, true);

    addToDisplay(new set<int>(AG.second.begin(), AG.second.end()));
    #endif

    // For each sliding windows
    for(size_t i = 0; i < nSl; ++i){

        string subSeq = cutSlides[i];
        vector<string> optimalStruct;
        std::pair<double, double> affs = T3.affinity(subSeq, false, &optimalStruct);
        size_t nS = optimalStruct.size();

        // only says what is happening for the top 100 sequences.
        cout << CDR3 << "\t" << subSeq << "\t" << affs.first << "\t" << nS;

        // For each optimal binding structure
        for(size_t j = 0; j < nS; ++j){
            cout << "\t" << optimalStruct[j];
            std::pair<int,string> parsed = retrieveStructureFromPosAndStructure( optimalStruct[j], '-');
            struct3D s = struct3D(parsed.second, UnDefined, parsed.first);
            superProtein* SasSuperProt = new superProtein(s);
            #ifdef ALLOW_GRAPHICS
            addToDisplay(SasSuperProt, false);
            #endif
            showParatopeEpitope(AG.first, SasSuperProt);
        }


        cout << endl;
    }
    #ifdef ALLOW_GRAPHICS
    glutMainLoop();
    #endif
}
// ==================== Option 5: Extract features: raw binding dataset => analyzed binding dataset

//possible format 1 (simplified, only one interaction code per sequence - do not include headers):\n";
//11AAseq,     Energy, interactionCode,                             nOptStr, (position, structure)*\n";
//RSRISRCHDMV  -84.15  iGjIkCeLgLdYkQhLjLkYjCbFdFkVbLcYaVkYjTiA     1        120675 UUDUSLSSLS\n";
//VSTAGIIGFTE  -82.8   cGbIaCaQbLaYbCdCfCgNaViRgFiFkFfFaYbTfTcAeA   2        129058 BSDRDLSRSR   129058 BSDRDLSRUR\n";

//possible format 2: (complete: slices of a CDR3, and all possible equally optimal interactions - do not include headers) \n";
//IDCDR3  CDR3(>=11AAs)       11AAseq,     Energy, nOptBind,  (pos,     structure,  interCode)*\n";
//36      CARDYYGSSYYFDYW	  RDYYGSSYYFD  -77.18  1          141157	BDDUSRSUUS	bKaDaIdSdTiWiVhDhViQjVkAjKcVcTeTkTjVeLgLhNiYadfkgj\n";
//36      CARDYYGSSYYFDYW	  DYYGSSYYFDY  -70.95  2          132896	LSLLDLDUUL	cViVkVcDhVdQfQhQhQiViTjKkPaEkEaQbYcRkViVbebgbkchgj 	132896	LDLLULUDLD  gViVkVgDhVdQfQhQhQiViTjKkPaEkEaQbYgRkViVbebgbkchcj\n";
//37      CARGETTGSSYDPYFDVW  CARGETTGSSY  -69.99  1          128996	RLLRLUDUUS	kLkRkPaWaVdDfDdVaQbVbKbVcLeQgQkQdNgWjLiDhGaYiGadgj\n";

void option5(string ID_antigen, string bindingDatasetFile, string outputFeaturesFile, int minDegree, bool includeDegreeInMotifs){
    //int degree = 1;
    //bool includeDegreeInMotifs = true;

    string AntigenName = "";
    if(ID_antigen.size() < 4){
        AntigenName = IDshortcut(std::stoi(ID_antigen));
    } else {
        AntigenName = IDshortcut(ID_antigen);
    }

    cout << "   ... loading antigen " << AntigenName << " from the library" << endl;
    std::pair<superProtein*, vector<int> > AG = getAntigen(AntigenName);
    {
        dataset<binding> data = dataset<binding>(bindingDatasetFile);
        dataset<analyzedBinding> res = dataset<analyzedBinding>();
        res.setNameAntigen(data.nameAntigen);


//        string originalFileName;
//        vector<string> listIDs;
//        vector<string> sequences;
//        vector<bool> bestForThisCDR3;       // just to make it easy to extract the best binding per CDR3
//        vector<T*> storage;
//        size_t nLines();
//        and binding is:
//        string AAseq;
//        double bindEnergy;
//        string structureID;
        features f(AntigenName, minDegree, includeDegreeInMotifs);
        size_t N = data.nLines();
        for(size_t i = 0; i < N; ++i){
            vector<string> vres = f.getProperties(*(data.storage[i]), data.sequences[i]);
            analyzedBinding* toAdd = new analyzedBinding(data.storage[i]->AAseq, data.storage[i]->bindEnergy, data.storage[i]->structureID, vres[interCodeWithIDpos], stringToSet(vres[positionsBound]), vres[hotspot_ID], vres);
            res.addLine(data.listIDs[i], data.sequences[i], toAdd, data.bestForThisCDR3[i]);
            // res << ID << sep << largeStr << sep << AAseq << sep << E << sep << interCode << sep << vres[interCodeWithIDpos] << sep << startPos << sep << structure; //<< endl;
            //res << "\t" << vres[seqAGEpitope] << "\t" << vres[seqABParatope] << "\t" << vres[motifAGEpitope] << "\t" << vres[motifABParatope] << "\t" << vres[motifsSizeGapsLigand] << "\t" << vres[motifsSizeGapsRec] << "\t" << vres[motifsChemicalLig] << "\t" << vres[motifsChemicalRec] << "\t" << vres[agregatesAGEpitope] << "\t" << vres[agregatesABParatope] << "\t" << vres[chemicalAGEpitope] << "\t" << vres[chemicalABParatope] << "\t" << vres[positionsBound] << endl;
            //cout << vres[interCodeWithIDpos] << endl;
            //cout << AAseq << sep << E << sep << concat(getProperties(AAseq, startPos, structure), sep) << endl;
        }
        res.write(outputFeaturesFile);
    }
    // note: the destructor of dataset will delete the analyzedBindings from memory
}









// In case we aim to show the antibody together with the bindings or hotspots, we first get the antibody
// into the correct coordinates, display it before option 7.
// If you wish to use your own PDB file, please also provide the discretized (difficult to rediscretize from a custom PDB by calling the PDB class.)
vector<vector<double> >* optionPre7(string ID_antigen){ //, string originalPDBonComputer = "", string discretizedPDBonComputer = ""){

    vector<string> cutName = split(ID_antigen, "_");
    if(cutName.size() != 2) {cerr << "ERR: optionPre7, could not cut the antigen ID (" << ID_antigen << " into antigen and chains." << endl; return nullptr;}
    cout << "Antigen PDB is " << cutName[0] << endl;
    cout << "Chains are " << cutName[1] << endl;

    // First, needs to rediscretize the original antigen to get the good coordinates,
#ifndef NOQT
    char* argv[] = {(char*) "Nothing!"};
    int argc = 1;
    QApplication appl(argc, argv);
    PDB* a2 = new PDB(cutName[0], cutName[1], false, 5.25, "FuC");

    // some antigens need nKeep=40
    if((!ID_antigen.compare("3WD5_A")) || (!ID_antigen.compare("4PP1_A")) || (!ID_antigen.compare("5JW4_A"))){
        a2->setnKeep(50);
    }
    a2->show();


    cerr << "First step OK" << endl;
    // note: this function will make sure that the PDB is downloaded or available in the running folder.

    string discretizedPDB = a2->getDiscretizedFileName();
    cout << "Discretized latfit output is in " << discretizedPDB << endl;
    //superProtein* discr = a2->inLattice;

    // test 2:
    // read a PDB from latfit, then read the same PDB but with another chain, include it in visualization
    latFitToLattice a = latFitToLattice();
    a.parseLatFitPDB(discretizedPDB);
    a.transform();
    superProtein* P3 = a.asSuperProtein();

    antigenInfo AI = getAntigenInfos(ID_antigen);
    string antibodyChains = AI.antibodyChains;
    cout << "The chains of binding antibodies in the original PDB are: " << antibodyChains << endl;

    vector<vector<double> > positionsAB = getPDBChainCoarseGrainedPositions(cutName[0] + ".pdb", antibodyChains, "FuC");
    vector<vector<double> >* transformed = new vector<vector<double> >(pooledPDBtoLattice(positionsAB, a.initXAxis, a.initYAxis, a.listPositions[0]));


    return transformed;
#else
    return nullptr;
#endif

//#ifdef ALLOW_GRAPHICS
//    glDisplay(); // this clears everything,
//    //addToDisplay(P3, false);
//    addToDisplay(transformed);
//    //glutMainLoop();
//#endif


}




//void option7(string ID_antigen, bool generateHotspots = false, string bindingDatasetFile = "", int maxStructures = 1e5){
void option7(string ID_antigen, bool generateHotspots, string bindingDatasetFile, double maxE, int sizeK, int degree){
    bool includeAntibodyFromPDB = defaultIncludePDBantibody;

    #ifndef NOQT
    setPictureFileNamePrefix(string("Hot") + ID_antigen);
    glDisplay(); // init glut
    #endif

    vector<vector<double> >* transformed = nullptr;
    #ifndef NOQT
    if(includeAntibodyFromPDB){
        transformed = optionPre7(ID_antigen);
    }
    #else
    if(includeAntibodyFromPDB) cerr << "WRN: for now, using option 7 (hotspots) with showing the antibody requires to use QT. Sorry. comment #define NO_QT" << endl;
    #endif

    cout << "Visualization of " << ID_antigen << "\t" << ((generateHotspots) ? "+Hotspots" : "") << "with maxE=" << maxE << ", K=" << sizeK << ", degree=" << degree << endl;
    double minE = -1e+6; // takes all
    //double maxE = 1e+6;

    string AntigenName = "";
    if(ID_antigen.size() < 4){
        AntigenName = IDshortcut(std::stoi(ID_antigen));
    } else {
        AntigenName = IDshortcut(ID_antigen);
    }

    cout << "   ... loading antigen " << AntigenName << " from the library" << endl;
    std::pair<superProtein*, vector<int> > AG = getAntigen(AntigenName);

    dataset<analyzedBinding> annotatedDataset;
    annotatedDataset.setNameAntigen(AntigenName);

    if(bindingDatasetFile.size() > 0){
        vector<struct3D*> goodFoldings = filterStructuresFromAbsolutRawData(bindingDatasetFile, minE, maxE);

        cout << "Got " << goodFoldings.size() << " structures for E within [" << minE << "," << maxE << "]" << endl;
        std::map<string, int> list = groupStructuresInClasses(goodFoldings);
        cout << "Got the following number of each non-redundant structures in space:" << endl;
        for(std::map<string, int>::iterator it = list.begin(); it != list.end(); ++it){

            std::pair<string, int> str = retrieveStructureFromID(it->first);
            struct3D* representativeStructClass = new struct3D(str.first, UnDefined, str.second);
            std::pair<string, int> infoOtherWay = oppositeEqualStructure(*representativeStructClass);

            cout << it->first << "\t" << it->second << "\t<=>\t" << infoOtherWay.first << "\t" << infoOtherWay.second << endl;

            if(generateHotspots){

                superProtein s2(*representativeStructClass);
                vector<string> analyzedFeatures = structuralFeatures(*(AG.first), s2, degree);
                set<int> IDresonAntigen = stringToSet(analyzedFeatures[positionsBound]);

                // Here, we cheat, the energy is replaced by number of sequences with this structure
                string writtenStructureID = getWrittenUniqueIDStructure(*representativeStructClass);
                analyzedBinding* bd = new analyzedBinding(writtenStructureID, it->second, writtenStructureID, "noInterCodeYet", IDresonAntigen, string("UnknownYet"), vector<string>());
                cerr << analyzedFeatures[positionsBound] << endl;

                // will use the structure name as ID
                annotatedDataset.addLine(writtenStructureID, writtenStructureID, bd, true);
            }
            #ifdef ALLOW_GRAPHICS
            //cout << "\t" << representativeStructClass->occupiedPositions.size() << "\t" <<
            addToDisplay(representativeStructClass, false);
            #else
            delete representativeStructClass;
            #endif
        }

//#ifndef NOQT
//        for(size_t i = 0; i < goodFoldings.size(); ++i){
//            addToDisplay(goodFoldings[i], false);
//        }
//#endif
    }

    // puts the protein
    #ifdef ALLOW_GRAPHICS
    addToDisplay(AG.first, true);
    if(transformed == nullptr){
        set<int>* s = new set<int>(AG.second.begin(), AG.second.end());
        addToDisplay(s);
    } else {
        addToDisplay(transformed);
    }

    if(generateHotspots){
        showBindingHotspots(annotatedDataset, AntigenName, sizeK);
    } //else {
      //  displayLigand(AG.first, AG.second, true, false); // always visible, yes do loop and do not show forbidden (that would remove the antibody if optionPre7)
    //}
    glutMainLoop();
    #endif
}

string generateBatchFLORIDA(string AntigenName, string fileCDR3s, bool test){
    stringstream res;

    res << ""
    "#!/bin/bash\n"
    "#SBATCH --job-name=phi_" << AntigenName << "\n"
    "#SBATCH --output=Res" << ((test)? "TestSlurmOutput" : "BatchOutput") << AntigenName << ".txt\n"
    "#SBATCH --error=Err" << ((test)? "TestSlurmOutput" : "BatchOutput") << AntigenName << ".txt\n"
    "#SBATCH --mail-type=BEGIN,END,FAIL\n"
    "#SBATCH --mail-user=pprobert@medisin.uio.no\n"
    "#SBATCH --time=" << ((test)? "00:05:0" : "72:00:0") << "\n"
    "\n"
    "#SBATCH --nodes=1\n"
    "#SBATCH --ntasks=1\n"
    "#SBATCH --cpus-per-task=32\n"
    "#SBATCH --mem=40000   #supposedly per node\n"
    "\n"
    "#SBATCH --account=brusko\n"
    "#SBATCH --qos=brusko-b\n"
    "\n"
    "./2020-07-23/AbsolutNoLib repertoire " << AntigenName << " " << fileCDR3s << " 16 /orange/brusko/philippe/\n";
    return res.str();
}

string generateBatchSAGA(string AntigenName, int nbNodes, string fileCDR3s, bool test){
    stringstream res;
    res << ""
    "#!/bin/bash\n"
    "#SBATCH --account=nn9603k\n"
    "#SBATCH --job-name=phi_" << AntigenName << "\n"
    "#SBATCH --output=Res" << ((test)? "TestSlurmOutput" : "BatchOutput") << AntigenName << ".txt\n"
    "#SBATCH --error=Err" << ((test)? "TestSlurmOutput" : "BatchOutput") << AntigenName << ".txt\n"
    "#SBATCH --time=" << ((test)? "00:05:0" : "120:00:0") << "\n"
    "#SBATCH --nodes=" << nbNodes << "\n"
    "#SBATCH --ntasks-per-node=2\n"
    "#SBATCH --cpus-per-task=16\n"
    "#SBATCH --mem=40000   #supposedly per node\n"
    << ((test)? "#SBATCH --qos=devel" : "") <<
    "\n"
    "set -o errexit # Make bash exit on any error\n"
    "set -o nounset # Treat unset variables as errors\n"
    "\n"
    "## Software modules\n"
    "module restore system\n"
    "module load intel/2018a\n"
    "\n"
    "mpiexec -n " << 2 * nbNodes << " ./2020-07-23/AbominationMPI repertoire " << AntigenName << " " << fileCDR3s << " 16 /cluster/work/users/pprobert/\n";
    return res.str();
}



string generateBatchFRAM(string AntigenName, int nbNodes, string fileCDR3s, bool test){

    stringstream res;
    res << ""
    "#!/bin/bash\n"
    "#SBATCH --account=nn9603k\n"
    "#SBATCH --job-name=phi_" << AntigenName << "\n"
    "#SBATCH --time=" << ((test)? "00:05:0" : "120:00:0") << "\n"
    "#SBATCH --nodes=" << nbNodes << "\n"
    "## Number of tasks to start on each node:\n"
    "##SBATCH --ntasks-per-node=2\n"
    "## Set OMP_NUM_THREADS\n"
    "#SBATCH --cpus-per-task=16\n"
    "\n"
    "## running ten replicates with jobarray and mpi scattering (121 param pairs)\n"
    "## Recommended safety settings:\n"
    "set -o errexit # Make bash exit on any error\n"
    "set -o nounset # Treat unset variables as errors\n"
    "\n"
    "## Software modules\n"
    "module restore system   # Restore loaded modules to the default\n"
    "module load intel/2017a\n"
    "\n"
    "## Make sure output is copied back after job finishes\n"
    "\n"
    "savefile *Binding*.txt\n"
    "\n"
    "mpiexec -n " << 2 * nbNodes << " ./09-April/AbominationMPI repertoire " << AntigenName << " " << fileCDR3s << " 16 &> "<< ((test)? "TestSlurmOutput" : "BatchOutput") << AntigenName << ".txt\n";
    return res.str();
}

// generate batch files
void option8(string fileCDR3s, int nbNodes, string shortTestFile, string architecture){
    enum clusters {FRAM, SAGA, FLORIDA, UNKNOWN};
    clusters foundCluster = UNKNOWN;
    if(!architecture.compare("FRAM"))foundCluster = FRAM;
    if(!architecture.compare("SAGA"))foundCluster = SAGA;
    if(!architecture.compare("FLORIDA"))foundCluster = FLORIDA;
    if(foundCluster == UNKNOWN){
        cerr << "ERR: Option Batch (option 8), the architecture " << architecture << " is unknown. Please use FRAM, SAGA or FLORIDA" << endl;
        return;
    }

    vector<string> AGs = listIDs();

    ofstream listBatchs("allBatches.sh");
    for(size_t i = 0; i < AGs.size(); ++i){

        string AntigenName = "";
        if(AGs[i].size() < 4){
            AntigenName = IDshortcut(std::stoi(AGs[i]));
        } else {
            AntigenName = IDshortcut(AGs[i]);
        }

        stringstream autoFileName;
        autoFileName << "Batch" << AntigenName << "For" << fileCDR3s << ".sh";
        ofstream f1(autoFileName.str());
        if(foundCluster == FRAM) f1 << generateBatchFRAM(AntigenName, nbNodes, fileCDR3s, false); // normal
        if(foundCluster == SAGA) f1 << generateBatchSAGA(AntigenName, nbNodes, fileCDR3s, false); // normal
        if(foundCluster == FLORIDA) f1 << generateBatchFLORIDA(AntigenName, fileCDR3s, false); // normal
        f1.close();

        stringstream autoFileNameTest;
        autoFileNameTest << "SlurmTest" << AntigenName << "For" << shortTestFile << ".sh";
        ofstream f2(autoFileNameTest.str());
        if(foundCluster == FRAM) f2 << generateBatchFRAM(AntigenName, nbNodes, shortTestFile, true); // testFile
        if(foundCluster == SAGA) f2 << generateBatchSAGA(AntigenName, nbNodes, shortTestFile, true); // testFile
        if(foundCluster == FLORIDA) f2 << generateBatchFLORIDA(AntigenName, shortTestFile, true); // testFile
        f2.close();

//        stringstream autoFileNameDirect;
//        autoFileNameDirect << "DirectTest" << AntigenName << "For" << shortTestFile << ".sh";
//        ofstream f3(autoFileNameDirect.str());
//        f3 << "module load intel/2017a\nmpiexec -n " << 4 << " ./09-April/AbominationMPI repertoire " << AGs[i] << " " << shortTestFile << " 10\n";
//        f3.close();

        listBatchs << "sbatch " << autoFileName.str() << "\n";

        stringstream getStructFile;
        std::pair<superProtein*, vector<int> > AG = getAntigen(AntigenName);
        getStructFile << "Get" << AntigenName << "_Structures.sh";
        ofstream f4(getStructFile.str());
        f4 << "wget philippe-robert.com/Absolut/Structures/" << fnameStructures(AG.first, 10, 11, AG.second) << "\n";
        f4.close();
    }
    listBatchs.close();
}

//void option9(string ID_antigen, string outputFolder = "", string folderWithBindings = "", string folderWithStructures = ""){
void option9(string ID_antigen, string outputFolder, string folderWithBindings, string folderWithStructures){
    //res << " " << programName << "html ID_antigen [outputFolfer] [folderWithBindings]\n";

    string AntigenName = "";
    if(ID_antigen.size() < 4){
        AntigenName = IDshortcut(std::stoi(ID_antigen));
    } else {
        AntigenName = IDshortcut(ID_antigen);
    }

    cout << "   ... loading antigen " << AntigenName << " from the library" << endl;
    std::pair<superProtein*, vector<int> > AG = getAntigen(AntigenName);

    string htmlPage = outputFolder + "/" + ID_antigen + "_infos.html";
    ofstream hout(htmlPage.c_str());
    if(!hout){
        cerr << "ERR: failing to create the html file in folder " << outputFolder << endl;
        return;
    }
    hout << html::getHtmlHeader();
    vector<string> AGs = listIDs();
    for(size_t i = 0; i < AGs.size(); ++i){
        if(!ID_antigen.compare(AGs[i])){
            html a(AGs[i]);
            hout << a.analyze(outputFolder, folderWithBindings, folderWithStructures);
        }
    }
    hout << html::getHtmlFooter();

    hout.close();
}

void info_fileNames(string ID_antigen){
    string AntigenName = "";
    if(ID_antigen.size() < 4){
        AntigenName = IDshortcut(std::stoi(ID_antigen));
    } else {
        AntigenName = IDshortcut(ID_antigen);
    }
    int receptorSize = DefaultReceptorSizeBonds; //10;
    int minInteract = DefaultContactPoints; //11;
    cout << "   ... loading antigen " << AntigenName << " from the library" << endl;
    std::pair<superProtein*, vector<int> > AG = getAntigen(AntigenName);
    cout << "Pre-calculated structures are in " << fnameStructures(AG.first, receptorSize, minInteract, AG.second) << endl;
    cout << "use:\nwget http://philippe-robert.com/Absolut/Structures/" << fnameStructures(AG.first, receptorSize, minInteract, AG.second)<< endl;
}
