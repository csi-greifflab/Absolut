#ifdef hideTasks
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

// for sleep
#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif

#include <cstdlib>
using namespace std;


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

struct oneTask{
   oneTask(int a, string b, string c, string d, string e) : IDsim(a), antigenID(b), origin(c), sequence(d), toProcess(e) {}
   int IDsim;
   string antigenID;
   string origin;
   string sequence;
   string toProcess;
};

#define MAX_ALLOWED_THREADS 50

// Each process will decide which ID task to take

void processTasks(string fileName, int nThreads, string prefix, int startIDsim, int endIDsim){

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
    if(nThreads > MAX_ALLOWED_THREADS) {cerr << "ERR: option tasks, Does not allow more than " << MAX_ALLOWED_THREADS << " threads (requested: " << nThreads << "), will take 50" << endl; nThreads = MAX_ALLOWED_THREADS;}
    if(startIDsim > endIDsim){cerr << "ERR: option tasks, the ending IDsim " << endIDsim << " is before the starting IDsim " << startIDsim << " -> ending" << endl; return;}

    // 1 - Reading the list of sequences to process:
    cout << myID << "   ... loading task file " << fileName << "\n";
    ifstream f(fileName);
    if(!f){cerr << "ERR: processTasks, file " << fileName << " not found " << endl; return;}

    int IDsim = -1, mutations = 0;
    string antigenID, origin, hotspot, location, sequence, toProcess;
    int cpt = 0;

    vector<oneTask*> pileOfTasks;
    stringstream plannedTasks;
    while((IDsim > 0) && (cpt < 10000000)) {
        IDsim = -1;
        f >> IDsim >> antigenID >> origin >> mutations >> hotspot >> location >> sequence >> toProcess;
        cpt++;

        if((IDsim % nJobs) == 0){
            plannedTasks << myID << "\t" << IDsim << "\t" << antigenID << "\t" << toProcess << "\n";
            pileOfTasks.push_back(new oneTask(IDsim, antigenID, origin, sequence, toProcess));
        }
    }
    f.close();
    cout << "Process " << myID << " will process the following tasks (on " << nThreads << " threads per task)" << endl;
    cout << plannedTasks.str() << endl;



    //<< "\n       expected to have 2 or 3 columns: ID , CDR3 sequence , [optional tag]" << endl;
    size_t nT = pileOfTasks.size();
    cout << myID << "   ... Found " << nT << " tasks to perform " << endl;
    if(nT < 1) return;

    for(size_t i = 0; i < nT; ++i){
        oneTask* currentTask = pileOfTasks[i];
        if(!currentTask) return;

        // 2a - Loading the antigen from the library. Either call with a number, or with a substring of its name.
        string AntigenName = currentTask->origin;
        cout << myID << "   ... loading origin antigen " << AntigenName << " from the library " << endl;
        std::pair<superProtein*, vector<int> > AG = getAntigen(AntigenName);
        if(!AG.first) continue;

        cout << myID << "   ... loading mutant sequence for antigen " << currentTask->antigenID << endl;
        AG.first->setAAs(currentTask->sequence);

        cout << myID << "   ... charging the sequences to process with this antigen " << endl;
        repertoire rep = repertoire(currentTask->toProcess);
        if(rep.nLines() < 1) {
            cerr << "ERR: No sequences to process, trying next task" << endl;
            continue;
        }

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
                    cout << "     the 'tasks' option, which is made to treat lots of sequences in multithreads, and the " << endl;
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
                continue;
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
        int minLine = 0;
        int maxLine = static_cast<int>(rep.nLines())-1;
        int nToProcess = maxLine - minLine + 1;


        // from now on, identical code as before
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


    }


    // This mutex control the access to shared memory, can be cleared now
    pthread_mutex_destroy(&lockAccessPrecompAffinities);
    pthread_mutex_destroy(&lockSaveCommonDataset);

    cout << myID << ", =============== process finished! ==============" << endl;


    #ifdef USE_MPI
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

#endif
