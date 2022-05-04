#ifndef PDB_H
#define PDB_H

#include <string>
#include "../Ymir/plot3d.h"
#include "../Ymir/ymir.h"
#include "discretize.h"
#include <cstdlib> // for system
#include <fstream>
#ifndef NOQT
#include <QThread>
#include <QtWidgets/QWidget>
#include <QMainWindow>
#include <QCoreApplication>
#include <QDir>
#endif
using namespace std;


// This uses wget to download a PDB file from the PDB database/webpage (if it is not already in the folder)
bool downloadPDB(string PDB_ID);

//returns the name of a prepared PDB file according to insertion and merging multiple chains (or empty if failed)
// note: the name is different if deletions or not, and if merging chains or not)
string prepareChainsIntoOneFile(string PDB_ID, string chains, bool deleteInsertions);


// small tool functions

// tells if a file exists in the current folder
inline bool exists(const std::string& name) ;

// Reads the full content of a text file and returns it as a (long) string
string readFile(string filename);



// A class to handle the discretizations without graphical interface (for command line use, and then no need of QT library)
struct discretization {

    discretization(){
        resolution = 0;
        initialize();
    }

    discretization(string _PDB_ID, string _chains, double _resolution, string _typeDiscrete) :
    PDB_ID(_PDB_ID), chains(_chains), resolution(_resolution), typeDiscrete(_typeDiscrete) {
        initialize();
    }

    string preparePDB(){
        // Parsing the inputs... This calss is made for the script language => Exit at error. When using PDB() GUI class, no exit(-1)
        if((resolution <= 0) || (resolution > 1000)){
            cerr << "ERR:PDB(), improper resolution: " << resolution << endl;
            return string("");
            //exit(-1);
        }
        if(typeDiscrete.compare("FuC") && typeDiscrete.compare("CA") && typeDiscrete.compare("CB")&& typeDiscrete.compare("CoM")){
            cerr << "ERR:PDB(), the only discretizations methods are : FuC for center of masses of full AAs, /n"
                    "CA or CB for carbon alphas/betas, and CoM for centroid center of mass (side-chain, not backbone).\n"
                    "You gave:" << typeDiscrete << endl;
            return string("");
            //exit(-1);
        }

        if(!downloadPDB(PDB_ID)) {
            return string("");
            //exit(-1);
        }

        string preparedFile = prepareChainsIntoOneFile(PDB_ID, chains, true);
        // if(preparedFile.size() == 0) exit(-1);

        cout << "Look for glycans" << endl;

        glycans = new readGlycans(PDB_ID + ".pdb", chains);
        cout << "Glycan search finished " << endl;

        return preparedFile;
    }

    string PDB_ID;
    string chains;
    double resolution;
    string typeDiscrete;
    readGlycans* glycans;

    // Refreshes the lastCRMSD, lastDRMSD and inLattice
    void initialize();

    // this function uses the information PDBfile, chains, resoltion, typeDiscrete, and outputs lastCMRSD and lastDRMSD
    // returns the name of the generated PDB and empty if failed
    string discretizeIntoFile(string preparedPDB, bool allow_jumps, int nKeep, bool silent);

    double lastCRMSD;
    double lastDRMSD;

    // This functions parses the discretized PDB (output from latfit) and transforms into a superProtein (lattice definition of this structure)
    string transformDiscretizedFileToLattice(string discretizedPDB);

    superProtein* inLattice;
    set<int> holes;
    set<int> positionGlycans;
};

void testDiscretization();

#ifndef NOQT
// A GUI class as a layer on top of the discretization class.
// The GUI contains all user options inside the window (fields ui->XXX), and uses the mother class
// discretization to perform all tasks.

namespace Ui {
class PDB;
}

class PDB : public QWidget, public discretization
{
    Q_OBJECT

public:
    string commandForRasmol;

    // First option: launches the GUI without any user options (blank)
    explicit PDB(QWidget *parent = nullptr);

    // Second option: launches the GUI and precalculates a discretization according to the user optpions (can be used by command line and
    // still the user can play with the GUI to visualize, for instance)
    PDB(string PDB_ID, string chains, bool pipeline = false, double resolution = 5.25, string typeDiscrete = "FuC", QWidget *parent = nullptr);


    void initialize();

    ~PDB();






    QStringList getChainsFromFasta();
    string getDiscretizedFileName();
    void setnKeep(int newValue);

public slots:
    void getPDB();
    void viewPDBin3D();
    void getFasta();

    void mergeChainsIntoFile();
    string discretizeIntoFile(bool silent = false);
    void readDiscretizedFileToLattice();

    void seePDB();
    void seeLatFitOutput();
    void calculateReceptors();
    string iterateBonds(string typeDiscr = "FuC", bool addHeader = false);
    void viewDiscretized(bool append = true);

private:
    Ui::PDB *ui;
};


/// impossible to have openGL on a thread
//class Worker : public QObject
//{
//    Q_OBJECT

//public slots:
//    void doWork(const QString &parameter) {
//        cout << "Launching" << endl;
//        QString result;
//        glDisplay();

//        struct3D* test = new struct3D(string("SSSSUUL"));
//        addToDisplay(test);

//        glutSwapBuffers();
//        //glutMainLoop();
//        cout << "Finished" << endl;
//        emit resultReady(result);
//    }

//signals:
//    void resultReady(const QString &result);
//};


//class Controller : public QObject
//{
//    Q_OBJECT
//    QThread workerThread;

//public:
//    Controller() {
//        Worker *worker = new Worker;
//        worker->moveToThread(&workerThread);
//        connect(&workerThread, &QThread::finished, worker, &QObject::deleteLater);
//        connect(this, &Controller::operate, worker, &Worker::doWork);
//        connect(worker, &Worker::resultReady, this, &Controller::handleResults);
//        workerThread.start();
//        worker->doWork(QString("test"));
//        cerr << "Free from slave" << endl;
//    }
//    ~Controller() {
//        workerThread.quit();
//        workerThread.wait();
//    }
//public slots:
//    void handleResults(const QString &){cout << "QThread done" << endl;}
//signals:
//    void operate(const QString &);
//};



//struct PDBmanips{
//    PDBmanips(){
//        Controller* test1 = new Controller();
//        cout << "test1 finished" << endl;
//    }

//    string getFasta(string PDB_ID_online);
//    void test2(){
//        Controller* test2 = new Controller();
//        cout << "test2 finished" << endl;
//    }



//};
#endif

#endif // PDB_H
