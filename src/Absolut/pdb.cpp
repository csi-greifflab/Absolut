#include "pdb.h"
#include "ui_pdb.h"
#include "discretize.h"
#include "antigenLib.h" // for testing
#include "ymir.h" // for testing

#include <fstream>
#include "../latFit/latFit.h"
#include "selfEvo.h"

inline bool exists(const std::string& name) {
    ifstream f(name.c_str());
    if(f.good()){f.close(); return true;}
    return false;
}

string readFile(string filename){
    ifstream fr(filename.c_str());
    stringstream res;
    if(fr){
        res << fr.rdbuf();
        fr.close();
    } else return filename + string(", file not found!\n");
    return res.str();
}

// This uses wget to download a PDB file from the PDB database/webpage (if it is not already in the folder)
bool downloadPDB(string PDB_ID){
    if(PDB_ID.size() != 4){
        return false;
    }

    if(exists(PDB_ID + string(".pdb"))){
        cout << "File " << PDB_ID << ".pdb was found in the folder, no need to download!" << endl;
        return true;
    }

    //string command = string("wget --no-check-certificate \"http://www.pdb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=") + PDB_ID + string("\" -O ") + PDB_ID + string(".pdb");
    //string command = string(" curl -O -J \"http://files.rcsb.org/download/") + PDB_ID + string(".pdb\" -O ") + PDB_ID + string(".pdb");
    string command = string(" curl \"http://files.rcsb.org/download/") + PDB_ID + string(".pdb\" > ") + PDB_ID + ".pdb";
    
    cout << "Executing >> " << command << endl;
    system(command.c_str());
    if(exists(PDB_ID + string(".pdb"))){
        return true;
    }


    string lowerCase = PDB_ID;
    for(size_t i = 0; i < lowerCase.size(); ++i){
        lowerCase[i] = lowerCase[i] + 32; // the trick to get lowercase
    }
    string command2 = string(" wget --no-check-certificate \"http://files.rcsb.org/download/") + lowerCase + string(".pdb\" -O ") + lowerCase + string(".pdb");
    cout << "Executing >> " << command2 << endl;
    system(command2.c_str());

    cerr << "ERR: for whatever reason, could not download the file " << PDB_ID << ".pdb" << endl;
    return false;
}

string getChainList(string chains){
    string chainList;
    for(unsigned int i = 0; i < chains.size(); ++i){
        char X = chains[i];
        if(!(((X <= 'Z') && (X >= 'A')) || (((X <= 'z') && (X >= 'a'))))) {
            cerr << "ERR: getChainList(), chain character: '" << X << "' from '" << chains << "',only letters are allowed for chains, do not put space or any punctuation." << endl;
            return string("");}
        if(chains.size() == 0) chainList = string(1, X);
        else chainList += string(",") + string(1, X);
    }
    return chainList;
}

//returns the name of a prepared PDB file according to insertion and merging multiple chains (or empty if failed)
// note: the name is different if deletions or not, and if merging chains or not)
string prepareChainsIntoOneFile(string PDB_ID, string chains, bool deleteInsertions){
    if(chains.size() == 0) return string("");

    if(!exists(PDB_ID + ".pdb")){
        cerr << "File " << PDB_ID << ".pdb not found in the current folder" << endl;
        return string("");
    }

    // 1 - checks the name of the chains [a-zA-Z] and put them separated with comas, to be run by
    string chainList = getChainList(chains);   // Chains will be reformatted to be: A,B,E ...
    if(chainList.size() == 0) return string("");

    if(deleteInsertions){
        stringstream command;
        command << "python \"../pdb-tools/pdbtools/pdb_delinsertion.py\" " << PDB_ID << ".pdb > " << PDB_ID << "deIns.pdb";
        cout << "Executing >> " << command.str() << endl;
        system(command.str().c_str());
    }

    // this calls pdb_reres to combine all chains of interest into a new PDB file, where the new chain is named 'A'
    if(chains.size() > 1){
        stringstream command;
        command << "python \"../pdb-tools/pdbtools/pdb_selchain.py\" -" << chainList << " " << PDB_ID << ((deleteInsertions) ? "deIns" : "") << ".pdb | "
                << "python \"../pdb-tools/pdbtools/pdb_chain.py\" -A | "
                << "python \"../pdb-tools/pdbtools/pdb_reres.py\" -1 > " << PDB_ID << "_" << chains << "prepared.pdb";
        cout << "Executing >> " << command.str() << endl;
        system(command.str().c_str());

        if(exists(PDB_ID + "_" + chains + "prepared.pdb")){
            return PDB_ID + "_" + chains + "prepared.pdb";
        } else return string("");
    }

    if(deleteInsertions) {
        if(exists(PDB_ID + "deIns.pdb")){
            return PDB_ID + "deIns.pdb";
        } else return string("");
    }

    // If only one chain and no insertion, no need to do anything => Will process the original PDB file
    return PDB_ID + ".pdb";
}




void discretization::initialize(){
    lastCRMSD = NAN;
    lastDRMSD = NAN;
    inLattice = nullptr;
}


string discretization::discretizeIntoFile(string preparedPDB, bool allow_jumps, int nKeep, bool silent){

    // Just checks that the chain names are correct (probably already tested before, just a re-test)
    string chainList = getChainList(chains);
    if(chainList.size() == 0) return string("");

    // generates the name of the output file for latfit
    stringstream niceName;
    niceName << PDB_ID << "_" << chains << "discretized" << typeDiscrete << resolution << ".pdb";

    // generates the command line for calling latfit
    stringstream cmd;
    // if there were multiple chains, they were merged into a chain called 'A'
    if(chains.size() > 1){
        cmd << "latFit.exe -pdbFile=" <<  preparedPDB << " -pdbAtom=" << typeDiscrete << " -pdbChain=A "
             << ((allow_jumps) ? "-pdbChainGaps " : "") << "-lat=CUB -outMode=PDB -outFile=" << niceName.str() << " -opt=D -outOrigData -pdbAtomAlt=A -nKeep=" << nKeep << " -s"; // -s is for silent
    } else { // if not, then give the original ID of the chain in the non-merged PDB
        cmd << "latFit.exe -pdbFile=" <<  preparedPDB << " -pdbAtom=" << typeDiscrete << " -pdbChain=" << chains << " "
             << ((allow_jumps) ? "-pdbChainGaps " : "") << "-lat=CUB -outMode=PDB -outFile=" << niceName.str()  << " -opt=D -outOrigData -pdbAtomAlt=A -nKeep=" << nKeep << " -s";
    }

    // if no resolution is given, will use latfit default resolution.
    if(resolution > 0.01) cmd << " -bondLength=" << resolution;

    //if(!silent)
    cout << "Executing >> " << cmd.str();

    // intead of using system(), transforms the command into in arguments and directly calls the main function of latfit in C++ with the same arguments as in cmd
    // Trick taken from: https://stackoverflow.com/questions/1706551/parse-string-into-argv-argc

//   Bad way of doing it: myStringStream << "MyProgram arg1 arg2";
//    char* cmdline = (char*) myStringStream.str().c_str();
//    leads to segfault (seems to point to the stringstream content, that is going to be deleted, then cmdline becomes random content)
 // good way of doing it:
    enum { kMaxArgs = 64 };
    int argc = 0;
    char *argv[kMaxArgs];
    string copyCommand = cmd.str(); // if using char* cmdline = (char*) cmd.str().c_str();
    char* cmdline = (char*) copyCommand.c_str();
    char *p2 = strtok(cmdline, (const char *) " ");
    while (p2 && argc < kMaxArgs-1){
        argv[argc++] = p2;
        p2 = strtok(nullptr, " ");
    }
    argv[argc] = nullptr;

    // Now calling latfit!!
    std::pair<double, double> resFit = mainLatFit(argc, argv);
    //cout << argc << "\t" << argv[0] << "\n";
    lastCRMSD = resFit.first;
    lastDRMSD = resFit.second;

    if(!exists(niceName.str())){
        cerr << "ERR: Discretization failed with latfit." << endl;
        return string("");
    }

    if(!silent) cout << "-> Discretization successful!" << endl;
    if(!silent) cout << "   ... output saved in " << niceName.str() << endl;

    return niceName.str();

    // Note: converting cmd into arguments with this solution didnt work
    //    char *args[] = {
    //        (char*)(string("LatFIT.exe")).c_str(),
    //        (char*)(string("-pdbFile=\"") + ui->lineEditFileMergedChains->text().toStdString() + string("\"")).c_str(),
    //        (char*)(string("-pdbAtom=") + ui->lineEditAtom->text().toStdString()).c_str(),
    //        (char*)(string("-pdbChain=A")).c_str(),
    //        (char*)(string("-lat=CUB")).c_str(),
    //        (char*)(string("-outOrigData")).c_str(),
    //        (char*)(string("-opt=D")).c_str(),
    //        (char*)(string("-lat=CUB")).c_str(),
    //        (char*)(string("-outFile=") + PDB_ID + string("discretized.pdb")).c_str(),
    //        (char*)(string("-outMode=PDB")).c_str(),
    //        (char*)(string("")).c_str(),
    //        NULL
    //    };
    //    if(ui->checkBoxAllowJumps->isChecked()){
    //        nargs++;
    //        string jmp = "-pdbChainGaps";
    //        args[9] = (char*)(jmp).c_str();
    //    }
    //    string command = string("latFit.exe ");
    //    for(int i = 0; i < nargs; ++i){
    //        cout << args[i] << endl;
    //        command += string(" ") + string(args[i]);
    //    }
}

string discretization::transformDiscretizedFileToLattice(string discretizedPDB){

     if(inLattice) delete inLattice;
     holes.clear();

     // latFitToLattice is a class that parses PDBs, either from latfit, or to extract glycans
     latFitToLattice a = latFitToLattice();
     a.parseLatFitPDB(discretizedPDB);

     // this function rotates the coordinates to the first turn (Ox, Oy), and rescales the structure
     // by the lattice resolution (the distance between two next AAs).
     // Then, each AA is separated by a distance of 1, and each next AA position is converted into
     // a move => This will be the lattice description of the structure.
     a.transform();



     // => The results are stored into the class a, and need to be written out
     // writes in parallel into a file (fwr) and into a string (display), to be displayed
     string fileNameOut = PDB_ID + "_" + chains + string("inlattice.txt");
     ofstream fwr(fileNameOut.c_str());

     fwr << "Number of subchains:\n";
     fwr << a.structures.size() << "\n";

     fwr << "Structure of each subchain:\n";
     // this is now redundant with function 'asSuperProtein' from the LatfitToLattice PDB. I kept the code here because it was used to generate the database.
     // don't want to take the risk to change it
     // Creates a superProtein, and adds the chains one by one
     inLattice = nullptr;                                           //new superProtein("",-1);
     for(unsigned int i = 0; i < a.structures.size(); ++i){

         struct3D toAdd = struct3D(a.structures[i], UnDefined,  a.positions[i]);
         if(i == 0){ inLattice = new superProtein(toAdd);
            //cout << "Starting protein: " << endl;
            //cout << print(*inLattice) << endl;
         }
         else {inLattice = new superProtein(insert(*inLattice, toAdd, a.IDeachStartingResidue[i]));
             //cout << "Adding the structure: " << endl;
             //cout << print(toAdd) << endl;
             //cout << "Resulting protein: " << endl;
             //cout << print(*inLattice) << endl;
         }
         fwr << "\t" << a.positions[i] << "\t" << a.structures[i] << endl;
     }

     if(!inLattice){
         cerr << "ERR: Could not transform discretized PDB structure into the lattice (from file: " << discretizedPDB << ")" << endl;
         return string("");
     }

     inLattice->setAAs(a.sequence);


     // Now adding glycans!!
     readGlycans GL(PDB_ID + ".pdb", chains);
     cout << "Found " << GL.nGlycans << " Glycans in the original PDB. Will patch them now " << endl;
     set<int> alreadyTaken = getOccupiedPositions(inLattice);
     positionGlycans = a.addGlycans(GL, alreadyTaken);



     fwr << "AA sequence (concatenated, all subchains):\n";
     fwr << inLattice->getAAseq() << endl; // couls also use a.sequence, but just wanted to check the AAs are put correctly into the superProtein.

     // Find the list of donut holes or topologically impossible points that should not be available to the antibodies
     fwr << "=> List additional inaccessible positions\n";
     holes = listEmbarrasingPoints(inLattice, positionGlycans);
     for(set<int>::iterator it = holes.begin(); it != holes.end(); ++it){
        fwr << *it << "\t" << printVector(lattice::positionFromID(*it)) << "\n";
     }
     fwr << "/n";

     fwr << "=> Among which the glycans are:\n";
     for(set<int>::iterator it = positionGlycans.begin(); it != positionGlycans.end(); ++it){
        fwr << *it << "\t" << printVector(lattice::positionFromID(*it)) << "\n";
     }
     fwr << "/n";


     // 0- defines the code for this antigen
     stringstream codeAG;
     codeAG << PDB_ID << "_" << chains << "_" << typeDiscrete << "_" << resolution;

     // 1-, position blocked into a vector,
     stringstream posBlocked;
     posBlocked << "blockV = {";
     int i = 0;
     for(set<int>::iterator it = holes.begin(); it != holes.end(); ++it){
         if(i > 0) posBlocked << ", ";
         posBlocked << *it;
         i++;
     }
     posBlocked << "};" << endl;


     // second, the C++ code

     fwr << "\n\n"
     "C++ code to add this structure in the library:\n"
     "string agStruct; \n"
     "string agSeq; \n"
     "int startPos = -1; \n"
     "vector<int> blockV; \n"
     " \n"
     "// This part can be copied into AntigenLib.cpp  \n"
     "if(codeID == str2int(\"" << codeAG.str() << "\")){ \n"
            "\t agStruct = string(\"" << (inLattice->structure->sequence) << "\"); \n"
            "\t agSeq = string(\"" << inLattice->getAAseq()  << "\"); \n"
            "\t startPos = " << inLattice->structure->startingPosition << "; \n"
            "\t " << posBlocked.str() << "\n"
     "} \n"
     " \n"
     "// This part can be used to run in any C++ code (then remove the 'if(codeID == str2int())' line) \n"
     "superProtein* sp = new superProtein(agStruct, startPos); \n"
     "sp->setAAs(agSeq); \n"
     "std::pair<superProtein*, vector<int> > Antigen(sp, blockV); \n";

     fwr.close();

     return fileNameOut;
}


void testDiscretization(){
//    //discretization a("5TZ2","C",5.25, "FuC");
    discretization a("3ULU","A",5.25, "FuC");
//    string file1 = a.preparePDB();
//    // //string file1 = "5TZ2deIns.pdb";
//    a.initialize();
//    string file2 = a.discretizeIntoFile(file1,true,25,true);
//    cout << file2 << endl;
    string file2 = "3ULU_AdiscretizedFuC5.25.pdb";
    a.transformDiscretizedFileToLattice(file2);
    cout << "List of glycans: " << print(a.positionGlycans) << endl;



    vector< std::pair<int, double >> heatmapGlycans;
    for(set<int>::iterator it = a.positionGlycans.begin(); it != a.positionGlycans.end(); ++it){
        heatmapGlycans.push_back(std::pair<int, double>(*it, 0.5));
        cout << " add " << *it << endl;
    }
    addToDisplay(new vector< std::pair<int, double> >(heatmapGlycans));

    //std::pair< superProtein*, vector<int> > AG = getAntigen("5TZ2_C");
    vector<int> blocked = setToVector(union_sets(a.holes, a.positionGlycans));
    displayLigand(a.inLattice, blocked);
    glutMainLoop();
}









PDB::PDB(QWidget *parent) :
    QWidget(parent), discretization(),
    ui(new Ui::PDB)
{
    ui->setupUi(this);
    initialize();
}

// This one is doing the full pipeline: Need the PDB and the chains of interest.
PDB::PDB(string _PDBfile, string _chains, bool pipeline, double _resolution, string _typeDiscrete, QWidget *parent) :
    QWidget(parent), discretization(_PDBfile, _chains, _resolution, _typeDiscrete),
    ui(new Ui::PDB)
{
    ui->setupUi(this);

    commandForRasmol = string("\"rasmol\" ");
    // the initialize() from the mother class is already called.

    // Sets up the graphical windows and assign buttons to functions
    PDB::initialize();

    // Filter the input parameters and assign the values to the graphical interface
    if((_resolution <= 0) || (_resolution > 1000)){
        cerr << "ERR:PDB(), improper resolution: " << _resolution << endl;
        return;
    } else {
        ui->doubleSpinBoxLatticeDistance->setValue(_resolution);
    }
    if(_typeDiscrete.compare("FuC") && _typeDiscrete.compare("CA") && _typeDiscrete.compare("CB")&& _typeDiscrete.compare("CoM")){
        cerr << "ERR:PDB(), the only discretizations methods are : FuC for center of masses of full AAs, /n"
                "CA or CB for carbon alphas/betas, and CoM for centroid center of mass (side-chain, not backbone).\n"
                "You gave:" << _typeDiscrete << endl;
        return;
    } else {
        ui->lineEditAtom->setText(QString(_typeDiscrete.c_str()));
    }
    ui->lineEditPDB_ID->setText(QString(_PDBfile.c_str()));
    ui->lineEditChainsOfInterest->setText(QString(_chains.c_str()));

    // Processes all the steps if the arguments were given from beginning
    if(_PDBfile.size() > 0){
        getPDB();
    }
    this->show();
    if(chains.size() > 0){
        mergeChainsIntoFile();
    }
    if(pipeline) iterateBonds();
}

void PDB::initialize(){
    discretization::initialize();
    QDir dir(".");
    //cout << dir.absolutePath().toStdString() << endl;
    ui->lineEditWorkingDirectory->setText(dir.absolutePath());
    ui->checkBoxDeInsert->setChecked(true);

    ui->checkBoxAllowJumps->setChecked(true);
    ui->doubleSpinBoxAverageDistancePDB->setValue(3.9);
    ui->lineEditAtom->setText(QString("CA"));
    ui->spinBoxSizeReceptors->setValue(7);
    ui->spinBoxMinInteract->setValue(10);
    ui->spinBoxnKeep->setValue(25);
    QObject::connect(ui->lineEditPDB_ID, SIGNAL(returnPressed()), this, SLOT(getPDB()));
    QObject::connect(ui->pushButtonDownloadPDB, SIGNAL(released()), this, SLOT(getPDB()));
    QObject::connect(ui->pushButtonMergeChains, SIGNAL(released()), this, SLOT(mergeChainsIntoFile()));
    QObject::connect(ui->pushButtonDiscretize, SIGNAL(released()), this, SLOT(discretizeIntoFile()));
    QObject::connect(ui->pushButtonConvertLattice, SIGNAL(released()), this, SLOT(readDiscretizedFileToLattice()));
    QObject::connect(ui->pushButtonViewPDB3D, SIGNAL(released()), this, SLOT(seePDB()));
    QObject::connect(ui->pushButtonViewDiscretized, SIGNAL(released()), this, SLOT(viewDiscretized()));
    QObject::connect(ui->pushButtonView3DLatfitOutput, SIGNAL(released()), this, SLOT(seeLatFitOutput()));
    QObject::connect(ui->pushButtonViewReceptorStructures, SIGNAL(released()), this, SLOT(calculateReceptors()));
    QObject::connect(ui->pushButtonIterate, SIGNAL(released()), this, SLOT(iterateBonds()));
}

PDB::~PDB()
{
    delete ui;
}


string PDB::getDiscretizedFileName(){
    return ui->lineEditFileOutputLatfit->text().toStdString();
}

void PDB::setnKeep(int newValue){
    ui->spinBoxnKeep->setValue(newValue);
}

void PDB::getPDB(){
    PDB_ID = ui->lineEditPDB_ID->text().toStdString();
    ui->comboBoxListChains->clear();
    ui->textEditFasta->clear();
    ui->textEditDiscretizedPDB->clear();
    ui->textEditMultiChainStructure->clear();
    if(PDB_ID.size() != 4) return;
    if(downloadPDB(PDB_ID)){
        ui->lineEditPDB_filename->setText(QString((PDB_ID + string(".pdb")).c_str()));
    } else {
        cerr << "PDB file not found!" << endl;
    }
    getFasta();

    ui->comboBoxListChains->addItems(getChainsFromFasta());
}

void PDB::viewPDBin3D(){
    //system(PDBopenGL);
}

QStringList PDB::getChainsFromFasta(){
    QStringList res;
    PDB_ID = ui->lineEditPDB_ID->text().toStdString();
    if(PDB_ID.size() != 4) return res;

    ifstream f;
    f.open((PDB_ID + string(".pdb")).c_str());
    if(!f) {
        cerr << "ERR: PDB::getChainsFromFasta(), file not found: " << PDB_ID + string(".pdb") << endl; return res;
    }
    //cout << "reading " << PDB_ID + string(".pdb") << endl;

    string lineType;
    string currentChain;
    bool gotMolecule = false;
    int line = 0;
    while(f.good()){
        f >> lineType;
        //cout << lineType << endl;
        bool lineAlreadyRead = false;
        if(!lineType.compare(string("COMPND"))){
            line++;
            int lineID = 0;
            string typeText;
            string nameChain;
            char IDchain;
            if(line > 1) f >> lineID; // for whatever reasonm, the first COMPND line has no number
            f >> typeText;
            //cout << typeText << endl;
            if(!typeText.compare(string("MOLECULE:"))){
                if(gotMolecule) cerr << "ERR: PDB::getChains(), last molecule without chain character name" << endl;
                gotMolecule = true;
                getline(f, currentChain);
                lineAlreadyRead = true;
            }
            if(!typeText.compare(string("CHAIN:"))){
                if(!gotMolecule) cerr << "ERR: PDB::getChains(), got molecule without chain character name, for: " << nameChain << endl;
                f >> IDchain;
                QString fullChain = QString(IDchain) + QString(currentChain.c_str());
                //cout << "Read chain : " << fullChain.toStdString() << endl;
                res << fullChain;
                currentChain.clear();
                gotMolecule = false;
            }
        }
        if(!lineAlreadyRead){
            char buf[100001];
            f.getline(buf, 100000);
        }
    }
    f.close();
    return res;
}

void PDB::seePDB(){
    PDB_ID = ui->lineEditPDB_ID->text().toStdString();
    string command = commandForRasmol + PDB_ID + string(".pdb");
    cout << "Executing >> " << command;
    system(command.c_str());
}

void PDB::seeLatFitOutput(){
    PDB_ID = ui->lineEditPDB_ID->text().toStdString();
    stringstream command; command << commandForRasmol + PDB_ID << "discretized" << ui->lineEditAtom->text().toStdString() << ui->doubleSpinBoxLatticeDistance->value() << ".pdb" ;
    cout << "Executing >> " << command.str();
    system(command.str().c_str());
}

void PDB::getFasta(){
    PDB_ID = ui->lineEditPDB_ID->text().toStdString();
    if(!exists(PDB_ID + string(".fa")) || readFile(PDB_ID + string(".fa")).size() < 1){
        stringstream command;
        // command << "wget -q --no-check-certificate -O " << PDB_ID << ".fa \"https://www.rcsb.org/pdb/download/viewFastaFiles.do?structureIdList=" << PDB_ID << "&compressionType=uncompressed\" | tar xfv -;\n";
        // suggested by popucui, but still gets some error on my system: command << "wget -q --no-check-certificate -O " << PDB_ID << ".fa \"https://www.rcsb.org/fasta/entry/" << PDB_ID << "/download\"";
        command << "curl \"https://www.rcsb.org/fasta/entry/" << PDB_ID << "/download\" > " << PDB_ID << ".fa";

        cout << "Executing >> " << command.str() << endl;
        system(command.str().c_str());
    }
    if(exists(PDB_ID + string(".fa"))){
        ui->textEditFasta->setText(QString((readFile(PDB_ID + string(".fa"))).c_str()));
    } else {
        ui->textEditFasta->setText(QString("Obtaining FASTA online failed"));
    }
}

void PDB::mergeChainsIntoFile(){
    PDB_ID = ui->lineEditPDB_ID->text().toStdString();
    chains = ui->lineEditChainsOfInterest->text().toStdString();
    string preparedFile = prepareChainsIntoOneFile(PDB_ID, chains, ui->checkBoxDeInsert->isChecked());

    if(preparedFile.size() > 0){
        ui->lineEditFileMergedChains->setText(QString(preparedFile.c_str()));
    } else {
        ui->lineEditFileMergedChains->setText(QString("File couldn't be generated"));
    }

    if(true){
        //cout << "Look for glycans" << endl;
        // Note: pdb-tools destroys the glycans when removing insertions => Work on the original PDB
        string processedPDB = PDB_ID + ".pdb"; //ui->lineEditFileMergedChains->text().toStdString();

        //string parse = ui->lineEditChainsOfInterest->text().toStdString();
        readGlycans a(processedPDB, chains);
        //cout << "Glycan search finished " << endl;
    }

    discretizeIntoFile(true);
    readDiscretizedFileToLattice();
}



string PDB::discretizeIntoFile(bool silent){
    // fills the fields of discretization before calling its function discretize
    PDB_ID = ui->lineEditPDB_ID->text().toStdString();
    chains = ui->lineEditChainsOfInterest->text().toStdString();

    typeDiscrete = ui->lineEditAtom->text().toStdString();
    resolution = ui->doubleSpinBoxLatticeDistance->value();

    int nKeep = ui->spinBoxnKeep->value();
    string preparedFile = ui->lineEditFileMergedChains->text().toStdString();
    bool allow_jumps = ui->checkBoxAllowJumps->isChecked();

    string discretizedName = discretization::discretizeIntoFile(preparedFile, allow_jumps, nKeep, silent);

    if(discretizedName.size() > 0){
        ui->lineEditFileOutputLatfit->setText(QString(discretizedName.c_str()));
        ui->textEditDiscretizedPDB->setText(QString(readFile(discretizedName).c_str()));
        return discretizedName;
    } else {
        ui->lineEditFileOutputLatfit->setText("Discretization with LatFit failed");
        ui->textEditDiscretizedPDB->setText("File not found!\n");
    }
    return string("");
}

void PDB::readDiscretizedFileToLattice(){
    // fills the fields of discretization before calling its function discretize
    PDB_ID = ui->lineEditPDB_ID->text().toStdString();
    chains = ui->lineEditChainsOfInterest->text().toStdString();
    typeDiscrete = ui->lineEditAtom->text().toStdString();
    resolution = ui->doubleSpinBoxLatticeDistance->value();

    string discretizedPDB =  ui->lineEditFileOutputLatfit->text().toStdString();

    string fileInLattice = discretization::transformDiscretizedFileToLattice(discretizedPDB);

    if(fileInLattice.size() > 0){
        ui->lineEditFileConvertedLattice->setText(QString(fileInLattice.c_str()));
        ui->textEditMultiChainStructure->setText(QString(readFile(fileInLattice).c_str()));
    } else {
        ui->lineEditFileConvertedLattice->setText("Could not read/convert the discretized PDB from latfit, sorry");
        ui->textEditMultiChainStructure->setText("File not found");
    }
}

void loopOnce(){
    static bool started = false;
#ifdef ALLOW_GRAPHICS
    if(!started) glutMainLoop();
    started = true;
#endif
}

void PDB::viewDiscretized(bool append){
    //cerr << "0" << endl;
    readDiscretizedFileToLattice();

    #ifdef ALLOW_GRAPHICS

    //cerr << "1" << endl;
    glDisplay(append);
    //cerr << "2" << endl;
    addToDisplay(new superProtein(*inLattice), false); // wrong wrong, when inLattice changes, then the previous one mught be deleted, be careful
    //cerr << "3" << endl;
    set<int>* holes = new set<int>(listEmbarrasingPoints(inLattice));
    addToDisplay(holes);

    vector<std::pair<int,double>> heatmapGlycans;
    for(set<int>::iterator it = positionGlycans.begin(); it != positionGlycans.end(); ++it){
        heatmapGlycans.push_back(std::pair<int,double>(*it, 0.5));
    }
    addToDisplay(new struct3D("S"),false);
    addToDisplay(new vector< std::pair<int, double> >(heatmapGlycans));
    //cerr << "4" << endl;
    if(append) loopOnce();
    //cerr << "5" << endl;
    return;


    // Old code: this code is putting a different color per chain.
    string fileName = ui->lineEditFileOutputLatfit->text().toStdString();
    latFitToLattice a = latFitToLattice();
    a.parseLatFitPDB(fileName);
    a.transform();

//#ifdef ALLOW_GRAPHICS
     glDisplay();
     if(a.structures.size() != a.positions.size()) cerr << "ERR: incompatible positions and structure numbers" << endl;
     for(unsigned int i = 0; i < a.structures.size(); ++i){
         cout << a.positions[i] << "\t" << a.structures[i] << endl;
         struct3D* test = new struct3D(a.structures[i], UnDefined, a.positions[i]);
         addToDisplay(test, false);
         //addDissected(test, UnDefined, a.positions[i]);
     }
     addToDisplay(&(a.listIntPositions));
     addToDisplay(&(a.listRotatedOriginPositions));


     loopOnce();
//#endif
    #endif
}

void PDB::calculateReceptors(){
    int sizeReceptors = ui->spinBoxSizeReceptors->value();
    int minInteract = ui->spinBoxMinInteract->value();
    affinityOneLigand a = affinityOneLigand(inLattice, sizeReceptors, minInteract, -1, 1.0);
    //protein insert(protein& existingProt, struct3D& toAdd, int IDfirstResidue);
    for(int i = 0; i < 1; ++i){
        string rec = randomProt(sizeReceptors+1);
        a.affinity(rec, true);
    }
}


string PDB::iterateBonds(string typeDiscr, bool addHeader){
    stringstream fLog;
    if(addHeader) fLog << "PDB_ID\tchains\ttypeDiscr\tbondDistance\tcRMSD\tdRMSD\n";

    PDB_ID = ui->lineEditPDB_ID->text().toStdString();
    chains = ui->lineEditChainsOfInterest->text().toStdString();

    if((chains.size() > 0) && (PDB_ID.size() > 0)){
        ui->lineEditAtom->setText(QString(typeDiscr.c_str()));
        for(double d = 3.0; d < 7; d = d + 0.05){
            ui->doubleSpinBoxLatticeDistance->setValue(d);
            string createdFile = discretizeIntoFile(true);
            if(createdFile.size() > 0){
                readDiscretizedFileToLattice();
                viewDiscretized(true);
                cout << PDB_ID << "_" << chains << "\t" << typeDiscr << "\tbond\t" << d << "\tcRMSD\t" << lastCRMSD << "\tdRMSD\t" << lastDRMSD;
                fLog << PDB_ID << "\t" << chains << "\t" << typeDiscr << "\t" << d << "\t" << lastCRMSD << "\t" << lastDRMSD << "\n";
            } else {
                cout << PDB_ID << "_" << chains << "\t" << typeDiscr << "\tbond\t" << d << "\tcRMSD\tNaN\tdRMSD\tNAN";
                fLog << PDB_ID << "\t" << chains << "\t" << typeDiscr << "\t" << d << "\tNaN\tNaN\n";
            }
        }
        // loopOnce(); //this would trigger glutMain

        //    ofstream fLogWrite = ofstream("logAllDiscrete", ios::app);
        //    fLogWrite << fLog.str();
        //    fLogWrite.close();
    }

    return fLog.str();
}
