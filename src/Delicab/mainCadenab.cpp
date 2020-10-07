
#include <iostream>
using namespace std;

// There are two .pro file. Zapotec.pro allows to use QT interface,
// but if you use ZapotecNoQt.pro, then 'NOQT' will be defined from the makefile and it is possible to run scripts on a computer without QT.
#include "Ymir/Ymir.h"
#ifndef NOQT
#include <QApplication>
#include "pdb.h"
#endif

int main(int argc, char** argv){

#ifndef NOQT
    QApplication appl(argc, argv);
    PDB* a2 = new PDB();

// You can also call the GUI with the PDB and chains of interests from the beginning:
//    PDB* a2 = new PDB("1CZ8", "VW", true);
//    PDB* a2 = new PDB("3I40", "AB", true);

    a2->show();
    appl.exec();
    return -1;
#else
    cerr << "You are intending to run the graphical interface of Zapotec. Please make sure to use Zapotec.pro and that NOQT is not defined." << endl;
    return -1;
#endif

//    {
//        string s2 = string("SURRDDLLURUSDSSLDSL");
//        int testsizeR = 6; //133152
//        affinityOneLigand a = affinityOneLigand(s2, "AGRHGRAHARGRAGRTWHSW" , -1, testsizeR, 7, -1, 1.0);
//        cout << a.affinity("AGNTWPL",false).first << endl;
//        // answer should be -29.36
//        //sanityCheckEasyRotate();
//        return 0;
//    }
//	  {
//    superProtein v0 = superProtein();
//    struct3D S01 = struct3D("SURRDRDDSDRLLDDRRUUL", UnDefined, 133152);
//    superProtein v1 = insert(v0, S01, 0);
//    struct3D S02 = struct3D("LUSDDURDRDDLLSURRUULDDUDSRLRR", UnDefined, 128867);
//    superProtein v2 = insert(v1, S02, 1000);
//    v2.setAAs(string("GIVEQCCTSICSLYQLENYCNFVNQHLCGSHLVEALYLVCGERGFFYTPKA"));
//    }



return 0;
}







