#include "html.h"

#include "../Ymir/ymir.h"
#include "../Tools/md5.h"

#ifndef NO_LIBS
#include "discretize.h"
#endif

#include "motifFeatures.h"
#include "quality.h"
#include "selfEvo.h"
#include "../Tools/zaprandom.h"
#include "../Tools/dirent.h"
#include <regex>
#include "antigenLib.h"
#include "importrepertoire.h"
#include "poolstructs.h"
#include "fileformats.h"

#include <set>

#ifndef NOQT
#include "pdb.h"
#include <QApplication>
#endif

#include "plot3d.h"

#include <iostream>
using namespace std;

string getPDB(string AGname){
    if(AGname.size() < 4) return string("");
    return AGname.substr(0,4);
}

string getChain(string AGname){
    if(AGname.size() < 6) return string("");
    string res = string(1, AGname[5]);
    for(size_t i = 6; i < AGname.size(); ++i){
        if(AGname[i] != '_') res = res + string(1, AGname[i]);
        else break;
    }
    return res;
}

html::html(std::string _AGname){
    AGname = _AGname;
}

string html::getHtmlHeader(){
    stringstream hout;
    hout << "<!DOCTYPE html>\n"
        "<html>\n"
        "<head>\n"
        "    <meta charset=\"utf-8\" />\n"
        "    <title>Philippe Robert's home page</title>\n"
        "    <link rel=\"stylesheet\" href=\"style3.css\" />\n"
        "</head>\n"
    "<body>\n"
    "<ul id=\"Menu\">\n"
    "    <li><a href=\"index.html\">Home</a></li>\n"
    "    <li><a href=\"#\">Docs</a><ul>\n"
    "                    <a href=\"#\">Docs</a>\n"
    "                    <li><a href=\"CV/A.pdf\">CV & Contact</a></li>\n"
    "                    <li><a href=\"#\">Teaching</a></li>\n"
    "<li><a href=\"Recherche/Modelisation/\">Research</a></li>\n"
    "    </ul></li>\n"
    "    <li><a href=\"#\">Programs</a><ul>	\n"
    "                    <li><a href=\"#\">Programs</a></li>\n"
    "                    <li><a href=\"Limonade/LimonadeV6.zip\">Limonade</a></li>\n"
    "                    <li><a href=\"cleozone.html\">Cleozone</a></li>\n"
    "    </ul></li>\n"
    "    <li><a href=\"#\">Links</a></li>\n"
    "</ul>\n";
    return hout.str();
}

string html::getHtmlFooter(){
    stringstream hout;
    hout << "</body>\n"
            "</html>\n";
    return hout.str();
}

// folderWithBindings will be where the bindings of THIS antigen are located.
string html::analyze(string outputFolder, string folderWithBindings, string folderWithStructures){


    // 0 - makes sure the folder name finished by '/'
    if((outputFolder.size() > 0) && (outputFolder.back() != '/')) outputFolder.push_back('/');
    if((folderWithBindings.size() > 0) && (folderWithBindings.back() != '/')) folderWithBindings.push_back('/');
    if((folderWithStructures.size() > 0) && (folderWithStructures.back() != '/')) folderWithStructures.push_back('/');

    cout << "Will generate HTML output for antigen " << AGname << ", output folder " << outputFolder << endl;
    if(folderWithBindings.size() > 0) cout << "   + analysis of bindings from folder " << folderWithBindings << endl;
    if(folderWithBindings.size() > 0) cout << "   + link to pre-computed structures from folder " << folderWithStructures << endl;

    // 1- check it is possible to write in the output folder.
    ofstream testFolder(outputFolder + "infos.txt", ios::out);
    if(!testFolder){
       cerr << "ERR: can not write in the output folder " << outputFolder << endl;
       return string("");
    }


    // 2- creates a html page. Will create files (pictures or links) and add them to the html page
    // Note: the page doesn't have header/footer, this is just content for the body.
    stringstream htmlout;

    string PDBname = getPDB(AGname);
    string chains = getChain(AGname);
    cout << "Generate analysis of antigen " << PDBname << ", chain(s) " << chains << endl;


    // 2a - checks the antigen exists and writes basic info.
    std::pair<superProtein*, vector<int> > FullAntigen = getAntigen(AGname);
    superProtein* AG = FullAntigen.first;
    set<int> AGholes = set<int>(FullAntigen.second.begin(), FullAntigen.second.end());
    //vector<int> startPos = startingPositions(AG);
    vector<std::pair< int, string> > subCh = getSubChains(AG);
    string originalAAseq = AG->getAAseq();
    if(!AG->structure) {cerr << "ERR: Empty antigen!" << endl; return string("");}
    string originalStructure = AG->structure->sequence;




    // 2b - re-discretizes the antigen and checks it returns the same positions,
    //      the forbidden positions should be included in the library antigen, where new forbidden positions might have been added manually.
    #ifndef NOQT
    // Step 1 - check that discretization is OK if would do it again.
    char* argv[] = {(char*) "a"};
    int motherFargc = 1;
    QApplication appl(motherFargc, argv);
    PDB* a2 = new PDB(PDBname, chains, false, 5.25, "FuC");

    superProtein* discr = a2->inLattice;
    set<int> holes = listEmbarrasingPoints(discr);
    vector<std::pair< int, string> > subChDiscr = getSubChains(discr);
    //vector<int> discrStartPos = startingPositions(discr);
    string AAs = discr->getAAseq();
    string structure = discr->structure->sequence;


    if(AAs.compare(originalAAseq)) {
        cerr << "ERR: the rediscretized antigen AA sequence is different: Registered in library is: \n"
                                        << originalAAseq << " while rediscretized is \n" << AAs << endl;
        return string("");
    }
    if(structure.compare(originalStructure)) {
        cerr << "ERR: the rediscretized antigen structure is different: Registered in library is: \n"
                                        << originalStructure << " while rediscretized is \n" << structure << endl;
        return string("");
    }

    if(subCh.size() != subChDiscr.size()) {
        cerr << "ERR: the rediscretized antigen structure has different number of chains" << endl;
        return string("");
    }

    htmlout << "<h1>" << AGname << "</h1>\n";


    a2->viewDiscretized(false);
    savePicture(outputFolder + AGname + string(".bmp"));

    htmlout << "<img src=\"" << AGname + string(".bmp") << "\" alt=\"AntigenPicture\" class=\"center\" style=\"filter:FlipH.\">\n";

    htmlout << "<h3>" << subChDiscr.size() << " Subchain(s): </h3>" << endl;
    htmlout << "<p>\n";
    for(size_t i = 0; i < subChDiscr.size(); ++i){
        if((subChDiscr[i].first) != (subCh[i].first)) {
            cerr << "ERR: Different starting positions between chains" << endl;
            return string("");
        }
        if(subChDiscr[i].second.compare(subCh[i].second)){
            cerr << "ERR: Different substructures between chains" << endl;
            return string("");
        }
        htmlout << subChDiscr[i].first << "\t" << subChDiscr[i].second << endl;
    }
    htmlout << AAs << "\n";
    htmlout << "</p>\n";

    // finding positions that were added,


    std::vector<int> common_points;
    set_intersection(AGholes.begin(),AGholes.end(),holes.begin(),holes.end(), std::back_inserter(common_points));
    if(common_points.size() != holes.size()){
        cerr << "ERR: the auto discretize has some positions that the library AG does not have" << endl;
    }
    cout << " ----------- Quality control passed successfully!  ----------- " << endl;


    htmlout << "<h3>Automatically identified aberrant positions, that will be be blocked:</h3>\n";
    htmlout << "<p>\n";
    for(size_t i = 0; i < common_points.size(); ++i){
         htmlout << common_points[i];
         if(i < common_points.size() - 1) htmlout << ", "; //"\t" << printVector(lattice::positionFromID(common_points[i])) << "\n";
    }
    htmlout << "\n";
    htmlout << "</p>\n";

    // these are the points that are added manyally
    std::vector<int> pointsOnlyInOne;
    set_symmetric_difference(AGholes.begin(),AGholes.end(),holes.begin(),holes.end(), std::back_inserter(pointsOnlyInOne));
    if(pointsOnlyInOne.size() > 0){
        htmlout << "<h3>Manually filled positions:</h3>\n";
        htmlout << "<p>\n";
        for(size_t i = 0; i < pointsOnlyInOne.size(); ++i){
             htmlout << pointsOnlyInOne[i]; // << "\t" << printVector(lattice::positionFromID(pointsOnlyInOne[i])) << "\n";
            if(i < pointsOnlyInOne.size() - 1) htmlout << ", ";
        }
        htmlout << "</p>\n";
        htmlout << "\n";
    }


    // 2c - Generates data for different modes of discretization
    ofstream latticeResolutions(outputFolder + AGname + "discretizations.txt", ios::out);
    latticeResolutions << a2->iterateBonds("FuC", true); // the true/false tells to add a header
    latticeResolutions << a2->iterateBonds("CoM", false);
    latticeResolutions << a2->iterateBonds("CA", false);
    latticeResolutions.close();

    // there should be a R script that later generates discretizations.png
    htmlout << "<p>\n";
    htmlout << "<a href=\"" << AGname << "discretizations.txt\">\n";
    htmlout << " <img src=\"" << AGname << "discretizations.png\" alt=\"Discretization quality for different lattice resolutions, for antigen " << AGname << "\" style=\"width:500px;height:600px;\"> \n";
    htmlout << "</a>\n";
    htmlout << "</p>\n";

    #endif

    // 2d - Now liks to available files in folderWithBindings and shows some statistics,

    // look for files:
    // - the structure files (are they available?)


    string fStruct = fnameStructures(FullAntigen.first, 10, 11, FullAntigen.second);
    htmlout << "<p>\n";
    htmlout << "The precomputed structure files for this antigen is named " << fStruct <<"\n";

    ifstream ftestStr(folderWithStructures + fStruct);
    if(ftestStr){
        htmlout << " and is available for download ";
        htmlout << "<a href=\"" << fStruct << "\"> here </a> \n";
    } else {
        htmlout << " and is not available for download (yet) - request it from us here <a href=\"probertmontp@gmail.com\">Link text</a>";
    }
    htmlout << "</p>\n";

    string fStructAll = fnameStructuresAndCompactForAASeqLigand(FullAntigen.first, 10, 11, FullAntigen.second);
    string fStructCompact = fileNameCompactForAASeqLigand(FullAntigen.first, 10, 11, FullAntigen.second);
    htmlout << "Additional files that will be generated automatically:\n" << fStructAll <<"\n" << fStructCompact << "\n";


    // - the bindings (whole data, SuperHeroes, Heroes, Mascotte, Loosers)
    // Will do epitope clustering for each?


    string bindingDatasetFile = folderWithBindings + AGname + "_HeroesforAbsolut.txt";
    vector<struct3D*> goodFoldings = filterStructuresFromAbsolutRawData(bindingDatasetFile, -1e6, +1e6);

    cout << "Got " << goodFoldings.size() << " structures for heroes of " << AGname << endl;
    std::map<string, int> list = groupStructuresInClasses(goodFoldings);
    cout << "Got the following number of each non-redundant structures in space:" << endl;

    for(std::map<string, int>::iterator it = list.begin(); it != list.end(); ++it){
        cout << it->first << "\t" << it->second << endl;
        #ifndef NOQT
        std::pair<string, int> str = retrieveStructureFromID(it->first);
        struct3D* representativeStructClass = new struct3D(str.first, UnDefined, str.second);
        //cout << "\t" << representativeStructClass->occupiedPositions.size() << "\t" <<
        addToDisplay(representativeStructClass, false);
        #endif
    }

    //#ifndef NOQT
    //        for(size_t i = 0; i < goodFoldings.size(); ++i){
    //            addToDisplay(goodFoldings[i], false);
    //        }
    //#endif

    displayLigand(FullAntigen.first, FullAntigen.second, true, true);







    return htmlout.str();
}
