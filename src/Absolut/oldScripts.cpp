
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




string listAvailableStructures(string fileWithList = ""){
    stringstream res;

    vector<string> listAGs = listIDs();
    std::map<string, string> IDtoFile;
    for(size_t i = 0; i < listAGs.size(); ++i){
        std::pair<superProtein*, vector<int>> AG = getAntigen(listAGs[i]);
        string fStruct = fnameStructures(AG.first, 10, 11, AG.second);
        IDtoFile[fStruct] = listAGs[i];
    }
    if(fileWithList.size() > 0){
        ifstream f(fileWithList);
        while(f.good()){
            string nextFile;
            f >> nextFile;
            if(IDtoFile.find(nextFile) != IDtoFile.end()){
                res << IDtoFile.find(nextFile)->second << "\t" << IDtoFile.find(nextFile)->first << "\tavailable\n";
                IDtoFile.erase(nextFile);
            }
        }
    }
    for(std::map<string,string>::iterator it = IDtoFile.begin(); it != IDtoFile.end(); ++it){
        res << it->second << "\t" << it->first << "\tneeds to be computed\n";
    }
    return res.str();
}

//template<class I, class O>
//function<O(I)> memoize(function<O(I)> f) {
//    // Keep copies in the map of parameters passed as references
//    map<typename std::decay<I>::type, O> memos;
//    return [=](I i) mutable -> O {
//        auto it = memos.lower_bound(i);
//        if (it == memos.end() || it->first != i)
//            it = memos.insert(it, make_pair(i, f(i)));
//        return it->second;
//    };
//};

template<class... I, class O>
function<O(I...)> memoize(function<O(I...)> f) {
    map<std::tuple<typename std::decay<I>::type...>, O> memos;
    return [=](I... i) mutable -> O {
        auto args(std::tie(i...));
        auto it(memos.lower_bound(args));
        if (it == memos.end() || it->first != args)
            it = memos.insert(it, make_pair(args, f(i...)));
        return it->second;
    };
};

string getRandomCDR3s(int sizeAAs, int number, string IDprefixTag = ""){
    if((sizeAAs < 2) || (sizeAAs > 1000)) return string("");
    if(number < 1) return string("");
    stringstream res;
    for(int i = 0; i < number; ++i){
        res << IDprefixTag << i+1 << "\t" << randomProt(sizeAAs) << "\n";
    }
    return res.str();
}

string getRandomSampling(string repertoireFile, int nbSeqToSample, string IDprefixTag = ""){
    stringstream res;
    repertoire a(repertoireFile);
    vector<pair<string,string> > lines = a.getRandomLines(nbSeqToSample);
    if(lines.size() < nbSeqToSample) cout << "WRN: getRandomSampling, you requested more sequences than the repertoire contains => you got less sequences than requested." << endl;
    for(size_t i = 0; i < lines.size(); ++i){
        res << IDprefixTag << lines[i].first << "\t" << lines[i].second << "\n";
    }
    return res.str();
}

void testRandomSeqGeneration(){
    cout << " A list of 100 random CDR3s of size 100" << endl;
    string list1 = getRandomCDR3s(15, 50000, "Ra");
    {ofstream f1("50000RandomCDR3sSize15.txt"); f1 << list1; f1.close();}
    cout << list1.substr(0, 10000) << endl;
    cout << " A list of 100 random CDR3s sampled from murine CDR3s" << endl;
    string list2 = getRandomSampling("C:/Users/pprobert/Desktop/Main/RepertoireVictor/UniqueCDR3s.txt", 50000, "Mu");
    {ofstream f2("50000MurineCDR3sSize15.txt"); f2 << list2; f2.close();}
    cout << list2.substr(0, 10000) << endl;

}
void testGlycans(){
    cerr << "Start" << endl;
    vector<string> AGs = listIDs();

#ifndef NO_LIBS
    for(size_t i = 0; i < min(size_t(500),AGs.size()); ++i){
        string PDB = getPDB(AGs[i]);
        string chain = getChain(AGs[i]);
        cout << "Glycans for " << AGs[i] << endl;
        readGlycans a(PDB+string(".pdb"),chain);
    }
    cerr << "End" << endl;
#endif
}


#include <sys/types.h>
vector<string> list_dir(const char *path) {
    vector<string> res;
    struct dirent *entry;
    DIR *dir = opendir(path);

    if (dir == nullptr) {
        return res;
    }
    while ((entry = readdir(dir)) != nullptr) {
        res.push_back(string(entry->d_name));
    }
    closedir(dir);
    return res;
}


vector<string> listFilesWithPattern(string directory, string filePatterns){

    vector<string> listFilesInFolder = list_dir(directory.c_str());
    vector<string> res;
    //std::string pattern("A|D");         // Regex expression
    std::regex rx(filePatterns);             // Getting the regex object

    for(size_t i = 0; i < listFilesInFolder.size(); ++i){
        string s = listFilesInFolder[i];
        std::ptrdiff_t number_of_matches = std::distance(  // Count the number of matches inside the iterator
            std::sregex_iterator(s.begin(), s.end(), rx),
            std::sregex_iterator());

        if(number_of_matches > 0){
            cout << "   ... Adding file to list: " << s << endl;
            res.push_back(s);


        }
    }
    return res;
}


void listdir(const char *name, int indent)
{
    DIR *dir;
    struct dirent *entry;

    if (!(dir = opendir(name)))
        return;

    while ((entry = readdir(dir)) != nullptr) {
        if (entry->d_type == DT_DIR) {
            char path[1024];
            if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0)
                continue;
            snprintf(path, sizeof(path), "%s/%s", name, entry->d_name);
            cout << indent << " " << entry->d_name;
            //printf("%*s[%s]\n", indent, "", entry->d_name);
            listdir(path, indent + 2);
        } else {
            cout << indent << " " << entry->d_name;
        }
    }
    closedir(dir);
}


// vector<string> listFilesWithPattern(string directory, string filePatterns);
// reads the raw binding. For each sequence with energy over a threshold, then analyze features and uses the binding positions on antigen to do epitope clustering
//void option7b(string ID_antigen, string bindingDatasetFile){

//    string AntigenName = "";
//    if(ID_antigen.size() < 4){
//        AntigenName = IDshortcut(std::stoi(ID_antigen));
//    } else {
//        AntigenName = IDshortcut(ID_antigen);
//    }

//    cout << "   ... loading antigen " << AntigenName << " from the library" << endl;
//    std::pair<superProtein*, vector<int> > AG = getAntigen(AntigenName);
//    {
//        dataset<binding> data = dataset<binding>(bindingDatasetFile);
//        dataset<analyzedBinding> res = dataset<analyzedBinding>();


////        string originalFileName;
////        vector<string> listIDs;
////        vector<string> sequences;
////        vector<bool> bestForThisCDR3;       // just to make it easy to extract the best binding per CDR3
////        vector<T*> storage;
////        size_t nLines();
////        and binding is:
////        string AAseq;
////        double bindEnergy;
////        string structureID;
//        features f(AG.first, 2);
//        size_t N = data.nLines();
//        for(size_t i = 0; i < N; ++i){
//            vector<string> vres = f.getProperties(*(data.storage[i]));
//            analyzedBinding* toAdd = new analyzedBinding(data.storage[i]->AAseq, data.storage[i]->bindEnergy, data.storage[i]->structureID, vres[interCodeWithIDpos], stringToSet(vres[positionsBound]), vres[hotspot_ID], vres);
//            res.addLine(data.listIDs[i], data.sequences[i], toAdd, data.bestForThisCDR3[i]);
//            // res << ID << sep << largeStr << sep << AAseq << sep << E << sep << interCode << sep << vres[interCodeWithIDpos] << sep << startPos << sep << structure; //<< endl;
//            //res << "\t" << vres[seqAGEpitope] << "\t" << vres[seqABParatope] << "\t" << vres[motifAGEpitope] << "\t" << vres[motifABParatope] << "\t" << vres[motifsSizeGapsLigand] << "\t" << vres[motifsSizeGapsRec] << "\t" << vres[motifsChemicalLig] << "\t" << vres[motifsChemicalRec] << "\t" << vres[agregatesAGEpitope] << "\t" << vres[agregatesABParatope] << "\t" << vres[chemicalAGEpitope] << "\t" << vres[chemicalABParatope] << "\t" << vres[positionsBound] << endl;
//            //cout << vres[interCodeWithIDpos] << endl;
//            //cout << AAseq << sep << E << sep << concat(getProperties(AAseq, startPos, structure), sep) << endl;
//        }
//        res.write(outputFeaturesFile);
//    }

//}


int oldMain(){

    if(0){
        pair<superProtein*, vector<int> > agTest = getAntigen("1FNS_A");
        int receptorSize = 10;
        affinityOneLigand T1 = affinityOneLigand(agTest.first, receptorSize, 11, -1, 1, agTest.second);
        for(int i = 0; i < 1; ++i){
            string Px = randomProt(receptorSize+1);
            std::pair<double, double> res = T1.affinity(Px);
            cout << "affinity(" << Px << ") = \t" << res.first << "\t" << res.second << endl;
        }

        for(int i = 0; i < 1; ++i){
            string Px = randomProt(receptorSize+1);
            std::pair<double, double> res = T1.affinity(Px, true);
            cout << "affinity(" << Px << ") = \t" << res.first << "\t" << res.second << endl;
        }


        return 0;
    }
    //pair<superProtein*, vector<int> > ag1N8Z = getAntigen("1N8Z_C");
    //cout << print(*ag1N8Z.first);
    //return 0;

//    function<long(const int&)> fib;
//        fib = [&](const int& n) -> long {
//            return n < 2 ? n : fib(n - 1) + fib(n - 2);
//        };
//        fib = memoize(fib);
//        for (int i = 0; i < 100; ++i) {
//            cout << i << ": " << fib(i) << endl;
//        }

//        return 0;

    // testRawBindingsRepertoire();
    // return 0;

    //testSetCoveringStructures();
    //return 0;
//    testGenerateSubsets();
//    return 0;
//    cerr << "Start" << endl;
//    vector<string> AGs = listIDs();


//    for(size_t i = 0; i < min(size_t(5000),AGs.size()); ++i){
//        string PDB = getPDB(AGs[i]);
//        string chain = getChain(AGs[i]);
//        cout << "Glycans for " << AGs[i] << "\t";
//
//        readGlycans a(PDB+string(".pdb"),chain);
//    }
//    cerr << "End" << endl;
//    //void testGlycans();
//    return 0;

    //makeHTML();
    // html a("1FNS_A");
    // html a("2J88_A");

    //return 0;

//    showListIDs();
//    vector<string> res = listIDs();
//    for(int i = 0; i < res.size(); ++i){
//        std::pair<superProtein*, vector<int>> AG = getAntigen(res[i]);
//        cout << AG.first->structure->sequence.size() << endl;
//    }
//    return 0;



    // getFeatures 1OB1 "C:/Users/pprobert/Desktop/Successful Bindings/ShortCDR3s/1OB1_CFinalBindings_Process_1_Of_1.txt"  testOutput.txt
    // testDiscretization();
    // return 0;
    //testSurface();
    //return 0;
    //showListIDs();
    //return 0;
    //cout << listAvailableStructures("d:/pprobert/Structures/listFiles.txt");
    //cout << listAvailableStructures("C:/Users/pprobert/Desktop/StructuresInGithub/Absolutions/Structures/ListStructureFilesFRAM.txt");
    //return 0;
    //testRawBindingsRepertoire();
    //return 0;


    //testQuality();
    //testSlide();
    //testFeatures();
    //showListIDs();
    //testImportRepertoire();
    //std::pair<superProtein*, vector<int> > s = getAntigen("1CZ8_P_ShufEpiOut1", true);
    //string agSeq = "VVKFMDVYQRSYCHPIETLVDIFQEYPDEIEYIFKPSCVPLMRCGGCCNDEGLECVPTEESNITMQIMRIKPHQGQHIGEMSFLQHNKCECRPKVVKFMDVYQRSYCHPIETLVDIFQEYPDEIEYIFKPSCVPLMRCGGCCNDEGLECVPTEESNITMQIMRIKPHQGQHIGEMSFLQHNKCECRPK";
    //cout << random::shuffle(agSeq) << endl;
    //return 0;
    //testIncreasingMutations();
    //return 0;
    // PDB 3I40, CA discretized 5.2 ångstrøms
    //    133152	SURRDRDDSDRLLDDRRUUL
    //    128867	LUSDDURDRDDLLSURRUULDDUDSRLRR
    //            >3I40:A|PDBID|CHAIN|SEQUENCE
    //            GIVEQCCTSICSLYQLENYCN
    //            >3I40:B|PDBID|CHAIN|SEQUENCE
    //            FVNQHLCGSHLVEALYLVCGERGFFYTPKA
    //    return -1;

    if(0){
        superProtein v0 = superProtein();
        struct3D S01 = struct3D("SURRDRDDSDRLLDDRRUUL", UnDefined, 133152);
        superProtein v1 = insert(v0, S01, 0);
        struct3D S02 = struct3D("LUSDDURDRDDLLSURRUULDDUDSRLRR", UnDefined, 128867);
        superProtein v2 = insert(v1, S02, 1000);
        v2.setAAs(string("GIVEQCCTSICSLYQLENYCNFVNQHLCGSHLVEALYLVCGERGFFYTPKA"));

//        {
//            features f(&v2, 1);
//            ofstream fwr("C:/Users/pprobert/Desktop/Main/B-CurrentZapotect/Zapotec/Analyzed120000.txt");
//            fwr << f.getPropertiesDeprecated("C:/Users/pprobert/Desktop/Main/B-CurrentZapotect/Zapotec/Dataset120000.txt");
//            fwr.close();
//        }
//        {
//            features f(&v2, 2);
//            ofstream fwr("C:/Users/pprobert/Desktop/Main/B-CurrentZapotect/Zapotec/Analyzed120000Degree2.txt");
//            fwr << f.getPropertiesDeprecated("C:/Users/pprobert/Desktop/Main/B-CurrentZapotect/Zapotec/Dataset120000.txt");
//            fwr.close();
//        }
//        {
//            features f(&v2, 3);
//            ofstream fwr("C:/Users/pprobert/Desktop/Main/B-CurrentZapotect/Zapotec/Analyzed120000Degree3.txt");
//            fwr << f.getPropertiesDeprecated("C:/Users/pprobert/Desktop/Main/B-CurrentZapotect/Zapotec/Dataset120000.txt");
//            fwr.close();
//        }
        return -1;
    }

    //testIncreasingMutations();
    //testEmbarrassing();
    // Sanity check, do not touch
    {
        string s2 = string("SURRDDLLURUSDSSLDSL");
        int testsizeR = 6; //133152
        affinityOneLigand a = affinityOneLigand(s2, "AGRHGRAHARGRAGRTWHSW" , -1, testsizeR, 7, -1, 1.0);
        cout << a.affinity("AGNTWPL",false).first << endl;
        // answer should be -29.36
        //sanityCheckEasyRotate();
        return 0;
    }

    //    //struct3D start = struct3D("");
    //    //protein empty = protein(start);
    //    superProtein v0 = superProtein();
    //    struct3D S01 = struct3D("SUSSDUSSDSSDRRUSSLUDDLRDLSSDSUSDDRUULDSDLUDRDLRSSSDUSRDLUDLDUSDLLURSLUUDDRSSURRSSULRDLUSSRUSURRDDLUURRLRSRUDUSDSUSSSRLDSDDUSLSDRUDDLLURLLRRUULDDUDSRLRDRRLDRSRDSLURSUDUUSSSULDSDLSSDUDURLUUDSUSSDRSRLLUSRSUDULSRLUDUSRLDUDUDSDRLSDLSURLLUUSULSSSSLRRSRURLSDDLLUURSLDDUULRUDULRLRDURSRLUDLULRURDLLDUSSDUDSDRRDLLDDRDLLURDRUDURRDLURRDD", UnDefined, 133152);
    //    superProtein v1 = insert(v0, S01, 0);
    //    //protein v1 = insert(empty, S01, 1000);
    //    struct3D S02 = struct3D("DRUDLRUDUSRSDLRUSSULUSLRSUSSRLRLURLSURRLURDDLDLRUDUSDRSSSULDRRSDSLUURDLUURRDLUSUSLRLRLRSSUURDLURDLLUUSUSRLSLUDRLURDUDULRUUDRULRDDLLURDLUUDRRLSURULDRSSURUSSUDLLUDRUSLSLSLURRDLSDSLSSRLRRUDRLRLSSSDLRLULRRSDUSLRSUUDSSSRLSSSSDDUDLRSLRDLLULUSSLDUDURSLLRLUDSLLUURRDLURRDRDSUSSRLSSDRDSLRUSDRLDURDSRSDUSSLUSDDLLUURSDLUUDLSDUDDLLRRLURR", UnDefined, 165399);
    //    superProtein v2 = insert(v1, S02, 1000);
    //    struct3D S03 =  struct3D("LUDSSSLRSDSDDLSLRSLDSLDUDLRDLSLRLDRDSLULSRUULDRLRUDSSRLRDLDLSRDDLLULRRLLDDRDUSDSDSRLRLUDSLUDDLLUURDDLLSSDULLDUDSRULSSSDULLSURRDDLUSDLUUSDDULRLLRLSDUSDLDSDULDSLRLDUSSLLSSSDULRRSRLSRLRLULLRUDURDLRUSSSLSDDRDURSDURLSRDUSSLRLSDLSSDLLUUDDUUSUSSSLRSSRDUURLSRULLDLURDLLUULRUSLRDULSDRDSSSRUDDLRSURDRSURDLSSDRDDLLUURRDSLRRDSUUSRDDUULSR", UnDefined, 189536);
    //    superProtein v3 = insert(v2, S03, 2000);
    //    v3.setAAs(randomProt(v3.size()));
    //    int sizeR = 3;
    //    //cout << print(v3) << endl;
    //    affinityOneLigand zz = affinityOneLigand(&v3, sizeR, 7, -1, 1.0);

    //    string s2 = string("SURRDDLLURUSDSSLDSLRURUSURLSULUDSUSSLLRLURRDDLLUURUSUDLLRSLSLURSDLLRRUDLLUURRDDLSDLLRSDLUUDRUURDUSRDDLLUURULURRDDULLDDUURRDSRURRDSLRRLLSUURDLSDDUURDDLUUDDLSDRLLUUR");
    //    int sizeR = 6; //133152
    //    affinityOneLigand a = affinityOneLigand(s2, randomProt(s2.size()+1) , -1, sizeR, 7, -1, 1.0);
    //    cout << a.affinity(randomProt(sizeR+1),true).first;
    //    return 0;

    //    // sanity test of the struct3D constructor is not copying pointers
    //    struct3D* v10 = new struct3D("SUUL");
    //    struct3D* v20 = new struct3D(*v10);
    //    delete v20;
    //    cout << *(v10->occupiedPositions.begin()) << endl;

    //    superProtein* v11 = new superProtein("SUUL");
    //    superProtein* v21 = new superProtein(*v11);
    //    delete v21;
    //    cout << *(v11->structure->occupiedPositions.begin()) << endl;
    //    return 0;

    //Figure1(1,8);
    //testMultichain();
    //return 0;

    //testFastAffinity();
    //testGenerateProteins();

    //Figure2(5,8);
    //    string s2 = string("SUSSDUSSDSSDLLUSSLSDDSLSSSDRLSUDDLURDSLDLUURRLUSRLSSSUDLSRDRLSRDDLUSLSRRDLLRSULRUSLRLRLRDURDUURDSURDDSLSLDSLRSSLRSSDSUUDRRSSLURRDDLSRLLUURLUURDSDSSSSRDRSRLRSURRLDUSSUUDUSLRLSDDSUSLSRDRRLRDLRLRLSDSSUURSRSSSLSUDLRULRLDUSDSURRSRULLRRDDLLSULLRSDUSSLSUUDSUULLRSDSRDLLULRLURDSSLURDDSURLRSRDSSUSSRSSUSDUUDLLSRDLLURRLLDDUDSDDUUDLSURDLSSSSSUSSRSSSSULRLRSSDULDSLRSDURURRSRDDLSLSRDDUDUDSRULDSSDRRDLRSRDDUULLDDRRUULRRLRSSUDUDDRRLLDSUURRDLUUDURLULRDULDRLUSDSRDDRLLUUSRDLUURRLSRRULLUULRUSLSDDSSDRURRLDURDUSURULRLSRLLDLUSSSRSRRULRLSSSUSLSRDSSLSSSUDDUDDUSSSUDSSSRRULDDUUDSRUULURLURLSDURSUSULDRDDLUULLDLURRSDUSRLRLRLDSDSSSURDDSRULDULSSSLRRLDSDDURDLLUUDDRSDRLRRDDLUSDLURSSSSUSSSSLRSSLRSSRUSUSUSRSLUSSRDUULSRDDUULDDLSRDUDUSDUSRSLLURDDSLUDDUURRDDLLRRDDSURLSDRLRSDLLDDLSRSDLLUDSLDSUDSSUDSDRURDRSRLLRRDDLDRUURDSSDLLRDUDSRDLURUUSSRLLUSSDLULDLSLRDULDDSDUSSLRSRSRSSLRDULRDSSLUURLRLRSDLDSRSLUDSRLURSLSUULDDUURRDDSLDSUDSSSSDDLUDSRRDDLUDSLUURRSRSLRSLURDURRULRUDRLDSLURRSSRLSUSSLSDLLURDLLUURDDLRSURRDLLUURRSSSSLSUDSSLSRSLRLRDSSLSSURSSSSSDUURLSRRDLDRRSSDUDRLRUDSSSRUSRLSRRLSRDDLLRRLLURDDUDSSRLSSLRDDLURDLSUUDDUSUSSULRDSUSSLSDSSRDDRLSURDDLLURDDLLUDUSLSUSLDSUSLSLSSLRDDUDUDLRSDRDSUDUDLLRLSRLSRDURUDRLSSSSDSLUUDRSDUSSUSSUDSUDUDRLSUDDSLRDDLLUSRDLDLSSSRSDSLSRDRSLRDDLLUURRDDLUUSSSRLUDSSRRSSLRLUSRSURDLULSSRLSDUDDRLURRDDLLUUDDLSRSRDDLUULLD");
    //    struct3D s1 = struct3D(string("SUSSDUSSDSSDLLUSSLSDDSLSSSDRLSUDDLURDSLDLUURRLUSRLSSSUDLSRDRLSRDDLUSLSRRDLLRSULRUSLRLRLRDURDUURDSURDDSLSLDSLRSSLRSSDSUUDRRSSLURRDDLSRLLUURLUURDSDSSSSRDRSRLRSURRLDUSSUUDUSLRLSDDSUSLSRDRRLRDLRLRLSDSSUURSRSSSLSUDLRULRLDUSDSURRSRULLRRDDLLSULLRSDUSSLSUUDSUULLRSDSRDLLULRLURDSSLURDDSURLRSRDSSUSSRSSUSDUUDLLSRDLLURRLLDDUDSDDUUDLSURDLSSSSSUSSRSSSSULRLRSSDULDSLRSDURURRSRDDLSLSRDDUDUDSRULDSSDRRDLRSRDDUULLDDRRUULRRLRSSUDUDDRRLLDSUURRDLUUDURLULRDULDRLUSDSRDDRLLUUSRDLUURRLSRRULLUULRUSLSDDSSDRURRLDURDUSURULRLSRLLDLUSSSRSRRULRLSSSUSLSRDSSLSSSUDDUDDUSSSUDSSSRRULDDUUDSRUULURLURLSDURSUSULDRDDLUULLDLURRSDUSRLRLRLDSDSSSURDDSRULDULSSSLRRLDSDDURDLLUUDDRSDRLRRDDLUSDLURSSSSUSSSSLRSSLRSSRUSUSUSRSLUSSRDUULSRDDUULDDLSRDUDUSDUSRSLLURDDSLUDDUURRDDLLRRDDSURLSDRLRSDLLDDLSRSDLLUDSLDSUDSSUDSDRURDRSRLLRRDDLDRUURDSSDLLRDUDSRDLURUUSSRLLUSSDLULDLSLRDULDDSDUSSLRSRSRSSLRDULRDSSLUURLRLRSDLDSRSLUDSRLURSLSUULDDUURRDDSLDSUDSSSSDDLUDSRRDDLUDSLUURRSRSLRSLURDURRULRUDRLDSLURRSSRLSUSSLSDLLURDLLUURDDLRSURRDLLUURRSSSSLSUDSSLSRSLRLRDSSLSSURSSSSSDUURLSRRDLDRRSSDUDRLRUDSSSRUSRLSRRLSRDDLLRRLLURDDUDSSRLSSLRDDLURDLSUUDDUSUSSULRDSUSSLSDSSRDDRLSURDDLLURDDLLUDUSLSUSLDSUSLSLSSLRDDUDUDLRSDRDSUDUDLLRLSRLSRDURUDRLSSSSDSLUUDRSDUSSUSSUDSUDUDRLSUDDSLRDDLLUSRDLDLSSSRSDSLSRDRSLRDDLLUURRDDLUUSSSRLUDSSRRSSLRLUSRSURDLULSSRLSDUDDRLURRDDLLUUDDLSRSRDDLUULLD"), UnDefined, 133152);
    //    if(s1.properlyFolded) cout << "Properly folded, size " << s1.sequence.size() << endl;
    //    else {cerr << "ERR: the structure entered is self-colliding!!\n   ... " << s1.sequence << endl; return(-1);}

    //TestPDBtoLat();
    //exit(-1);

    //    set<int> emptyset = set<int>();
    //    receptorLigand rl = receptorLigand(s1, 2, 1, emptyset);
    //    rl.generateReceptors();
    //testBigReceptorLigand(2, 1, s2, 133152);
    //testRotate3D();
    //testDiscretize();

    string simpleAccessible = string("DUUSURUSLDDUUDDUDDSUUSRDDRUDDUUDDUDDSUULS");
    string AAsimple = "YFHGCARRATLNTTISWEYVSVDMEKIRVGGNEWFNHTMYVT"; //randomProt(simpleAccessible.size()+1);
    int receptorSize = 8;

    //cerr << "Now displaying the found structures " << endl;
    //struct3D* ligand = new struct3D(simpleAccessible, UnDefined, lattice::idFromPosisition(31,30,34));

    vector<int> forbOptions = {POSITIONS_SQUARE_BLOCKED};
    //set<int>* forb =  generateForbidden(forbOptions);

    /*  glDisplay(argc, argv);
        addToDisplay(ligand, true);
        addToDisplay(forb);
    //  glutPostRedisplay();
    //  glutMainLoop(); */

    affinityOneLigand T1 = affinityOneLigand(simpleAccessible, AAsimple, lattice::idFromPosisition(31,30,34), receptorSize, 4, 4, 0.4, forbOptions);

    cout << "Details of the structures and affinities for " << simpleAccessible << " (" << AAsimple << "), receptors " << receptorSize << " minI=4" << endl;
    for(int i = 0; i < 1; ++i){
        string Px = randomProt(receptorSize+1);
        std::pair<double, double> res = T1.affinity(Px, true);
        cout << "affinity(" << Px << ") = \t" << res.first << "\t" << res.second << endl;
    }

    return 0;
}














//unsigned int fnv_hash (string s)
//{
//    int len = s.size();
//    unsigned int h = 2166136261;
//    int i;

//    for (i = 0; i < len; i++)
//        h = (h*16777619) ^ ((char) s[i]);

//    return h;
//}
/*struct residue {
                            /// \brief Basic constructor for a residue. By default, the AA is undefined, and only the position in space matters \ingroup Prot
    residue(int _IDposition, int _IDresidue = -1, AA _TypeResidue = UndefinedYet)
        : IDposition(_IDposition), IDresidue(_IDresidue), TypeResidue(_TypeResidue) {}

                            /// \brief Copies a Residue \ingroup Prot
    residue(const residue &toCopy) : IDposition(toCopy.IDposition), IDresidue(toCopy.IDresidue), TypeResidue(toCopy.TypeResidue) {}
                            /// \brief Position in space in the lattice \ingroup Prot
    int IDposition;
                            /// \brief An ID that can be given to the residue. By default it can be the position inside a protein : 1,2,3,4, ..., N, but can also be changed. This is not used for any algorithms, so can be user-defined.
    int IDresidue;
                            /// \brief The type of residue, see enum AA above (or UndefinedYet) \ingroup Prot
    AA TypeResidue;
};
Struct3D
    int startingPosition;
                            /// \brief Absolute sequence of moves
    string sequence;
                            /// \brief Decomposition of each step direction in space (i.e. not relative to the previous move, just as translation in 3D)
                            /// this notation is helpful when performing rotations of structures: the singlemoves can be rotated easily, then the sequence is updated.
    string separatedSingleMoves;
                            /// \brief List of observer coordinates (yAxis) after each move. Used mainly for checking consistency.
    string listYAxis;
                            /// \brief Position in space of the last point of the structure (integer encoding, see lattice)
    int endingPosition;
                            /// \brief Set of occupied positions in space by this structure. Useful for collision check
    set<int> occupiedPositions;
                            /// \brief During construction of the structure, this flag says whether the protein is self-colliding or not.
    bool properlyFolded;


prot.points
*/










void print(int argc, char** argv){
    cout << argc << " arguments:\n";
    for(int i = 0; i < argc; ++i){
        cout << "   " << argv[i] << endl;
    }
    cout << endl;
}




int main4(int argc, char** argv){


    //testDistribsFoldedFree();
    //return 0;


    //testVectorsDirections();

    //cerr << 8940696716. / 18446744073. << endl;

//printVector(generateRandomSeqSize(6));
//return 0;


/*
    // Epitope 1
    if(0) displayEpitope("SSRLLRLLRLLRS", lattice::idFromPosisition(32,32,33));
    testBigReceptorLigand(9, 4, "SSRLLRLLRLLRS", lattice::idFromPosisition(32,32,33));

    // Epitope 2
    if(0){
        set<int> block = struct3D("SSSSSUUSSSS", UnDefined, lattice::idFromPosisition(30,33,33)).occupiedPositions;
        vector<int> blockV = vector<int>(1, -2); // -2 inserts the flat square
        std::copy(block.begin(), block.end(), std::back_inserter(blockV));
        //for(int i = 0; i < blockV.size(); ++i){
        //    cout << blockV[i] << endl;
        //}
        displayEpitope("RRUUSURRDSSRSRSSDDSSRS", lattice::idFromPosisition(30,32,33), blockV);
    }

    // Epitope 3
    if(1){
        vector<int> blockV = vector<int>(1, -2); // -2 inserts the flat square
        set<int> block = struct3D("SSSSSUUSSSS", UnDefined, lattice::idFromPosisition(31,34,33)).occupiedPositions;
        std::copy(block.begin(), block.end(), std::back_inserter(blockV));
        set<int> block2 = struct3D("RSSSUUSSS", UnDefined, lattice::idFromPosisition(36,33,33)).occupiedPositions;
        std::copy(block2.begin(), block2.end(), std::back_inserter(blockV));
        //for(int i = 0; i < blockV.size(); ++i){
        //    cout << blockV[i] << endl;
        //}
        displayEpitope("DRRDSRSRSDSDSRSRRLLRSSRSSRSSDDSSRSSRS", lattice::idFromPosisition(32,31,35), blockV);
    }
    */



    // Using getopt to parse arguments.
    //print(argc, argv);
    // short version: https://stackoverflow.com/questions/9642732/parsing-command-line-arguments
    /*int mode = -1;
    int opt;
    while ((opt = getopt(argc, argv, "ilw")) != -1) {
         switch (opt) {
         case 'i': isCaseInsensitive = true; break;
         case 'l': mode = 1; break;
         case 'w': mode = 0; break;
         default:
             fprintf(stderr, "Usage: %s [-ilw] [file...]\n", argv[0]);
             exit(EXIT_FAILURE);
         }
     }*/

    // long version:See https://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Option-Example.html#Getopt-Long-Option-Example

    /*
    static int verbose_flag;        // Flags
    int c = 0;
    while (c != -1){
        static struct option long_options[] =
        {
          {"verbose", no_argument,       &verbose_flag, 1}, // Options to set a flag.
          {"brief",   no_argument,       &verbose_flag, 0},
          {"add",     no_argument,       0, 'a'},   // These options don’t set a flag. We distinguish them by their indices ??.
          {"append",  no_argument,       0, 'b'},
          {"delete",  required_argument, 0, 'd'},
          {"create",  required_argument, 0, 'c'},
          {"file",    required_argument, 0, 'f'},
          {0, 0, 0, 0}
        };

      int option_index = 0;     // getopt_long stores the option index here.
      c = getopt_long (argc, argv, "abc:d:f:", long_options, &option_index);
      switch (c){
        case 0:
          if (long_options[option_index].flag != 0)            // If this option set a flag, do nothing else now.
            break;
          //printf ("option %s", long_options[option_index].name);
          if (optarg)
            //printf (" with arg %s", optarg);
          //printf ("\n");
          break;

        case 'a':
          cout << "-a detected " << endl;
          //puts ("option -a\n");
          break;

        case 'b':
          cout << "-b detected " << endl;
          //puts ("option -b\n");
          break;

        case 'c':
          cout << "-c detected with " << optarg << endl;
          //printf ("option -c with value `%s'\n", optarg);
          break;

        case 'd':
          cout << "-d detected with " << optarg << endl;
          //printf ("option -d with value `%s'\n", optarg);
          break;

        case 'f':
          cout << "-e detected with " << optarg << endl;
          //printf ("option -f with value `%s'\n", optarg);
          break;

        case '?':

          // getopt_long already printed an error message.
          break;

        default:
          // no arguments
          cout << "no argument" << endl;
          //abort ();
          break;
        }
    }

  // Print any remaining command line arguments (not options).
  if (optind < argc)
    {
      cout << "Detected additional parameter: ";
      while (optind < argc){
        cout << argv[optind++] << "\t";
      }
    }
    //return 0;
*/


    //test();
    //cout << printVector(resized(generateProteinsSizeLAndLess(5),6));
    //cout << intToID(0,0,0,0,16384) << endl;

    //test2();
    //testCodingOnExampleSequences();
//    testVectorsDirections();
  // testLattice();
  //  testProteins();
    //testGenerateProteins();


    // test functions one by one

//testEncoding();
    //testRotations();
    //testBigReceptorLigand();
    //showStructures("../xiv/cubicStructuresL26.txt", false);//, 100000, 101000);
    //showStructures("../xiv/cubicStructuresL26NonRedundant.txt", true);//, 100000, 101000);
    //readCubes("../xiv/cubicStructuresL26.txt");
    //testCubicAffinities();
    //testCubicFoldings();
    //testEnsProts();
    //testFoldingStructures();

  //testBigReceptorLigand();
  //testFastAffinity();
  //testProteins();
    //testAAaffinities();

    //testFaces();

    cerr << "continued " << endl;
     //testRecepLigands();
    return 0;
}





















/*
// extending the protein with '-' on the left to reach size 30
char buf[30];
buf[29] = '\0';
for(int i = 0; i < 29; ++i){
    if(i < 29-L) buf[i] = '-';
    else buf[i] = prot[i-29+L];
}

char part1[6], part2[7], part3[7], part4[7], part5[7];
part1[5] = '\0';
part2[6] = '\0';
part3[6] = '\0';
part4[6] = '\0';
part5[6] = '\0';
for(int i = 0; i < 5; ++i){
    part1[i] = buf[i];
}
for(int i = 0; i < 6; ++i){

}
*/
/*
string print(protein& op){
    stringstream res;
    res << "Protein:(ID=" << op.ID << ", Size=" << op.size() << ")\n";
    res << "Residues : " << endl;
    for(unsigned int i = 0; i < op.points.size(); ++i){
        res << "\t" << print(op.points[i]) << endl;
    }
    res << print(op.structure);
    for(int i = 0; i < op.size(); ++i){
        res << AAname(op.points[i].TypeResidue);
    }
    res << "\n\tIDres\tResidue\tIDpos\tx\ty\tz\n";
    for(int i = 0; i < op.size(); ++i){
        res << "\t" << print(op.points[i]) << "\n";
    }
    res << "\nOccupiedPositions\n";
    std::set<int>::iterator it;
    for(it = op.occupiedPositions.begin(); it != op.occupiedPositions.end(); ++it){
        res << *it << " ";
    }

    return res.str();
}*/






// Note: initializing argsThread with malloc raises a segfault on linux compilers but not windows.
// this is likely because no initializer is called.
// one solution here: https://tia.mat.br/posts/2015/05/01/initializing_a_heap_allocated_structure_in_c.html
// but it's a fail because only works in C++, not C

//void *memdup(const void *src, size_t sz) {
//        void *mem = malloc(sz);
//        return mem ? memcpy(mem, src, sz) : nullptr;
//}

//#define ALLOC_INIT(type, ...)   \
//        (type *)memdup((type[]){ __VA_ARGS__ }, sizeof(type))

// to create an argsThread as void*:
//struct foobar *foobar = ALLOC_INIT(struct foobar, {
//        .field = value,
//        .other_field = other_value,
//        .yet_another_field = yet_another_value
//});

//struct argsThread *arguments = ALLOC_INIT(struct argsThread, {
//        .resultFile = resultFile.str(),
//        .receptorSize = receptorSize,
//        .T3 = &T3,
//        .listToProcess = listToProcess,
//        .identificationThread = threadName.str()
//});



//struct3D S = struct3D(structLigand, UnDefined, posStart);
//superProtein Rec = superProtein(S);
//Rec.setAAs(AAligand);
////cout << "The antigen is " << endl;
////cout << print(*currentLigand);
////cout << "The receptor is " << endl;
////cout << print(Rec);
//return structuralFeatures(*currentLigand, Rec, minDegree);




void option5old(string ID_antigen, string bindingDatasetFile, string outputFeaturesFile){

    string AntigenName = "";
    if(ID_antigen.size() < 4){
        AntigenName = IDshortcut(std::stoi(ID_antigen));
    } else {
        AntigenName = IDshortcut(ID_antigen);
    }

    cout << "   ... loading antigen " << AntigenName << " from the library" << endl;
    //std::pair<superProtein*, vector<int> > AG = getAntigen(AntigenName);

    {
        features f(AntigenName, 2);
        ofstream fwr(outputFeaturesFile.c_str());
        string tryFormat1 = f.getPropertiesDeprecated(bindingDatasetFile);
        if(tryFormat1.size() > 0){
            fwr << tryFormat1;
        } else {
            fwr << f.getPropertiesFormat2Deprecated(bindingDatasetFile);
        }
        fwr.close();
    }

}



//#define testPooling
//void newOption5b(string ID_antigen, string folder, int NbProcesses, string outputFeaturesFile){
//    string AntigenName = "";
//    if(ID_antigen.size() < 4){
//        cerr << "Antigen name too short, look afo a number : " << ID_antigen << endl;
//        AntigenName = IDshortcut(std::stoi(ID_antigen));
//    } else {
//        AntigenName = IDshortcut(ID_antigen);
//    }

//    cout << "   ... loading antigen " << AntigenName << " from the library" << endl;
//    std::pair<superProtein*, vector<int> > AG = getAntigen(AntigenName);

//    #ifdef testPooling
//    affinityOneLigand test(AG.first, 10, 11, -1, 1, AG.second);
//    test.setUltraFast(true);
//    #endif

//    dataset<binding> pooledRes; // this one will be sorted and annotated true/false for best binder

//    for(int i = 0; i < NbProcesses; ++i){
//        stringstream fName; fName << AntigenName << "FinalBindings_Process_" << i+1 << "_Of_" << NbProcesses << ".txt";
//        cout << " -> Treating: " << fName.str() << endl;

//        dataset<binding> individualFile = dataset<binding>(fName.str());


//            // file format to be read: AAsequence, Energy, InteractionProfile (for checking), nbStructures, StartPos-Structure1, StartPos-Structure2 ...
//            stringstream res;
//            string sep = "\t";
//            ifstream fr = ifstream(fileName);
//            if(!fr) cerr << "ERR: " << fileName << ", file not found." << endl;
//            int lineNr = 0;

//            while(fr.good()){
//                lineNr++;
//                string AAseq;

//                double E;
//                string interCode;
//                int nS;
//                int startPos;
//                string structure;
//                if(fr >> AAseq){

//                    rep.addSequence(AAseq, -1, true, true);

//                    if(is_number(AAseq)){
//                        //cerr << AAseq << " is number " << endl;
//                        int testFormat = stoi(AAseq);
//                        if((lineNr == 1) && (testFormat > -1) && (testFormat < 1e10)){
//                            cout << "   ... Recognized alternate format 2 for binding datasets\n";
//                            return string();
//                        } else {
//                            cout << "   ... ERR in the formatting. got a weird integer as AA sequence in the first line instead of an AA sequence (format 1)" << endl;
//                            return string();
//                        }
//                    }

//                    fr >> E >> interCode >> nS;
//                    //cout << AAseq << "\t" << E << "\t" << interCode << "\t" << nS << endl;
//                    if((nS < 1) || (nS > 1000)){cerr << "ERR: features::getProperties, getting a too high number of structures (" << nS <<"), around line " << lineNr << ", from reading file " << fileName << endl; return res.str();}
//                    for(int i = 0; i < nS; ++i){
//                        fr >> startPos >> structure;

//                        // for last line, stupid iostream reads empty stuff
//                        if(AAseq.size() > 0){
//        #ifdef ALLOW_GRAPHICS
//                        if(display && (lineNr <= nToDisplay) && (lineNr < 10000)){
//                            addToDisplay(new struct3D(structure, UnDefined, startPos));
//                        }
//        #endif
//                        // normally, will just take the last structure, because the motifs will be the same.
//                        //res << AAseq << sep << E << sep << concat(getProperties(AAseq, startPos, structure), sep) << endl;
//                        // 4 Valentin:
//                        vector<string> vres = getProperties(AAseq, startPos, structure);

//                        if((lineNr == 1) && (i == 0)){
//                            exampleSeq = AAseq;
//                            exampleEnergy = E;
//                            exampleFirstInterCode = interCode;
//                            exampleFirstStruct = structure;
//                        }


//                        if(i == 0) res << AAseq << sep << E << sep << interCode << sep << vres[interCodeWithIDpos] << sep << startPos << sep << structure << endl;
//                        }
//                        //cout << vres[interCodeWithIDpos] << endl;
//                        //cout << AAseq << sep << E << sep << concat(getProperties(AAseq, startPos, structure), sep) << endl;
//                    }
//                }
//            }
//            return res.str();
//        }




//        string tryFormat1 = f.getProperties(bindingDatasetFile);
//        string readyToWrite;
//        if(tryFormat1.size() > 0){
//            readyToWrite = tryFormat1;
//        } else {
//            readyToWrite = f.getPropertiesFormat2(bindingDatasetFile);
//        }

//        #ifdef testPooling
//        // checks that the read file belongs to the correct antigen
//        cout << "   ... testing this file on AA sequence " << f.exampleSeq << " that should have E=" << f.exampleEnergy << endl;
//        double testE = test.affinity(f.exampleSeq).first;
//        if(fabs(testE - f.exampleEnergy) < 0.0001){
//            // then writes the features
//        #endif

//            wrFeatures << readyToWrite;
//            cout << "       test OK" << endl;
//        #ifdef testPooling
//        } else {
//            cerr << "ERR: A file has been generated with another antigen. Skipping this file and continueing:\n   -> " << bindingDatasetFile << endl;
//        }
//        #endif
//        wrFeatures.close();
//        f.rep.writeUpdatedRepertoire(outputListTreatedSequencesFile);


//}




// poolFeatures 1OB1 "C:/Users/pprobert/Desktop/15percentsBis/" "BindingsFor1OB1*.*" whatever1.txt whatever2.txt
// poolFeatures ID_antigen folderFiles patternFileNames outputFeaturesFile outputListTreatedSequencesFile
// Directory should finish by /, very important,
void option5b(string ID_antigen, string folder, string patternFileNames, string outputFeaturesFile, string outputListTreatedSequencesFile){

    string AntigenName = "";
    if(ID_antigen.size() < 4){
        cerr << "Antigen name too short, look afo a number : " << ID_antigen << endl;
        AntigenName = IDshortcut(std::stoi(ID_antigen));
    } else {
        AntigenName = IDshortcut(ID_antigen);
    }

    cout << "   ... loading antigen " << AntigenName << " from the library" << endl;
    std::pair<superProtein*, vector<int> > AG = getAntigen(AntigenName);

    vector<string> liste = listFilesWithPattern(folder, patternFileNames);
    if(liste.size() == 0){
        cout << " ... No file could be found with pattern " << patternFileNames << "\n      in folder " << folder << endl;
        return;
    }

    features f(AntigenName, 2);

    ofstream wrFeatures(outputFeaturesFile);
    if(!wrFeatures){
        cerr << "ERR: poolBindingsFiles, could not open features file output: " << outputFeaturesFile <<endl;
        return;
    }
    // just testing the file can be created before doing all the job
    ofstream wrTreated(outputListTreatedSequencesFile);
    if(!wrTreated){
        cerr << "ERR: poolBindingsFiles, could not open treated list output: " << outputListTreatedSequencesFile << endl;
        return;
    }
    wrTreated.close();

#ifdef testPooling
    affinityOneLigand test(AG.first, 10, 11, -1, 1, AG.second);
    test.setUltraFast(true);
#endif
    for(size_t i = 0; i < liste.size(); ++i){
        string bindingDatasetFile = folder + liste[i];
        cout << " -> Treating: " << liste[i] << endl;

        string tryFormat1 = f.getPropertiesDeprecated(bindingDatasetFile);
        string readyToWrite;
        if(tryFormat1.size() > 0){
            readyToWrite = tryFormat1;
        } else {
            readyToWrite = f.getPropertiesFormat2Deprecated(bindingDatasetFile);
        }

        #ifdef testPooling
        // checks that the read file belongs to the correct antigen
        cout << "   ... testing this file on AA sequence " << f.exampleSeq << " that should have E=" << f.exampleEnergy << endl;
        double testE = test.affinity(f.exampleSeq).first;
        if(fabs(testE - f.exampleEnergy) < 0.0001){
            // then writes the features
        #endif

            wrFeatures << readyToWrite;
            cout << "       test OK" << endl;
        #ifdef testPooling
        } else {
            cerr << "ERR: A file has been generated with another antigen. Skipping this file and continueing:\n   -> " << bindingDatasetFile << endl;
        }
        #endif
    }
    wrFeatures.close();
    f.rep.writeUpdatedRepertoire(outputListTreatedSequencesFile);
}
