#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>    // std::count
using namespace std;

int main(void){
    string folder = "C:/Users/pprobert/Desktop/Etulos/4_AnalyzeDataMLTask2/";

    // Infos: This code takes:
    // - a text file with a list of antigen IDs

    string fileAntigens = folder + "ListAntigens142.txt";

    // - a text file with slices and their binding information:
    //      columns: slice label duplicates nDiffAntigens className
    //      the class is "NonBinder" or one of the antigenIDs. The aim is to merge the antigens bound by
    //      the slies into a string with size nAntigens, with 1 at the position of bound antigens IDs.
    //      example: "00100..." means it binds antigen number 3.
    // Please, NO HEADERS in this file

    string fileWithSlices = folder + "Task2Annotated_142_nonredundant.txt";

    //      We manually insert "NonBinder" as first antigenID.

    bool optionAddNonBinderToAntigens = true;

    //      Afterwards, we remove the first column when outputting.

    bool removeFirstColumnNonBinders = true;

    // Output: slice -> "00001010111"

    string fileOutput = folder + "Treated142.txt";






    // First, will allocate an ID to each antigen. Stores them into a vector<string>
    vector<string> IDtoAntigen;
    std::map<string, size_t> AntigenToID;

    // Option: Manually put NonBinder
    if(optionAddNonBinderToAntigens){
        IDtoAntigen.push_back("NonBinder");
        AntigenToID.insert(std::pair<string, size_t>("NonBinder", 0));
    }

    // Read antigen names from the text file 1
    ifstream fAg;
    fAg.open(fileAntigens.c_str());
    if(!fAg) {cerr << fileAntigens << " not found" << endl; return 0;}
    string antigenID = "?";     // sometimes, >> continues after end of file. I parse by always putting "?". If still "?" then >> could not read anything new: end of file.
    while(fAg.good() && (fAg >> antigenID) && (antigenID.compare("?"))){
        // Inserts the antigenID into the vector and dictionary
        IDtoAntigen.push_back(antigenID);
        AntigenToID.insert(std::pair<string, size_t>(antigenID, IDtoAntigen.size()-1));
        antigenID = "?";
    }
    fAg.close();

    size_t nAntigens = IDtoAntigen.size();

    // Now creates a map, for each slice, will store which antigens are bound in a string of size nAntigens
    // database: slice -> "000010100000000"
    std::map<string, string> database;

    ifstream f;
    f.open(fileWithSlices.c_str());
    if(!f) {cerr << fileWithSlices << " not found" << endl; return 0;}

    string slice = "?"; // to check if file finished reading.
    while(f.good() && (f >> slice) && (slice.compare("?"))){

        if(!slice.compare("Slide")){
            char buf[10000];
            f.getline(buf, 9999);
            continue;
        }

        //&& (cout << slice << endl)
        //cout << slice << endl;

        // Other columns. Only className will matter.
        string label;
        int duplicates;
        int nDiffAntigens;
        string className;
        f >> label >> duplicates >> nDiffAntigens >> className;
        //cout << label << duplicates << nDiffAntigens << className << endl;

        // Inm case of header, disc
        if(slice.compare("Slide")) {

            // Creates a new "00000000..." if the slide is new.
            if(database.find(slice) == database.end()){
                database[slice] = string(nAntigens, '0');
            }

            // Now finds the ID of the antigen bound (className) and adds it into the string at position ID.
            if(AntigenToID.find(className) != AntigenToID.end()){
                size_t ID = AntigenToID[className];

                // paranoid error parsing
                if(ID >= nAntigens){cerr << "ERR: AntigenToID has an out of bounds antigen ID" << ID;}
                if(database[slice].size() < nAntigens){cerr << "ERR: Patate en raclette" << endl;}
                if((ID > 0) && (database[slice][ID] == '1')) cerr << "ERR: Duplicate slice " << slice << " ID " << ID << endl;

                database[slice][ID] = '1';
            } else {
                cerr << "ERR: Antigen/Class " << className << " not found" << endl;
            }
        }
        slice = "?";
    }

    // Now only output
    ofstream fwr(fileOutput.c_str());
    if(!fwr) {cerr << fileOutput << " could not be written" << endl; return 0;}
    for(std::map<string,string>::iterator it = database.begin(); it != database.end(); ++it){
        string bindingCode = it-> second;

        if(removeFirstColumnNonBinders){
            bindingCode.erase(bindingCode.begin(), bindingCode.begin()+1); // up to the second position that is not erased
        }

        int counted = std::count(bindingCode.begin(), bindingCode.end(), '1');
        if(counted > 0){
            fwr << it->first << "\t" << counted << "\t" << bindingCode << "\n";
        }
    }
    fwr.close();
}
