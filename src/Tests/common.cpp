
#include "common.h"
#include <algorithm> // for replace

#include <iostream>
using namespace std;

#ifdef WINDOWS
#include <windows.h>
#endif
#ifdef UNIX
#include <sys/stat.h>
#endif

void createFolder(string folderName){

#ifdef WINDOWS
// .................. Now this is all the shit required to create a folder in windows : .............. need to be in the wchar* type, not char*, so need to convert.
const char *p= folderName.c_str(); const WCHAR *pwcsName;
int nChars = MultiByteToWideChar(CP_ACP, 0, p, -1, NULL, 0);
pwcsName = new WCHAR[nChars];
MultiByteToWideChar(CP_ACP, 0, p, -1, (LPWSTR)pwcsName, nChars);
CreateDirectory(pwcsName, NULL); delete [] pwcsName;
#endif

#ifdef UNIX
const int dir_err = mkdir(folderName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
if ((-1 == dir_err) && (errno != EEXIST)){cerr << "Error creating directory : " << folderName << endl;}
#endif
}

string removeFolderFromFile(string file){
    int S = file.size();
    int cpt = 0;
    for(int i = 0; i < S; ++i){
        if((file[i] == '\\') || (file[i] == '/')){
            cpt = i+1;
        }
    }
    return file.substr(cpt, S - cpt);
}


// listing the subfolders in a folder :
// from http://stackoverflow.com/questions/6133647/how-do-i-list-subdirectories-in-windows-using-c
#ifdef WINDOWS
#include <io.h>
#endif
#ifdef UNIX
#include <unistd.h>
#include <dirent.h>
#endif
#include <iostream>
#include <stdio.h>

using namespace std;

#ifdef WINDOWS
// here, keep the whole path
vector<string> listSubDirectories(string dir)
{
    vector<string> res;
    char originalDirectory[_MAX_PATH];

    // Get the current directory so we can return to it
    _getcwd(originalDirectory, _MAX_PATH);

    //SetCurrentDirectory(dir.c_str())
    chdir(dir.c_str());  // Change to the working directory
    _finddata_t fileinfo;

    // This will grab the first file in the directory
    // "*" can be changed if you only want to look for specific files
    intptr_t handle = _findfirst("*", &fileinfo);

    if(handle == -1)  // No files or directories found
    {
        perror("Error searching for file");
        exit(1);
    }

    do
    {
        if(strcmp(fileinfo.name, ".") == 0 || strcmp(fileinfo.name, "..") == 0)
            continue;
        if(fileinfo.attrib & _A_SUBDIR){ // Use bitmask to see if this is a directory
            //cout << "This is a directory : " << fileinfo.name << endl;
            res.push_back(dir + string(fileinfo.name) + string("/"));
        }
        else
        {//cout << "This is a file : " << fileinfo.name << endl;
        }
    } while(_findnext(handle, &fileinfo) == 0);

    _findclose(handle); // Close the stream

    chdir(originalDirectory);
    return res;
}
#endif
#ifdef UNIX
// here, keep the whole path
vector<string> listSubDirectories(string namedir)
{
    vector<string> res;
    DIR* dir = opendir(namedir.c_str());
    struct dirent *entry;
    if (!(entry = readdir(dir)))  {return res;}
    do {
        if ((!string(entry->d_name).compare(string(".")))|| (!string(entry->d_name).compare(".."))) // didnt find strcomp ... don't like it anyways
           continue;
        if (entry->d_type == DT_DIR){
            res.push_back(namedir + string(entry->d_name) + string("/"));
        }
        else
        {//cout << "This is a file : " << fileinfo.name << endl;
        }
    } while ((entry = readdir(dir)) != NULL);
    closedir(dir);
    return res;
}
#endif


#ifdef WINDOWS
string currentDir(){
    char originalDirectory[_MAX_PATH];
    _getcwd(originalDirectory, _MAX_PATH);
    string res = string(originalDirectory);
    replace(res.begin(), res.end(), '\\', '/');
    return res;
}
#endif

#ifdef UNIX
string currentDir(){
    char originalDirectory[1024];
    getcwd(originalDirectory, 1024);
    string res = string(originalDirectory);
    replace(res.begin(), res.end(), '\\', '/');
    return res;
}
#endif


// from http://www.cplusplus.com/reference/string/string/find_last_of/
/*void SplitFilename (string& str)
{
  size_t found;
  cout << "Splitting: " << str << endl;
  found=str.find_last_of("/\\");
  cout << " folder: " << str.substr(0,found) << endl;
  cout << " file: " << str.substr(found+1) << endl;
  return
}*/

string getParentFolder(string dir){
    if(dir.size() == 0) return dir;
    if((dir[dir.size()-1] == '\\') || (dir[dir.size()-1] == '/')) {
        dir = dir.substr(0, dir.size()-1);
    }
    size_t found=dir.find_last_of("/\\"); //,dir.find_last_of("/\/")) ; /// §§§§ THIS SHOULD BE TESTED IN LINUX
    return dir.substr(0,found) + string("/");
}

vector<string> getAllResultSubFolders(string dir){
    vector<string> res;
    //if(dir.size() == 0) dir = folder;
    if(dir.size() == 0) return res;

    vector<string> foldersToRead;
    foldersToRead.push_back(dir);
    int cpt = 0;

    while((foldersToRead.size() > 0) && (cpt < 1000)){
        string nextFolder = foldersToRead.back();
        foldersToRead.pop_back();
        res.push_back(nextFolder);
        vector<string> newFolders = listSubDirectories(nextFolder);
        for(int i = 0; i < (int) newFolders.size();++i){
            foldersToRead.push_back(newFolders[i]);
        }
        cpt++;
    }
    if(cpt == 1000) cerr << "ERR:  getAllResultSubFolders(" << dir << "), too many subfolder, or infinite loop." << endl;
    return res;
}

#ifdef WINDOWS
vector<string> listFilesInDir(string dir, string containing)
{
    vector<string> res;
    //cerr << "Check Folder :" << dir << endl;
    char originalDirectory[_MAX_PATH];

    // Get the current directory so we can return to it
    _getcwd(originalDirectory, _MAX_PATH);

    //SetCurrentDirectory(dir.c_str())
    chdir(dir.c_str());  // Change to the working directory
    _finddata_t fileinfo;

    // This will grab the first file in the directory
    // "*" can be changed if you only want to look for specific files
    intptr_t handle = _findfirst("*", &fileinfo);

    if(handle == -1)  {
        //perror("Error searching for file");
        //exit(1);
        return res;
    }
    do {
        if(strcmp(fileinfo.name, ".") == 0 || strcmp(fileinfo.name, "..") == 0)
            continue;
        if(fileinfo.attrib & _A_SUBDIR){ // Use bitmask to see if this is a directory
            //cout << "This is a directory : " << fileinfo.name << endl;
            //res.push_back(string(fileinfo.name));
        }
        else {
            string tp = string(fileinfo.name);
            if(((containing.size() > 0) && (tp.find(containing) != std::string::npos)) || (containing.size() == 0)){
                res.push_back(tp);
                //cerr << "GOT IT !" << endl;
            }
            //cout << "This is a file : " << fileinfo.name << endl;
        }
    } while(_findnext(handle, &fileinfo) == 0);
    _findclose(handle); // Close the stream
    chdir(originalDirectory);
    return res;
}
#endif

#ifdef UNIX
// here, keep the whole path
vector<string> listFilesInDir(string namedir, string containing)
{
    vector<string> res;
    DIR* dir = opendir(namedir.c_str());
    struct dirent *entry;
    if (!(entry = readdir(dir)))  {return res;}
    do {
        if ((!string(entry->d_name).compare(string(".")))|| (!string(entry->d_name).compare(".."))) // didnt find strcomp ... don't like it anyways
           continue;
        if (entry->d_type ==  DT_REG){ // regular file
            string tp =  string(entry->d_name);
            if(((containing.size() > 0) && (tp.find(containing) != std::string::npos)) || (containing.size() == 0)){
                res.push_back(tp);
            }
        }
        else
        {//cout << "This is a file : " << fileinfo.name << endl;
        }
    } while ((entry = readdir(dir)) != NULL);
    closedir(dir);
    return res;
}
#endif



// I don't use qtcreator/qt routines here because it should be able to work in a Qt free computer as well.
// looks for the THdiff.pro folder.
// first solution : look, from the exe folder, ../Sources/ and look for THDiff.pro
// second solution : look inside the exe folder if there is THDiff.pro
// third solution : look ../ and look for THdiff.pro

string locateProjectDirectory(){
    string exeFolder = currentDir();
    if(listFilesInDir(getParentFolder(exeFolder) + string("Sources/"), string("THDiff.pro")).size() > 0) return getParentFolder(exeFolder) + string("Sources/");
    if(listFilesInDir(exeFolder, string("THDiff.pro")).size() > 0) return exeFolder;
    if(listFilesInDir(getParentFolder(exeFolder), string("THDiff.pro")).size() > 0) return getParentFolder(exeFolder);
    return string("NotFound!!!");
}

vector<string> findAllResultFolders(string dir){
    vector<string> liste1 = getAllResultSubFolders(dir);
    vector<string> res;
    for(int i = 0; i < (int) liste1.size(); ++i){
        if(listFilesInDir(liste1[i],string("History.txt")).size() > 0){
            res.push_back(liste1[i]);
        }
    }
    return res;
}

/*void printVector(vector<string> l){
    cout << "[" << l.size() << "] ";
    for(int i = 0; i < (int) l.size(); ++i){
        cout << l[i] << ((i < (int) l.size()) ? "\n" : "");
    }
}*/

void testDirectoryFunctions(){
  /*  string exeFolder = currentDir();
    printVector(getAllResultSubFolders(exeFolder));
    cout << getParentFolder(exeFolder) << endl;
    printVector(listFilesInDir(exeFolder));*/
}


/*void printVector(vector<string> v){
    for(int i = 0; i < (int) v.size(); ++i){
        cout << v[i] << endl;
    }
}*/


