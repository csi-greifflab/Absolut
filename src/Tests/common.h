#ifndef COMMON_H
#define COMMON_H




#ifdef _WIN32
   //define something for Windows (32-bit and 64-bit, this part is common)
    #define WINDOWS
#ifdef _WIN64
      //define something for Windows (64-bit only)
   #endif
#elif __APPLE__
    #define MAC
    #include "TargetConditionals.h"
    #if TARGET_IPHONE_SIMULATOR
         // iOS Simulator
    #elif TARGET_OS_IPHONE
        // iOS device
    #elif TARGET_OS_MAC
        // Other kinds of Mac OS
    #else
    #   error "Unknown Apple platform"
    #endif
#elif __linux__
    #define UNIX
    // linux
#elif __unix__ // all unices not caught above
    #define UNIX
    // Unix
#elif defined(_POSIX_VERSION)
    // POSIX
#else
#   error "Unknown compiler"
#endif

#include <string>
#include <vector>
using namespace std;

void createFolder(string folderName);       // because this sh** is OS dependent, better in a function

vector<string> listSubDirectories(string dir);

string currentDir();
string getParentFolder(string dir);
vector<string> getAllResultSubFolders(string dir); // folders containing a 'History' file
vector<string> listFilesInDir(string dir, string containing = string(""));
string locateProjectDirectory();
vector<string> findAllResultFolders(string dir);
//void printVector(vector<string> v);
void testDirectoryFunctions();

#endif // COMMON_H
