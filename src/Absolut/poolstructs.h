#ifndef POOLSTRUCTS_H
#define POOLSTRUCTS_H

#include "ymir.h"
#include <set>
#include <string>

bool areStructEqual(std::string struc1, std::string struc2, std::set<int> & pos1, std::set<int> & pos2);

bool areStructEqual(int pos1, std::string struc1, int pos2, std::string struc2);

bool areStructEqual(struct3D& s1, struct3D& s2);

std::string testEqual(struct3D& s1, struct3D& s2);

std::map<string, int> groupStructuresInClasses(vector<struct3D*> toParse);

std::pair<string, int> retrieveStructureFromID(string ID, char sep = '\t');

std::pair<int, string> retrieveStructureFromPosAndStructure(string ID, char sep = '\t');


string getUniqueIDStructure(struct3D& s1);

string getWrittenUniqueIDStructure(struct3D& s1);

std::pair<string, int> uniqueStructure(struct3D& s1);
std::pair<string, int> oppositeEqualStructure(struct3D& s1);

void testAreStructEqual();

#endif // POOLSTRUCTS_H
