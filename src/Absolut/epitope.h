#ifndef EPITOPE_H
#define EPITOPE_H

// This class aims at defining binding hotspots (epitopes) on the antigen,

#include <set>
#include <string>
#include <iostream>
#include <vector>
#include "fileformats.h"
#include "../Ymir/compact.h" // to get print(set)
using namespace std;

void testGenerateSubsets();


void testSetCoveringStructures();
void showBindingHotspots(dataset<analyzedBinding>& annotatedDataset, string antigenID, int sizeSet = 4);
void showParatopeEpitope(superProtein* prot1, superProtein* prot2, vector<int> forbiddenPos = vector<int>());
//class epitope
//{
//public:
//    epitope();
//    vector< binding > annotatedDataset;



//};

#endif // EPITOPE_H
