#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include "ymir.h"

set<int> getAtmosphere(int onePos, set<int>& totOccupPos, set<int>& alreadyTaken);
set<int> getAtmosphere(superProtein* S);
set<int> getAtmosphere(set<int> &occupiedPositions);
set<int> getSurface(superProtein* S);
string getSurfaceAAs(superProtein* S);
void testSurface();

#endif // TOPOLOGY_H
