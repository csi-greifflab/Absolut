#include "topology.h"
#include "../Ymir/ymir.h"
#include "../Ymir/plot3d.h"
#include "fileformats.h"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include "antigenLib.h" // for testing only
using namespace std;


std::pair<int,int> oneAirAndOneSurfacePositions(set<int> occupiedPositions){
    std::pair<int,int> res(-1,-1);

    // No point should touch the border of the lattice, if not it would fail
    for(set<int>::iterator it = occupiedPositions.begin(); it!= occupiedPositions.end(); ++it){
        vector<int> IDneighb = lattice::idNeighbors(*it);
        for(size_t j = 0; j < IDneighb.size(); ++j){
            if(!lattice::testPos(lattice::positionFromID(IDneighb[j]))){
                 cerr << "ERR: surfaceIDResidues, the protein touches the border of the lattice, the algorithm can not get surface residues when it happens" << endl;
                 return res;
            }
        }
    }
    if(occupiedPositions.size() == 0) return res;

    // Finds a first residue outside, starting from an extreme of the matrix
    int startPos = *(occupiedPositions.begin());
    vector<int> vecStartPos = lattice::positionFromID(startPos);

    // Will go from X=0 to reach the starting position of the protein. The first point on this line that is in the protein will be
    // a 'surface' residue and will be the start of the search.
    vector<int> borderPoint = vecStartPos;
    borderPoint[0] = 0; // starts X = 0;

    int aSurfaceResidue = -1;
    int i = 0;
    int previous = -1;
    int testPos = -1;
    while((aSurfaceResidue < 0) && (i <= vecStartPos[0])) {
        previous = testPos;
        testPos = lattice::idFromPosisition({i, vecStartPos[1], vecStartPos[2]});
        if(occupiedPositions.find(testPos) != occupiedPositions.end()){
            aSurfaceResidue = testPos;
        }
        ++i;
    }
    int anAtmosphericPoint = previous;
    if(aSurfaceResidue == -1){
        cerr << "ERR: surfaceIDResidues, could not find a first surface point starting from the projection of starting position with X=0. Don't know why" << endl;
        return res;
    }
    if(occupiedPositions.find(anAtmosphericPoint) != occupiedPositions.end()){
        cerr << "The last point before entering the structure(anAtmosphericPoint) seems to be in the structure" << endl;
    }
    return std::pair<int,int>(anAtmosphericPoint, aSurfaceResidue);

}



set<int> getAtmosphere(int onePos, set<int>& totOccupPos, set<int>& alreadyTaken);
set<int> getAtmosphere(set<int>& occupiedPositions){
    std::pair<int,int> twoPos = oneAirAndOneSurfacePositions(occupiedPositions);

    set<int> buffer = set<int>();

    return getAtmosphere(twoPos.first, occupiedPositions, buffer);
}

set<int> getAtmosphere(superProtein* S){

    set<int> occupiedPositions = getOccupiedPositions(S);
    return getAtmosphere(occupiedPositions);
}


set<int> getSurface(superProtein* S){
    set<int> res;
    set<int> occup = getOccupiedPositions(S);
    set<int> atmos = getAtmosphere(occup);
    size_t N = S->points.size();
    for(size_t i = 0; i < N; ++i){
        int pos = S->points[i].IDposition;
        set<int> neigh = vectorToSet(lattice::idNeighbors(pos));
        if(intersection_sets<int>(neigh, atmos).size() > 0){
             res.insert(pos);
        }
    }
    return res;
}

string getSurfaceAAs(superProtein* S){
    string res;
    set<int> occup = getOccupiedPositions(S);
    set<int> atmos = getAtmosphere(occup);
    size_t N = S->points.size();
    for(size_t i = 0; i < N; ++i){
        int pos = S->points[i].IDposition;
        set<int> neigh = vectorToSet(lattice::idNeighbors(pos));
        if(intersection_sets<int>(neigh, atmos).size() > 0){
             res.push_back(AAname(S->points[i].TypeResidue));
        } else {
            res.push_back('-');
        }
    }
    return res;
}




// BAD
// should start with a position touching the structure.
// already taken will be modified
//set<int> getAtmosphere(int onePos, set<int>& totOccupPos, set<int>& alreadyTaken){
//    // search only air positions that touch each-other and touch the antigen...
//    set<int> res;
//    // treats also the non-nearest neighbors, but will put them in the list only if they touch
//    vector<int> toTreat = lattice::idFreeLargeNeighbors(onePos, totOccupPos); //lattice::idFreeNeighbors(onePos, totOccupPos); // lattice::idFreeLargeNeighbors(onePos, totOccupPos);
//    for(size_t i = 0; i < toTreat.size(); ++i){
//        int neigh = toTreat[i];
//        if(alreadyTaken.find(neigh) == alreadyTaken.end()){
//            set<int> aroundSet = vectorToSet(lattice::idNeighbors(neigh));
//            if(intersection_sets<int>(aroundSet, totOccupPos).size() > 0){
//                alreadyTaken.insert(neigh);
//                res.insert(neigh);
//                set<int> retrieved = getAtmosphere(neigh, totOccupPos, alreadyTaken);
//                res = union_sets(res, retrieved);
//            }
//        }
//    }
//    return res;
//}

// THIS IS THE GOOD oNE
// should start with a position touching the structure.
// already taken will be modified
set<int> getAtmosphere(int onePos, set<int>& totOccupPos, set<int>& alreadyTaken){
    // search only air positions that touch each-other and touch the antigen...
    set<int> res;
    // treats also the non-nearest neighbors, but will put them in the list only if they touch
    vector<int> toTreat = lattice::idFreeNeighbors(onePos, totOccupPos); //lattice::idFreeLargeNeighbors(onePos, totOccupPos); //lattice::idFreeNeighbors(onePos, totOccupPos); // lattice::idFreeLargeNeighbors(onePos, totOccupPos);
    for(size_t i = 0; i < toTreat.size(); ++i){
        int neigh = toTreat[i];
        if(alreadyTaken.find(neigh) == alreadyTaken.end()){
            set<int> aroundSet = vectorToSet(lattice::idLargeNeighbors(neigh));
            if(intersection_sets<int>(aroundSet, totOccupPos).size() > 0){
                alreadyTaken.insert(neigh);
                res.insert(neigh);
                set<int> retrieved = getAtmosphere(neigh, totOccupPos, alreadyTaken);
                res = union_sets(res, retrieved);
            }
        }
    }
    return res;
}

///// This algo is not good, it will include inside residues touching holes as atmosphere (in particular if a hole is half-buried (45 dgree=
/// => Do not use!!
// surface propagation:
// deep tree search when only considering neighboring points that are touching a free point, and starting from an outside position

// recursive function that, starting from a surface residues, deep searches not-yet-treated neighbor residues with a free neighbor position (in contact with the air)
// Very important, remainingResiduesToTreat is a GLOBAL variable to all recursive calls, that will be MODIFIED, if not the recursion would treat multiple times the same point from different branches!

// note: if you start inside a pocket, all neighbors are inside => Need to get the 8 neeighbors...
// also, if hole is only two steps down, should not take the residue under you => the free neighbors should intersect... (same air)
/// Bad, do not use.
set<int> recAddSurfaceNeighbors(int thisPos, set<int>& totOccupPos, set<int> &remainingResiduesToTreat){

    set<int> res;


    // which neighbor points are remaining to treat (inside the protein)
    vector<int> neighborPoints = lattice::idLargeNeighbors(thisPos);
    set<int> neighborsAsSet = set<int>(neighborPoints.begin(), neighborPoints.end());
    set<int> nonTreatedNeighbors = intersection_sets<int>(neighborsAsSet, remainingResiduesToTreat);

    // Set of outside positions seen by this point

    // for each neighbor, check if it touches air (number of neighbors outside occupied+blocked)
    vector<int> posToTreatNext;
    for(set<int>::iterator it = nonTreatedNeighbors.begin(); it != nonTreatedNeighbors.end(); ++it){
        vector<int> neighborOutsidePoints = lattice::idFreeNeighbors(*it, totOccupPos);
        // if this position 'touches' the outside, adds it to the treated list
        if(neighborOutsidePoints.size() > 0){
            if(res.find(*it) != res.end()){
                cerr << "ERR: addSurfaceNeighbors, position treated twice: " << *it <<  endl;
            }
            //cout << "Adding position " << *it << endl;
            res.insert(*it);
            if(remainingResiduesToTreat.find(*it) == remainingResiduesToTreat.end()){
                cerr << "ERR: addSurfaceNeighbors, a position is removed twice? " << *it << endl;
            } else {
                remainingResiduesToTreat.erase(*it);
            }
            posToTreatNext.push_back(*it);
        }
    }
    for(size_t i = 0; i < posToTreatNext.size(); ++i){
        set<int> recursiveRes = recAddSurfaceNeighbors(posToTreatNext[i],  totOccupPos, remainingResiduesToTreat );
        res = union_sets(res, recursiveRes);
    }
    return res;
}

/// Bad, do not use.
set<int> addSurfaceNeighbors(int thisPos, set<int>& totalBlocked, set<int> &occupiedPositions){
    set<int> remainingToTreat = set<int>(occupiedPositions);
    set<int> res = recAddSurfaceNeighbors(thisPos, totalBlocked, remainingToTreat);
    if(lattice::idFreeNeighbors(thisPos, totalBlocked).size() > 0){
        res.insert(thisPos);
    }
    return res;
}


/// compares different methods. Only the atmosphere has to be trusted, nothing else
set<int> surfaceIDResidues(superProtein * S, vector<int>& forbiddenPoints){
    if((!S) | (S->size() == 0)) return set<int>();
    set<int> OP = getOccupiedPositions(S);
    set<int> forb = set<int>(forbiddenPoints.begin(), forbiddenPoints.end());
    set<int> totalBlocked = union_sets(OP, forb);
    if(intersection_sets<int>(OP, forb).size() > 0) cerr << "ERR: an occupied position is also forbidden" << endl;

    std::pair<int,int> twoPos = oneAirAndOneSurfacePositions(OP);
    int anAtmosphericPoint = twoPos.first;
    int aSurfaceResidue = twoPos.second;

    // First, naive test with all positions that do touch free positions
    set<int> naivePoints;
    for(set<int>::iterator it = OP.begin(); it != OP.end(); ++it){
        vector<int> neighborOutsidePoints = lattice::idFreeNeighbors(*it, totalBlocked);
        if(neighborOutsidePoints.size() > 0) naivePoints.insert(*it);
    }

    // Third, naive test with all positions that do touch the atmosphere
    set<int> atmosphericPoints;
    set<int> buffer = set<int>();
    set<int> atmosphere = getAtmosphere(anAtmosphericPoint, totalBlocked, buffer);
    cout << "Starging atmos point " << anAtmosphericPoint << ", Got atmosphere of size " << atmosphere.size() << endl;
    for(set<int>::iterator it = OP.begin(); it != OP.end(); ++it){
        vector<int> directNeighbors = lattice::idNeighbors(*it);
        set<int> neighbSet = set<int>(directNeighbors.begin(), directNeighbors.end());
        if(intersection_sets<int>(neighbSet, atmosphere).size() > 0) atmosphericPoints.insert(*it);
    }

    // Now, surface search (not naive algorithm) - careful, OP will be modified => I think this algo doesn't work for half hidden positions!
    set<int> res = addSurfaceNeighbors(aSurfaceResidue, totalBlocked, OP);

    // Now compares:
    //res = atmosphericPoints;
    if((atmosphericPoints.size() > 0) && (union_sets(naivePoints, atmosphericPoints).size() == naivePoints.size())){
        cout << "For this case, the naive and search approach raised the same results with " << res.size() << " elements " << endl;
    } else {
        cout << "Naive approach raised " << naivePoints.size() << " points while search raised " << res.size() << "points. details:" << endl;
        cout << " -> Naive: " << print(naivePoints) << endl;
        cout << " -> Search:" << print(res) << endl;
        cout << "Additional positions that are in the Naive: " << endl;
        for(set<int>::iterator it=naivePoints.begin(); it != naivePoints.end(); ++it){
            if(res.find(*it) == res.end()){
                cout << *it << "\t";
            }
        }
        cout << "Additional positions that are in the Search: " << endl;
        for(set<int>::iterator it=res.begin(); it != res.end(); ++it){
            if(naivePoints.find(*it) == naivePoints.end()){
                cout << *it << "\t";
            }
        }
        cout << endl;
    }

    //cout << "Surface is " << endl;
    vector<std::pair<int, double> > heatmap;
    for(set<int>::iterator it = res.begin(); it != res.end(); ++it){
         heatmap.push_back(std::pair<int, double>(*it, 0.5));
         //cout << *it << "\t";
    }

    //cout << "Atmosphere is " << endl;
    vector<std::pair<int, double> > heatmapAtmosphere;
    for(set<int>::iterator it = atmosphere.begin(); it != atmosphere.end(); ++it){
         heatmapAtmosphere.push_back(std::pair<int, double>(*it, 0.9));
         //cout << *it << "\t";
    }

#ifdef ALLOW_GRAPHICS
    glDisplay();
    addToDisplay(S, false);
    //addToDisplay(heatmap);
    addToDisplay(new vector< std::pair<int, double> >(heatmapAtmosphere));
    set<int>* forbidden = new set<int>(vectorToSet(forbiddenPoints));
    addToDisplay(forbidden);
    addToDisplay(new struct3D("S", UnDefined, 137633 ));
    glutMainLoop();
#endif

    return res;
}

// how do we do for glycans? => can be helpful to study the effect of the glycans
void testSurface(){
    if(0){
        // test 1: cube with empty spot in the middle
        struct3D v1 = struct3D("SRRSRSRSDDSRSRSRUURSRSRS"); // cube with 2 holes,
        vector<int> startVec = lattice::positionFromID(lattice::centralPosition());
        vector<int> blockPosCloseCube = {lattice::idFromPosisition({startVec[0], startVec[1], startVec[2]-2})};
        //displayLigand(v1.sequence, lattice::centralPosition(), {blockPosCloseCube});

        cout << "Case 1: Cube with central position empty, blocked by a forbidden position => Should raise 25 positions" << endl;
        struct3D v2 = struct3D("SSDRRSRSRSDDSRSRSRUURSRSRS"); // full cube with 2 holes,
        //displayLigand(v2.sequence, lattice::centralPosition());
        //glutMainLoop();
        superProtein* S1 = new superProtein(v1);
        set<int> res2 = surfaceIDResidues(S1, blockPosCloseCube  );
        cout << print(res2) << endl;

        cout << "Case 2: Completely full cube => Should raise 26 positions. Note: StartPos is" << v2.startingPosition << endl;
        vector<int> emptySet = vector<int>();
        superProtein* S2 = new superProtein(v2);
        set<int> res3 = surfaceIDResidues(S2, emptySet );
        cout << print(res3) << endl;
    }
    cout << "Case 3: All antigens from the library" << endl;
//    std::pair<superProtein* , vector<int> > AG = getAntigen("1FNS_A");
//    set<int> res = surfaceIDResidues(AG.first, AG.second);


    //vector<string> toDo = listIDs();
    //for(size_t i = 0; i < toDo.size(); ++i){
        //string AGID = toDo[i];
        string AGID = "3ULU_A";
        std::pair<superProtein* , vector<int> > AG = getAntigen(AGID);
        //displayLigand(AG.first, AG.second);
        //return;
        set<int> res = surfaceIDResidues(AG.first, AG.second);
        size_t occup = getOccupiedPositions(AG.first).size();
        double ratio = (double) res.size() / (double) occup;
        cout << AGID << ", " << getOccupiedPositions(AG.first).size() << " Occup pos, " << res.size() << " surface => " << 100*ratio  << "%" << endl;
        if(ratio < 0.5) {
            cout << print(*AG.first);
            cout << print(res);
        }
    //}
}



