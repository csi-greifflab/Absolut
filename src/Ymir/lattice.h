#ifndef LATTICE_H
#define LATTICE_H

#include <vector>
#include <string>
#include <iostream>
using namespace std;

#include "compact.h"


/// \file
/// \brief Encoding of 3D positions in a (euclidian/cubic lattice) into an integer code
/// \date 9th October 2019 \author Philippe A. Robert
/// \defgroup Lattice Encoding of 3D positions in a (euclidian/cubic lattice) into an integer code (lattice.h/cpp)

                        /// \brief Predefined dimensions of the cubic lattice  \ingroup Lattice
#define XWidth 64
                        /// \brief Predefined dimensions of the cubic lattice  \ingroup Lattice
#define YWidth 64
                        /// \brief Predefined dimensions of the cubic lattice  \ingroup Lattice
#define ZWidth 64

struct superProtein;

                        /// \brief (static) Class to convert positions into an integer ID \ingroup Lattice
struct lattice{
                        /// \brief Position integer ID to coordinates \ingroup Lattice
    static vector<int> positionFromID(int IDcoord);
                        /// \brief Coordinates to position integer ID \ingroup Lattice
    static int idFromPosisition(vector<int> position);
                        /// \brief Coordinates to position integer ID \ingroup Lattice
    static int idFromPosisition(double x, double y, double z);
    static int idFromPosisition(int x, int y, int z);

                        /// \brief integer ID of all neighbors of a position (from its ID) \ingroup Lattice
    static vector<int> idNeighbors(int IDcoord);
                        /// \brief integer ID of neighbors of a position that are not within a list of positions
                        /// (like the occupied positions of a structD3 for instance) \ingroup Lattice
    static vector<int> idFreeNeighbors(int IDcoord, set<int> &listPositions);
                        /// \brief Position reached after moving from an initial position following a certain direction \ingroup Lattice

    // This takes the 6 nearest neighbors (neighbors) + the 12 non-nearest but NOT the 8 complete diagonals.
    static vector<int> idLargeNeighbors(int IDcoord);
    static vector<int> idFreeLargeNeighbors(int IDcoord, set<int> &listPositions);

    static int getIdFromAbsoluteMove(int IDcurrent, moveDirection MoveToNext);

                        /// \brief tests whether a position is within the lattice boundaries \ingroup Lattice
    static bool testPos(int x, int y, int z);
                        /// \brief tests whether a position is within the lattice boundaries \ingroup Lattice
    static bool testPos(vector<int> v);
                        /// \brief Returns the ID of the most central position in the lattice \ingroup Lattice
    static int centralPosition();
                        /// \brief Tests whether two positions in space are nearest neighbors \ingroup Lattice
    static bool areNeighbors(int pos1, int pos2);
};

                        /// \brief example on how to use the lattice (see code) \ingroup Lattice
void testLattice();


#endif // LATTICE_H
