#ifndef DISCRETIZE_H
#define DISCRETIZE_H

/// This file handles the conversion from and to lattice coordinates (tools used by the PDB class)
///
/// In particular, it can:
/// - read a PDB and store the position of AAs (when reading a latfit output, HETATM are the discretized positions)
///     parseLatFitPDB myLatticeProt();
///     myLatticeProt.parseLatFitPDB(discretizedPDB);       //can omit second argument
///
/// - convert latfit PDB files into lattice coordinates and description with multi-chain SLURD
///     myLatticeProt.transform() // note: it only works if the positions are following 90 degrees angles (latfit output).
///
/// - extract the lattice protein as a 'superProtein' class
///     superProtein* inLatticeProtein = myLatticeProt.asSuperProtein()
///
/// (optional: in order to add glycans to the lattice protein:)
/// - read the original PDB file and extract the glycan positions on selected chains
///     readGlycans glycansPos("OriginalPDB.pdb", "AB")
///
/// - add the glycans onto a lattice
///     occupiedPositions = getOccupiedPositions(inLatticeProtein)
///     set<int> listGlycanPositions = myLatticeProt.addGlycans(glycansPos, occupiedPositions)
///
/// All these steps are normally handled by the PDB class, that can read a PDB, call latfit, and calls the functions
/// of this file to read the PDB output of latfit..
///
/// See TestPDBtoLat() function for example how to use.

#include <vector>
#include <cmath>
#include <string>
#include <map>
//Ideally, the 'using namespace std;' should be placed after the system #include <...>s and before the user #include ""s.
using namespace std;

#include "../Ymir/proteins.h"

// Basic definitions for rotating elements
typedef vector<double> vector3D;
bool isVector3D(vector3D &v);
typedef vector< vector<double>> mat3D;
bool isMat3D(mat3D & M);
typedef std::pair<double, double> polarCoord;   // (theta, phi), theta around z, within 0-360 and phi around x

vector3D operator *(mat3D M, vector3D v);
mat3D emptyMat();
mat3D Rotx(double a);
mat3D Roty(double a);
mat3D Rotz(double a);
vector3D rotate (vector3D v, double theta, double phi);


///@brief Parsing PDB files for glycans, and extracts:
/// - the position of the NAN of NAG that start each glycan
/// - the AA and ID of the AA bound by each glycan
/// - the 3D position of these bound AAs in the same order.
///       This will help to retrieve these AAs in the latfit generated PDB,
///       because AA IDs are shifted in the pipeline => we use the position to retrieve the AAs.
struct readGlycans {
    readGlycans(string filename, string chain);
    size_t nGlycans;
    vector< vector3D > listPositionGlycans;
    vector< std::pair<char, int> > AAboundByGlycans;
    vector< vector<double> > PositionsAAboundByGlycans; // position of the CA actually,
    string print();

    // Position of all CAs of the PDB, as side informations of the PDB file in case for later
    std::map<string, vector<double> > positionsAllCAs;
};

vector< vector3D > getPDBChainCoarseGrainedPositions(string PDB, string chains, string typeDiscretization);
string generatePDBfromLattice(vector3D xAxis, vector3D yAxis, vector3D originPosition, superProtein* p1, superProtein* p2=NULL, superProtein* p3=NULL);
void testInOutLatticePDB();

///@brief Parsing the PDB file outputed by Latfit. It should contain floating position of
/// a discretized chain where each bond is of same distance, and angles are 90 degrees.
/// We just need to rescale and rotate to fit it in integer coordinates.
struct latFitToLattice {
    latFitToLattice() {}

    /// @brief reads the PDB and extracts the positions of the AAs of a certain type.
    /// This function is only for coarse grained PDBs where an HETATM is actually a full AA - like in latfit output
    /// Note: the PDB class contains other functions to read a PDB with full atoms description.
    void parseLatFitPDB(string filename, string typeAtoms = "HETATM");
    vector< vector3D > listPositions;
    vector<int> IDresidues;
    string sequence;
    vector<int> IDeachStartingResidue;

    /// @brief takes the first bond as X axis and the first next non-colinear bond as Y axix.
    /// then transforms all points into the lattice: the distance of each bond should be identical
    /// and the angles should be 90 degrees (i.e. should always fall integer in the (Xaxis, Yaxis, Xaxis^Yaxis) coordinates,
    /// if not it will raise an error.
    void transform();
    vector< vector3D > listRotatedOriginPositions;
    vector< vector3D > listIntPositions;
    vector< vector3D > listOriginPositions;
    vector<string> structures;
    vector<int> positions;

    /// @brief Vector in (real floating space), defining the X axis of the lattice.
    vector3D initXAxis;
    /// @brief Vector in (real floating space), defining the Y axis of the lattice.
    vector3D initYAxis;

    /// @brief Export the lattice chains as a superProtein
    superProtein* asSuperProtein();

    /// @brief Takes a list of glycans (read from a PDB file, see class readGlycans)
    /// This function will use the occupied positions of the protein to define its outer
    /// free positions (atmosphere) and assign glycans only within these atmospheric positions.
    set<int> addGlycans(readGlycans& glycansFromOriginalPDB, set<int> &occupiedPositions);

    /// @brief Converts a floating position into the lattice position. The positions do not need
    /// to be integer in the lattice (for instance bystander proteins that were not discretized in the original PDB)
    vector3D PDBtoLattice(vector3D realPos);

    /// @brief Converts an integer lattice position into the real floating position from the original PDB.
    /// Useful to generate a PDB from lattice proteins.
    vector3D LatticetoPDB(vector3D realPos);

    /// @brief test function for PDBtoLattice and LatticetoPDB
    void testConversions();
};


/// @brief Main test function for the latFitToLattice class.
void TestPDBtoLat();

/// @brief test function for an old piece of code that samples N equally distributed directions in space. Not used in ABsolut package.
void testEquallyDistributedPositions();

/// @brief test function for rotation functions
void testRotate3D();

/// tool functions for the transformations. The functions within the latFitToLattice just do the same but with the axis from the class.
vector3D toolPDBtoLattice(vector3D realPos, vector3D initXAxis, vector3D initYAxis, vector3D OriginPos);
vector< vector3D> pooledPDBtoLattice(vector<vector3D> &toTransform, vector3D initXAxis, vector3D initYAxis, vector3D OriginPos);





//class discretize
//{
//public:
//    discretize();
//};

#endif // DISCRETIZE_H
