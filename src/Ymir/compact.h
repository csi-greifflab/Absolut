#ifndef COMPACT_H
#define COMPACT_H

#include <vector>
#include <string>
#include <iostream>
#include <set>
#include <algorithm>
using namespace std;

/// \file
/// \brief Functions for encoding and manipulating lattice structures using the alphabet SUDLR
/// \date 9th October 2019 \author Philippe A. Robert



/// \defgroup BasicMoves Definitions of moves and basic functions (compact.h/cpp).

/// Possible moves are Straight, Right, Left, Up, Down, Backwards, encoded as 'S', 'R', 'L', 'U', 'D', 'B'


                        /// \brief Directions defined as integers. A direction is always relative to the previous move. However, the first one only can be 'Backwards'.
                        /// Absolute moves include 'B'. Relative moves exclude 'B'. \ingroup BasicMoves
enum moveDirection{Straight, Right, Left, Up, Down, Backwards, UnDefined};
#define Nb_Moves_Relative 5
#define Nb_Moves_Absolute 6

                        /// \brief Main function: moving the observer, from a current direction (Ox) and a current yaxis (Oy) according to the wished RELATIVE next move.
                        /// \param previous Current direction of observer (previous X axis)
                        /// \param yaxis Current Y axis of observer.
                        /// \param nextToTranslate Next move (relative to observer)
                        /// \return the new x axis and y axis after moving, as a pair. \image html dox/Compact1-ListMoves.png  \ingroup BasicMoves
std::pair<moveDirection, moveDirection> nextAbsoluteMove(moveDirection previous, moveDirection yaxis, moveDirection nextToTranslate);


                        /// \brief Default predefined initial Y axis for the each first move of the observer  \ingroup BasicMoves
moveDirection initialYaxis(moveDirection dir);

                        /// \brief Gets a random direction, excluding Backwards (Relative move)  \ingroup BasicMoves
moveDirection randRelMove();
                        /// \brief Gets a random direction, including Backwards (Absolute move)  \ingroup BasicMoves
moveDirection randAbsMove();

                        /// \brief Conversions for printing: Char code for each move.  \ingroup BasicMoves
char intToMoveChar(moveDirection m);
char intToMoveChar(int m);
                        /// \brief Conversions moveDirection associated with a Char 'S' 'U' 'D' 'L' 'S'.  \ingroup BasicMoves
int charMoveToInt(char c);  // same as outputing moveDirection. Casts might be necessary.

                        /// \brief Transforms a direction into a 3D vector  \ingroup BasicMoves
vector<int> moveVector(moveDirection dir);
                        /// \brief Transforms a unit integer 3D vector into a direction.  \ingroup BasicMoves
moveDirection moveID(vector<int> v);

                        /// \brief vectorial product between vectors  \ingroup BasicMoves
vector<int> vectorialProduct(vector<int> x, vector<int> y);
                        /// \brief vectorial product between directions  \ingroup BasicMoves
moveDirection vectorialProduct(moveDirection x, moveDirection y);

                        /// \brief Example / Test function for moves  \ingroup BasicMoves
void testVectorsDirections();






/// \defgroup StructManip Functions to manipulate structures in space (compact.h/cpp).

/// A structure will be described as a set of moves. Inside a structure, the next move is defined relative to the
/// position of the observer (like a plane pilot would do). No self-collision is allowed, therefore 'B' is only possible for the first move.
/// We will separate absolute structure (any string of expression '[BSRLUD][SRLUD]*') to relative structures that are rotated
/// to start by S and have 'U' as first turn (expression 'S(S* | S*U[SRLUD]*)' ).
/// Therefore, all absolute structures that are identical by rotation can be described by the same relative structure.
/// The main class to manipulate structures is called struct3D. Most functions to work on struct3D are defined outside the class
/// to minimize space consumption when handling millions of struct3D.


                            /// \brief Compact structures IDs. Compact stuctures should start by S(S*)U..  \ingroup StructManip
typedef long long unsigned int CompactStructure;
                            /// \brief Absolute structures IDs. Absolute stuctures can be whatever sequence (URLDS), and might as well start by B  \ingroup StructManip
typedef long long unsigned int AbsoluteStructure;
                            /// \brief Space structures description: a position in space and an absolute structure. Described by 2 IDs: position and absolute ID \ingroup StructManip
struct SpaceStructure {
    int gridStartPosition;
    CompactStructure AbsoluteStructure;
};
                            /// \brief Transforms an absolute ID structure into the compact version ID (rotated to start by SS*U \ingroup StructManip
CompactStructure compacte(AbsoluteStructure protID);

// tool functions, output info on a structure as string
string printCompact(CompactStructure cs);
string printAbs(AbsoluteStructure as);
string print(SpaceStructure &ss);

/// \brief Class to generate 3D structures in memory from a string of moves in space, and to manipulate them (rotation, fusion) \ingroup StructManip
struct struct3D{

                            /// \brief Main function: Building a 3D structure from an absolute sequence of moves, and from an initial position in space. \ingroup StructManip
                            /// \param AbsoluteSequence   String of moves that may start with B
                            /// \param IDinitposition     Starting position in space of the structure (as a single integer number).
                            ///                           Using -1 by default leads to use the center of the lattice. Use lattice::idFromPosisition(x,y,z) to specify a position
                            /// \param initYAxis          By default (use Undefined), the observer coordinates is predefined to be (Ox, Oy) before the first move, and the coordinate
                            /// after the first move is given by the initialYaxis() function. It is possible to define an alternative observer coordinates
                            /// by giving another yAxis after first move (i.e. after the first move, the xAxis of the observer will still be the direction of first move).
    struct3D(string AbsoluteSequence, moveDirection initYAxis = UnDefined, int IDinitposition = -1);

                            /// \brief Create a structure from an existing one \ingroup StructManip
    struct3D(const struct3D &toCopy);
                            // do not recommend to use this constructor, it puts a bad starting position by default.
    struct3D();
    ~struct3D(){}//{std::cout << "Delete" << &occupiedPositions << endl;}

    // Containers, created automatically:
                            /// \brief Lattice integer code for the starting position of the structure in 3D
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

                            /// \brief Function to elongate a structure by one absolute direction in space (not relative to anything that happened to the structure before)
                            /// the function will find what was the equivalent 'relative move' to add to the tail of the structure to extend in the same direction. \ingroup StructManip
    bool pushBackAbsoluteMove(moveDirection d);
};


                            /// \brief Takes the absolute sequence of a structure in space, and rotates into the relative description of the same structure,
                            /// i.e., starting from S and with U as first turn. \ingroup StructManip
string normalizeAbsolute(string seq);


                            /// \brief Performs a rotation around the firt step/bond/direction of the structure.  \ingroup StructManip
string relativeRotate(string seq, bool clockwise = true);


                            /// \brief Rotate a structure around one axis (defined as a direction), clockwise or anticlockwise. \ingroup StructManip
                            /// To achieve this, it takes the separated list of step directions in space, rotate them and reconstructs it as a sequence of moves.
                            /// \param seq absolute sequence of a structure,
                            /// \param absAxis axis to rotate around (S,U,L,D,R,B), (S being Ox, L being Oy and U being Oz)
                            /// \param clockwise turn clockwise ? Note: call this function twice to rotate 180 degrees.
                            /// \return the absolute sequence of the rotated structure (result.first)
                            /// \return also, the list of separated moves in space (result.second)
string easyRotate(string seq, moveDirection absAxis, bool clockwise = true);


                            /// \brief Concatenates (fuse) two absolute structures. The second structure is rewritten according
                            /// to the position of the observer at the end of the first structure in order to describe the same shape. \ingroup StructManip
                            /// \param seq1 \param seq2 absolute sequences to concatenate
string fuse(string seq1, string seq2);


                            /// \brief describes the same structure from the other side, with the *exact same* configuration in space
                            /// i.e. the absolute moves are the same but reverted and in opposite order. \ingroup StructManip
string revert(string seq1);

// §§§ Philippe: this function should be checked properly
                            /// \brief Generate all possible rotations of an absolute structure around its starting point. By default, all rotations around an axis (Ox,Oy and Oz), clockwise or not,
                            /// including 180 degree rotations, are returned \ingroup StructManip
vector<string> allRotationsStruct(string base, moveDirection ForcingStart = UnDefined, moveDirection firstTurn = UnDefined);
                            /// \brief Generate all possible rotations of a compact structure \ingroup StructManip
vector<AbsoluteStructure> allRotationsStruct(CompactStructure base, moveDirection ForcingStart = UnDefined, moveDirection firstTurn = UnDefined);

void testRotate();

// tool functions
bool contains(set<int> s, int v);
set<int> union_sets(set<int>& s1, set<int>& s2);
//set<int> intersection_sets(set<int> &s1, set<int> &s2);

template <typename T>
set<T> intersection_sets(set<T> &s1, set<T> &s2){
    vector<T> v1 = vector<T>(s1.begin(), s1.end());
    vector<T> v2 = vector<T>(s2.begin(), s2.end());

    // Precondition: make sure they're sorted
    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());

    std::vector<T> v_intersection;

    std::set_intersection(v1.begin(), v1.end(),
                          v2.begin(), v2.end(),
                          std::back_inserter(v_intersection));

    set<T> res = set<T>(v_intersection.begin(), v_intersection.end());
    return res;
}

                            /// \brief Return the intersection of points occupied by two structures. \ingroup StructManip
vector<int> intersection(struct3D& s1, struct3D& s2);
                            /// \brief Says if two structures overlap in space (shared occupied position). \ingroup StructManip
bool collide(struct3D& s1, struct3D& s2);
                            /// \brief Says if a position is a nearest neighbor of any point in a structure. \ingroup StructManip
bool touch(struct3D& s1, int pos);
                            /// \brief Number of positions that are nearest neighbors of a position, inside a structure. \ingroup StructManip
int nbTouchPoints(struct3D& s1, int pos);
                            /// \brief Number of positions that are nearest neighbors of a position, inside a set of positions.
int nbTouchPoints(set<int>& occupPos, int pos);
                            /// \brief Number of points that are nearest neighbors between two structures (they are the interacting residues).  \ingroup StructManip
int nbTouchPoints(struct3D& s1, struct3D & s2);
                            /// \brief List all the neighbor positions of a structure, with possibility to exclude 'forbidden positions'. \param s1 structure \param toExclude list of positions to be excluded. \ingroup StructManip
set<int> neighborPositions(struct3D& s1, set<int> toExclude = set<int>());
                            /// \brief same function, but taking a set<int> of occupied positions. \ingroup StructManip
set<int> neighborPositions(set<int> occupPos, set<int> toExclude = set<int>());

                            /// \brief Retrieves the direction in space of a structure last move. \ingroup StructManip
moveDirection lastAbsDirection(struct3D & s);




/// \defgroup compactEncoding Compact memory encoding of absolute and relative sequences into 64 bits (compact.h/cpp).

                            /// \brief Maximum size of relative sequences that can be represented inside 64 bits. \ingroup compactEncoding
#define MAXL 29

                            /// \brief Sequences in memory will be stored like a string of size exactly MAXL = 29 (filled with '-'s before the sequence starts).
                            /// tool functions that remove the '-' (cut) \ingroup compactEncoding
string cut(string); // removes the initial '-'s. Note: will not detect if there are '-' inside the sequence, such as SU---SSU
                            /// \brief tool function to complete a sequence with '-' at the befinning to reach size MAXL again. \ingroup compactEncoding
string fill(string);

                            /// \brief Checks if a sequence has a correct relative structure description. (any typo or problem will lead to false).
                            /// Self-colliding is not checked yet. \ingroup compactEncoding
bool checkSyntaxSequence(string toCheck);


// 1 - Enumerations

                            /// \brief: enumerates all possible sequences of alphabet [SUDLR] of size exactly L \ingroup compactEncoding
vector<string> generateRandomSeqSize(int L);
                            /// \brief: enumerates all possible *** relative sequences *** from rule S* | S(S*)U[SUDLR]* of size exactly L \ingroup compactEncoding
vector<string> generateRelativeSeqSizeL(int L);
                            /// \brief: enumerates possible relative sequences of size 0 to L included. \ingroup compactEncoding
vector<string> generateRelativeSeqSizeLAndLess(int L);


// 2 - Compact / Relative representations and encoding

                            /// \brief ERROR 64-bits ID code for improperly formatted sequences. (16000 is actually an ID for a sequence of size 6 that is not feasible) \ingroup compactEncoding
#define ERRORSEQ 16000
                            /// \brief Converts a relative structure into a 64 bits ID. Note: Input size can be anything up to 29 characters, and will
                            /// automatically be completed by '-' if necessary.
                            /// \param relSeq relative sequence describing a structure.
                            /// \return a unique 64-bits (long long int) ID for this structure.
                            /// \return the value ERRORSEQ (code '16000', an impossible ID) if there is an error on the sequence. \ingroup compactEncoding
unsigned long long int relativeToInt(string relSeq);
                            /// \brief Retrieves the relative sequence from a 64-bits ID
                            /// \param a unique 64-bits (long long int) ID for this structure.
                            /// \return the absolute sequence matching this ID
                            /// \return the returned string will be 'ERROR' if the ID is not a correct one. \ingroup compactEncoding
string intToRelative(unsigned long long protID);


// 3 - Absolute representations and encoding

                            /// \brief converts an absolute structure into a 64-bits ID. Maximum allowed size is 26 characters.
                            /// returns the value ERRORSEQ (16000, an impossible ID) if there is an error on the sequence. \ingroup compactEncoding
                            /// \param absSeq absolute sequence describing a structure. (note that the starting position has to be stored separately)
                            /// \return a unique 64-bits (long long int) ID for this structure.
                            /// Note: The absolute is first converted into a fake relative sequence by adding the prefix: "SSU" for BACK or "SU" if not starting by Backwards.
                            /// and then computing the ID this absolute sequence would have. Therefore, some IDs can describe both a relative or an absolute sequence,
                            /// and one need to remember what type of sequence was encoded in this ID.
                            /// Note : We only describe absolute structures starting with the default implicit Yaxis plane vector (see initialYaxis()).
unsigned long long absoluteToInt(string absSeq);
                            /// \brief retrieves the absolute structure from its ID. 'ERROR' will be returned if the ID is not valid. \ingroup compactEncoding
                            /// \param a unique 64-bits (long long int) ID for this absolute structure
string intToAbsolute(unsigned long long protID);
// Beware, with absolute sequences, the first movement can be 'Backwards' (B).
// In order to give an ID to absolute structures, they are first transformed into the syntax of relative structures, by
// adding the prefix: "SSU" for Backwards or "SU" if not starting by Backwards (then the next letter is the first move).
// Exemple: BULRDS => SSUULRDS, and RUDDS => SURUDDS. therefore, the maximum length of an absolute move is 26.


// 4 - Tool functions to compute an ID, cutting strings into 5-6-6-6-6 characters (= 8-14-14-14-14 bits)

// Gives the ID (0-15824) of a string of size 6. Returns ERRORSEQ if incorrect. \ingroup compactEncoding
int String6ToInt(string S);
// Merges the ID of each substrings into a large ID (64 bits) \ingroup compactEncoding
unsigned long long int intToID(int part1, int part2, int part3, int part4, int part5);

// 5 - Test/Example functions

                            /// \brief Shows examples of encodings and how to use the functions.  \ingroup compactEncoding
void testCodingOnExampleSequences();
                            /// \brief Tests whether the coding is working: From a string, converts into an ID and back to sequence and tests it matches the input.  \ingroup compactEncoding
void testCoding(string totest);
                            /// \brief Efficiency test: shows the numbers of possible sequences and efficiency of encoding, depending on the number of bytes used. \ingroup compactEncoding
void testEncoding();


// basic tool functions
string print(set<int> &s);
string print(set<string> &s);
string print(struct3D & s);
string printVector(vector<int> v);
string printVector(vector<double> v);
string printVector(vector<size_t> v);
string printVector(vector<string> v);
int norm2(vector<int> toTest);
int norm2(vector<int> v1, vector<int> v2);
double norm2(vector<double> toTest);
double norm2(vector<double> v1, vector<double> v2);

set<int> stringToSet(string text, char sep = '\t');
string setToString(set<int> mySet, char sep = ' ');
bool isIncluded(set<int>& thisSet, set<int>& intoThatOne);
set<int> vectorToSet(vector<int> v);
vector<int> setToVector(set<int> v);
// Remaining TODO make an error in the copy constructor of complex structures, to make sure there is no copy !!

#endif // COMPACT_H
