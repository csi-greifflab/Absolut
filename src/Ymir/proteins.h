#ifndef PROTEINS_H
#define PROTEINS_H

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <set>
#include <iterator>
using namespace std;

#include "compact.h"

/// \file
/// \brief Encoding of proteins = Structures with AAs inside. Note that structures are define by bonds (moves), so there is one more AA than the number of moves.
/// \date 10th October 2019 \author Philippe A. Robert
/// \defgroup Prot Encoding of proteins = Structures with AAs inside (proteins.h/cpp)

//Note: This enum should follow the same order as the values of interaction in the matrix !!!
                            /// \brief Each AA has an integer ID, or UndefinedYet \ingroup Prot
enum AA{Cys, Met, Phe, Ile, Leu, Val, Trp, Tyr, Ala, Gly, Thr, Ser, Asn, Gln, Asp, Glu, His, Arg, Lys, Pro, UndefinedYet, NB_AAs};
                            /// \brief Conversion ID to char \ingroup Prot
char AAname(AA aa);
                            /// \brief Conversion char to ID \ingroup Prot
AA AA_ID(char AA_name);
                            /// \brief Generates random AA sequence \ingroup Prot
string randomProt(int size);
                            /// \brief Checks if an AA exist \ingroup Prot
bool correctAA(char AA_name);

bool checkAAseq(string s);
                            /// \brief Returns a random AA \ingroup Prot
char randomAA();


                            /// \brief Experimental interaction strength (KbT units) between any pair of neighboring AAs, from  Miyazawa 1996 J Mol Biol \ingroup Prot
double AAaffinity(AA a1, AA a2);
                            /// \brief Example how to use the affinity matrix \ingroup Prot
void testAAaffinities();

// from compact
//std::pair<moveDirection, moveDirection> nextAbsoluteMove(moveDirection previous, moveDirection yaxis, moveDirection nextToTranslate);
//moveDirection initialYaxis(moveDirection dir);

                            /// \brief Class storing the information on one residue in a protein. \ingroup Prot
struct residue {
                            /// \brief Basic constructor for a residue. By default, the AA is undefined, and only the position in space matters \ingroup Prot
    residue(int _IDposition, int _IDresidue = -1, AA _TypeResidue = UndefinedYet)
        : IDposition(_IDposition), IDresidue(_IDresidue), TypeResidue(_TypeResidue) {}

                            /// \brief Copies a Residue \ingroup Prot
    residue(const residue &toCopy) : IDposition(toCopy.IDposition), IDresidue(toCopy.IDresidue), TypeResidue(toCopy.TypeResidue) {}
                            /// \brief Position in space in the lattice \ingroup Prot
    int IDposition;
                            /// \brief An ID that can be given to the residue. By default it can be the position inside a protein : 1,2,3,4, ..., N, but can also be changed. This is not used for any algorithms, so can be user-defined.
    int IDresidue;
                            /// \brief The type of residue, see enum AA above (or UndefinedYet) \ingroup Prot
    AA TypeResidue;
};

// Note: functions working with proteins are written outside the classes (i.e. not as methods) to save memory space when handling millions of proteins.

                            /// \brief Printing info on a residue as string \ingroup Prot
string print(residue & a);


                            /// \brief Class for Proteins: Structure in 3D + additional list of residues and their positions in space \ingroup Prot
                            /// Philisophy: A struct3D should be continuous. We want to use discontinuous proteins as well. So wm
                            /// Note that the struct3D class doesn't know/keep the position of each successive point in space, only all occupied positions,
                            /// therefore, the storage of residues (vector 'points') will be used to store each residue position instead.
                            /// There is no restriction on the size of a protein. It also becomes possible to use this class without conventional 'structure'
                            /// (defined by a sequence of moves), but just a list of residues in space (start with an empty structure and add AAs one by one)
struct superProtein {
                            /// \brief A protein is mainly a list of residues, but it can also store a structure together.
    struct3D* structure;

    superProtein();
    ~superProtein();
                            /// \brief Creates a struct3D from an absolute sequence/starting point, and prepare the list of residues for each position as type 'Undefined Yet' with ID their position in the structure \ingroup Prot
    superProtein(string absoluteSequence, int initialPos = -1);
                            /// \brief Copy a struct3D and extends it by creating the list of residues at each position, a 'UndefinedYet' with ID their position in the structure  \ingroup Prot
    superProtein(const struct3D& source);
                            /// \brief Copies an existing protein \ingroup Prot
    superProtein (const superProtein& toCopy);

                            /// \brief A unique new ID is generated each time a protein is created. \ingroup Prot
    int ID;
                            /// \brief Storage of each residue in memory, (on top of inherited struct3D fields) \ingroup Prot
    vector<residue> points;



                            /// \brief Number of stored residues. This is one more than the size of the structure (defined as number of bonds) \ingroup Prot
    int size();
                            /// \brief Says if a position is free (not occupied) \ingroup Prot
    bool isFree(int idPosition);

                            /// \brief  \ingroup Prot
    void push_back(AA newAA, int IDnewPos);

                            /// \brief  \ingroup Prot

    void push_back(residue &toAddToTail);

                            /// \brief Access a residue \ingroup Prot
    residue operator[](int _IDresidue);
                            /// \brief First residue \ingroup Prot
    residue* begin();
                            /// \brief Last residue \ingroup Prot
    residue* end();
                            /// \brief Browse through the residues to get AA sequence in that order. Adds + when jumps.  \ingroup Prot
    string getAAseq();
                            /// \brief Function to set the type of each AA in a protein \ingroup Prot
    void setAAs(string seqAAs);
                            /// \brief Checks that the residues IDs are contiguous and that there is no gap in the structure sequence.
    bool contiguous();

protected:
    // Tool function, called by the constructors: recreates the list of residues and their position from the underlying struct3D structure
    void createPoints();
};

vector<int> startingPositions(superProtein* P);
vector<std::pair<int, string> > getSubChains(superProtein* P);
set<int> getOccupiedPositions(superProtein *P);
                            /// \brief Print informations on a Protein as string \ingroup Prot
//string print(protein& op);

string print(superProtein& P);

// exact same function (needed a function called different than print() - see in cpp file
string printProtein(superProtein& P);

// Tool function to generate new IDs independently per 'Channel'
int getNewID(int channel = 0);



                            /// \brief Structure to store ensemble of proteins as pointers, without copying them, can be useful for
                            /// more comnplex enumerations than struct3D . Right now it is defined identical to a vector<proteins>
                            /// but can be customized later. \ingroup Prot
struct ensProts {
    ensProts(){}
                            /// \brief Storage. \ingroup Prot
    vector<superProtein*> list;
                            /// \brief Add pointer \ingroup Prot
    void add(superProtein* toAdd){
        list.push_back(toAdd);
    }
                            /// \brief Get pointer to ith protein \ingroup Prot
    superProtein* operator[](int i) {
       if((i < 0) || (i > (int) list.size())) {cerr << "ERR: ensProts[" << i << "], out of bounds \n"; return nullptr;}
          return list[i];
    }
                            /// \brief Number of protein stored \ingroup Prot
    int size(){return list.size();}
};
                            /// \brief Prints information stored in a ensProt as string \ingroup Prot
string print(ensProts ep);



                            /// \brief Test function that shows how every function can be used and check they do well \ingroup Prot
void testProteins();
                            /// \brief Test function that shows how every function can be used and check they do well \ingroup Prot
void testEnsProts();

// extending the function collide(struct3D, struct3D) to the type protein
bool collide(superProtein& s1, superProtein& s2);
bool collide(superProtein& s1, struct3D& s2);



superProtein insert(superProtein* existingProt, struct3D* toAdd, int IDfirstResidue);
superProtein insert(superProtein& existingProt, struct3D& toAdd, int IDfirstResidue);
superProtein insert(superProtein* existingProt, string sequence, int startPos, int IDfirstResidue);


#endif // PROTEINS_H
