#ifndef RECEPTORLIGAND_H
#define RECEPTORLIGAND_H

#include <vector>
#include <string>
#include "compact.h"
#include "lattice.h"
#include "proteins.h"
using namespace std;

// Option:
// You can define a specific location where the program will generate the two processed structure files
// example: #define manualStructureFilesLocation "$SCRATCH/"
// #define manualStructureFilesLocation "/cluster/work/users/pprobert/Temp/"
// if you wish to use the current folder, just use ""
// of note, the input structure files will still be searched in the local folder (or its parent)
#define manualStructureFilesLocation ""

/// \file
/// \brief Enumeration of (receptor) structures around a predefined protein (ligand: structure + AAs)
/// \date 11th October 2019 \author Philippe A. Robert
/// \defgroup EnumRL Enumeration of (receptor) structures around a predefined protein (ligand: structure + AAs) (receptorligand.h/cpp)


                            /// \brief Class to define a receptor-ligand problem, and to enumerate all possible structure of specified length (in bonds), around the predefined ligand
                            ///  structure, with a minimum number of contacts and possibly excluding a list of 'forbidden/inaccessible' positions \ingroup EnumRL
struct receptorLigand{
                            /// \brief Constructor, to define all constraints/parameters for a receptor-ligand problem \ingroup EnumRL
                            /// \param prot an existing 3D structure of the ligand (it already contains a starting position and starting direction)
                            /// \param sizeReceptor length of receptor structures to generate, in number of bonds/moves. This is the number of AAs minus 1
                            /// \param minimalNInteract Minimum contacts required to the ligand. Minumum 1.
                            /// \param forbiddenVolume list of predefined inaccessible positions
                            /// \return Nothing yet, just storing all the options inside the class.
    receptorLigand(struct3D & _ligand, int _sizeReceptor, int _minimalNInteract, set<int> & _forbiddenVolume) :
        ligand(_ligand), sizeReceptor(_sizeReceptor), minimalNInteract(_minimalNInteract), forbiddenVolume(_forbiddenVolume) {}

                            /// \brief MAIN FUNCTION Generates the list of all possible structures around the ligand according to the predefined options \ingroup EnumRL
                            /// This is performed by recursively calling the tool function generateProteins() (see below) from each possible nearest neighbor of the ligand
                            /// and in all possible directions.
                            /// \return All the structures are stored inside the field vector<struct3D*> possibleReceptors;
    void generateReceptors();

    // inputs:
    struct3D ligand;
    int sizeReceptor;
    int minimalNInteract;
    set<int> forbiddenVolume;

    /// \brief To manually define a set of forbidden/inaccessible positions. note: the codes '-1 or -2 are not recognized here for predefined sets of inaccessible positions.
    /// Please use the function generateForbidden in that case.
    void setForbiddenVolume(vector<int> positions);
    void setForbiddenVolume(set<int> positions);

    // output stored in:
    vector<struct3D*> possibleReceptors;

    /// \brief Additional input: Enumerating structures is independent of their AA content. However, this class can also output 'interaction codes'
    /// that describe how a receptor structure interacts with the ligand, and this incorporates which AA of the ligand is bound (see below for interaction codes)
    /// Therefore, one can specify the AA sequence of the ligand with this function, before calling printToFile(), the function that generate the codes.
    void putSequenceLigand(string _AAsequence);
    string AAsequence;


                            /// \brief MAIN FUNCTION: Writes all the possible structures inside files \ingroup EnumRL
                            /// \param fnameOnlyStruct file where only the list of receptor structures: starting position and absolute sequence, will be written
                            /// \param fnameAllInfo file where the list of receptor structures is written, together with their 'interaction code' to the ligand:
                            /// starting position, absolute sequence, interaction code will be written. This helps in the future retrive all the structures
                            /// with the same interaction code.
                            /// \param fnameCompressed file where only the list of interaction codes is written, together with the number of structures with this code.
    void printToFile(string fnameOnlyStruct);
    int printToFile(string fnameOnlyStruct, string fnameAllInfo, string fnameCompressed);
};
// Put this to 1 or 2 for more text output (debugging) inside the receptorligand class function
#define DBGgenRecept 0



                            /// \brief MAIN TOOL FUNCTION: Recursive function to enumerate all possible structure *tails* of specified length (in bonds), around a predefined ligand
                            ///  structure, with a minimum number of contacts \ingroup EnumRL
                            /// \param length length of receptor structures to generate, in number of bonds/moves. This is the number of AAs minus 1
                            /// \param IDposStart starting point of the tail
                            /// \param absStartDir First move in space of the wished tails
                            /// \param minNbInteract Minimum contacts required for the tail (to the ligand), excluding the contacts already formed by the starting position. Use 0 or -1 if shouldNotTouchProt is true.
                            /// \param shouldNotTouchProt Request no contact to the ligand.
                            /// \param prot structure of the ligand
                            /// \param forbiddenVolume list of predefined inaccessible positions
                            /// \return a list of all possible tails following those constraints (pointers)
                            /// Note that the starting position should have been checked to be feasible (not inside the ligand nor forbidden) before calling this function (not checked inside)
vector<struct3D*> generateProteins(int length, int IDposStart, moveDirection absStartDir, int minNbInteract, bool shouldNotTouchProt, struct3D & prot, set<int> & forbiddenVolume);
// Put this to 1 or 2 for more text output (debugging) inside the generateProteins function
#define DBGgenProts 0

                            /// \brief function that generate all possible structures of a certain length with a minimum number of self-itneractions. \ingroup EnumRL
vector<struct3D*> generateSelfFoldings(int length, int minNbInteract);
                            /// \brief Tool function that generate all possible subtails of a certain length, starting from a position with a certain direction, and with a minimum number of
                            /// self-itneractions to the already generated head of a structure (and without self-colliding nor colliding to the head). \ingroup EnumRL
vector<struct3D*> generateSelfFoldings(int length, int IDposStart, moveDirection absStartDir, int minNbInteract, struct3D & alreadyFoldedTail);


                            /// \brief Transform a vector of 'forbidden positions' into a set. Additionally, the first value can instead define usual sets of
                            /// forbidden positions like a square in the (Oxy) plane. \ingroup EnumRL
set<int>* generateForbidden(vector<int> listForbiddenPositions);
                            /// \brief Code to say no forbidden position. Put it as first value in the vector of positions, the following ones will be excluded. \ingroup EnumRL
#define POSITIONS_ALL_FREE -1
                            /// \brief A square plane in (Oxy) of 10x10 is defined as forbidden positions. Put it as first value in the vector of forbidden positions, and the square will be added. \ingroup EnumRL
#define POSITIONS_SQUARE_BLOCKED -2

                            /// \brief Display a structure and forbidden positions \ingroup EnumRL
void displayLigand(string structureSeq, int startPosition, vector<int> listForbiddenPositions  = vector<int>());
                            /// \brief Display a protein and forbidden positions \ingroup EnumRL
void displayLigand(superProtein* P, vector<int> listForbiddenPositions, bool alwaysVisible = false, bool doNotLoop = false);
                            /// \brief Tutorial/test function on how to use the receptorLigand class. \ingroup EnumRL
void testRecepLigands();
                            /// \brief Tutorial/test function on how to use the GenerateProtein tails tool function. \ingroup EnumRL
void testGenerateProteins();
                            /// \brief Runs the receptorLigand class for any parameter values. \ingroup EnumRL
void testBigReceptorLigand(int sizeReceptors = 6, int minimalNInteract = 4, string structSeq = "ULSRDSLRDUS", int startPos = -1,  vector<int> listForbidden = vector<int>());
                            /// \brief Test function on how to use self-foldings \ingroup CodInter
void testFoldingStructures();




/// \defgroup CodInter Encoding receptor-ligand interactions and saving enumerated receptor structures into files (receptorLigand.h/cpp).

                            /// \brief MAIN FUNCTION: returns an interaction code between a ligand and the receptor. \ingroup CodInter
                            /// an interaction code describes all the binding pairs of AAs, either between the receptor and ligand, or inside the receptor.
                            /// The AAs inside the receptors are coded by their position as a character: a,b,c,d ... (because they are not yet known)
                            /// so receptors can not exceed size 26 (which anyway is too big).
                            /// The AAs inside the ligand are coded by their type of AA, as a capital letter: Leu as L, etc ...
                            /// Therefore, there is no size constraint on the size of the antigen.
                            /// Interactions inside the receptor are named by the two characters (small letters) for the position of the two AAs in the receptor [a-z][b-z],
                            /// Interactions between the receptor and ligand are named as the position in the receptor and the type of AA on the ligand [a-z][A-z]
string interactions(superProtein &ligand, superProtein &receptor);

                            /// \brief Same function as interactions() but from two structures (the ligand and the receptor) and the AA sequence of the ligand, \ingroup CodInter
                            /// instead of using two 'proteins'. This function just creates two proteins classes from both structures, fills the ligand with AAs and calls interactions(protein, protein)
string interactions(struct3D& ligand, struct3D & s2, string AAsequenceLigand);

                            /// \brief Code of self-interactions only inside a structure. See description of the interactions() function for explanations. \ingroup CodInter
string codeSelfInteractions(struct3D& s);

                            /// \brief Generate a unique name for a file with all the receptor structures (position, absolute sequence) of receptors around a ligand. All parameters are put in the filename
                            /// so re-calculations will be avoided by looking if such a file already exists. \ingroup CodInter
string fnameStructures(superProtein* ligand, int sizeReceptors, int minimalNInteract, vector<int> forbiddenPos = vector<int>());
                            /// \brief Generate a unique name for a file with all the receptor structures (position, absolute sequence) and interaction codes, of receptors around a ligand. All parameters are put in the filename
                            /// so re-calculations will be avoided by looking if such a file already exists. This file can be recalculated for a new AA sequence of the ligand by reading the structure file and recomputing the interaction codes. \ingroup CodInter
string fnameStructuresAndCompactForAASeqLigand(superProtein* ligand, int sizeReceptors, int minimalNInteract, vector<int> forbiddenPos = vector<int>());
                            /// \brief Generate a unique name for a file with all the interaction codes of receptors around a ligand. All parameters are put in the filename
                            /// so re-calculations will be avoided by looking if such a file already exists. This file can be recalculated for a new AA sequence of the ligand by reading the structure file and recomputing the interaction codes. \ingroup CodInter
string fileNameCompactForAASeqLigand(superProtein* ligand, int sizeReceptors, int minimalNInteract, vector<int> forbiddenPos = vector<int>());
                            /// \brief Tool function to recreate the files with interaction codes of receptors around a ligand when a ligand has a different AA sequence but the
                            /// file with the list of structures has already been computed \ingroup CodInter
int reGenerateCompressedStructures(string fnameAllStructures, string newSequenceLigand, string fnameOutAll, string fnameOutCompressed, superProtein* forceLigand = nullptr);
                            /// \brief Generate a unique name for a file (library) with receptor sequences and their affinity to a ligand. If only good receptor are kept
                            /// in the file, the energy threshold is included in its name. \ingroup CodInter
string fnameLibrary(string ligandStructSeq, string ligandSeq, int sizeReceptors, int minimalNInteract, double threshold, vector<int> forbiddenPos  = vector<int>());
                            /// \brief Displays the top N structures in a file containing list of structures (containsThirdColumn = false) or list of structures and interaction codes (containsThirdColumn = true). \ingroup CodInter
void showStructures(string fileCompact, bool containsThirdColumn = false, int fromNr=0, int toNr = 1000);
                            /// \brief Tool function: Adds a string to a .txt file name just before an extension. The return file will be .txt whatever the input file type was. \ingroup CodInter
string addBeforeFileType(string fname, string toAdd);
                            /// \brief Tool function: Generates a 'code' summing up a list of positions in space. Used to encode the list of forbidden/inaccessible positions into a filename \ingroup CodInter
string codePos(vector<int> pos);
                            /// \brief Writes the list of self-interactions of a structure in a file. \ingroup CodInter
int exportSelfInteractions(vector<struct3D*> &structures, string fname);
                            /// \brief Looks for an existing file with all self-folding structures with at least k interactions. If not found, recomputes them.
vector<struct3D*> loadOrGenerateSelfFoldings(int minInteractForSelfFold, string folderToLookFor = string(""));

// Potential TODO : replace vector<struct3D*> by vector<string> as output of generateProteins
// Potential TODO idea: for print to file inside generateReceptors, make 'append to file' instead, to be called during generateReceptors and free memory...

#endif // RECEPTORLIGAND_H
