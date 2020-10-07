
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <time.h>
#include <iostream>
#include <iomanip>
#include <limits.h>
#include <climits>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>

#include <biu/Point.hh>
#include <biu/Matrix.hh>
#include "biu/Rotator3D.hh"
#include <biu/OptionParser.hh>
#include <biu/LatticeDescriptorSQR.hh>
#include <biu/LatticeDescriptorCUB.hh>
#include <biu/LatticeDescriptorFCC.hh>
#include <biu/LatticeDescriptorCKW.hh>
#include <biu/LatticeModel.hh>

#include <biu/LatticeProteinUtil.hh>
#include <biu/SuperPos_Kabsch.hh>

#include "version.hh"


/* Implements and extends a greedy CRMSD fitting method as used in

	 @Article{ Miao.et.al:JMB:H-collapse:07,
	  author =   {Jiangbo Miao and Judith Klein-Seetharaman and Hagai Meirovitch},
	  title =    {The Optimal Fraction of Hydrophobic Residues Required to Ensure Protein Collapse},
	  journal =  {Journal of Molecular Biology},
	  year =     {2003},
	  volume =   {344},
	  pages =    {797-811},
	  doi  =     {doi:10.1016/j.jmb.2004.09.061}
	 }

 */

/* TODO
 *
 * - print neigh vectors of used lattices in verbose mode
 * - absolute move string output
 * - cRMSD : gapped sequences : selfavoidingness check in each subchain only currently!!!
 *
 */


//////////////////////////////////////////////////////////////////////////


#define PI 3.14159265

  //! maximal length of a chain gap that is bridged by a selfavoiding walk
#define MAX_SAW_LENGTH_TO_CONSIDER 5

enum OutMode { CML, XYZ, PDB };
enum Coordinates { X, Y, Z };

typedef std::vector< std::string > StrVec;

 // standard length of base vector
const double STANDARD_BASE_LENGTH = 3.8;

 // info for Jmol output
const std::string JMOL_SCRIPT = "select all; color bonds grey;";


  // the globally used rotation handler
biu::Rotator3D curRot;

biu::DblPoint operator*(const double d, const biu::IntPoint& p) {
	return biu::DblPoint(p)*d;
}

biu::DblPoint operator*(const biu::IntPoint& p, const double d) {
	return biu::DblPoint(p)*d;
}



//////////////////////////////////////////////////////////////////////////

/**
 * A container class to store information of a subchain
 */
class Subchain {
public:
	 //! index of first amino acid of this sub chain in original PDB file
	int idBegin;
	 //! which atom represents this sub chain
	std::string atom;
	 //! the coordinates of this sub chain
	biu::DPointVec points;
	 //! the sequence information of this subchain (3 letter code)
	StrVec seq;

	Subchain()
	 :	idBegin(1),
	 	atom(""),
	 	points(),
	 	seq()
	{}
};

//////////////////////////////////////////////////////////////////////////

/**
 * A container class to store information of a subchain of the lattice model.
 */
class SubchainLat {
public:
	 //! the coordinates of this sub chain
	biu::IPointVec points;
	 //! the absolute move string representation of the sub chain
	biu::MoveSequence moves;

	SubchainLat()
	 :	points(),
	 	moves()
	{}
};



//////////////////////////////////////////////////////////////////////////


 /*! Represents an amino acid within a protein chain 
  * and all assiciated data needed for fitting
  */
template < class POINT >
class AminoAcidT {
public:
	
	  //! the position with protein chain
	int id;
	  //! the amino acid name (3 letter code)
	std::string type;
	  //! position of the C_alpha atom
	POINT posCA;
	  //! position of the C_beta atom
	POINT posCB;
	  //! whether or not position of the C_beta atom is set
	bool posCBset;
	
	AminoAcidT()
	 :	id(-1)
	 	, type("")
	 	, posCA()
	 	, posCB()
	 	, posCBset(false)
	{}
};

typedef AminoAcidT< biu::DblPoint > AminoAcid_D;
typedef AminoAcidT< biu::IntPoint > AminoAcid_I;


//////////////////////////////////////////////////////////////////////////

  //! represents a protein chain to fit which might include chain gaps
typedef std::vector< AminoAcid_D > AAchain_D;
typedef std::vector< AminoAcid_D >::const_iterator AAcit_D;

  //! represents a protein chain to fit which might include chain gaps
typedef std::vector< AminoAcid_I > AAchain_I;
typedef std::vector< AminoAcid_I >::const_iterator AAcit_I;

//////////////////////////////////////////////////////////////////////////


  /*! Derives the AminoAcid vector for given sequence coordinate data that 
   * ignores chain gaps (just encoded within the amino acid positions)
   * 
   * @param chainCA the C_alpha atom chain information 
   * @param chainCB the C_beta atom chain information
   * 
   * @return the encoded vector of the according AminoAcid objects 
   */
AAchain_D
deriveAAvec(	const std::vector<Subchain> & chainCA
				, const std::vector<Subchain> * chainCB = NULL
			)
{
	if (chainCB != NULL) {
		assertbiu(chainCA.size()==chainCB->size(),"different length of chain information for CA and CB");
		for (size_t c=0; c<chainCA.size(); c++) {
			assertbiu(chainCA.at(c).points.size() == chainCB->at(c).points.size(),"some chain fragements differ in length");
		}
	}
	  // the data structure to fill
	AAchain_D aa;
	size_t posAA = 0;
	
	for (size_t c=0; c<chainCA.size(); c++) {
		  // add new amino acids
		aa.resize(aa.size()+chainCA.at(c).seq.size());
		  // set amino acid data
		for (size_t i=0; i<chainCA.at(c).points.size(); i++) {
			  // position in chain
			aa[posAA].id = chainCA.at(c).idBegin + int(i);
			  // amino acid type
			aa[posAA].type = chainCA.at(c).seq.at(i);
			  // C_alpha atom position
			aa[posAA].posCA = chainCA.at(c).points.at(i);
			  // C_beta atom position
			if (chainCB != NULL) {
				aa[posAA].posCB = chainCB->at(c).points.at(i);
				aa[posAA].posCBset = true;
			}
			  // increase access to writing position 
			posAA++;
		}
	}
	
	return aa;
}


//////////////////////////////////////////////////////////////////////////


  /*! Derives the AminoAcid vector for given sequence coordinate data that 
   * ignores chain gaps (just encoded within the amino acid positions)
   * 
   * @param chain the atom chain information 
   * @param sideChainPresent whether or not the chain data contains the 
   *          sidechain position after each backbone position 
   * 
   * @return the encoded vector of the according AminoAcid objects 
   */
AAchain_D
deriveAAvec(	const std::vector<Subchain> & chain
				, const bool sideChainPresent = false
			)
{
	for (size_t c=0; c<chain.size(); c++) {
		if (sideChainPresent) {
			assertbiu(chain.at(c).points.size()%2 == 0,"some chain fragements have incomplete backbone/sidechain data");
			assertbiu(chain.at(c).points.size()/2 == chain.at(c).seq.size(),"some chain fragements have different length of position and sequence");
		} else {
			assertbiu(chain.at(c).points.size() == chain.at(c).seq.size(),"some chain fragements have different length of position and sequence");
		}
	}
	
	  // the data structure to fill
	AAchain_D aa;
	size_t posAA = 0;
	
	for (size_t c=0; c<chain.size(); c++) {
		  // add new amino acids
		aa.resize(aa.size()+chain.at(c).seq.size());
		  // set amino acid data
		for (size_t i=0; i<chain.at(c).seq.size(); i++) {
			  // position in chain
			aa[posAA].id = chain.at(c).idBegin + int(i);
			  // amino acid type
			aa[posAA].type = chain.at(c).seq.at(i);
			  // atom position
			if (sideChainPresent) {
				  // C_alpha atom position
				aa[posAA].posCA = chain.at(c).points.at(i*2);
				  // C_beta atom position
				aa[posAA].posCB = chain.at(c).points.at((i*2)+1);
				aa[posAA].posCBset = true;
			} else {
				  // C_alpha atom position
				aa[posAA].posCA = chain.at(c).points.at(i);
				  // no sidechain position set
				aa[posAA].posCBset = false;
			}
			  // increase access to writing position 
			posAA++;
		}
	}
	
	return aa;
}


//////////////////////////////////////////////////////////////////////////


 /*! Checks whether or not the given position is part of the AminoAcidT iterator
  * range.
  * 
  * @param p the point to check for
  * @param begin the iterator pointing to the beginning of the range to check
  * @param end the iterator marking the end of the range to check (not checked)
  * 
  * @return true if p is NOT contained in the range; false otherwise 
  */
template <class CHAIN_IT>
bool
isFreePosition(const biu::IntPoint& p, CHAIN_IT begin, CHAIN_IT end ) {
	while (begin != end) {
		  // check backbone position
		if (begin->posCA == p) {
			return false;
		}
		  // check sidechain position if present
		if ( begin->posCBset && (begin->posCB == p)) {
			return false;
		}
		begin++;
	}
	return true;
}


//////////////////////////////////////////////////////////////////////////

 /*! Odering function for STL container or functions working on the first 
  * element of two given pairs only. It checks if the first element is smaller
  * than the one of the other pair.
  * 
  * @param p1 the first pair to compare 
  * @param p2 the second pair to compare 
  * 
  * @return (p1.first < p2.second)
  */
template <class PAIR>
bool
pairFirstLess(const PAIR& p1, const PAIR& p2 ) {
	return p1.first < p2.first;
}

 /*! Odering function for STL container or functions working on the first 
  * element of two given pairs only. It checks if the first element is equal to
  * the one of the other pair.
  * 
  * @param p1 the first pair to compare 
  * @param p2 the second pair to compare 
  * 
  * @return (p1.first == p2.second)
  */
template <class PAIR>
bool
pairFirstEqual(const PAIR& p1, const PAIR& p2 ) {
	return p1.first == p2.first;
}
#include <utility>

 /*! Odering function for STL container or functions working on the first 
  * element of two given pairs only. It checks if the first element is smaller
  * than the one of the other pair.
  */
template <class PAIR>
struct PairFirstLess  : public std::binary_function< PAIR, PAIR, bool> {
	 /*! Performs the check.
	  * 
	  * @param p1 the first pair to compare 
	  * @param p2 the second pair to compare 
	  * 
	  * @return (p1.first < p2.second)
	  */
	bool
	operator()(const PAIR& p1, const PAIR& p2 ) {
		return p1.first < p2.first;
	}
};

//#define isFreePosition isFreePositionT< AAchain_I::const_reverse_iterator >

//////////////////////////////////////////////////////////////////////////

#include <list>

//	typedef std::set< std::pair< double, AAchain_I >, PairFirstLess<std::pair< double, AAchain_I > > > STORAGE; 
//	typedef std::vector< std::pair< double, AAchain_I > > STORAGE; 
/*! containers to store extended structures in;
 * in each pair : FIRST = square deviation, SECOND = structure
 */
typedef std::list< std::pair< double, AAchain_I > > STORAGE; 
typedef STORAGE::value_type VT;
typedef STORAGE::iterator STit;

 /*! Inserts a new element to the storage. The insertion is done if
  * the storage size is less than best2store or the element to insert is smaller
  * than the maximal element in the storage.
  * 
  * @param toInsert the element to insert
  * @param newFits the storage to add to
  * @param best2store the maximal storage size allowed
  * 
  * @return true if the element was inserted; false otherwise
  */
bool 
insertToStorage( const VT & toInsert, STORAGE & newFits, const size_t best2store ) 
{
	  // copy structure and corresponding distance sum to newFits
	if (newFits.size() < best2store) {
		  // find location where it should be inserted
		STit insertPos = std::lower_bound(newFits.begin(), newFits.end(), toInsert, pairFirstLess<VT>);
		  // handle location possibilities
		if ( insertPos == newFits.end() ) {
			if (newFits.rbegin()->first < toInsert.first ) {
				  // add at the end
				newFits.push_back( toInsert );
				return true;
			}
		} else if ( insertPos == newFits.begin() ) {
			  // check if tailing position is not equal
			if ( toInsert.first < newFits.begin()->first ) {
				  // add to the front
				newFits.push_front( toInsert );
				return true;
			}
		} else {
			  // check if tailing position is not equal
			if ( toInsert.first < insertPos->first ) {
				insertPos--;
				  // check if leading position is really small and NOT equal
				if ( insertPos->first < toInsert.first ) {
					insertPos++;
					newFits.insert( insertPos, toInsert );
					return true;
				}
			}
		}
		
	} else if (toInsert.first < newFits.rbegin()->first) {

		  // find location where it should be inserted
		STit insertPos = std::lower_bound(newFits.begin(), newFits.end(), toInsert, pairFirstLess<VT>);
		
		if ( insertPos == newFits.begin() ) {
			  // check if tailing position is not equal
			if ( toInsert.first < insertPos->first) {
				  // insert at the front
				newFits.push_front( toInsert );
				  // erase last that had greater DSD than current
				newFits.pop_back();
				return true;
			}
		} else {
			  // check if tailing position is not equal
			if ( toInsert.first < insertPos->first ) {
				insertPos--;
				  // check if leading position is really small and NOT equal
				if ( insertPos->first < toInsert.first ) {
					insertPos++;
					newFits.insert( insertPos, toInsert );
					  // erase last that had greater DSD than current
					newFits.pop_back();
					return true;
				}
			}
		}
	}
	return false;
}


//////////////////////////////////////////////////////////////////////////


 /*! Recursive enumeration of all selfavoiding walks of a certain length. All
  * end points of the generated selfavoiding walks are stored in the provided
  * container.
  * 
  * @param saw the selfavoiding walk to be filled
  * @param curToSet the current position in saw that is to be set within this
  *              recursion iteration (if > saw.size() the recursion aborts)
  * @param sawEnds the container of all selfavoiding walk ends to be filled
  * @param lattice the lattice in which the selfavoiding walk is to be placed
  */
void
getSAWendsRecursion( 
			biu::IPointVec& saw,
			const size_t curToSet,
			biu::IPointSet& sawEnds,
			const biu::LatticeModel& lattice )
{
	if ( curToSet == saw.size() ) {
		  // add last position to result set
		sawEnds.insert( *(saw.rbegin()) );
		  // break recursion
		return;
	}
	  // check for each neighbor in the lattice
	for (biu::LatticeNeighborhood::const_iterator neigh = lattice.getNeighborhood().begin();
			neigh != lattice.getNeighborhood().end(); neigh++)
	{
		  // calculate next position to try
		saw[curToSet] = saw[curToSet-1] + (*neigh);
		bool isSA = true;
		  // check SAW condition
		for (size_t i=curToSet-1; isSA && i>0; i--) {
			isSA = (saw[i-1] != saw[curToSet]);
		}
		  // if selfavoiding proceed recursion
		if (isSA) {
			  // recursive call
			getSAWendsRecursion( saw, curToSet+1, sawEnds, lattice );
		}
	}	
}


//////////////////////////////////////////////////////////////////////////


#include <map>

typedef std::map< size_t, biu::IPointSet > SawEndStorage; 

 /*! Container to store the selfavoiding walk ends of the required lengths to
  * enable a fast lookup and to avoid unneccessary recomputation.
  */
SawEndStorage sawEndStorage;

 /*! Access to the end positions of selfavoiding walks of a given length.
  * 
  * @param sawLength the length of the selfavoiding walks required
  * @param lattice the lattice in which the selfavoiding walk is to be placed
  * 
  * @return the set of selfavoiding walk end points of this length originating
  *         in (0,0,0)
  */
const biu::IPointSet&
getSAWends( 
			const size_t sawLength,
			const biu::LatticeModel& lattice )
{
	SawEndStorage::iterator sawEnds = sawEndStorage.find( sawLength );
	  // check if already computed
	if (sawEnds != sawEndStorage.end()) {
		  // access to SAW ends
		return sawEnds->second;
	}
	
	  // selfavoiding walk to fill
	biu::IPointVec saw(sawLength+1);
	  // set start of SAW
	saw[0] = biu::IntPoint(0,0,0);
	
	biu::IPointSet curSawEnds;
	  // call recursive SAW enumeration
	getSAWendsRecursion( saw, 1, curSawEnds, lattice );
	
	  // add to storage
	sawEnds = sawEndStorage.insert(SawEndStorage::value_type(sawLength, curSawEnds)).first;
	  // access to SAW ends
	return sawEnds->second;
}


//////////////////////////////////////////////////////////////////////////


 /*! Calculates the best lattice fit of a given amino acid chain within a 
  * lattice. The fitting is done utilizing a greedy dRMSD optimizing approach.
  * 
  * @param orig the original amino acid data to fit
  * @param lp the lattice protein approximation to fill
  * @param fitSideChain whether or not backbone and sidechain monomer should be
  *            placed (true) or if a backbone-only model should be derived 
  *            (false)
  * @param lattice the lattice in which the lattice model is to be placed
  * @param best2store holds how many best structures should be maintained along
  *            the greedy elongation for the next iteration
  * @param baseLength the factor to be multiplied with lattice vectors to derive
  *            well scaled distances in real space in accordance to the original
  *            data
  * @param verbose if true a progress bar output is done to std::cout
  * 
  * @return the final dRMSD of the returned lp and orig chains
  */
double
fit2lattice_dRMSD(	
				const AAchain_D& orig,
				AAchain_I& lp,
				const bool fitSideChain,
				const biu::LatticeModel& lattice,
				const size_t best2store = 1,
				const double baseLength = STANDARD_BASE_LENGTH,
				const bool verbose = false )
{
	assertbiu(orig.size() > 1, "given original data has to be at least 2");
	
	STORAGE fits1, fits2;
	  // access to the storage containers
	STORAGE& bestFits = fits1, newFits = fits2, tmpFits = fits2;
	
	  // access to neighboring vectors
	const biu::LatticeNeighborhood& neighbors = lattice.getNeighborhood();
	biu::LatticeNeighborhood::const_iterator neigh;

	  // number of summation events used to get the square deviation
	double curDiff = 0.0;

	  // initialise bestFits with a single uninitialized amino acid
//	bestFits.insert( VT(0.0,AAchain_I(1)) );
	bestFits.push_back( VT(0.0,AAchain_I(1)) );

	  // set first monomer
	size_t i=0;
	bestFits.begin()->second[i].id = orig[i].id;
	bestFits.begin()->second[i].type = orig[i].type;
	bestFits.begin()->second[i].posCA = biu::IntPoint(0,0,0);
	
	if (fitSideChain) {
		  // set sidechain position
		bestFits.begin()->second[i].posCB = bestFits.begin()->second[i].posCA + *(neighbors.begin());
		  // tag that side chain is present
		bestFits.begin()->second[i].posCBset = true;
		  // update square deviation sum for first BB-SC link
		curDiff = ((biu::DblPoint(bestFits.begin()->second[i].posCA)*baseLength).distance(biu::DblPoint(bestFits.begin()->second[i].posCB)*baseLength)) 
				- (orig.at(i).posCA.distance(orig.at(i).posCB));
		bestFits.begin()->first += pow (curDiff, 2);
	} else
	  // (no sidechain) set second backbone monomer if successive IDs in original chain
	if ((orig[i].id +1)==orig[i+1].id) {
		i++;
		  // resize structure
		bestFits.begin()->second.resize(i+1);
		  // copy amino acid data
		bestFits.begin()->second[i].id = orig[i].id;
		bestFits.begin()->second[i].type = orig[i].type;
		  // set second to an arbitrary neighbored position 
		  // (since we are optimizing dRMSD)
		bestFits.begin()->second[i].posCA = bestFits.begin()->second[i-1].posCA + *(neighbors.begin());
		  // update square deviation sum for first link
		curDiff = ((biu::DblPoint(bestFits.begin()->second[i-1].posCA)*baseLength).distance(biu::DblPoint(bestFits.begin()->second[i].posCA)*baseLength)) 
				- (orig.at(i-1).posCA.distance(orig.at(i).posCA));
		bestFits.begin()->first += pow (curDiff, 2);
	}
	
	  // width of the progress bar to print
	const size_t PROGRESS_BAR_WIDTH = 64;
	  // size of the already printed progress bar
	size_t progressBarSize = 0;

	if (verbose) {
		  // print progress bar limits
		std::cout 	<<"\n"
					<<"  Fitting PDB structure onto lattice :\n"
					<<"\n"
					<<"    |0%" <<std::setw(PROGRESS_BAR_WIDTH-3) <<100 <<"%|\n"
					<<"    |";
		std::cout.flush();
	}


	
	// do chain elongations
	for (i++; i<orig.size(); i++) {
		  // clear container where new structures are added to
		newFits.clear();
		
		  // extend all of best of last round
		for (STORAGE::const_iterator it = bestFits.begin();
				it != bestFits.end(); it++)
		{
			  // copy and extend current best to extend
			AAchain_I next = it->second;
			assertbiu( next.size() == i, "chain length differs from expectation...");
			
			  // add new amino acid 
			next.resize(i+1);
			  // copy amino acid data
			next[i].id = orig[i].id;
			next[i].type = orig[i].type;
			
				
			  // get iterator to position INFRONT of last position
			AAchain_I::const_reverse_iterator nextRbegin = next.rbegin(); 
			nextRbegin++;
			AAchain_I::const_reverse_iterator nextRend = next.rend(); 

            // Philippe: here, could open a new chain
            bool detectedJump = false;
            double distance = orig.at(i).posCA.distance(orig.at(i-1).posCA);
            if(distance > 10) {
                 detectedJump = true;
                 //std::cout << "Got distance " << distance << std::endl;
            }

            // check if no chain gap present
            if ( (!detectedJump) && (next[i].id == (it->second.at(i-1).id +1)) ) {
				
				  // check for each neighbor in the lattice
				for (biu::LatticeNeighborhood::const_iterator neigh = lattice.getNeighborhood().begin();
						neigh != lattice.getNeighborhood().end(); neigh++)
				{
					  // SET BACKBONE POSITION
					  // set new position of the monomer to extend
					next[i].posCA = next[i-1].posCA + (*neigh);
					  // check if selfavoiding
					if ( isFreePosition( next[i].posCA, nextRbegin, nextRend) ) 
					{
						if (fitSideChain) {
							  // SET SIDECHAIN POSITION
							  // check for each neighbor in the lattice
							for (biu::LatticeNeighborhood::const_iterator neighSC = lattice.getNeighborhood().begin();
									neighSC != lattice.getNeighborhood().end(); neighSC++)
							{
								  // set new position of the monomer to extend
								next[i].posCB = next[i].posCA + (*neighSC);
								next[i].posCBset = true;
								  // check if selfavoiding
								if ( isFreePosition( next[i].posCB, nextRbegin, nextRend) ) 
								{
									  // update dRMSD according to added position
									double curDSD = it->first; // == current distance square deviation
									double curDiff = 0.0; 
									  // CA-CB distance of this amino acid
									curDiff = ((biu::DblPoint(next[i].posCA)*baseLength).distance(biu::DblPoint(next[i].posCB)*baseLength)) 
											- (orig.at(i).posCA.distance(orig.at(i).posCB));
									curDSD += pow (curDiff, 2);
									  // distances to all other amino acids
									for (size_t j=0; j<i; j++) {
										curDiff = ((biu::DblPoint(next[i].posCA)*baseLength).distance(biu::DblPoint(next[j].posCA)*baseLength)) 
												- (orig.at(i).posCA.distance(orig.at(j).posCA));
										curDSD += pow (curDiff, 2);
										curDiff = ((biu::DblPoint(next[i].posCA)*baseLength).distance(biu::DblPoint(next[j].posCB)*baseLength)) 
												- (orig.at(i).posCA.distance(orig.at(j).posCB));
										curDSD += pow (curDiff, 2);
										curDiff = ((biu::DblPoint(next[i].posCB)*baseLength).distance(biu::DblPoint(next[j].posCB)*baseLength)) 
												- (orig.at(i).posCB.distance(orig.at(j).posCB));
										curDSD += pow (curDiff, 2);
										curDiff = ((biu::DblPoint(next[i].posCB)*baseLength).distance(biu::DblPoint(next[j].posCA)*baseLength)) 
												- (orig.at(i).posCB.distance(orig.at(j).posCA));
										curDSD += pow (curDiff, 2);
									}
									  // prepare data object to insert
									VT toInsert( curDSD, next );
									  // add to storage
									insertToStorage( toInsert, newFits, best2store );
								}
							}
						}  // end sidechain handling
						  // HANDLE BACKBONE ONLY CASE
						else {
							  // update dRMSD according to added position
							double curDSD = it->first;
							double curDiff = 0.0; 
							for (size_t j=0; j<i; j++) {
								curDiff = ((biu::DblPoint(next[i].posCA)*baseLength).distance(biu::DblPoint(next[j].posCA)*baseLength)) 
										- (orig.at(i).posCA.distance(orig.at(j).posCA));
								curDSD += pow (curDiff, 2);
							}
							  // prepare data object to insert
							VT toInsert( curDSD, next );
							  // add to storage
							insertToStorage( toInsert, newFits, best2store );
						}  // end backbone only handling
					}
				} // end try all neighbored positions
			} // end of no chain gap
			
			  // ELSE : CHAIN GAP PRESENT
			else {
				  // the length of the gap to bridge
				const size_t gapLength = next[i].id - it->second.at(i-1).id +1;
				
				  // container that can be filled to hold the new fragment start
				  // positions to try
				biu::IPointSet fragmentStarts;
				  // access to the new fragment starts to try on
				const biu::IPointSet* fragStarts = &fragmentStarts;

                  // decide on gap handling depending on gap length // Philippe modified here
                if ((!detectedJump) && (gapLength <= MAX_SAW_LENGTH_TO_CONSIDER)) { // short gap handling
					
					  // selfavoiding walk ends of given length to bridge chain gap
					fragStarts = &getSAWends( gapLength, lattice );
					assertbiu( fragStarts->size() > 0, "no SAW ends computed");
	
				
				} // end short gap handling
				else { // long gap handling
					
					  // copy position data of fit and original data so far processed
					biu::DPointVec soFarP((i-1)*(fitSideChain?2:1));
					biu::DPointVec soFarP_((i-1)*(fitSideChain?2:1));
					biu::DPointVec soFarL((i-1)*(fitSideChain?2:1));
					bool atLeastOneXnonZero = false;
					for (size_t j=0; j<(i-1); j++) {
						  // copy C_alpha coordinates
						soFarP[j] = orig.at(j).posCA;
						soFarP_[j] = soFarP.at(j);
						atLeastOneXnonZero = atLeastOneXnonZero || (soFarP.at(j).getX() != 0.0); 
						soFarL[j] = next.at(j).posCA * baseLength;
						if (fitSideChain) { 
							  // copy side chain coordinates
							soFarP[j+i-1] = orig.at(j).posCB;
							soFarP_[j+i-1] = soFarP.at(j+i-1);
							atLeastOneXnonZero = atLeastOneXnonZero || (soFarP.at(j+i-1).getX() != 0.0); 
							soFarL[j+i-1] = next.at(j).posCB * baseLength;
						}
					}

					  // create reflection of lattice data
					if (atLeastOneXnonZero) {
						  // flip X value
						for (size_t j=0; j<soFarL.size(); j++)
							soFarP_[j].setX( -(soFarP_.at(j).getX()) );
					} else {
						  // flip Y value
						for (size_t j=0; j<soFarL.size(); j++)
							soFarP_[j].setY( -(soFarP_.at(j).getY()) );
					}
					
					  // get best superpositioning of chains so far
					const biu::SuperPos_Kabsch::Transformation trans
						= biu::SuperPos_Kabsch::superposition( soFarP, soFarL, false );
					  // get best superpositioning of chains so far reflected
					const biu::SuperPos_Kabsch::Transformation trans_
						= biu::SuperPos_Kabsch::superposition( soFarP_, soFarL, false );

					  // calculate new point to fit according to best 
					  // superpositioning
					biu::DblPoint newCAi;
					if (	biu::LatticeProteinUtil::cRMSD( soFarP, soFarL )
						< biu::LatticeProteinUtil::cRMSD( soFarP_, soFarL ) )
					{
						  // calculate the C_alpha point to fit according to 
						  // the superpositioning of the already fitted chain
                        // Philippe: I guess they translate the previous chain's end (the +1 AA) to join the new position.
						newCAi = trans.translation + (biu::SuperPos_Kabsch::rotate(trans.rotation, orig.at(i).posCA)*trans.scale);
						  // substract end to get relative distance vector
						newCAi -= soFarP.at(i-2);
					} else {
						  // calculate the C_alpha point to fit according to 
						  // the superpositioning of the reflected fitted chain
						biu::DblPoint origAtI = orig.at(i).posCA;
						if (atLeastOneXnonZero) {
							origAtI.setX(-(origAtI.getX()));
						} else {
							origAtI.setY(-(origAtI.getY()));
						}
						newCAi = trans_.translation + (biu::SuperPos_Kabsch::rotate(trans_.rotation, origAtI)*trans_.scale);
						  // substract end to get relative distance vector
						newCAi -= soFarP_.at(i-2);
					}
					  // scale to lattice length
					newCAi /= baseLength;
					
					  // fill possible fragment starts (C_alpha positions)
					  // RELATIVE to the end of the current chain
					
                    // all possible combinations of floor and ceil
					const biu::LatticeDescriptor * const ld = lattice.getDescriptor();
					biu::IntPoint p = biu::IntPoint((int)floor(newCAi.getX()), (int)floor(newCAi.getY()), (int)floor(newCAi.getZ()));
						if (ld->isLatticeNode(p)) {fragmentStarts.insert(p);}
					p = biu::IntPoint((int)ceil(newCAi.getX()), (int)floor(newCAi.getY()), (int)floor(newCAi.getZ()));
						if (ld->isLatticeNode(p)) {fragmentStarts.insert(p);}
					p = biu::IntPoint((int)floor(newCAi.getX()), (int)ceil(newCAi.getY()), (int)floor(newCAi.getZ()));
						if (ld->isLatticeNode(p)) {fragmentStarts.insert(p);}
					p = biu::IntPoint((int)ceil(newCAi.getX()), (int)ceil(newCAi.getY()), (int)floor(newCAi.getZ()));
						if (ld->isLatticeNode(p)) {fragmentStarts.insert(p);}
					p = biu::IntPoint((int)floor(newCAi.getX()), (int)floor(newCAi.getY()), (int)ceil(newCAi.getZ()));
						if (ld->isLatticeNode(p)) {fragmentStarts.insert(p);}
					p = biu::IntPoint((int)ceil(newCAi.getX()), (int)floor(newCAi.getY()), (int)ceil(newCAi.getZ()));
						if (ld->isLatticeNode(p)) {fragmentStarts.insert(p);}
					p = biu::IntPoint((int)floor(newCAi.getX()), (int)ceil(newCAi.getY()), (int)ceil(newCAi.getZ()));
						if (ld->isLatticeNode(p)) {fragmentStarts.insert(p);}
					p = biu::IntPoint((int)ceil(newCAi.getX()), (int)ceil(newCAi.getY()), (int)ceil(newCAi.getZ()));
						if (ld->isLatticeNode(p)) {fragmentStarts.insert(p);}

					  // generate all neighbors of these rounded ends and add
					  // them to possible fragment start positions
					biu::IPointSet tmp = fragmentStarts;
					for (biu::IPointSet::const_iterator roundedStart = tmp.begin();
							roundedStart != tmp.end(); roundedStart++) {
						  // add each neighbor of the fragment start in the lattice
						for (biu::LatticeNeighborhood::const_iterator neigh = lattice.getNeighborhood().begin();
								neigh != lattice.getNeighborhood().end(); neigh++)
						{
							fragmentStarts.insert( (*roundedStart) + (*neigh));
						}	
					}
					
				} // end long gap handling
				
				assertbiu( fragStarts != NULL, "fragment starts container not set (==NULL)");
				assertbiu( fragStarts->size() != 0, "no possible fragment starts determined");
				
				  // try all fragment starts
				for (biu::IPointSet::const_iterator fStart = fragStarts->begin(); fStart != fragStarts->end(); fStart++) {
					  // set possible next fragment start
					next[i].posCA = next.at(i-1).posCA + (*fStart);
					  // check if selfavoiding
					if ( isFreePosition( next[i].posCA, nextRbegin, nextRend) ) 
					{
						if (fitSideChain) {
							  // SET SIDECHAIN POSITION
							  // check for each neighbor in the lattice
							for (biu::LatticeNeighborhood::const_iterator neighSC = lattice.getNeighborhood().begin();
									neighSC != lattice.getNeighborhood().end(); neighSC++)
							{
								  // set new position of the monomer to extend
								next[i].posCB = next[i].posCA + (*neighSC);
								next[i].posCBset = true;
								  // check if selfavoiding
								if ( isFreePosition( next[i].posCB, nextRbegin, nextRend) ) 
								{
									  // update dRMSD according to added position
									double curDSD = it->first;
									double curDiff = 0.0; 
									  // CA-CB distance of this amino acid
									curDiff = ((biu::DblPoint(next[i].posCA)*baseLength).distance(biu::DblPoint(next[i].posCB)*baseLength)) 
											- (orig.at(i).posCA.distance(orig.at(i).posCB));
									curDSD += pow (curDiff, 2);
									  // distances to all other amino acids
									for (size_t j=0; j<i; j++) {
										curDiff = ((biu::DblPoint(next[i].posCA)*baseLength).distance(biu::DblPoint(next[j].posCA)*baseLength)) 
												- (orig.at(i).posCA.distance(orig.at(j).posCA));
										curDSD += pow (curDiff, 2);
										curDiff = ((biu::DblPoint(next[i].posCA)*baseLength).distance(biu::DblPoint(next[j].posCB)*baseLength)) 
												- (orig.at(i).posCA.distance(orig.at(j).posCB));
										curDSD += pow (curDiff, 2);
										curDiff = ((biu::DblPoint(next[i].posCB)*baseLength).distance(biu::DblPoint(next[j].posCB)*baseLength)) 
												- (orig.at(i).posCB.distance(orig.at(j).posCB));
										curDSD += pow (curDiff, 2);
										curDiff = ((biu::DblPoint(next[i].posCB)*baseLength).distance(biu::DblPoint(next[j].posCA)*baseLength)) 
												- (orig.at(i).posCB.distance(orig.at(j).posCA));
										curDSD += pow (curDiff, 2);
									}
									  // prepare data object to insert
									VT toInsert( curDSD, next );
									  // add to storage
									insertToStorage( toInsert, newFits, best2store );
								}
							}
						}  // end sidechain handling
						  // HANDLE BACKBONE-ONLY CASE
						else {
							  // update dRMSD according to added position
							double curDSD = it->first;
							double curDiff = 0.0; 
							for (size_t j=0; j<i; j++) {
								curDiff = ((biu::DblPoint(next[i].posCA)*baseLength).distance(biu::DblPoint(next[j].posCA)*baseLength)) 
										- (orig.at(i).posCA.distance(orig.at(j).posCA));
								curDSD += pow (curDiff, 2);
							}
							  // prepare data object to insert
							VT toInsert( curDSD, next );
							  // add to storage
							insertToStorage( toInsert, newFits, best2store );
						} // end backbone-only handling
					}
				} // end try all SAW ends for fragment start
			} // END CHAIN GAP HANDLING
		}

		  // make the new fits the best ones for next iteration
		tmpFits = bestFits;
		bestFits = newFits;
		newFits = tmpFits;
		
		if (verbose) {
			  // do progress bar output
			  // calculate necessary progress bar extension
			size_t progExt = ((i * PROGRESS_BAR_WIDTH) / (orig.size()-1)) - progressBarSize;
			if (progExt > 0) {
				progressBarSize += progExt;
				  // print progress bar elements 
				while ( progExt-- ) {
					std::cout	<<'#';
				}
				std::cout.flush();
			}
		}

		  // check if at least one structure was extensible
		if (bestFits.size() == 0) {
			if (verbose) {
				std::cout <<"...\n" <<std::endl;
			}
			return -1;
		}
	}
	
	if (verbose) {
		while(progressBarSize < PROGRESS_BAR_WIDTH) {
			std::cout 	<<'#';
		}
		std::cout 	<<"|\n" <<std::endl;
	}
	
	  // copy best structure found
	lp = bestFits.begin()->second;
	  // calculate and return dRMSD of that structure
	const double n = bestFits.begin()->second.size() * ( fitSideChain ? 2 : 1 );
	const double sumCount = ((n*(n-1))/2.0);
	return sqrt( bestFits.begin()->first / sumCount);
}



//////////////////////////////////////////////////////////////////////////

 /*! Calculates the centroid position of a given set of coordinates, i.e. the
  * mean coordinate of all given ones.
  * 
  * @param sc the set of coordinates
  * 
  * @return the centroid of the given set of points
  */
biu::DblPoint
getCentroid( const biu::DPointVec& sc ) {
	  // if no coordinates are given return cAlpha atom
	  // calculate center of mass
	biu::DblPoint com(0.0,0.0,0.0);
    if(sc.size() == 0) {
        std::cerr << "There is probably a variant conformation (atoms name +A/B/C at the end). Use -pdbAtomAlt=X where X is the character of the conformer.";
        //return com;
    }
    assertbiu(sc.size() != 0, "no coordinate data given to derive centroid from");

	  // calc mean value of all coordinates
	for (biu::DPointVec::const_iterator i = sc.begin(); i!=sc.end(); i++) {
		com += *i;
	}
	com.setX( com.getX() / (double)sc.size() );
	com.setY( com.getY() / (double)sc.size() );
	com.setZ( com.getZ() / (double)sc.size() );
	  // return center of mass
	return com;
}

// Philippe: added the option to include the backbone atoms in the mass center calculation
int parsePDB_CoM(	std::istream& in, const char chainID, const char atomAlt,
				std::vector<Subchain> & chain, std::string& pdbID,
                const bool chainGaps, const int pdbModel, bool includeBackboneInMassCenter = false)
{

    // Philippe: Idea: browse line by line and build the sideChain. seqNum is the ID of the residue. When seqNum changes, i.e. completion of the residue, calls getCentroid to get the center.
    // note that the Calpha and chain are stored differently, because the centroid means only the side chain, not the backbone
	using namespace std;

	  // clean data storage
	chain.clear();

	  // add first element to fill
	chain.push_back(Subchain());

	biu::DPointVec sideChain;
	biu::DblPoint cAlpha(0.0,0.0,0.0);
	std::string aaName = "???";

	  // temp data structures
	char cline[201];
	string line;

	size_t lineCount = 0;
	int seqNumLast = INT_MAX;
	bool cAlphaFound = false;
	bool pdbModelFound = false;
	bool pdbModelPresent = false;  // assume no model information present
	  // read istream line by line
	while ( in.good() ) {
		lineCount++;
		  // get next line content
		in.getline( cline, 200, '\n' );
		line = string(cline);
		  // line too short
		if ( line.size() < 6)
			continue;

		  // get PDB ID out of HEADER line
		if (	line.compare( 0, 6, "HEADER") == 0 && line.size() > 66 ) {
			pdbID = line.substr(62,4);
		}
		
		  // a model was found !!!
		  // check if current model is the one to screen 
		if (	line.compare( 0, 5, "MODEL") == 0 ) {
			pdbModelPresent = true;
			  // check if this model is the one to screen
			pdbModelFound = atoi( line.substr( 6, 20).c_str() ) == pdbModel;
		}
		
		  // in case models are present:
		  // check if we found already the model to screen --> otherwise keep on reading
		if ( pdbModelPresent && !pdbModelFound ) {
			continue;
		}
		
		  // in case models are present:
		  // when here than we are already screening a model : check for end and abort reading
		if ( pdbModelPresent && line.compare( 0, 6, "ENDMDL") == 0 ) {
			  // end of model to read was reached : end loop
			break;
		}



		  // handle ATOM line of correct chain ID
		if (	line.compare( 0, 6, "ATOM  ") == 0
			 &&	line.at(21) == chainID
			 && (line.at(16) == ' ' || line.at(16) == atomAlt) // atom alternative ok
			)
		{
			  // handle alternative atom location
//			if ( line.at(16) != ' ' ) {
//				std::cerr <<"\n   PDB : line " <<lineCount <<" : alternate location indicator found --> not handled yet\n\n"
//							<<line <<"\n\n";
//				return -1;
//			}
			  // handle residue insertion
			if ( line.at(26) != ' ' ) {
				std::cerr <<"\n   PDB : line " <<lineCount <<" : residue insertion found --> not handled yet\n\n"
							<<line <<"\n\n";
				return -1;
			}

			  // handle and check sequence number
			int seqNum = atoi(line.substr(22,4).c_str());
			if (seqNumLast == INT_MAX) {
			  // handle first atom of the chain
				sideChain.clear();
				cAlphaFound = false;
				seqNumLast = seqNum;
				aaName = line.substr(17,3);
				cAlpha = biu::DblPoint(0.0,0.0,0.0);
				chain.rbegin()->idBegin = seqNum;
				chain.rbegin()->atom = "CoM";
			} else if (seqNum == seqNumLast+1) {
              // next amino acid reached / previous one is completed
				  // check if the cAlpha atom was found for the last amino acid
				if (cAlphaFound) {
					if (sideChain.size() == 0 && aaName == "GLY") {
						  // no side chain data available for glycin
						  // --> store backbone position again
						chain.rbegin()->points.push_back(cAlpha);
					} else  {
						  // data ok :
						  // calculate center of mass and add to points
//                        if(sideChain.size() == 0) {
//                            cerr << "Residue(completed) " <<  aaName << " seqNum " << seqNum-1;
//                            cerr << "Case 1" << endl;
//                        }
						chain.rbegin()->points.push_back(getCentroid(sideChain));
					}
					  // add amino acid name
					chain.rbegin()->seq.push_back(aaName);
				} else {
					  // no cAlpha atom found for last amino acid :
					  // check if the incomplete information is the start of this chain
					if ( chain.rbegin()->idBegin == seqNumLast) {
						  // update sequence start to this amino acid
						chain.rbegin()->idBegin = seqNum;
					} else if (chainGaps) {
						// handle new sub chain
						  // add next subchain to fill
						chain.push_back(Subchain());
						  // fill new subchain element
						chain.rbegin()->idBegin = seqNum;
						chain.rbegin()->atom = "CoM";
					} else {
						std::cerr <<"\n   PDB : line " <<lineCount 
								<<" : no CA coordinates for previous amino acid "<<seqNumLast 
								<<" and no chain gaps allowed\n\n"
									<<line <<"\n\n";
						return -1;
					}
				}
				  // clear side chain container
				sideChain.clear();
				seqNumLast = seqNum;
				aaName = line.substr(17,3);
				cAlpha = biu::DblPoint(0.0,0.0,0.0);
				cAlphaFound = false;
			} else if (seqNum != seqNumLast) {
			  // non-successive identifier found
				  // handle error case, if no sub-chain support
				if (chainGaps) {
					// handle subchain support --> end of subchain reached
					// handle last sub chain
					  // calculate center of mass and add to points
                    //cerr << "Case 2" << endl;
					chain.rbegin()->points.push_back(getCentroid(sideChain));
					  // add amino acid name
					chain.rbegin()->seq.push_back(aaName);
//std::cerr <<" read subchain : fistToRead=" <<chain.rbegin()->idBegin <<" : last=" <<chain.rbegin()->idBegin+chain.rbegin()->seq.size()-1 <<" length=" <<chain.rbegin()->seq.size()<<"\n";
					// handle new sub chain
					  // add next subchain to fill
					chain.push_back(Subchain());
					  // update begin of chain
					seqNumLast = seqNum;
					sideChain.clear();
					aaName = line.substr(17,3);
					cAlpha = biu::DblPoint(0.0,0.0,0.0);
					cAlphaFound = false;
					  // fill new subchain element
					chain.rbegin()->idBegin = seqNum;
					chain.rbegin()->atom = "CoM";
				} else {
					std::cerr <<"\n   PDB : line " <<lineCount <<" : next sequence number not consecutive and no chain gaps allowed !\n\n"
								<<line <<"\n\n";
					return -1;
				}
			}

			  // get cAlpha atom
			if (line.compare( 13, 3, "CA ") == 0) {
				  // coordinates
				double x = atof(line.substr(30,8).c_str());
				double y = atof(line.substr(38,8).c_str());
				double z = atof(line.substr(46,8).c_str());
				cAlpha = biu::DblPoint(x,y,z);
				cAlphaFound = true;
			}

			  // only add atoms not belonging to the backbone to the side chain
            if (includeBackboneInMassCenter ||
                      ( line.compare(13, 3, "CA ") != 0
                    &&	line.compare(13, 3, "N  ") != 0
                    &&	line.compare(13, 3, "C  ") != 0
                    &&	line.compare(13, 3, "O  ") != 0)
                )
			{
				  // coordinates
				double x = atof(line.substr(30,8).c_str());
				double y = atof(line.substr(38,8).c_str());
				double z = atof(line.substr(46,8).c_str());
                //cerr << "Residue " <<  aaName << " seqNum " << seqNum << " got atom " << line.substr(13,3) << endl;
				sideChain.push_back(biu::DblPoint(x,y,z));
			}

		} // if ATOM and chainID
	} // while
	
	  // check if a model was present and we have found the one to screen
	if (pdbModelPresent && !pdbModelFound) {
		std::cerr <<"\n   PDB : line " <<lineCount <<" : models are present but we have not found the model with id "
				<<pdbModel <<" !\n\n"
					<<"\n\n";
		return -1;
	}

	  // check for and add center of mass of last chain element if cAlpha coordinates present
	if (cAlphaFound && aaName.compare("???") != 0) {
		  // calculate center of mass and add to points
        //cerr << "Case 3" << endl;
		chain.rbegin()->points.push_back(getCentroid(sideChain));
		  // add amino acid name
		chain.rbegin()->seq.push_back(aaName);
	}

//std::cerr <<" read subchain : fistToRead=" <<chain.rbegin()->idBegin <<" : last=" <<chain.rbegin()->idBegin+chain.rbegin()->seq.size()-1 <<" length=" <<chain.rbegin()->seq.size()<<"\n";
	if (seqNumLast == INT_MAX) {
		std::cerr <<"\n   PDB : line " <<lineCount <<" : up to here no center of mass could be computed !\n\n"
					<<"\n\n";
		return -1;
	} else {
		return seqNumLast+(int)chain.rbegin()->seq.size()-1;
	}
}


//////////////////////////////////////////////////////////////////////////


/**
 * Reads a stream in PDB format and parses the given atoms and amino acid
 * sequence.
 *
 * @param in the stream to read from
 * @param chainID the chain ID of the atoms to read
 * @param atom the atom type to read
 * @param atomAlt the atom alternative identifier to use if one is present
 * @param chain OUT-parameter: the read 3D-coordinates and sequence of the atoms
 *        of each subchain read
 * @param pdbID OUT-parameter: the PDB ID of the file processed
 * @param chainGaps
 *        IF equal to FALSE, a gap in the numbering will lead to an error.
 *        IF equal to TRUE, a set of consecutive subchains is created.
 *
 * @return the last amino acid number in the file read or a negative error code
 *          0 : no amino acid atom was read from string, nothing found for the
 *              given arguments
 *         -1 : read error
 *
 */
int parsePDB(	std::istream& in, const char chainID, const std::string &atom,
				const char atomAlt,
				std::vector< Subchain > & chain, std::string& pdbID,
                const bool chainGaps, const int pdbModel, double distanceDisrupt = +1e16)
{
	if (atom.compare("CoM") == 0) {
		return parsePDB_CoM( in, chainID, atomAlt, chain, pdbID, chainGaps, pdbModel);
	}    
    if (atom.compare("FuC") == 0) {
        return parsePDB_CoM( in, chainID, atomAlt, chain, pdbID, chainGaps, pdbModel, true);
    }
	using namespace std;
    double lastx = NAN;
    double lasty = NAN;
    double lastz = NAN;

	  // clean data storage
	chain.clear();
	pdbID = "1???";

	  // add new element to fill by read
	chain.push_back(Subchain());

	  // temp data structures
	char cline[201];
	string line;

	size_t lineCount = 0;
	int seqNumFirst = INT_MAX;
	bool pdbModelFound = false;
	bool pdbModelPresent = false;  // assume no model information present
	  // read istream line by line
	while ( in.good() ) {
		lineCount++;
		  // get next line content
		in.getline( cline, 200, '\n' );
		line = string(cline);
		  // line too short
		if ( line.size() < 6)
			continue;

		  // get PDB ID out of HEADER line
		if (	line.compare( 0, 6, "HEADER") == 0 && line.size() > 66 ) {
			pdbID = line.substr(62,4);
		}

		
		  // a model was found !!!
		  // check if current model is the one to screen 
		if (	line.compare( 0, 5, "MODEL") == 0 ) {
			pdbModelPresent = true;
			  // check if this model is the one to screen
			pdbModelFound = atoi( line.substr( 6, 20).c_str() ) == pdbModel;
		}
		
		  // in case models are present:
		  // check if we found already the model to screen --> otherwise keep on reading
		if ( pdbModelPresent && !pdbModelFound ) {
			continue;
		}
		
		  // in case models are present:
		  // when here than we are already screening a model : check for end and abort reading
		if ( pdbModelPresent && line.compare( 0, 6, "ENDMDL") == 0 ) {
			  // end of model to read was reached : end loop
			break;
		}

		  // handle ATOM line of CA atoms of correct chain ID
		if (	line.compare( 0, 6, "ATOM  ") == 0
			&&	line.at(21) == chainID )
		{
			  // handle residue insertion
			if ( line.at(26) != ' ' ) {
				std::cerr <<"\n   PDB : line " <<lineCount <<" : residue insertion found --> not handled yet\n\n"
							<<line <<"\n\n";
				return -1;
			}

			  // check if correct line found
			  // check if alternative identifier is wanted if present
			  // take CA atom instead CB for glycin
			if ( 	(line.at(16) != ' ' && line.at(16) != atomAlt) // atom alternative not wanted
				|| ( line.compare( 13, atom.size(), atom) != 0  // wrong atom
					&& (line.substr(17,3).compare("GLY")!=0 || atom.compare("CB")!=0  || line.compare( 13, 2, "CA") != 0 )
					)
				)
			{
				continue;
			}

            // Philippe 2019-10-04 Moved the coordinates up so we can initiate new chain if atoms are too far away
            // store data
            // coordinates
            double x = atof(line.substr(30,8).c_str());
            double y = atof(line.substr(38,8).c_str());
            double z = atof(line.substr(46,8).c_str());
            double distance = sqrt((x - lastx)*(x - lastx) + (y - lasty)*(y - lasty) + (z - lastz)*(z - lastz));
            if(distance > distanceDisrupt) cout << "Got jump of distance " << distance << endl;

			  // handle and check sequence number
			int seqNum = atoi(line.substr(22,4).c_str());

			if (seqNumFirst == INT_MAX) {
				seqNumFirst = seqNum;
				chain.rbegin()->idBegin = seqNumFirst;
				chain.rbegin()->atom = atom;
			} else {
                if ((distance > distanceDisrupt) || (seqNum != (seqNumFirst + (int)chain.rbegin()->points.size()))) {
					  // handle error case, if no sub-chain support
					if (chainGaps) {
                        //cout << "Creating new sub chain" << endl;
						// handle subchain support --> end of subchain reached
//std::cerr <<" read subchain : fistToRead=" <<chain.rbegin()->idBegin <<" : last=" <<seqNumFirst+chain.rbegin()->seq.size()-1 <<" length=" <<chain.rbegin()->seq.size()<<"\n";
						  // update begin of chain
						seqNumFirst = seqNum;
						  // add next subchain to fill
						chain.push_back(Subchain());
						  // fill new subchain element
						chain.rbegin()->idBegin = seqNum;
						chain.rbegin()->atom = atom;
					} else {
						std::cerr <<"\n   PDB : line " <<lineCount <<" : next sequence number not consecutive and no chain gaps allowed !\n\n"
									<<line <<"\n\n";
						return -1;
					}
				}
			}


			chain.rbegin()->points.push_back(biu::DblPoint(x,y,z));
			  // amino acid
			chain.rbegin()->seq.push_back(line.substr(17,3));

            lastx = x;
            lasty = y;
            lastz = z;

		} // if ATOM and chainID
	} // while
	
	  // check if a model was present and we have found the one to screen
	if (pdbModelPresent && !pdbModelFound) {
		std::cerr <<"\n   PDB : line " <<lineCount <<" : models are present but we have not found the model with id "
				<<pdbModel <<" !\n\n"
					<<"\n\n";
		return -1;
	}

//std::cerr <<" read subchain : fistToRead=" <<chain.rbegin()->idBegin <<" : last=" <<seqNumFirst+chain.rbegin()->seq.size()-1 <<" length=" <<chain.rbegin()->seq.size()<<"\n";
	  // return the file ID of the last AA read from stream
	if (seqNumFirst == INT_MAX) {
		std::cerr <<"\n   PDB : line " <<lineCount <<" : up to here no '"
				<<atom
				<<"' atom coordinates were found !\n\n"
					<<"\n\n";
		return -1;
	} else {
		return seqNumFirst+(int)chain.rbegin()->points.size()-1;
	}
}

//////////////////////////////////////////////////////////////////////////


//double
//distP2 (	const biu::DblPoint& d1,
//			const biu::DblPoint& d2,
//			const double d2Stretch = 1.0,
//			const bool rotate = true )
//{
//	biu::DblPoint t2 = d2;
//	double x = d1.getX() - t2.getX();
//	double y = d1.getY() - t2.getY();
//	double z = d1.getZ() - t2.getZ();
//	return x*x + y*y + z*z;
//}

double
distP2_d (	const biu::DblPoint& d1,
			const biu::DblPoint& d2,
			const double d2Stretch,
			const bool rotate );
double
distP2 (	const biu::DblPoint& d1,
			const biu::IntPoint& d2,
			const double d2Stretch,
			const bool rotate )
{
	biu::DblPoint t2 = (biu::DblPoint)d2;
	return distP2_d(d1, t2, d2Stretch, rotate);
}
double
distP2_d (	const biu::DblPoint& d1,
			const biu::DblPoint& d2,
			const double d2Stretch ,
			const bool rotate )
{
	biu::DblPoint t2 = (biu::DblPoint)d2;
	t2.setX(t2.getX()*d2Stretch);
	t2.setY(t2.getY()*d2Stretch);
	t2.setZ(t2.getZ()*d2Stretch);

	if (rotate) {
		t2 = curRot.rotate(t2);
	}

	double x = d1.getX() - t2.getX();
	double y = d1.getY() - t2.getY();
	double z = d1.getZ() - t2.getZ();
	return x*x + y*y + z*z;
}

double
cRMSD ( const double distSum, const size_t sumNum ) {
	return sqrt( distSum / (double) sumNum);
}

double
dRMSD ( const biu::DPointVec& p1, const biu::DPointVec & p2) {
	if (p1.size() != p2.size())
		return -1.0;
	if (p1.size() < 2)
		return -2.0;
	double diff2Sum = 0.0;
	for (size_t i=0; (i+1)<p1.size(); i++) {
		for (size_t j=i+1; j<p1.size(); j++) {
			double diff =	(p1.at(i)-p1.at(j)).vectorLength()
							- (p2.at(i)-p2.at(j)).vectorLength();
			diff2Sum += diff*diff;
		}
	}
	return sqrt( diff2Sum / (double(p1.size()*(p1.size()-1)) / 2.0) );
}
double
dRMSD ( const biu::DPointVec& p1, const biu::IPointVec & p2, const double baseFactor) {
	if (p1.size() != p2.size())
		return -1.0;
	if (p1.size() < 2)
		return -2.0;
	double diff2Sum = 0.0;
	for (size_t i=0; (i+1)<p1.size(); i++) {
		for (size_t j=i+1; j<p1.size(); j++) {
			double diff =	(p1.at(i)-p1.at(j)).vectorLength()
							- ((baseFactor*p2.at(i))-(baseFactor*p2.at(j))).vectorLength();
			diff2Sum += diff*diff;
		}
	}
	return sqrt( diff2Sum / (double(p1.size()*(p1.size()-1)) / 2.0) );
}
double
dRMSD ( const std::vector<Subchain> &chain, const std::vector<SubchainLat> & fit, const double baseFactor) {
	double distSum = 0.0;
	if (chain.size() != fit.size())
		return -3.0;
	
	size_t num = 0;
	  // for all chains and all position pairs
	for (size_t ci=0; ci<chain.size(); ci++) {
		for (size_t ii=0; ii<chain.at(ci).points.size(); ii++) {
			
			for (size_t cj=ci; cj<chain.size(); cj++) {
				for (size_t ij=(ii+1); ij<chain.at(cj).points.size(); ij++) {

					  // calculate distance difference
					double diff 
						=	(chain.at(ci).points.at(ii)-chain.at(cj).points.at(ij)).vectorLength()
							- ((baseFactor*fit.at(ci).points.at(ii))-(baseFactor*fit.at(cj).points.at(ij))).vectorLength();
					  // add squared distance difference to sum
					distSum += diff*diff;
					  // count the add
					num++;
				}
			}
			
		}
	}
	
	  // calculate square root of distance difference mean
	return sqrt( distSum / double(num) );
}
double
dRMSD ( const std::vector<Subchain> &chain, const std::vector<biu::DPointVec> & fit ) {
	double distSum = 0.0;
	if (chain.size() != fit.size())
		return -3.0;
	
	size_t num = 0;
	  // for all chains and all position pairs
	for (size_t ci=0; ci<chain.size(); ci++) {
		for (size_t ii=0; ii<chain.at(ci).points.size(); ii++) {
			
			for (size_t cj=ci; cj<chain.size(); cj++) {
				for (size_t ij=(ii+1); ij<chain.at(cj).points.size(); ij++) {

					  // calculate distance difference
					double diff 
						=	(chain.at(ci).points.at(ii)-chain.at(cj).points.at(ij)).vectorLength()
							- ((fit.at(ci).at(ii))-(fit.at(cj).at(ij))).vectorLength();
					  // add squared distance difference to sum
					distSum += diff*diff;
					  // count the add
					num++;
				}
			}
			
		}
	}
	
	  // calculate square root of distance difference mean
	return sqrt( distSum / double(num) );
}

//////////////////////////////////////////////////////////////////////////

bool RMS_equal (std::pair< double, biu::IPointVec > i,std::pair< double, biu::IPointVec > j)
{ 
	return (i.first == j.first); 
}


bool RMS_ordering (std::pair< double, biu::IPointVec > i,std::pair< double, biu::IPointVec > j)
{ 
	return (i.first < j.first); 
}


//////////////////////////////////////////////////////////////////////////


double
fit2lattice(	const biu::DPointVec& p,
				biu::IPointVec& l_res,
				const biu::LatticeModel& lattice,
				const size_t best2store = 1,
				const double baseLength = STANDARD_BASE_LENGTH,
				const bool verbose = false )
{
	  // containers to store extended structures in
	std::vector< std::pair< double, biu::IPointVec > > fits1, fits2;
	  // access to the storage containers
	std::vector< std::pair< double, biu::IPointVec > >& bestFits = fits1, newFits = fits2, tmpFits = fits2;

	  // initialise bestFits with empty structure
	bestFits.push_back( std::pair< double, biu::IPointVec >(0.0,biu::IPointVec()));
	  // add first monomer and update sum of distances
	bestFits[0].second.push_back(biu::IntPoint(0,0,0));
	bestFits[0].first = distP2(p.at(0)-p.at(0),bestFits[0].second.at(0), baseLength, true);

	  // extend all in bestFits with one atom (if possible) until all added
	for (size_t i=1; i<p.size(); i++) {
		  // clear storage for extended structures
		newFits.clear();

		  // extend all of best of last round
		for (std::vector< std::pair< double, biu::IPointVec > >::const_iterator it = bestFits.begin();
				it != bestFits.end(); it++)
		{
			  // copy and extend current best to extend
			biu::IPointVec next = it->second;
			size_t lastID = next.size()-1;
			  // add dummy position to rewrite using the neighborhood of the lattice
			next.push_back(biu::IntPoint(0,0,0));
			  // check for each neighbor in the lattice
			for (biu::LatticeNeighborhood::const_iterator neigh = lattice.getNeighborhood().begin();
					neigh != lattice.getNeighborhood().end(); neigh++)
			{
				  // possible new position of the monomer to extend
				biu::IntPoint newAdd = next[lastID] + *neigh;
				  // check if selfavoiding
				if (std::find(next.begin(), next.end(), newAdd) == next.end()) {
					  // overwrite dummy position ad the end of the structure
					*(next.rbegin()) = newAdd;
					  // copy structure and corresponding distance sum to newFits
					newFits.push_back( std::pair< double, biu::IPointVec >
										( it->first + distP2(p.at(i)-p.at(0),*next.rbegin(), baseLength, true),
											next ));
				}
			}
		}

		  // check if selfavoidin extension was not possible
		if (newFits.size() == 0) {
			if (verbose)
				std::cout <<"\n  --> no selfavoiding extension of " <<bestFits.size() <<" structures of length " <<i <<" possible!" <<std::endl;
			  // go for loop break
			i = p.size();
		}

		  // reduce extended structures to the first best2store best
		if (newFits.size() > best2store) {
			std::sort(newFits.begin(),newFits.end(),RMS_ordering);
			newFits.resize(best2store);
		}
		  // switch bestFits and newFits
		tmpFits = bestFits;
		bestFits = newFits;
		newFits = tmpFits;
	}

	if (bestFits.size() > 0) {
		  // store best fit (== at the beginning of bestFits)
		l_res.resize( bestFits.begin()->second.size());
		std::copy(bestFits.begin()->second.begin(), bestFits.begin()->second.end(), l_res.begin());
	
		  // return minimal cRMS found (of l_res)
		return cRMSD(bestFits.begin()->first,bestFits.begin()->second.size());;
		
	} else {
		return -1.0;
	}
}

//////////////////////////////////////////////////////////////////////////


double
fit2lattice_sc(	const biu::DPointVec& p, // iterating cAlpha and side chain atoms
				biu::IPointVec& l_res,
				const biu::LatticeModel& lattice,
				const size_t best2store,
				const double baseLength,
				const biu::LatticeModel& scLattice,
				const double scContrib,
				const bool verbose )
{
	  // containers to store extended structures in
	std::vector< std::pair< double, biu::IPointVec > > fits1, fits2;
	  // access to the storage containers
	std::vector< std::pair< double, biu::IPointVec > >& bestFits = fits1, newFits = fits2, tmpFits = fits2;

	  // initialise bestFits with empty structure
	bestFits.push_back( std::pair< double, biu::IPointVec >(0.0,biu::IPointVec()));
	  // add first monomer and update sum of distances
	bestFits[0].second.push_back(biu::IntPoint(0,0,0));
	bestFits[0].first = distP2(p.at(0)-p.at(0),bestFits[0].second.at(0), baseLength, true);

	  // extend all in bestFits with one atom (if possible) until all added
	for (size_t i=1; i<p.size(); i++) {
		  // clear storage for extended structures
		newFits.clear();

		  // extend all of best of last round
		for (std::vector< std::pair< double, biu::IPointVec > >::const_iterator it = bestFits.begin();
				it != bestFits.end(); it++)
		{
			  // copy and extend current best to extend
			biu::IPointVec next = it->second;
			  // set lastID to last cAlpha index
			size_t lastID = next.size()-(next.size()%2==0?2:1);
			  // add dummy position to rewrite using the neighborhood of the lattice
			next.push_back(biu::IntPoint(0,0,0));
			  // get current lattice model
			const biu::LatticeModel& curLat = (i%2==0?lattice:scLattice);
			  // check for each neighbor in the current lattice
			for (biu::LatticeNeighborhood::const_iterator neigh = curLat.getNeighborhood().begin();
					neigh != curLat.getNeighborhood().end(); neigh++)
			{
				  // possible new position of the monomer to extend
				biu::IntPoint newAdd = next[lastID] + *neigh;
				  // check if selfavoiding
				if (std::find(next.begin(), next.end(), newAdd) == next.end()) {
					  // overwrite dummy position ad the end of the structure
					*(next.rbegin()) = newAdd;
					  // copy structure and corresponding distance sum to newFits
					newFits.push_back( std::pair< double, biu::IPointVec >
										( it->first + ((i%2==0?1.0:scContrib)*distP2(p.at(i)-p.at(0),*next.rbegin(), baseLength, true)),
											next ));
				}
			}
		}

		  // check if selfavoidin extension was not possible
		if (newFits.size() == 0) {
			if (verbose)
				std::cout <<"\n  --> no selfavoiding extension of " <<bestFits.size() <<" structures of length " <<((i/2)+1) <<" possible!" <<std::endl;
			  // go for loop break
			i = p.size();
		}

		  // reduce to best2store only if cAlpha AND side chain atom are added
		if ( i%2 != 0 ) {
			  // reduce extended structures to the first best2store best
			if (newFits.size() > best2store) {
				std::sort(newFits.begin(),newFits.end(),RMS_ordering);
				newFits.resize(best2store);
			}
		}
		  // switch bestFits and newFits
		tmpFits = bestFits;
		bestFits = newFits;
		newFits = tmpFits;
	}

	if (bestFits.size() > 0)  {

		  // store best fit (== at the beginning of bestFits)
		l_res.resize( bestFits.begin()->second.size());
		std::copy(bestFits.begin()->second.begin(), bestFits.begin()->second.end(), l_res.begin());
	
		  // return minimal cRMSD found (of l_res)
		return cRMSD(bestFits.begin()->first,bestFits.begin()->second.size());;
	} else {
		
		return -1.0;
	}
}

//////////////////////////////////////////////////////////////////////////


std::string getTimeString() {
	  // return string
	std::string strTime = "";
	  // get current time
	time_t tim = time(NULL);
	tm *curTime = localtime(&tim);
	  // generate time string
	strTime += char(48+curTime->tm_mday/10);
	strTime += char(48+curTime->tm_mday%10);
	strTime += '-';
	switch (curTime->tm_mon) {
		case 0  : strTime += "JAN"; break;
		case 1  : strTime += "FEB"; break;
		case 2  : strTime += "MAR"; break;
		case 3  : strTime += "APR"; break;
		case 4  : strTime += "MAY"; break;
		case 5  : strTime += "JUN"; break;
		case 6  : strTime += "JUL"; break;
		case 7  : strTime += "AUG"; break;
		case 8  : strTime += "SEP"; break;
		case 9  : strTime += "OCT"; break;
		case 10 : strTime += "NOV"; break;
		case 11 : strTime += "DEC"; break;
	}
	strTime += '-';
	int year10 = curTime->tm_year%10;
	int year100 = (curTime->tm_year-year10)%100;
	strTime += char(48+year100);
	strTime += char(48+year10);
	  // return time string
	return strTime;
}


void
writeXYZ(	const std::string & sourceInfo,
			const std::string & pdbID,
			const std::vector<biu::DPointVec> & points,
			const double cRMSD,
			const double dRMSD,
			const biu::LatticeModel & lattice,
			std::ostream & out)
{
	out	<<"# lattice protein in absolute positions " <<getTimeString()
		<<"\n# lattice model    : " <<lattice.getDescriptor()->getName()
		<<"\n# PDB entry fitted : " <<pdbID
		<<"\n# coordinate RMSD  : " <<cRMSD
		<<"\n# distance RMSD    : " <<dRMSD
		<<"\n# created by LatFit -- (c) Martin Mann 2008 - " <<getTimeString()
		<<"\n#   http://www.bioinf.uni-freiburg.de/Software/"
		<<"\n# parameters : " <<sourceInfo
		<<"\n#"
		<<std::endl;
	for(size_t i=0; i<points.size(); i++) {
		for (biu::DPointVec::const_iterator p = points[i].begin(); p!=points[i].end(); p++) {
			out <<p->getX() <<" " <<p->getY() <<" " <<p->getZ() <<"\n";
		}
		std::cout <<"\n";
	}
	out	<<"#\n# EOF #" <<std::endl;
	out.flush();
}

void
writeCML(	const std::string & sourceInfo,
			const std::string & pdbID,
			const std::vector<Subchain> & orig,
			const std::vector<biu::DPointVec> & fit,
			const double cRMSD,
			const double dRMSD,
			const biu::LatticeModel & lattice,
			std::ostream & out,
			bool hasSideChain,
			const biu::LatticeModel & sclattice)
{
	assertbiu(orig.size()==fit.size(), "original data and fit have different numbers of subchains!");
	// shift output to display with positive dimensions only
	double shift = 0.0;	
	// calculate the minimal value of all coordinates
	for (size_t c=0; c<fit.size(); c++) {
		for (size_t i=0; i<fit[c].size(); i++) {
			shift = std::min( shift, std::min( fit[c][i].getX(), std::min( fit[c][i].getY(), fit[c][i].getZ() )));
		}
	}
	  // calculate shift of all into first quadrant (ensure all coordinates are positive)
	shift = std::max( 0.0, -shift ) + 10.0;

	out	<<"<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>"
		<<"\n\n<!-- FOR BEST VIEWING e.g. IN JMOL USE AFTER LOADING THE SCRIPT -->"
		<<"\n<!-- "<<JMOL_SCRIPT<<" -->"
		<<"\n\n<!-- created with LatFit - (c) Martin Mann 2008 -->"
		<<"\n<!--                 "<<getTimeString()<<"                  -->"
		<<"\n"
		<<"\n<list dictRef=\"cdk:model\" xmlns=\"http://www.xml-cml.org/schema\">"
		<<"\n  <list>"
		<<"\n    <molecule id=\"lattice-protein\">"
		<<"\n    <!-- parameters : " <<sourceInfo <<" -->"
		<<"\n    <!-- lattice    : " <<lattice.getDescriptor()->getName() <<" -->";
	if (hasSideChain) {
		out	<<"\n    <!-- sc-lattice : " <<sclattice.getDescriptor()->getName() <<" -->";
	}
	out	<<"\n    <!-- PDB fitted : " <<pdbID <<" -->"
		<<"\n    <!-- side chain : " <<(hasSideChain?"":"not ") <<"represented -->"
		<<"\n    <!-- cRMSD      : " <<cRMSD <<" -->"
		<<"\n    <!-- dRMSD      : " <<dRMSD <<" -->"
		<<"\n    <!-- sequence : ";
	for (size_t c=0; c<orig.size(); c++) {
		std::copy(orig[c].seq.begin(), orig[c].seq.end(), std::ostream_iterator<std::string>(out, " "));
		out <<" ";
	}
	out	<<" -->"
		<<"\n      <atomArray>";

	for (size_t c=0; c<orig.size(); c++) {
		assertbiu(orig[c].points.size()==fit[c].size(), "original data and fit have different lengths of subchains!");
		for (size_t i=0; i<fit[c].size(); i++) {
			out <<"\n        <atom id=\""<<orig[c].idBegin+(int)i<<"_"<<orig[c].seq[hasSideChain?i/2:i]<<"\""
				<<" elementType=\"C"<<(hasSideChain&&(i%2==1)?"l":"") <<"\""
				<<" x3=\""<<(shift + fit[c][i].getX())<<"\""
				<<" y3=\""<<(shift + fit[c][i].getY())<<"\""
				<<" z3=\""<<(shift + fit[c][i].getZ())<<"\""
				<<" />";
		}
	}
	out	<<"\n      </atomArray>"
		<<"\n      <bondArray>";
	for (size_t c=0; c<orig.size(); c++) {
		for (size_t i=1; i<fit[c].size(); i++) {
			if (hasSideChain) {
				if ( i%2 == 0) {
					out	<<"\n        <bond id=\"b"<<orig[c].idBegin+(int)i<<"\" atomRefs2=\""<<orig[c].idBegin+(int)(i-2)<<"_"<<orig[c].seq[(i-1)/2]<<" "<<orig[c].idBegin+(int)(i)<<"_"<<orig[c].seq[i/2]<<"\" order=\"S\"/>";
				} else {
					if ( i != fit[c].size() ) {
						out	<<"\n        <bond id=\"b"<<orig[c].idBegin+(int)i<<"\" atomRefs2=\""<<orig[c].idBegin+(int)(i-1)<<"_"<<orig[c].seq[(i-1)/2]<<" "<<orig[c].idBegin+(int)(i)<<"_"<<orig[c].seq[i/2]<<"\" order=\"S\"/>";
					}
				}
			} else {
				out	<<"\n        <bond id=\"b"<<orig[c].idBegin+(int)i<<"\" atomRefs2=\""<<orig[c].idBegin+(int)(i-1)<<"_"<<orig[c].seq[i-1]<<" "<<orig[c].idBegin+(int)(i)<<"_"<<orig[c].seq[i]<<"\" order=\"S\"/>";
			}
		}
	}
	out	<<"\n      </bondArray>"
		<<"\n    </molecule>"
		<<"\n  </list>"
		<<"\n</list>"
		<<"\n\n<!-- FOR BEST VIEWING e.g. IN JMOL USE THE SCRIPT -->"
		<<"\n<!-- "<<JMOL_SCRIPT<<" -->"
		<<std::endl;
	out.flush();
}


void
writePDBseq( const std::vector<std::string> & seq,
			char chainID,
			std::ostream & out)
{
	  // print sequence information
	size_t seqResLine = 1;
	out <<"SEQRES "<<std::setw(3) <<seqResLine++ <<" " <<chainID <<" "<<std::setw(4) <<seq.size()<<" ";

	for (size_t i=1; i<=seq.size(); i++) {
		out <<" " <<seq[i-1];
		if (i%13==0 && (i+1) < seq.size()) {
			out <<"\nSEQRES "<<std::setw(3) <<seqResLine++ <<" " <<chainID <<" "<<std::setw(4) <<seq.size()<<" ";
		}
	}
	out <<"\n";
}

	/**
	 *
	 * @return last HETATM ID printed
	 */
int
writePDBhetatm( const std::vector<std::string> & seq,
			const biu::DPointVec & points,
			const char chainID,
			const bool hasSideChain,
			std::ostream & out,
			const int firstAANum,
			const int firstAtomID)
{
	int atomID = 0;
	for (size_t i=0; i<seq.size(); i++) {

		out <<"HETATM" <<std::setw(5)<<(atomID+firstAtomID) <<"  "<<"CA " <<" "<<seq[i]
			<<" " <<chainID <<std::setw(4) <<(firstAANum+(int)i)
			<<" " <<"   " <<std::fixed
			<<std::setprecision(3) <<std::setw(8) <<points[atomID].getX()
			<<std::setprecision(3) <<std::setw(8) <<points[atomID].getY()
			<<std::setprecision(3) <<std::setw(8) <<points[atomID].getZ()
			<<std::setprecision(2) <<std::setw(6) <<1.0
			<<std::setprecision(2) <<std::setw(6) <<0.0
			<<"           C  "
			<<"\n";
		atomID++;
		if (hasSideChain) {
			out <<"HETATM" <<std::setw(5)<<(atomID+firstAtomID) <<"  "<<"CB " <<" "<<seq[i]
				<<" " <<chainID <<std::setw(4) <<(firstAANum+(int)i)
				<<" " <<"   " <<std::fixed
				<<std::setw(8) <<std::setprecision(3) <<points[atomID].getX()
				<<std::setw(8) <<std::setprecision(3) <<points[atomID].getY()
				<<std::setw(8) <<std::setprecision(3) <<points[atomID].getZ()
				<<std::setw(6) <<std::setprecision(2) <<1.0
				<<std::setw(6) <<std::setprecision(2) <<0.0
				<<"          Cl  "
				<<"\n";
			atomID++;
		}
	}
	return atomID-1+firstAtomID;
}

void
writePDBconect( const StrVec & seq,
			const biu::DPointVec & points,
			const char chainID,
			const bool hasSideChain,
			std::ostream & out,
			const int firstBondID)
{
//                  1         2         3         4         5         6         7         8
//         12345678901234567890123456789012345678901234567890123456789012345678901234567890
	 // CA -> CB

	if (hasSideChain) {
		  // give connections explicitly
		int i = 0;
		out	<<"CONECT" <<std::setw(5) <<firstBondID+(2*i)
			<<std::setw(5) <<firstBondID+((2*i)+1) <<std::setw(5) <<firstBondID+(2*(i+1))
			<<"                                                           \n";
		out	<<"CONECT" <<std::setw(5) <<firstBondID+((2*i)+1)
			<<std::setw(5) <<firstBondID+(2*i)
			<<"                                                                \n";
		for (i++; (i+1)< int(seq.size()); i++) {
			out	<<"CONECT" <<std::setw(5) <<firstBondID+(2*i)
				<<std::setw(5) <<firstBondID+(2*(i-1)) <<std::setw(5) <<firstBondID+((2*i)+1) <<std::setw(5) <<firstBondID+(2*(i+1))
				<<"                                                      \n";
			out	<<"CONECT" <<std::setw(5) <<firstBondID+((2*i)+1)
				<<std::setw(5) <<firstBondID+(2*i)
				<<"                                                                \n";
		}
		out	<<"CONECT" <<std::setw(5) <<firstBondID+(2*i)
			<<std::setw(5) <<firstBondID+(2*(i-1)) <<std::setw(5) <<firstBondID+((2*i)+1)
			<<"                                                           \n";
		out	<<"CONECT" <<std::setw(5) <<firstBondID+((2*i)+1)
			<<std::setw(5) <<firstBondID+(2*i)
			<<"                                                                \n";
	} else {
		  // give connections explicitly
		int i = firstBondID-1;
		out	<<"CONECT" <<std::setw(5) <<i+1
			<<std::setw(5) <<(i+1)+1
			<<"                                                                \n";
		for (i++; (i+2)< int(firstBondID+seq.size()); i++) {
			out	<<"CONECT" <<std::setw(5) <<i+1
				<<std::setw(5) <<i <<std::setw(5) <<(i+1)+1
				<<"                                                           \n";
		}
		out	<<"CONECT" <<std::setw(5) <<i+1
			<<std::setw(5) <<i
			<<"                                                                \n";
	}
}
void
writePDB(	const std::string & sourceInfo,
			const std::string & pdbID,
			const std::vector<Subchain> & orig,
			const std::vector<biu::DPointVec> & fit,
			const double cRMSD,
			const double dRMSD,
			const biu::LatticeModel & lattice,
			std::ostream & out,
			bool hasSideChain,
			const biu::LatticeModel & sclattice,
			const bool writeOrigPoints )
{
	const char chainID = 'L';
	  // get lattice name
	std::string latName = lattice.getDescriptor()->getName();
	latName.resize(58,' ');

//                  1         2         3         4         5         6         7         8
//         12345678901234567890123456789012345678901234567890123456789012345678901234567890
	out	<<"HEADER    LATTICE PROTEIN STRUCTURE               "<<getTimeString()<<"   "<<pdbID<<"              \n"
		<<"TITLE     FIT OF THE PDB STRUCTURE "<<pdbID<<" ONTO A LATTICE                          \n"
		<<"COMPND    MOL_ID: 1;                                                            \n"
		<<"COMPND   2 MOLECULE: LATTICE PROTEIN;                                           \n"
		<<"COMPND   3 CHAIN: "<<chainID<<(writeOrigPoints?", P;":";   ")<<"                                                         \n"
		<<"COMPND   4 ENGINEERED: YES                                                      \n"
		<<"SOURCE    MOL_ID: 1                                                             \n"
		<<"KEYWDS    LATTICE FITTING                                                       \n"
		<<"EXPDTA    THEORETICAL MODEL                                                     \n"
		<<"AUTHOR    LATFIT SOFTWARE                                                       \n"
		<<"REMARK  40                                                                      \n"
		<<"REMARK  40 GENERATED WITH LATFIT (C) MARTIN MANN 2008                           \n";
	out	<<"REMARK  40  LATTICE = "<<latName <<"\n"
		<<"REMARK  40  PDB ID  = "<<pdbID<<"                                                      \n"
		<<"REMARK  40  CRMSD   = "<<std::setw(8) <<std::fixed <<std::setprecision(4)<<cRMSD<<"                                                  \n"
		<<"REMARK  40  DRMSD   = "<<std::setw(8) <<std::fixed <<std::setprecision(4)<<dRMSD<<"                                                  \n";
	if (hasSideChain) {
		std::string latName = sclattice.getDescriptor()->getName();
		latName.resize(47,' ');
	out	<<"REMARK  40  SIDE CHAIN MODEL                                                    \n"
		<<"REMARK  40  SIDE CHAIN LATTICE = " <<latName <<"\n";
	}
	out	<<"REMARK 220                                                                      \n"
		<<"REMARK 220 EXPERIMENTAL DETAILS                                                 \n"
		<<"REMARK 220 EXPERIMENT TYPE : THEORETICAL MODELLING                              \n"
		<<"REMARK 220                                                                      \n"
		<<"REMARK 220 PARAMETER SET APPLIED TO LATFIT :                                    \n"
		<<"REMARK 220                                                                      \n";
	  // write parameter set out of infoString
	size_t cut = 0;
	while(cut < sourceInfo.size() && sourceInfo[cut]==' ')
	{ cut++; }
	while(cut < sourceInfo.size() && cut != std::string::npos) {
		size_t start = cut;
		cut = sourceInfo.find(" ", start+1);
		std::string output = "REMARK 220  ";
		output += sourceInfo.substr(start,cut-start);
		output.resize(80,' ');
		out <<output <<"\n";
		if (cut!=std::string::npos)
			cut++;
	}

	out	<<"REMARK 220                                                                      \n"
		<<"REMARK 220 REMARK: THE LATTICE STRUCTURE WAS FITTED USING A GREEDY METHOD       \n"
		<<"REMARK 225                                                                      \n"
		<<"REMARK 225 THEORETICAL MODEL                                                    \n"
		<<"REMARK 225 THE COORDINATES IN THIS ENTRY REPRESENT A LATTICE MODEL STRUCTURE.   \n"
		<<"REMARK 225                                                                      \n"
		<<"REMARK 225 THE LATTICE PROTEIN STRUCTURE IS GIVEN AS CHAIN "<<chainID<<".                   \n"
		<<"REMARK 225 THE BACKBONE IS REPRESENTED BY A CA-ATOMS AS A 'C'-elements.         \n"
		<<(hasSideChain?"REMARK 225 THE SIDE CHAIN IS REPRESENTED BY A CB-ATOMS AS A 'Cl'-elements.      \n":"")
		<<(writeOrigPoints?"REMARK 225                                                                      \n":"")
		<<(writeOrigPoints?"REMARK 225 THE INITIAL PROTEIN POSITIONS THAT WERE FITTED ARE GIVEN AS CHAIN P. \n":"")
		<<"REMARK 225                                                                      \n"
		;

	  // handle sequence gaps
	if (orig.size() > 1) {
		out	<<"REMARK 465 MISSING RESIDUES                                                     \n"
			<<"REMARK 465 THE FOLLOWING RESIDUES WERE NOT LOCATED IN THE                       \n"
			<<"REMARK 465 EXPERIMENT. (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN               \n"
			<<"REMARK 465 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)                \n"
			<<"REMARK 465                                                                      \n"
			<<"REMARK 465   M RES C SSSEQI                                                     \n"
			;
		for (size_t c=1; c<orig.size(); c++) {
			for (int i=orig[c-1].idBegin+orig[c-1].seq.size(); i<orig[c].idBegin; i++) {
				out <<"REMARK 465     GLY "<<chainID <<std::setw(5) <<i <<std::setw(55) <<" " <<"\n";
			}
		}
		if (writeOrigPoints) {
			for (size_t c=1; c<orig.size(); c++) {
				for (int i=orig[c-1].idBegin+orig[c-1].seq.size(); i<orig[c].idBegin; i++) {
					out <<"REMARK 465     GLY P" <<std::setw(5) <<i <<std::setw(55) <<" " <<"\n";
				}
			}
		}
		out	<<"REMARK 465                                                                      \n";
	}

	  // print sequence information
	{
		StrVec seq = orig[0].seq;
		for (size_t c=1; c<orig.size(); c++) {
			  // insert dummy GLY for missing residues
			seq.insert(seq.end(), orig[c].idBegin-orig[c-1].idBegin-orig[c-1].seq.size(), "GLY");
			seq.insert(seq.end(), orig[c].seq.begin(), orig[c].seq.end());
		}
		writePDBseq( seq, chainID, out);
		if (writeOrigPoints) {
			writePDBseq( seq, 'P', out);
		}
	}

	  // print coordinates
	int hetAtomsWritten = 0;
	int lastHetAtomID = 0;
	for (size_t c=0; c<fit.size(); c++) {
		lastHetAtomID = writePDBhetatm( orig[c].seq, fit[c], chainID, hasSideChain, out, orig[c].idBegin, lastHetAtomID+1);
		hetAtomsWritten += fit[c].size();
	}
	if (writeOrigPoints) {
		for (size_t c=0; c<fit.size(); c++) {
			lastHetAtomID = writePDBhetatm( orig[c].seq, orig[c].points, 'P', hasSideChain, out, orig[c].idBegin, lastHetAtomID+1);
		}
	}
	  // mark end of chain
	out <<"TER                                                                             \n";

	  // print contacts
	int IDshift = 0;
	for (size_t c=0; c<fit.size(); c++) {
		writePDBconect( orig[c].seq, fit[c], chainID, hasSideChain, out, IDshift+1);
		IDshift += orig[c].seq.size() * (hasSideChain ? 2 : 1);
	}
	if (writeOrigPoints) {
		IDshift = 0;
		for (size_t c=0; c<fit.size(); c++) {
			writePDBconect( orig[c].seq, orig[c].points, 'P', hasSideChain, out, hetAtomsWritten+IDshift+1);
			IDshift += orig[c].seq.size() * (hasSideChain ? 2 : 1);
		}
	}

	  // mark EOF
	out <<"END   "
		<<std::endl;

	out.flush();
}




//////////////////////////////////////////////////////////////////////////

void
initAllowedArguments(biu::OptionMap & allowedArgs, std::string &infoText )
{
	allowedArgs.push_back(biu::COption(
							"pdbFile", true, biu::COption::STRING,
							"the pdb file to read or reading STDIN if not given","STDIN"));
	allowedArgs.push_back(biu::COption(
							"pdbAtom", true, biu::COption::STRING,
							"PDB : the atom to fit or 'CoM' for center of mass", "CA"));
	allowedArgs.push_back(biu::COption(
							"pdbAtomAlt", true, biu::COption::CHAR,
							"PDB : the identifier to use in case alternative atoms are present"));
	allowedArgs.push_back(biu::COption(
							"pdbChain", true, biu::COption::CHAR,
							"PDB : chain identifier (use '_' for white space)", "A"));
	allowedArgs.push_back(biu::COption(
							"pdbChainGaps", true, biu::COption::BOOL,
							"PDB : if given, gaps in the chain will be allowed and a subchain set will be fitted instead"));
	allowedArgs.push_back(biu::COption(
							"pdbModel", true, biu::COption::INT,
							"PDB : if the PDB file contains different models, the specified model is used only",
							"1"));
	allowedArgs.push_back(biu::COption(
							"lat", true, biu::COption::STRING,
							"lattice (SQR, CUB, FCC, 210)", "FCC"));
	allowedArgs.push_back(biu::COption(
							"bondLength", true, biu::COption::DOUBLE,
							"length of a bond == the base vectors in the lattice", "3.8"));
	allowedArgs.push_back(biu::COption(
							"nKeep", true, biu::COption::INT,
							"number of best greedy extended structures to keep", "5"));
	allowedArgs.push_back(biu::COption(
							"opt", true, biu::COption::CHAR,
							"the RMSD to optimize : (C)oordinate or (D)istance", "D"));
	allowedArgs.push_back(biu::COption(
							"rotSteps", true, biu::COption::INT,
							"number of steps each of the XYZ rotations (see rotMax) should be split into (>0)", "5"));
	allowedArgs.push_back(biu::COption(
							"rotMax", true, biu::COption::DOUBLE,
							"factor k that limits the rotation for XYZ to [0..k*PI] in radian measure", "0.5"));
	allowedArgs.push_back(biu::COption(
							"refRotSteps", true, biu::COption::INT,
							"refinement phase : number of steps each of the XYZ rotations (see refRotMax) should be split into [0 == no refinement]", "0"));
	allowedArgs.push_back(biu::COption(
							"refRotMax", true, biu::COption::DOUBLE,
							"refinement phase : factor k that limits the rotation for XYZ to +-(k*PI) in radian measure around the best rotation found", "0.1"));
	allowedArgs.push_back(biu::COption(
							"fitSideChain", true, biu::COption::BOOL,
							"fit 'pdbAtom' as additional side chain atom"));
	allowedArgs.push_back(biu::COption(
							"scContrib", true, biu::COption::DOUBLE,
							"contribution of the side chain deviation to the RMSD", "1.0"));
	allowedArgs.push_back(biu::COption(
							"scLat", true, biu::COption::STRING,
							"lattice for side chain placement or '-lat' used instead (SQR, CUB, FCC)"));
	allowedArgs.push_back(biu::COption(
							"fitDirVec", true, biu::COption::BOOL,
							"fits a point |k*(pdbAtom-cAlpha)|==dirVecLength as additional side chain atom"));
	allowedArgs.push_back(biu::COption(
							"dirVecLength", true, biu::COption::DOUBLE,
							"length of the direction vector to fit if 'fitDirVec' is given", "3.8"));
	allowedArgs.push_back(biu::COption(
							"outMode", true, biu::COption::STRING,
							"output mode (CML,XYZ,PDB)", "PDB"));
	allowedArgs.push_back(biu::COption(
							"outFile", true, biu::COption::STRING,
							"name of the output file or written to STDOUT", "STDOUT"));
	allowedArgs.push_back(biu::COption(
							"outAllBest", true, biu::COption::BOOL,
							"writes output to file 'outFile' each time a better fit was found"));
	allowedArgs.push_back(biu::COption(
							"outLatPnt", true, biu::COption::BOOL,
							"output the non rotated lattice fit with integer positions"));
	allowedArgs.push_back(biu::COption(
							"outOrigData", true, biu::COption::BOOL,
							"adds the points to fit (of original pdb structure) to output"));
	allowedArgs.push_back(biu::COption(
							"v", true, biu::COption::BOOL,
							"do verbose output"));
	allowedArgs.push_back(biu::COption(
							"s", true, biu::COption::BOOL,
							"silent : no output except results"));
	allowedArgs.push_back(biu::COption(
							"help", true, biu::COption::BOOL,
							"program parameters and help"));
	allowedArgs.push_back(biu::COption(
							"version", true, biu::COption::BOOL,
							"version information of this program"));

	infoText =	std::string("LatFit does a fitting of a proteins given in PDB format to lattice positions. This can be done either for a backbone model or for a side chain model.\n")
			+	std::string("\n")
			+	std::string("In backbone mode, the given amino acid atoms ('pdbAtom') are represented as close as possible by lattice positions.\n")
			+	std::string("\n")
			+	std::string("If side chains are enabled, the given 'pdbAtom' is fitted as the side chain atom and the cAlpha atoms are used for the fitting of the backbone.\n")
			+	std::string("To fit the side chains a different lattice for atom placement can be given.\n")
			+	std::string("\n")
			+	std::string("Note: Glycine is represented with a side chain as well in this lattice model even it has none in real proteins.\n")
			+	std::string("\n") ;

} // initArguments


//////////////////////////////////////////////////////////////////////////


int printOutput(	const OutMode& outMode,
					const biu::LatticeModel& lattice,
					const biu::LatticeModel& scLattice,
					const std::string& fileName,
					const std::string& paramInfo,
					const std::vector<SubchainLat>& latFit,
					const double baseFactor,
					const double cRMSD,
					const double dRMSD_,
					const std::vector<Subchain>& origChain,
					const std::string& pdbID,
					const bool outLatPnt,
					const bool outOrigData,
					const bool hasSideChain,
					const bool noRotation
					)
{
	assertbiu(latFit.size()==origChain.size(), "fit and original have different number of subchains");
	std::vector<biu::DPointVec> fitPoints;

	if (noRotation) {
		  // just copy scaled data
		for (size_t c=0; c<origChain.size(); c++) {
			  // create new fit point chain vector
			fitPoints.push_back(biu::DPointVec(origChain.at(c).points.size()));
			  // copy scaled data
			for (size_t i=0; i<origChain.at(c).points.size(); i++) {
				fitPoints.at(c)[i] = biu::DblPoint(latFit.at(c).points.at(i)) * baseFactor;
			}
		}
		  // calculate superpositioning
		if (outOrigData) {
			  // container of all points of all subchains
			biu::DPointVec allO, allF;
			  // copy coordinate data to container
			size_t j=0;
			for (size_t c=0; c<origChain.size(); c++) {
				allO.resize( allO.size() + origChain.at(c).points.size() );
				allF.resize( allF.size() + fitPoints.at(c).size() );
				for (size_t i=0; i<fitPoints.at(c).size(); i++) {
					allO[j] = origChain.at(c).points.at(i);
					allF[j] = fitPoints.at(c).at(i);
					j++;
				}
			}
			  // get superpositioning of all points
			biu::SuperPos_Kabsch::bestsuperposition( 
								allF, allO
								, lattice.getDescriptor()->getAutomorphisms());
			
			const biu::DblPoint shiftVec = origChain.at(0).points.at(0) - allO.at(0);
			  // shift superpositioned data and copy back
			j=0;
			for (size_t c=0; c<origChain.size(); c++) {
				for (size_t i=0; i<origChain.at(c).points.size(); i++) {
					fitPoints.at(c)[i]  = allF.at(j) + shiftVec;
					j++;
				}
			}
		}
		
	} else {
		if (hasSideChain) {
			
			  // side chain model
			for (size_t i = 0; i < latFit.size(); i++) {
				assertbiu(latFit.at(i).points.size()==origChain.at(i).points.size(), "fit and original subchain have different number of elements");
				fitPoints.push_back(biu::DPointVec());
				if (outLatPnt) {
					  // do unrotated integer lattice output
					fitPoints[i].push_back(biu::DblPoint(latFit[i].points[0]));
				} else { 
					  // stretch vectors and rotate and finally shift by first original position
					fitPoints[i].push_back(origChain[i].points[0]+(curRot.rotate(baseFactor*latFit[i].points[0])));
				}
				for (size_t j = 1; j < latFit[i].points.size(); j++) {
					if (outLatPnt) {
						  // do unrotated integer lattice output
						fitPoints[i].push_back(biu::DblPoint(latFit[i].points[j]));
					} else { 
						  // stretch vectors and rotate and finally shift by first original position
	
						const size_t last = (j%2==0 ? j-2 : j-1); 
						  // get distance vector between this and the last position
						biu::DblPoint p(latFit[i].points[j]-latFit[i].points[last]);
						  // stretch distance vector
						p *= baseFactor;
						  // rotate new position relative to last position
						p = curRot.rotate(p);
						  // calculate new shifted position
						p = p + fitPoints[i][last];
						  // add to output vector
						fitPoints[i].push_back(p);
					}
				}
			}
		} else  {
			  // backbone model
			for (size_t i = 0; i < latFit.size(); i++) {
				assertbiu(latFit[i].points.size()==origChain[i].points.size(), "fit and original subchain have different number of elements");
				fitPoints.push_back(biu::DPointVec());
				if (outLatPnt) {
					  // do unrotated integer lattice output
					fitPoints[i].push_back(biu::DblPoint(latFit[i].points[0]));
				} else { 
					  // stretch vectors and rotate and finally shift by first original position
					fitPoints[i].push_back(origChain[i].points[0]+(curRot.rotate(baseFactor*latFit[i].points[0])));
				}
				for (size_t j = 1; j < latFit[i].points.size(); j++) {
					if (outLatPnt) {
						  // do unrotated integer lattice output
						fitPoints[i].push_back(biu::DblPoint(latFit[i].points[j]));
					} else { 
						  // stretch vectors and rotate and finally shift by first original position
						
						  // get distance vector between this and the last position
						biu::DblPoint p(latFit[i].points[j]-latFit[i].points[j-1]);
						  // stretch distance vector
						p *= baseFactor;
						  // rotate new position relative to last position
						p = curRot.rotate(p);
						  // calculate new shifted position
						p = p + fitPoints[i][j-1];
						  // add to output vector
						fitPoints[i].push_back(p);
					}
				}
			}
		}
	}

	std::ostream* out = &std::cout;

	std::ofstream* fout = NULL;
	if (fileName.compare("STDOUT") != 0) {
		fout = new std::ofstream(fileName.c_str());
		if (!fout->is_open()) {
			std::cerr <<"\n  ERROR : can not open output file '" <<fileName <<"' !\n";
			return -2;
		}
		out = fout;
	}


	switch (outMode) {
		case CML :
			writeCML(paramInfo,pdbID,origChain,fitPoints,cRMSD,dRMSD_,lattice,*out,hasSideChain, scLattice);
			break;
		case XYZ :
			writeXYZ(paramInfo,pdbID,fitPoints,cRMSD,dRMSD_,lattice,*out);
			break;
		case PDB :
			writePDB(paramInfo,pdbID,origChain,fitPoints,cRMSD,dRMSD_,lattice,*out,hasSideChain,scLattice,outOrigData);
			break;
		default :
			std::cerr <<"\n  ERROR : outMode '" <<outMode <<"' not known !\n";
			return -2;
	}
	if (fileName.compare("STDOUT") != 0) {
		fout->close();
		delete fout; fout = NULL;
		out = &std::cout;
	}

	return 0;
}



//////////////////////////////////////////////////////////////////////////
// Usage:
// -pdbFile="C:/Qt/Zapotec/PDB/3ECA/newFromPDBtools.pdb" -pdbAtom=CA -pdbChain=A -pdbChainGaps -lat=CUB -outMode=PDB -outFile="c:/Qt/out.pdb" -opt=D
// returns now best cRMSD and best dRMSD
std::pair<double, double> mainLatFit( int argc, char** argv ) {

    std::cout << "Latfit received " << argc << " arguments " << std::endl;
    for(int i = 0; i < argc; ++i){
        std::cout << argv[i] << "\n";
    }
    // since now I return a pair, each time there is an error
    std::pair<double, double> returnError = std::pair<double, double>(NAN,NAN);
    std::pair<double, double> returnNotFound = std::pair<double, double>(-1,-1);
	using namespace std;

	typedef std::vector< Subchain > ChainVec;

	//////////////////////////////////////////////////////////////
	// data to fill
	//////////////////////////////////////////////////////////////


	char chainID = ' ';
	ChainVec chain;
	std::string pdbID;
	string pdbAtom = "CA";
	biu::LatticeDescriptor* latDescr = NULL;
	biu::LatticeDescriptor* latDescrSC = NULL;
	std::vector<SubchainLat> latFit;
	size_t nKeep = 1;
	int pdbModel = 1;
	double baseFactor = 1.0;
	double bondLength = STANDARD_BASE_LENGTH;
	double scBondLength = -1.0; // has to be calculated
	bool verbose = false;
	OutMode outMode = CML;
	std::string outFile = "tmp.out";
	bool fitSideChain = false;
	double scContrib = 1.0;
	bool fitDirVec = false;
	double dirVecLength = STANDARD_BASE_LENGTH;
	size_t rotSteps = 1;
	double rotStepWidth = PI;
	std::vector< double > rotMax( 3, 0.5*PI ); // maximal rotation angle
	bool doRefinement = false; // whether or not to do refinement step
	double refRotInterval = 0.1*PI; // +- interval around best rotation angles to do refinement within
	size_t refRotSteps = 0; // number of refinement steps the interval is cut into
	const double precisionAdd = 0.001; // used to allow for < comparison instead of == 0.0
	bool outOrigData = false;
	bool outLatPnt = false;
	bool outAllBest = false;
	double curDRMSD = -1.0;
	bool silent = false;
	char pdbAtomAlt = ' ';
	bool chainGaps = false;	// whether or not subchaining due to PDB numbering gaps are allowed or not
	enum OPT_MODE { OPT_CRMSD, OPT_DRMSD };
	OPT_MODE optMode = OPT_CRMSD;

	//////////////////////////////////////////////////////////////
	// parameter parsing and checking
	//////////////////////////////////////////////////////////////

	biu::OptionMap allowedArgs;
	std::string infoText;
	initAllowedArguments(allowedArgs,infoText);	// init

		// parse programm arguments
	biu::COptionParser opts = biu::COptionParser(	allowedArgs, argc,
													argv, infoText);
		// arguments parseable and all mandatory arguments given
	if (opts.noErrors()) {
			// help output
		if (opts.getBoolVal("help")) {
			opts.coutUsage();
            return returnError;
		}
		if (opts.getBoolVal("version")) {
			giveVersion();
            return returnError;
		}
	} else {
        return returnError;
	}

	fitSideChain = opts.argExist("fitSideChain");
	fitDirVec = opts.argExist("fitDirVec");
	outOrigData = opts.argExist("outOrigData");
	outLatPnt = opts.argExist("outLatPnt");
	outAllBest = opts.argExist("outAllBest");
	silent = opts.argExist("s");
	verbose = opts.argExist("v");
	
	pdbModel = opts.getIntVal("pdbModel");

	if (silent && verbose) {
		cerr <<"\n   Error: can not do silent (-s) and verbose (-v) output at once !\n";
        return returnError;
	}

	if (outOrigData && outLatPnt) {
		cerr <<"\n   Error: can not do lattice point output (-outLatPnt) combined with original data (-outOrigData) at once !\n";
        return returnError;
	}
	
	switch( opts.getCharVal("opt")) {
	case 'c':
	case 'C': optMode = OPT_CRMSD; break;
	case 'd':
	case 'D': optMode = OPT_DRMSD; break;
	default : 
		cerr <<"\n   Error: optimization mode '" <<opts.getCharVal("opt") <<"' not supported !\n";
        return returnError;
	}

	string pdbFile = opts.getStrVal("pdbFile");
	if ( pdbFile.size() == 0 ) {
		cerr <<"\n   Error: no pdb file given !\n";
        return returnError;
	}

	chainID = opts.getCharVal("pdbChain");
	  // check if user is asking for empty chain ID
	if (chainID == '_') {
		chainID = ' ';
	}

	chainGaps = opts.argExist("pdbChainGaps");

	pdbAtom = opts.getStrVal("pdbAtom");
	if (pdbAtom.size()<1 || pdbAtom.size() > 3) {
		cerr <<"\n   Error: pdbAtom '"<<pdbAtom<<"' has to be 1-3 characters long !\n";
        return returnError;
	}
	if (fitSideChain && ( pdbAtom.compare(" ") == 0 || pdbAtom.compare("CA") == 0 ) ) {
		cerr <<"\n   Error: pdbAtom '"<<pdbAtom<<"' is not allowed for side chain atom !\n";
        return returnError;
	}

	if (opts.argExist("pdbAtomAlt")) {
		if (opts.getStrVal("pdbAtomAlt").size() != 1) {
			cerr <<"\n   Error: pdbAtomAlt has to be ONE character !\n";
            return returnError;
		}
		pdbAtomAlt = opts.getCharVal("pdbAtomAlt");
	}

	double dtmp = opts.getDoubleVal("bondLength");
	if ((int)(dtmp*100000) <= 0) {
		cerr <<"\n   Error: bond length = '"<<dtmp<<"' has to be > 0 !\n";
        return returnError;
	}
	bondLength = dtmp;

	dtmp = opts.getDoubleVal("scContrib");
	if ((int)(dtmp*100000) <= 0) {
		cerr <<"\n   Error: side chain contribution = '"<<dtmp<<"' has to be > 0 !\n";
        return returnError;
	}
	scContrib = dtmp;

	dtmp = opts.getDoubleVal("dirVecLength");
	if ((int)(dtmp*100000) <= 0) {
		cerr <<"\n   Error: side chain direction vector length = '"<<dtmp<<"' has to be > 0 !\n";
        return returnError;
	}
	dirVecLength = dtmp;

	dtmp = opts.getDoubleVal("rotMax");
	if ((int)(dtmp*100000) <= 0) {
		cerr <<"\n   Error: maximal rotation factor k for k*PI = '"<<dtmp<<"' has to be > 0 !\n";
        return returnError;
	}
	rotMax[X] = dtmp*PI;
	rotMax[Y] = dtmp*PI;
	rotMax[Z] = dtmp*PI;

	dtmp = opts.getDoubleVal("refRotMax");
	if ((int)(dtmp*100000) <= 0) {
		cerr <<"\n   Error: refRotMax = '"<<dtmp<<"' has to be > 0 !\n";
        return returnError;
	}
	refRotInterval = dtmp*PI;

	int itmp = opts.getIntVal("nKeep");
	if (itmp <= 0) {
		cerr <<"\n   Error: nKeep = '"<<itmp<<"' has to be > 0 !\n";
        return returnError;
	}
	nKeep = (size_t) itmp;

	itmp = opts.getIntVal("rotSteps");
	if (itmp <= 0) {
		cerr <<"\n   Error: rotSteps = '"<<itmp<<"' has to be > 0 !\n";
        return returnError;
	}
	rotSteps = (size_t)itmp;
	rotStepWidth = rotMax[X] / double(itmp);

	itmp = opts.getIntVal("refRotSteps");
	if (itmp < 0) {
		cerr <<"\n   Error: refRotSteps = '"<<itmp<<"' has to be >= 0 !\n";
        return returnError;
	}
	refRotSteps = (size_t)itmp;
	doRefinement = refRotSteps != 0;

	if (opts.getStrVal("outMode").compare("CML") == 0)
		outMode = CML;
	else if (opts.getStrVal("outMode").compare("XYZ") == 0)
		outMode = XYZ;
	else if (opts.getStrVal("outMode").compare("PDB") == 0)
		outMode = PDB;
	else {
		cerr <<"\n   ERROR : Unknown output mode '"+opts.getStrVal("outMode")+"'\n\n";
        return returnError;
	}
	if (fitSideChain && outMode==XYZ) {
		cerr <<"\n   ERROR : output mode '" <<outMode <<"' not possible for side chain fitting !\n\n";
        return returnError;
	}
	if (outOrigData && outMode!=PDB) {
		cerr <<"\n   ERROR : output mode '" <<outMode <<"' currently not supported for output of fitted and original chain !\n\n";
        return returnError;
	}

	outFile = opts.getStrVal("outFile");
	if (outFile.size() == 0) {
		cerr <<"\n   Error: no output file given !\n";
        return returnError;
	}
	if (outAllBest && outFile.compare("STDOUT") == 0) {
		cerr <<"\n   Error: output file has to be given in 'outAllBest' mode !\n";
        return returnError;
	}

	std::string latStr = opts.getStrVal("lat");
	if (latStr.compare("SQR") == 0)
		latDescr = new biu::LatticeDescriptorSQR();
	else if (latStr.compare("CUB") == 0)
		latDescr = new biu::LatticeDescriptorCUB();
	else if (latStr.compare("FCC") == 0)
		latDescr = new biu::LatticeDescriptorFCC();
	else if (latStr.compare("210") == 0)
		latDescr = new biu::LatticeDescriptorCKW();
	else {
		cerr <<"\n   ERROR : Unknown lattice type '"+latStr+"'\n\n";
        return returnError;
	}

	 // create lattices for cAlpha and side chain neighboring
	biu::LatticeModel lattice(latDescr);
	if (opts.argExist("scLat")) { // own lattice for side chain neighboring set
		std::string latStr = opts.getStrVal("scLat");
		if (latStr.compare("SQR") == 0)
			latDescrSC = new biu::LatticeDescriptorSQR();
		else if (latStr.compare("CUB") == 0)
			latDescrSC = new biu::LatticeDescriptorCUB();
		else if (latStr.compare("FCC") == 0)
			latDescrSC = new biu::LatticeDescriptorFCC();
		else {
			cerr <<"\n   ERROR : Unknown lattice type '"+latStr+"' for side chain (scLat)\n\n";
            return returnError;
		}
	} else { // use cAlpha lattice for neighboring
		latDescrSC = latDescr;
	}
	biu::LatticeModel scLattice(latDescrSC);


	{ // calculate base multiplier to fit a neighbor vector to the  given bond length
	biu::LatticeNeighborhood::const_iterator neigh = lattice.getNeighborhood().begin();
	baseFactor = sqrt( (bondLength*bondLength) / (neigh->getX()*neigh->getX() + neigh->getY()*neigh->getY() + neigh->getZ()*neigh->getZ()) );
	}

	{ // calculate the length of a side chain neighbor vector which is depending on the backbone lattice
	biu::LatticeNeighborhood::const_iterator neigh = scLattice.getNeighborhood().begin();
	scBondLength = (baseFactor*(*neigh)).vectorLength();
//	scBaseFactor = sqrt( (scBondLength*scBondLength) / (neigh->getX()*neigh->getX() + neigh->getY()*neigh->getY() + neigh->getZ()*neigh->getZ()) );
	}


	//////////////////////////////////////////////////////////////
	// input parameter output
	//////////////////////////////////////////////////////////////

	std::ostringstream infostream(std::ostringstream::out);

	if (!silent) std::cout <<"\n  pdb-file       = " <<pdbFile;
	infostream <<" pdb-file=" <<pdbFile;
	if (!silent) std::cout <<"\n  pdb-chain      = " <<chainID;
	infostream <<" pdb-chain=" <<chainID;
	if (!silent) std::cout <<"\n  pdb-chain-gaps = " <<(chainGaps?"allowed":"forbidden");
	infostream <<" pdb-chain-gaps=" <<(chainGaps?"true":"false");
	if (!silent) std::cout <<"\n  pdb-atom       = " <<pdbAtom;
	infostream <<" pdb-atom=" <<pdbAtom;
	if (opts.argExist("pdbAtomAlt")) {
		if (!silent) std::cout <<"\n  pdb-atom-alt   = " <<pdbAtomAlt;
		infostream <<" pdb-atom-alt=" <<pdbAtomAlt;
	}
	if (!silent) std::cout <<"\n  pdb-model      = " <<pdbModel;
	infostream <<" pdb-model=" <<pdbModel;
	if (!silent) std::cout <<"\n  lattice        = " <<latDescr->getName();
	infostream <<" lattice=" <<latDescr->getName();
	if (!silent) std::cout <<"\n  bond-len       = " <<bondLength;
	if (!silent) std::cout <<"\n  base-factor    = " <<baseFactor;
	infostream <<" bond-len=" <<bondLength;
	if (!silent) std::cout <<"\n  nKeep          = " <<nKeep;
	infostream <<" nKeep=" <<nKeep;
	if (optMode == OPT_CRMSD) {
		if (!silent) std::cout <<"\n  optimizing     = cRMSD";
		infostream <<" optimizing=cRMSD";
		if (!silent) std::cout <<"\n  rot-steps      = " <<opts.getIntVal("rotSteps");
		infostream <<" rot-steps=" <<opts.getIntVal("rotSteps");
		if (!silent) std::cout <<"\n  rot-max        = " <<opts.getDoubleVal("rotMax") <<"*PI = "<<rotMax[X];
		infostream <<" rot-max=" <<opts.getDoubleVal("rotMax")<<"*PI";
		if (!silent) std::cout <<"\n  do refinement  = " <<(doRefinement?"true":"false");
		infostream <<" do-refinement=" <<(doRefinement?"true":"false");
		if (!silent) std::cout <<"\n  ref-rot-steps  = " <<opts.getIntVal("refRotSteps");
		infostream <<" ref-rot-steps=" <<opts.getIntVal("refRotSteps");
		if (!silent) std::cout <<"\n  ref-rot-max    = " <<opts.getDoubleVal("refRotMax") <<"*PI = "<<refRotInterval*PI;
		infostream <<" ref-rot-max=" <<opts.getDoubleVal("refRotMax")<<"*PI";
	}
	if (optMode == OPT_DRMSD) {
		if (!silent) std::cout <<"\n  optimizing     = dRMSD";
		infostream <<" optimizing=dRMSD";
	}
	if (!silent) std::cout <<"\n  sidechain      = " <<(fitSideChain?"true":"false");
	infostream <<" side-chain=" <<(fitSideChain?"true":"false");
	if (!silent) std::cout <<"\n  sc-lattice     = " <<latDescrSC->getName();
	infostream <<" sc-lattice=" <<latDescrSC->getName();
	if (!silent) std::cout <<"\n  sc-bond-len    = " <<scBondLength;
	if (!silent) std::cout <<"\n  sc-contrib     = " <<scContrib;
	infostream <<" sc-contrib=" <<scContrib;
	if (!silent) std::cout <<"\n  fit dir-vec    = " <<(fitDirVec?"true":"false");
	infostream <<" fit-dir-vec=" <<(fitDirVec?"true":"false");
	if (!silent) std::cout <<"\n  dirVecLeng     = " <<dirVecLength;
	infostream <<" dirVecLength=" <<dirVecLength;
	if (!silent) std::cout <<"\n  outOrigData    = " <<(outOrigData?"true":"false");
	infostream <<" outOrigData=" <<(outOrigData?"true":"false");
	if (!silent) std::cout <<"\n  outMode        = ";
	infostream <<" outMode=";
	switch (outMode) {
		case CML :
			if (!silent) std::cout <<"CML";
			infostream <<"CML";
			break;
		case XYZ :
			if (!silent) std::cout <<"XYZ";
			infostream <<"XYZ";
			break;
		case PDB :
			if (!silent) std::cout <<"PDB";
			infostream <<"PDB";
			break;
		default :
			std::cerr <<"\n  ERROR : outMode '" <<outMode <<"' not known !\n";
            return returnError; // was -2
	}
	if (!silent) std::cout <<"\n\n";

	//////////////////////////////////////////////////////////////
	// input parsing
	//////////////////////////////////////////////////////////////

	ifstream *inFile = NULL;
	istream* in = &std::cin;
	  // open stream if file given
	if (pdbFile.compare("STDIN") != 0) {
		inFile = new ifstream( pdbFile.c_str() );
		if (!inFile->is_open()) {
			cerr <<"\n   ERROR : can not open pdb file '"+pdbFile+"' !\n\n";
            return returnError;
		}
		in = inFile;
	} else if (fitSideChain) {
		cerr <<"\n   ERROR : side chain fit not possible for PDB stream from STDIN due to necessary double read of the stream !\n\n";
        return returnError;
	}

	  // add next subchain to fill
	chain.push_back(Subchain());
	chain[0].atom = pdbAtom;

	  // read the position of the pdbAtoms
    int retVal = parsePDB( *in, chainID, pdbAtom, pdbAtomAlt, chain, pdbID, chainGaps, pdbModel, 5.0);

	  // print the read chain AA number flanks
	if (verbose) {
		for (size_t i=0; i<chain.size(); i++)
			std::cout <<"  read " <<(chain.size()>1?"sub":"")<<"chain from " <<chain[i].idBegin <<" to " <<chain[i].idBegin+(int)chain[i].seq.size()-1 <<"\n";
		std::cout <<std::endl;
	}
	  // number of amino acids in the whole chain (sum over all subchain lengths)
	size_t aaNumber = 0;
	for (size_t i=0; i<chain.size(); i++)
		aaNumber += chain[i].seq.size();

	  // close stream if necessary
	if (pdbFile.compare("STDIN") != 0) {
		in = &std::cin;
		inFile->close();
		delete inFile;
		inFile = NULL;
	}
	
//	std::cout <<"\n read information : \n";
//	for (size_t i=0; i<chain.size();i++) {
//		std::cout <<" chain " <<i <<" : idBegin \t" <<chain[i].idBegin
//				<<"\t atom \t" <<chain[i].atom
//				<<std::endl <<" -> points : "<< chain[i].points.size()<<" : ";
//		for (size_t j=0; j<chain[i].points.size(); j++) {
//			std::cout <<" "<< chain[i].points[j]; std::cout.flush();
//		}
//		std::cout <<"\n -> seq   : "<< chain[i].seq.size()<<" : "; std::cout.flush();
//		for (size_t j=0; j<chain[i].seq.size(); j++) {
//			std::cout <<" "<< chain[i].seq[j]; std::cout.flush();
//		}
//		std::cout <<std::endl;
//	}

	  // get cAlpha atoms for side chain mode
	if (fitSideChain && retVal > 0) {
		ifstream inF( pdbFile.c_str() );
		if (!inF.is_open()) {
			cerr <<"\n   ERROR : can not open pdb file '"+pdbFile+"' for second read !\n\n";
            return returnError;
		}
		ChainVec cAlpha;
        retVal = parsePDB( inF, chainID, "CA", pdbAtomAlt, cAlpha, pdbID, chainGaps, pdbModel, 5.0);
		inF.close();
		assertbiu(cAlpha.size()>0, "no C_alpha atoms found");
		
//		std::cout <<"\n\n read information : \n";
//		for (size_t i=0; i<cAlpha.size();i++) {
//			std::cout <<" cAlpha " <<i <<" : idBegin \t" <<cAlpha[i].idBegin
//					<<"\t atom \t" <<cAlpha[i].atom
//					<<std::endl <<" -> points : "<< cAlpha[i].points.size()<<" : ";
//			for (size_t j=0; j<cAlpha[i].points.size(); j++) {
//				std::cout <<" "<< cAlpha[i].points[j]; std::cout.flush();
//			}
//			std::cout <<"\n -> seq   : "<< cAlpha[i].seq.size()<<" : "; std::cout.flush();
//			for (size_t j=0; j<cAlpha[i].seq.size(); j++) {
//				std::cout <<" "<< cAlpha[i].seq[j]; std::cout.flush();
//			}
//			std::cout <<std::endl;
//		}
		
		bool equalNumber = retVal > 0 && chain.size() == cAlpha.size();
		for (size_t i=0; equalNumber && i<cAlpha.size(); i++) {
			equalNumber =	chain[i].seq.size() == cAlpha[i].seq.size()
							&& chain[i].points.size() == cAlpha[i].points.size();
		}

		if (equalNumber && chain[0].seq.size() > 0) {

			// check if only direction vector instead of side chain points have to be used
			if (fitDirVec) {
				size_t aaNum = 0;
				double dirVecLengthSum = 0.0;
				for (size_t c=0; c<chain.size(); c++) {
					assertbiu(chain[c].seq.size()==cAlpha[c].seq.size(), "different sequence lengths of C_alpha and side chain atoms");
					for (size_t i = 0; i < chain[c].points.size(); i++) {
						  // get direction vector
						biu::DblPoint dir = chain[c].points[i] - cAlpha[c].points[i];
						dirVecLengthSum += dir.vectorLength();
						  // calculate length correction factor
						  // f = sqrt( l^2 / (x^2 + y^2 + z^2) )
						double lengthFix = 0.0;
						if (dir != biu::DblPoint(0.0,0.0,0.0))
							lengthFix = sqrt( (dirVecLength*dirVecLength) / (dir.getX()*dir.getX() + dir.getY()*dir.getY() + dir.getZ()*dir.getZ()) );
						  // set side chain point to a cAlpha+(f*dir)
						chain[c].points[i] = cAlpha[c].points[i] + (dir*lengthFix);
					}
				}
				if (!silent) {
					std::cout	<<"\n  avg dir vec length  = " <<(dirVecLengthSum/double(aaNum))
								<<"\n  new dir vec length  = " <<dirVecLength <<"\n"<<endl;
				}
			}

			for (size_t c=0; c<chain.size(); c++) {
				  // blow up to double size
				cAlpha[c].points.resize(cAlpha[c].points.size()*2, biu::DblPoint(0.0,0.0,0.0));
				  // insert side chain atoms
				for (size_t i=chain[c].points.size(); i>0; i--) {
					  // shift cAlpha to new position
					cAlpha[c].points[(i-1)*2] = cAlpha[c].points[i-1];
					  // insert side chain atom
					cAlpha[c].points[((i-1)*2)+1] = chain[c].points[i-1];
				}
				  // overwrite old points data
				chain[c].points = cAlpha[c].points;
			}
		} else {
//			for (size_t c=0; c<cAlpha.size(); c++) {
//				std::cerr <<" cAlpha : " <<c <<" : seq.size() = " <<cAlpha[c].seq.size() <<"  points.size() = " <<cAlpha[c].points.size() <<std::endl;
//				std::cerr <<" chain : " <<c <<" : seq.size() = " <<chain[c].seq.size() <<"  points.size() = " <<chain[c].points.size() <<std::endl;
//			}
			if (!equalNumber) {
				cerr <<"\n   ERROR : number of CA positions is different to "<<pdbAtom<<" positions !";
				if (!opts.argExist("pdbAtomAlt"))
					cerr <<"\n       --> Maybe due to alternative atom positions ?"
						 <<"\n       --> Check '-pdbAtomAlt' parameter !"
						 <<"\n       --> Or incomplete amino acid atom information in PDB file ?";
				cerr <<"\n\n";
			}
			if (chain[0].seq.size() == 0)
				cerr <<"\n   ERROR : no "<<pdbAtom<<" positions in the first subchain !\n\n";
            return returnError;
		}
	}

//	for (size_t c=0; c<chain.size(); c++) {
//		std::cerr <<" chain : " <<c <<" : seq.size() = " <<chain[c].seq.size() <<"  points.size() = " <<chain[c].points.size() <<std::endl;
//	}

//	std::ofstream tmpFile("tmp.out");
//	writeCML("",pdbID,seq,points,0.0,0.0,lattice,tmpFile,fitSideChain, scLattice);
//	tmpFile.close();

	if (retVal < 0) {
		std::cerr <<"\n\n   ERROR : could not parse the PDB file !\n" <<endl;
        return returnError; // was retVal to propagate error type;
	}
	if (aaNumber < 2) {
		std::cerr <<"\n\n   ERROR : could only read "<<aaNumber<<" atoms from file '"<<outFile<<"' !\n" <<endl;
        return returnError; // was -2;
	}

	//////////////////////////////////////////////////////////////
	// lattice fitting
	//////////////////////////////////////////////////////////////

	  // print the positions to fit (currently non-side chain only)
	if (verbose && !fitSideChain) {
		for (size_t k=0; k < chain.size(); k++) {
			assertbiu(chain[k].points.size() == chain[k].seq.size(), "points and sequence data differ in length");
			for (size_t i = 0; i < chain[k].points.size(); i++) {
				std::cout <<"  " <<std::setw(6)<<chain[k].idBegin+(int)i <<" = " <<chain[k].seq[i] <<" : " <<chain[k].points[i] <<"\n";
			}
			std::cout <<std::endl;
		}
	}

	
	  // stores best cRMSD found
	double bestCRMSD = (double)INT_MAX;
	
	switch (optMode) {
	
	case OPT_CRMSD : {
		//////////////////////////////////////////////////////////////
		// CRMSD OPTIMIZATION
		//////////////////////////////////////////////////////////////
		  // number of amino acids of the last fit
		size_t fitAAnumber = 0;
		  // stores the rotation information for the best fit found (init using curRot)
		biu::Rotator3D bestRotation = curRot;
	
	//	if (!silent) std::cout <<"\n  cRMSD values :";
	//	std::cout.flush();
	
		  // overall number of rotation steps to do
		size_t numOfRotSteps = (rotSteps+1)*(rotSteps+1)*(rotSteps+1);
		  // width of the progress bar to print
		const size_t PROGRESS_BAR_WIDTH = 64;
		  // size of the already printed progress bar
		size_t progressBarSize = 0;
		  // number of rotation steps done so far
		size_t curRotSteps = 0;
	
		  // flag for refinement iteration
		bool doOneMoreFit = true;
		  // flag that becomes true if a better solution was found
		  // used for the progress bar output to mark when a new solution was found
		bool foundNewBest = false;
	
		if (!silent) {
			  // print progress bar limits
			std::cout 	<<"\n"
						<<"  Fitting PDB structure to different lattice rotations :\n"
						<<"\n"
						<<"    |0%" <<std::setw(PROGRESS_BAR_WIDTH-3) <<100 <<"%|\n"
						<<"    |";
			std::cout.flush();
		}
	
		  // the rotation start angle (changed in case of refinement)
		std::vector< double > startRot( 3, 0.0 );
	
		while (doOneMoreFit)
		{
			  // rotate around X axis --> angle increase
			for (double xRot = startRot[X]; xRot < rotMax[X]+precisionAdd; xRot += rotStepWidth) {
	
			  // rotate around Y axis --> angle increase
			for (double yRot = startRot[Y]; yRot < rotMax[Y]+precisionAdd; yRot += rotStepWidth) {
	
			  // rotate around Z axis --> angle increase
			for (double zRot = startRot[Z]; zRot < rotMax[Z]+precisionAdd; zRot += rotStepWidth) {
				
				  // set angles in rotation handler
                // This is changing all the space computations in the BIU namespace. Hidden jobs there
				curRot.setRotation( xRot, yRot, zRot );
	
				  // find best fit for current lattice rotation
	
				std::vector<SubchainLat> lastFitted; // best fit for this rotation
				  // add first fitted chain to fill
				lastFitted.push_back(SubchainLat());
				fitAAnumber = 0;
				double cRMSD = 0.0;
				if (fitSideChain) {
					double cSquareDist = 0.0;
					for (size_t i=0; i<chain.size(); i++) {
						cRMSD = fit2lattice_sc(	chain[i].points,
												lastFitted[i].points,
												lattice,
												nKeep,
												baseFactor,
												scLattice,
												scContrib,
												verbose);
						  // check if a fit was found
						if (cRMSD < 0.0) {
							cRMSD = (double)INT_MAX;
							 // no extension was possible -> stop for this rotation
							break;
						}
						  // update number of fitted amino acids
						fitAAnumber += lastFitted[i].points.size()/2;
						  // recalculate cSquareDist and add
						cSquareDist += cRMSD*cRMSD*double(lastFitted[i].points.size());
						if (i+1<chain.size()) {
							  // add next fitted chain to fill
							lastFitted.push_back(SubchainLat());
						}
					}
					  // calculate overall cRMSD of all subchains
					cRMSD = sqrt(cSquareDist/double(fitAAnumber*2));
				} else {
					double cSquareDist = 0.0;
					for (size_t i=0; i<chain.size(); i++) {
                        //cout << "CMRSD-Chain " << i << endl;
						cRMSD = fit2lattice(	chain[i].points,
												lastFitted[i].points,
												lattice,
												nKeep,
												baseFactor,
												verbose);
						  // check if a fit was found
						if (cRMSD < 0.0) {
							cRMSD = (double)INT_MAX;
							 // no extension was possible -> stop for this rotation
							break;
						}
						  // recalculate cSquareDist and add
						fitAAnumber += lastFitted[i].points.size();
						cSquareDist += cRMSD*cRMSD*double(chain[i].points.size());
						if (i+1<chain.size()) {
							  // add next fitted chain to fill
							lastFitted.push_back(SubchainLat());
						}
					}
					  // calculate overall cRMSD of all subchains
					cRMSD = sqrt(cSquareDist/double(aaNumber));
				}
	
				  // check if new fit is complete and better than the old best one
				if (	lastFitted.size()==chain.size()
						&& fitAAnumber == aaNumber
						&& cRMSD < bestCRMSD)
				{
					  // remind we have found something better
					foundNewBest = true;
					  // store new best solution
					latFit = lastFitted;
					  // store new best cRMSD
					bestCRMSD = cRMSD;
					  // store new best rotation matrices
					bestRotation = curRot;
	//				  // print new best cRMSD so far
	//				if (!silent && verbose) {
	//					cout <<"\n";
	//				}
	//				if (!silent) cout <<" > "<<std::setw(7) <<std::fixed <<std::setprecision(4) <<cRMSD;
	//				cout.flush();
					if (outAllBest) {
						curDRMSD = dRMSD( chain, latFit, baseFactor );
						printOutput(	outMode,
										lattice,
										scLattice,
										outFile,
										infostream.str(),
										latFit,
										baseFactor,
										bestCRMSD,
										curDRMSD,
										chain,
										pdbID,
										outLatPnt,
										outOrigData,
										fitSideChain,
										(optMode != OPT_CRMSD)
									);
					}
				}
				  // increase number of done rotation steps;
				curRotSteps++;
	
				if (!silent) {
					  // do progress bar output
					  // calculate necessary progress bar extension
					size_t progExt = ((curRotSteps * PROGRESS_BAR_WIDTH) / numOfRotSteps) - progressBarSize;
					if (progExt > 0) {
						progressBarSize += progExt;
						if (progExt--) {
							std::cout	<<(foundNewBest?'+':'#');
						}
						  // mark that we know that a better one was found
						foundNewBest = false;
						  // print additional bar elements if necessary (for low rot steps)
						while ( progExt-- ) {
							std::cout	<<'#';
						}
						std::cout.flush();
					}
				}
	
			} // z rotation
			} // y rotation
			} // x rotation
	
			  // handle refinement
			if (doRefinement) {
	
				if (!silent) {
					std::cout 	<<"|\n\n"
								<<"  Refined structure fitting :\n"
								<<"\n"
								<<"    |0%" <<std::setw(PROGRESS_BAR_WIDTH-3) <<100 <<"%|\n"
								<<"    |";
					std::cout.flush();
				}
				  // reinit progress variables
				progressBarSize = 0;
				curRotSteps = 0;
				foundNewBest = false;
				numOfRotSteps = (2*refRotSteps+1)*(2*refRotSteps+1)*(2*refRotSteps+1);
	
	
	
				  // store new best rotation angle
				rotMax[X] = bestRotation.getRotationX() + refRotInterval;
				rotMax[Y] = bestRotation.getRotationY() + refRotInterval;
				rotMax[Z] = bestRotation.getRotationZ() + refRotInterval;
	
				  // shift to begin of refinement interval around best rotation angle
				startRot[X] = bestRotation.getRotationX() - refRotInterval;
				startRot[Y] = bestRotation.getRotationY() - refRotInterval;
				startRot[Z] = bestRotation.getRotationZ() - refRotInterval;
	
				rotStepWidth = (rotMax[X] - bestRotation.getRotationX()) / double(refRotSteps);
	
				  // avoid another iteration of refinement after the coming one
				doRefinement = false;
	
			} else { // no further refinement needed
				  // avoid next iteration
				doOneMoreFit = false;
			}
	
		} // end while
	
		std::cout <<"|" <<std::endl;
		

		if (!silent) std::cout <<"\n\n  best cRMSD  = " << bestCRMSD <<" Angstroms\n";

		  // set rotation matrix to the one corresponding to the best fit
		curRot = bestRotation;
		  // calculate fitted points for 
		std::vector<biu::DPointVec> fitPoints;
		double cAlphaDist = 0.0, scDist = 0.0;
		size_t atomNumber = 0;
		  // set rotation matrices to best found
		curRot = bestRotation;
		for (size_t c = 0; c < latFit.size(); c++) {
			fitPoints.push_back(biu::DPointVec());
			for (size_t j=0; j<latFit[c].points.size();j++) {
				atomNumber++;
				fitPoints[c].push_back(chain[c].points[0]+(curRot.rotate(baseFactor*latFit[c].points[j])));
				  // calc independent cRMSDs for cAlphas and side chains
				if (fitSideChain) {
					if (j%2==0) {
						double dist = (chain[c].points[j]-fitPoints[c][j]).vectorLength();
						cAlphaDist += dist*dist;
					} else {
						double dist = (chain[c].points[j]-fitPoints[c][j]).vectorLength(); 
						scDist += dist*dist;
					}
				}
			}
		}
		  // give cRMSD for backbone and side chain atoms
		if (!silent && fitSideChain) {
			std::cout <<"  cAlpha      = " << cRMSD(cAlphaDist,atomNumber/2) <<" Angstroms\n";
			std::cout <<"  side chain  = " << cRMSD(scDist,atomNumber/2) <<" Angstroms\n";
		}
		  // calculate global dRMSD for all subchains
		curDRMSD = dRMSD( chain, fitPoints );

		break;
	/////////////////////////////////////////////////////////////////////////
	} // END CRMSD OPTIMIZATION
	/////////////////////////////////////////////////////////////////////////

	case OPT_DRMSD : {
		//////////////////////////////////////////////////////////////
		// DRMSD OPTIMIZATION
		//////////////////////////////////////////////////////////////

		std::vector<SubchainLat> lastFitted;
		  // add first fitted chain to fill
		lastFitted.push_back(SubchainLat());

		latFit.resize(chain.size());
				
		  // derive data in chain format
		AAchain_D aaOrig = deriveAAvec( chain, fitSideChain );
		  // the lattice protein data to fill
		AAchain_I aaLP;
		
		  // perform fitting
		curDRMSD = fit2lattice_dRMSD( aaOrig, aaLP, fitSideChain, lattice, nKeep, baseFactor, !silent);
			
// TODO : +function : dRMSD( AAchain, AAchain ) 
// TODO : +function : cRMSD( AAchain, AAchain )
// TODO : +function : superposition( AAchain, AAchain )
// TODO : +function : pdb output of AAchain
		
// TODO : cRMSD optimization : AAchain based version ...
// TODO : make latFit variable obsolete
		
		
		  // check if a fit was found
		if (curDRMSD < 0.0) {
			curDRMSD = (double)INT_MAX;
			bestCRMSD = (double)INT_MAX;
			latFit.clear();
		} else {
			  // copy data to "old" representation and get CRMSD
			  // (currently needed for structure output)
			size_t i=0;
			latFit.resize(chain.size());
			biu::DPointVec allO(aaOrig.size()*(fitSideChain?2:1));
			biu::DPointVec allL(aaOrig.size()*(fitSideChain?2:1));
			for (size_t c=0; c<chain.size(); c++) {
                //cout << "DMRSD-Chain " << i << endl;
				latFit[c].points.resize(chain.at(c).seq.size()*(fitSideChain?2:1));
				for (size_t j=0; j<chain.at(c).seq.size(); j++) {
					if (fitSideChain) {
						latFit[c].points[j*2] = aaLP[i].posCA;
						latFit[c].points[(j*2)+1] = aaLP[i].posCB;
					} else {
						latFit[c].points[j] = aaLP[i].posCA;
					}
					allO[i] = aaOrig[i].posCA;
					allL[i] = biu::DblPoint(aaLP[i].posCA)*baseFactor;
					if (fitSideChain) {
						allO[aaOrig.size()+i] = aaOrig[i].posCB;
						allL[aaOrig.size()+i] = biu::DblPoint(aaLP[i].posCB)*baseFactor;
					}
					i++;
				}
			}
			  // calculate final coordinate RMSD
			bestCRMSD = biu::SuperPos_Kabsch::bestsuperposition( allL, allO, lattice.getDescriptor()->getAutomorphisms());
			
		}

		
	/////////////////////////////////////////////////////////////////////////
	} // END DRMSD OPTIMIZATION
	/////////////////////////////////////////////////////////////////////////
	} // end switch for optimization mode
	
	
	  // check if no approximation was found
	if (latFit.size() == 0) {
		std::cout	<<"\n\n  No selfavoiding lattice approximation found !\n\n"
					<<"  --> maybe try with increased nKeep or rotSteps parameter.\n"
					<<std::endl;
        return returnNotFound; // was 0
		
	} else if (!silent) {
		std::cout	<<"\n\n  best cRMSD  = " << bestCRMSD <<" Angstroms\n"
					<<"\n  dRMSD       = " << curDRMSD <<" Angstroms\n";
	}

	 // print final structure output
	printOutput(	outMode,
					lattice,
					scLattice,
					outFile,
					infostream.str(),
					latFit,
					baseFactor,
					bestCRMSD,
					curDRMSD,
					chain,
					pdbID,
					outLatPnt,
					outOrigData,
					fitSideChain,
					(optMode != OPT_CRMSD)
				);
	
	if (verbose) {
		std::cout <<"\n absolute move strings per subchain:\n\n";
		for (size_t i=0; i<latFit.size(); i++) {
			std::cout <<std::setw(3) <<(i+1) <<" : "
					<<biu::LatticeProteinUtil::toString( 
							biu::LatticeProteinUtil::toMoveSequence(
									latFit.at(i).points, *latDescr, fitSideChain )
							, *latDescr, fitSideChain )
					<<"\n";
		}
		std::cout <<std::endl;
	}

	if (latDescrSC != NULL && latDescrSC != latDescr) {
		delete latDescrSC; 
	}
	latDescrSC = NULL;
	if (latDescr != NULL) {
		delete latDescr;
		latDescr = NULL;
	}

	  // final output flush
	std::cout <<std::endl;

    return std::pair<double,double>(bestCRMSD, curDRMSD);
}
