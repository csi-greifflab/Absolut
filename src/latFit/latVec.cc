
#include <sstream>
#include <fstream>
#include <map>
#include <iomanip>
#include <algorithm>

#include <biu/util/Util_String.h>

#include <biu/Point.hh>
#include <biu/Matrix.hh>
#include <biu/OptionParser.hh>
#include <biu/LatticeDescriptorSQR.hh>
#include <biu/LatticeDescriptorCUB.hh>
#include <biu/LatticeDescriptorFCC.hh>
#include <biu/LatticeModel.hh>
#include <biu/DistanceEnergyFunction.hh>

#include "energyFileSupport.hh"
#include "PDBsupport.hh"
#include "version.hh"


///////////////////////////  TYPEDEFS  /////////////////////////////////////////

	  // store energy and COMPRESSED structure representation in absolute moves
	typedef std::pair< double, biu::MoveAlphabet::CSequence > ES_pair;
	  // type of structure container with their corresponding energy
	typedef std::set< ES_pair > StructStore;
	
	enum OutMode { OUT_SILENT, OUT_NORMAL, OUT_VERBOSE }; 

////////////////////////  GLOBAL CONSTANTS  ////////////////////////////////////

	const size_t OUT_PRECISION = 3;		//!< precision of floating number output
	const size_t MIN_SEQ_LENGTH = 4;	//!< minimal sequence length to process

////////////////////////////////////////////////////////////////////////////////


/**
 * Decodes the move string representation of a side chain protein into 3D
 * coordinates.
 * 
 * @param absMoves the internal absolute move string representation
 * @param points (IN-OUT) the 3D coordinate container to fill
 * @param lattice the lattice to use
 * @return reference to the changed input coordinate container
 */
biu::IPointVec &
abs2pointsSC(	const biu::MoveSequence & absMoves,
				biu::IPointVec & points,
				const biu::LatticeModel & lattice
				)
{
	  // resize container
	points.resize(absMoves.size()+1, biu::IntPoint(0,0,0));
	
	points[0] = biu::IntPoint(0,0,0);
	
	for (size_t i = 0; i< absMoves.size(); i++) {
		  // add side chain monomer
		if (i%2 == 0) {
			points[i+1] = points[i] + lattice.getNeighborhood().getElement(absMoves[i]);
		} else { // add next backbone
			points[i+1] = points[i-1] + lattice.getNeighborhood().getElement(absMoves[i]);
		}
	}
	
	return points;
}

////////////////////////////////////////////////////////////////////////////////


/**
 * Converts the internal used absolute move string representation into the user
 * readable format with brackets around side chain moves
 * 
 * @param absMoves the absolute move string to convert
 * @param lattice  the lattice used
 * @return the converted string
 */
std::string
abs2absSC( const std::string absMoves, const biu::LatticeModel & lattice ) 
{
	  // no move
	if (absMoves.size() == 0)
		return "";

	  // get length of a single move string in this lattice
	const size_t elemLength = lattice.getDescriptor()->getAlphabet()->getElementLength();

	  // only one side chain move 
	if (absMoves.size() == elemLength)
		return "("+absMoves+")";
	
	  // add first side chain move
	size_t pos = 0;
	std::string ret = "(" + absMoves.substr(pos, elemLength) + ")";
	pos += elemLength;
	
	  // add remaining backbone and side chain moves
	while (pos < absMoves.size()) {
		ret += absMoves.substr(pos, elemLength) + "(";
		ret += absMoves.substr(pos+elemLength, elemLength) + ")";
		pos += 2*elemLength;
	}
	  // final string	
	return ret;
}

////////////////////////////////////////////////////////////////////////////////


/**
 * Checks whether or not position with index idx of the given structure has a 
 * free neighbored position in the lattice.
 * 
 * @param str the structure in 3D coordinates to check
 * @param idx the index to check
 * @param lattice the lattice in use
 * @return true if the last position has at least one free neighbor, false otherwise
 */
bool
hasFreeNeighbor(	const biu::IPointVec & str,
					const size_t idx,
					const biu::LatticeModel & lattice)
{
	if (str.size() < 3)
		return true;
	
	assertbiu(idx < str.size(), "given index exceeds structure length");
	
	// check if the end position has a free neighbored position in the lattice
	
	biu::IntPoint actCheck = biu::IntPoint(0,0,0);

	// try all possible extensions
	for (	biu::LatticeNeighborhood::const_iterator actNeigh = lattice.getNeighborhood().begin(); 
			actNeigh != lattice.getNeighborhood().end(); actNeigh++ )
	{
		actCheck = str[idx] + *actNeigh;
		  // check if this neighbor is already part of the chain
		bool isNotFree = false;
		for (biu::IPointVec::const_iterator it = str.begin(); !isNotFree && it != str.end(); it++) {
			isNotFree = (*it == actCheck);
		}
		  // this neighbor is not part of the chain --> at least one free neighbor of idx !
		if (!isNotFree)
			return true;
	}
	
	return false;
}


////////////////////////////////////////////////////////////////////////////////


/**
 * Implements a chain-growth algorithm that extends all structures within a 
 * given energy interval deltaE from the minimal energy reachable so far. All
 * other structures are rejected. This can lead to a non-extensible chain. 
 * Therefore, a extensibility check can be done that validates if the last 
 * monomer position has at least one free neighbored position. This is 
 * sufficient in most cases but might still be not enough for long sequences.
 * 
 * @param seq the sequence of the lattice protein
 * @param structStore (IN-OUT) container that store structures from the last 
 *                    extension and the current one
 * @param curIdx (IN-OUT) the monomer index to append in the next iteration
 * @param lattice the lattice model to use
 * @param energyC the contact energy function to use if not NULL 
 * @param energyD the distance based energy function used if 
 *                energyC is equal to NULL
 * @param deltaE the energy interval used to reject the extension of structures
 *               that have an energy higher than the minimal energy so far plus
 *               deltaE
 * @param minExtCheckLength the minimal chain length that makes an extensibility
 *                          check necesary. If == 0 than NO check is done.
 * @param outMode the output level
 * @param maxStoreSize maximal number of structures to extend
 * @return TRUE if the chain was extended to full length, FALSE if the algorithm
 *         got stucked.  
 * 
 */
bool
chainGrowth(	const biu::Alphabet::Sequence & seq,
				std::vector<StructStore*> & structStore,
				size_t & curIdx,
				const biu::LatticeModel & lattice,
				const biu::DistanceEnergyFunction * const energy,
				const double deltaE,
				const size_t minExtCheckLength, // if == 0 --> no check
				const OutMode& outMode,
				const size_t maxStoreSize 
			)
{
	assertbiu(energy!=NULL , "no energy function given");
	
	  // 3D points of structure that is currently processed
	biu::IPointVec actStruct;
	  // relative move sequence of structure that is currently processed
	biu::MoveSequence actMoves;
	  // energy of structure that is currently processed
	double actE;

	// width of the progress bar to print
	const size_t PROGRESS_BAR_WIDTH = 64;
	  // size of the already printed progress bar
	size_t progressBarSize = 0;
	
	if (outMode == OUT_NORMAL) {
		  // print progress bar limits
		std::cout 	<<"\n"
					<<"  Simulating co-translational folding via chain-growth :\n"
					<<"\n"
					<<"    |0%" <<std::setw(PROGRESS_BAR_WIDTH-3) <<100 <<"%|\n"
					<<"    |";
		std::cout.flush();
	}

	// do elongation
	for (; curIdx < seq.size(); curIdx++) {
		  // get structure store access
		size_t toFill = curIdx%2;
		size_t toRead = (curIdx-1)%2;
		  // clear structure store that will be filled now
		structStore[toFill]->clear();
		
		assertbiu(structStore[toRead]->size() > 0, "nothing in the structure store to extend");
		  // generate all extend versions of structure from last iteration
		for (	StructStore::const_iterator act=structStore[toRead]->begin();
				act != structStore[toRead]->end(); act++)
		{
			  // get energy of current sub structure
			actE = act->first;
			  // get relative move string
			actMoves = lattice.getDescriptor()->getAlphabet()->decompress(act->second, curIdx-1);
			  // get 3D points of current relative move string from store
			actStruct = lattice.absMovesToPoints( actMoves );
			  // extend structure with dummy monomer
			actStruct.push_back(biu::IntPoint(0,0,0));
			actMoves.push_back(lattice.getDescriptor()->getAlphabet()->getElement(0));
			  // try all possible extensions
			for (	biu::LatticeNeighborhood::const_iterator actNeigh = lattice.getNeighborhood().begin(); 
					actNeigh != lattice.getNeighborhood().end(); actNeigh++ )
			{
				bool isSelfavoiding = true;
				double energyGain = 0.0;
				  // get correct last position according to current neighbor vector
				actStruct[curIdx] = actStruct[curIdx-1] + *actNeigh;
				
				// do selfavoiding check and calculate energy gain
				for (size_t i = 0; i < curIdx-1; i++) {
					  // selfavoidingness check
					if ( actStruct[i] == actStruct[curIdx] ) {
						isSelfavoiding = false;
						break;
					}
					  // get distance dependent energy contribution
					energyGain += energy->getEnergy(	seq[i], 
														seq[curIdx],
														actStruct[i],
														actStruct[curIdx]);
				}
				double tmpE = actE + energyGain;
				
				  // store new structure if selfavoiding and energy is in energy band to store
				if (	isSelfavoiding 
						&&	(	structStore[toFill]->empty()
								|| (structStore[toFill]->begin()->first + deltaE) >= tmpE)
					) 
				{
					  // check if a new lower energy bound has been found
					bool shrinkStore = (!structStore[toFill]->empty())
										&& (tmpE < structStore[toFill]->begin()->first);
					bool isExtensible = true;
					
					  // check extensibility if
					  // - first is added or new lower energy bound has been found
					  // - we have to do
					  // - minimal sequence length has been reached
					  // - we are not adding the last monomer
					if (	shrinkStore
							&& 	minExtCheckLength != 0  
							&&	curIdx >= minExtCheckLength  
							&&	(curIdx+1) < seq.size()
						)
					{
						isExtensible = hasFreeNeighbor( actStruct, curIdx, lattice);
					}

					if (isExtensible) {
						// set correct move
						actMoves[actMoves.size()-1] = actNeigh->getMove();
						// add to storage
						structStore[toFill]->insert( ES_pair(actE + energyGain, 
										lattice.getDescriptor()->getAlphabet()->compress(
												lattice.getDescriptor()->normalizeSequence(actMoves)) 
										));
						
						  // check if new energy is smaller than best so far
						if (shrinkStore) {
							  // remove all above maximal allowed energy
							double maxAllowedE = tmpE + deltaE;
							StructStore::iterator delPos = structStore[toFill]->begin();
							  // find pos from where to delete from the back 
							while (delPos != structStore[toFill]->end()) {
								if (delPos->first <= maxAllowedE) {
									delPos++;
								} else {
									break;
								}
							}
							  // do deletion beginning from delPos up to the end
							structStore[toFill]->erase(delPos, structStore[toFill]->end());
						} else {
							if (structStore[toFill]->size() > maxStoreSize) {
								std::cout	<<"\n\n --> Maximal number of structures to extend exceeded for monomer " <<(curIdx+1)
											<<"! Chain growth stopped.\n" <<std::endl;
								return false;
							}
						}
					}
				}
			}
		}
		if (structStore[toFill]->size() == 0) {	// appending failed
			std::cout	<<"\n\n --> Appending of monomer " <<(curIdx+1)
						<<" failed! Chain growth stopped.\n" <<std::endl;
			return false;
		}
		if (outMode == OUT_VERBOSE) {
			if (curIdx+1 < seq.size())
				std::cout <<std::setw(5)<<(curIdx+2)
					<<". monomer appending to " <<std::setw(12) <<structStore[toFill]->size()
					<<" structures with minimal E = "
					<<std::fixed <<std::setprecision(OUT_PRECISION)
					<<std::setw(10+OUT_PRECISION) <<structStore[toFill]->begin()->first 
					<<std::endl;
			else
				std::cout 
					<<"\n       ==>   ending up with " <<std::setw(12) <<structStore[toFill]->size()
					<<" structures with minimal E = "
					<<std::fixed <<std::setprecision(OUT_PRECISION)
					<<std::setw(10+OUT_PRECISION) <<structStore[toFill]->begin()->first 
					<<"\n                                        "
					<<"             and maximal E = "
					<<std::fixed <<std::setprecision(OUT_PRECISION)
					<<std::setw(10+OUT_PRECISION) <<structStore[toFill]->rbegin()->first 
					<<std::endl;
		} else 	if (outMode == OUT_NORMAL) {
			  // do progress bar output
			  // calculate necessary progress bar extension
			size_t progExt = ((curIdx * PROGRESS_BAR_WIDTH) / (seq.size()-1)) - progressBarSize;
			if (progExt > 0) {
				progressBarSize += progExt;
				  // print bar elements 
				while ( progExt-- ) {
					std::cout	<<'#';
				}
				std::cout.flush();
			}
		}		
	}


	if (outMode == OUT_NORMAL) {
		std::cout <<"|\n" <<std::endl;
	}

	return true;
}


//////////////////////////////////////////////////////////////////////////


/**
 * Implements a chain-growth algorithm for side chain lattice models that 
 * extends all structures within a 
 * given energy interval deltaE from the minimal energy reachable so far. All
 * other structures are rejected. This can lead to a non-extensible chain. 
 * Therefore, a extensibility check can be done that validates if the last 
 * monomer position has at least one free neighbored position. This is 
 * sufficient in most cases but might still be not enough for long sequences.
 * 
 * @param seq the sequence of the lattice protein
 * @param structStore (IN-OUT) container that store structures from the last 
 *                    extension and the current one
 * @param curIdx (IN-OUT) the monomer index to append in the next iteration
 * @param lattice the lattice model to use
 * @param energyC the contact energy function to use if not NULL 
 * @param energyD the distance based energy function used if 
 *                energyC is equal to NULL
 * @param deltaE the energy interval used to reject the extension of structures
 *               that have an energy higher than the minimal energy so far plus
 *               deltaE
 * @param minExtCheckLength the minimal chain length that makes an extensibility
 *                          check necesary. If == 0 than NO check is done.
 * @param outMode the output level
 * @param maxStoreSize maximal number of structures to extend
 * @return TRUE if the chain was extended to full length, FALSE if the algorithm
 *         got stucked.  
 * 
 */
bool
chainGrowthSC(	const biu::Alphabet::Sequence & seq,
				std::vector<StructStore*> & structStore,
				size_t & curIdx,
				const biu::LatticeModel & lattice,
				const biu::DistanceEnergyFunction * const energy,
				const double deltaE,
				const size_t minExtCheckLength, // if == 0 --> no check
				const OutMode& outMode,
				const size_t maxStoreSize
			)
{
	assertbiu(energy!=NULL, "no energy function given");
	
	  // 3D points of structure that is currently processed
	biu::IPointVec actStruct;
	  // relative move sequence of structure that is currently processed
	biu::MoveSequence actMoves;
	  // energy of structure that is currently processed
	double actE;

	// width of the progress bar to print
	const size_t PROGRESS_BAR_WIDTH = 64;
	  // size of the already printed progress bar
	size_t progressBarSize = 0;
	
	if (outMode == OUT_NORMAL) {
		  // print progress bar limits
		std::cout 	<<"\n"
					<<"  Simulating co-translational folding via chain-growth :\n"
					<<"\n"
					<<"    |0%" <<std::setw(PROGRESS_BAR_WIDTH-3) <<100 <<"%|\n"
					<<"    |";
		std::cout.flush();
	}

	// do elongation
	for (; curIdx < seq.size(); curIdx++) {
		  // get structure store access
		size_t toFill = curIdx%2;
		size_t toRead = (curIdx-1)%2;
		  // clear structure store that will be filled now
		structStore[toFill]->clear();
		
		assertbiu(structStore[toRead]->size() > 0, "nothing in the structure store to extend");
		  // generate all extend versions of structure from last iteration
		for (	StructStore::const_iterator act=structStore[toRead]->begin();
				act != structStore[toRead]->end(); act++)
		{
			  // current index in the side chain representation
			size_t actIdx = curIdx*2;
			  // get energy of current sub structure
			actE = act->first;
			  // get relative move string
			actMoves = lattice.getDescriptor()->getAlphabet()->decompress(act->second, actIdx-1);
			// get 3D points of current relative move string from store
			actStruct = abs2pointsSC( actMoves, actStruct, lattice );
			// extend structure with dummy monomer for the next backbone position
			actStruct.push_back(biu::IntPoint(0,0,0));
			actMoves.push_back(lattice.getDescriptor()->getAlphabet()->getElement(0));
			  // extend structure with dummy monomer for the next side chain position
			actStruct.push_back(biu::IntPoint(0,0,0));
			actMoves.push_back(lattice.getDescriptor()->getAlphabet()->getElement(0));
			  // try all possible placements of the next backbone monomer
			for (	biu::LatticeNeighborhood::const_iterator actNeighB = lattice.getNeighborhood().begin(); 
					actNeighB != lattice.getNeighborhood().end(); actNeighB++ )
			{
				bool isSelfavoidingB = true;
				  // get correct last position according to current neighbor vector
				actStruct[actIdx] = actStruct[actIdx-2] + *actNeighB;
				
				// do selfavoiding check and calculate energy gain
				for (size_t i = 0; i < actIdx; i++) {
					  // selfavoidingness check
					if ( actStruct[i] == actStruct[actIdx] ) {
						isSelfavoidingB = false;
						break;
					}
				}
				  // try next neighboring position if this one is not selfavoiding
				if (!isSelfavoidingB)
					continue;
				
				  // set correct move
				actMoves[actIdx-1] = actNeighB->getMove();
				
				  // try all possible placements of the next side chain monomer 
				for (	biu::LatticeNeighborhood::const_iterator actNeighSC = lattice.getNeighborhood().begin(); 
						actNeighSC != lattice.getNeighborhood().end(); actNeighSC++ )
				{
					size_t isSelfavoidingSC = true;
					double energyGain = 0.0;
					  // get correct last position according to current neighbor vector
					actStruct[actIdx+1] = actStruct[actIdx] + *actNeighSC;
					
					// do selfavoiding check and calculate energy gain
					for (size_t i = 0; i < actIdx; i++) {
						  // selfavoidingness check
						if ( actStruct[i] == actStruct[actIdx+1] ) {
							isSelfavoidingSC = false;
							break;
						}
						  // if i is sidechain position
						if (i%2==1) {
							  // get distance dependent energy contribution
							energyGain += energy->getEnergy(	seq[i/2], 
																seq[curIdx],
																actStruct[i],
																actStruct[actIdx+1]);
						}
					}
					double tmpE = actE + energyGain;
					
					// set correct move
					actMoves[actIdx] = actNeighSC->getMove();
					
					  // store new structure if selfavoiding and energy is in energy band to store
					if (	isSelfavoidingB && isSelfavoidingSC 
							&&	(	structStore[toFill]->empty()
									|| (structStore[toFill]->begin()->first + deltaE) >= tmpE)
						) 
					{
						  // check if a new lower energy bound has been found
						bool shrinkStore = (!structStore[toFill]->empty())
											&& (tmpE < structStore[toFill]->begin()->first);
						bool isExtensible = true;
						
						  // check extensibility if
						  // - first is added or new lower energy bound has been found
						  // - we have to do
						  // - minimal sequence length has been reached
						  // - we are not adding the last monomer
						if (	shrinkStore 
								&& 	minExtCheckLength != 0  
								&&	actIdx >= minExtCheckLength  
								&&	(curIdx+1) < seq.size()
							)
						{
							isExtensible = hasFreeNeighbor( actStruct, actIdx, lattice);
						}
	
						if (isExtensible) {
							// add to storage
							structStore[toFill]->insert( ES_pair(actE + energyGain, 
											lattice.getDescriptor()->getAlphabet()->compress(
													lattice.getDescriptor()->normalizeSequence(actMoves)) 
											));
							
							  // check if new energy is smaller than best so far
							if (shrinkStore) {
								  // remove all above maximal allowed energy
								double maxAllowedE = tmpE + deltaE;
								StructStore::iterator delPos = structStore[toFill]->begin();
								  // find pos from where to delete from the back 
								while (delPos != structStore[toFill]->end()) {
									if (delPos->first <= maxAllowedE) {
										delPos++;
									} else {
										break;
									}
								}
								  // do deletion beginning from delPos up to the end
								structStore[toFill]->erase(delPos, structStore[toFill]->end());
							} else {
								if (structStore[toFill]->size() > maxStoreSize) {
									std::cout	<<"\n\n --> Maximal number of structures to extend exceeded for monomer " <<(curIdx+1)
												<<"! Chain growth stopped.\n" <<std::endl;
									return false;
								}
							}
						}
					}
				}
			}
		}
		if (structStore[toFill]->size() == 0) {	// appending failed
			std::cout	<<"\n\n --> Appending of monomer " <<(curIdx+1)
						<<" failed! Chain growth stopped.\n" <<std::endl;
			return false;
		}
		if (outMode == OUT_VERBOSE) {
			if (curIdx+1 < seq.size())
				std::cout <<std::setw(5)<<(curIdx+2)
					<<". monomer appending to " <<std::setw(12) <<structStore[toFill]->size()
					<<" structures with minimal E = "
					<<std::fixed <<std::setprecision(OUT_PRECISION)
					<<std::setw(10+OUT_PRECISION) <<structStore[toFill]->begin()->first 
					<<std::endl;
			else
				std::cout 
					<<"\n       ==>   ending up with " <<std::setw(12) <<structStore[toFill]->size()
					<<" structures with minimal E = "
					<<std::fixed <<std::setprecision(OUT_PRECISION)
					<<std::setw(10+OUT_PRECISION) <<structStore[toFill]->begin()->first 
					<<"\n                                        "
					<<"             and maximal E = "
					<<std::fixed <<std::setprecision(OUT_PRECISION)
					<<std::setw(10+OUT_PRECISION) <<structStore[toFill]->rbegin()->first 
					<<std::endl;
		} else 	if (outMode == OUT_NORMAL) {
			  // do progress bar output
			  // calculate necessary progress bar extension
			size_t progExt = ((curIdx * PROGRESS_BAR_WIDTH) / (seq.size()-1)) - progressBarSize;
			if (progExt > 0) {
				progressBarSize += progExt;
				  // print bar elements 
				while ( progExt-- ) {
					std::cout	<<'#';
				}
				std::cout.flush();
			}
		}		
	}

	if (outMode == OUT_NORMAL) {
		std::cout <<"|\n" <<std::endl;
	}

	return true;
}


//////////////////////////////////////////////////////////////////////////

void
initAllowedArguments(biu::OptionMap & allowedArgs, std::string &infoText )  
{
	allowedArgs.push_back(biu::COption(	
							"seq", false, biu::COption::STRING, 
							"sequence of the protein conform to the alphabet from the 'energyFile'"));
	allowedArgs.push_back(biu::COption(	
							"energyFile", false, biu::COption::STRING, 
							"the energy function file that contains also the alphabet"));
	allowedArgs.push_back(biu::COption(	
							"energyForDist", true, biu::COption::BOOL, 
							"if present a distance based energy function is used, otherwise a contact energy function is applied"));
	allowedArgs.push_back(biu::COption(	
							"energyCalphaDist", true, biu::COption::DOUBLE, 
							"if '-energyForDist' is present, this value is used to scale the C_alpha distances of the given energy function to the C_alpha monomer distances in the used lattice model",
							"3.8"));
	allowedArgs.push_back(biu::COption(	
							"sideChain", true, biu::COption::BOOL, 
							"use a side chain model and calculate energy for side chain monomers only"));
	allowedArgs.push_back(biu::COption(	
							"lat", true, biu::COption::STRING, 
							"lattice model to use (SQR, CUB, FCC)", "FCC"));
	allowedArgs.push_back(biu::COption(	
							"deltaE", true, biu::COption::DOUBLE, 
							"energy band of structure to keep from minimal energy possible (>=0)", "0.0"));
	allowedArgs.push_back(biu::COption(	
							"noExtCheck", true, biu::COption::BOOL, 
							"if NO extensibility check of the last monomer should be done (faster)"));
	allowedArgs.push_back(biu::COption(	
							"extMax", true, biu::COption::INT, 
							"extend maximally X structures per iteration", "1000000"));
	allowedArgs.push_back(biu::COption(	
							"print", true, biu::COption::INT, 
							"print X structures in absolute moves and increasing energy", "10"));
	allowedArgs.push_back(biu::COption(	
							"printBest", true, biu::COption::BOOL, 
							"print only structures with minimal energy"));
	allowedArgs.push_back(biu::COption(	
							"best2PDB", true, biu::COption::STRING, 
							"prints the best structure in PDB format to the given file name"));
	allowedArgs.push_back(biu::COption(	
							"v", true, biu::COption::BOOL, 
							"verbose output"));
	allowedArgs.push_back(biu::COption(	
							"s", true, biu::COption::BOOL, 
							"minimal output only"));
	allowedArgs.push_back(biu::COption(	
							"help", true, biu::COption::BOOL, 
							"program parameters and help"));
	allowedArgs.push_back(biu::COption(
							"version", true, biu::COption::BOOL,
							"version information of this program"));

	infoText =	std::string("LatSeF implements a chain-growth algorithm that simulates 'Sequential Folding' of lattice proteins.\n")
			+	std::string("\n")
			+	std::string("The application on side chain models is supported.")
			+	std::string(" Here, the backbone monomers have no energy contribution")
			+	std::string(" and only the contacts of side chain monomers are considered.\n")
			+	std::string("\n")
			+	std::string("The energy function is contact based and has to be")
			+	std::string(" provided as a text file that contains the alphabet")
			+	std::string(" and the energy contributions in matrix form.\n")
			+	std::string("\n")
			+	std::string("An example energy file : \n")
			+	std::string("------------------------------\n")
			+	std::string("HPNX\n")
			+	std::string("-4.0  0.0  0.0  0.0\n")
			+	std::string(" 0.0 +1.0 -1.0  0.0\n")
			+	std::string(" 0.0 -1.0 +1.0  0.0\n")
			+	std::string(" 0.0  0.0  0.0  0.0\n")
			+	std::string("------------------------------\n")
			+	std::string("\n") ;

} // initArguments


//////////////////////////////////////////////////////////////////////////

/**
 * main procedure 
 */
int mainVec( int argc, char** argv ) {
	using namespace std;
	
	//////////////////////////////////////////////////////////////
	// data to fill
	//////////////////////////////////////////////////////////////
	
	biu::Alphabet * alph = NULL;				//!< alphabet used
	std::string alphString = "";				//!< allowed characters of the alphabet
	biu::EnergyMatrix * energyMatrix = NULL;	//!< energy matrix of contact energy function
	biu::DistanceEnergyFunction * energy = NULL;	//!< distance based energy function to use
	biu::LatticeDescriptor* latDescr = NULL;	//!< descriptor of the lattice model
	biu::Alphabet::Sequence seq;			//!< sequence to fold sequentially
	double deltaE = 0.0;		//!< the energy band used to remove structures from storage
	size_t minExtCheckLength = 3;	//!< minimal structure length that makes a extensionable check necessary
	size_t toPrint = 5;			//!< maximal number of best structures to print in the end
	bool printBestOnly = false;	//!< print only structures with minimal energy
	OutMode outMode = OUT_NORMAL;	//!< output mode
	bool sideChain = false;			//!< use a side chain model
	std::string pdbFile = "";		//!< PDB file to write 
	size_t maxStoreSize = UINT_MAX;	//!< maximal number of structure to extend in next iteration
	double cAlphaDist = 3.8;		//!< the C_alpha distance used to scale the distances of the energy function
	
	//////////////////////////////////////////////////////////////
	// parameter parsing and checking
	//////////////////////////////////////////////////////////////
	
	biu::OptionMap allowedArgs;	//!< map of allowed arguments
	std::string infoText;		//!< info string of the program
	initAllowedArguments(allowedArgs,infoText);	// init
	
		// parse programm arguments	
	biu::COptionParser opts = biu::COptionParser(	allowedArgs, argc, 
													argv, infoText);
		// arguments parseable and all mandatory arguments given
	if (opts.noErrors()) {
			// help output
		if (opts.getBoolVal("help")) {
			opts.coutUsage();
			return 0;
		}
		if (opts.getBoolVal("version")) {
			giveVersion();
			return 0;
		}
	} else {
		return -1;
	}
	
	
	if (opts.argExist("energyCalphaDist")) {
		double tmp = opts.getDoubleVal("energyCalphaDist");
		if (tmp <= 0.0) {
			cerr	<<"\n   ERROR : Given C_alpha distance '" <<tmp
					<<"' is smaller than or equal to ZERO\n\n";
			return -1;
		}
		cAlphaDist = tmp; 
	}
	
	  // init lattice model
	std::string latStr = opts.getStrVal("lat");
	if (latStr.compare("SQR") == 0) {
		latDescr = new biu::LatticeDescriptorSQR();
		minExtCheckLength = 7;
	} else if (latStr.compare("CUB") == 0) {
		latDescr = new biu::LatticeDescriptorCUB();
		minExtCheckLength = 11;
	} else if (latStr.compare("FCC") == 0) {
		latDescr = new biu::LatticeDescriptorFCC();
		minExtCheckLength = 11;
	} else {
		cerr <<"\n   ERROR : Unknown lattice type '"+latStr+"'\n\n";
		return -1;
	}
	biu::LatticeModel lattice(latDescr);
	
		// init energy function
	std::string energyFile = opts.getStrVal("energyFile");
	if ( energyFile.size() == 0 ) {
		cerr <<"\n   Error: no energy file given ('-energyFile=XXX') !\n";
		return -1;
	}
	{
		  // temporary data structures
		std::ifstream *inFile = NULL;
		std::istream* in = &std::cin;
		  // open stream if file given
		if (energyFile.compare("STDIN") != 0) {
			inFile = new std::ifstream( energyFile.c_str() );
			if (!inFile->is_open()) {
				std::cerr <<"\n   ERROR : can not open energy file '"+energyFile+"' !\n\n";
				return -1;
			}
			in = inFile;
		}
		if (in->bad()) {
			cerr <<"\n   Error: cannot read energy function from '" <<energyFile <<"' !\n";
			return -1;
		}
		if (opts.argExist("energyForDist")) {
		/////// DISTANCE BASED ENERGY FUNCTION ///////////////////////////////////
	
			double cAlphaDistScale = (lattice.getNeighborhood().getElementByIndex(0).distance(biu::IntPoint(0,0,0)))
										/ cAlphaDist;
			
			// do parsing
			if (initIntervalEnergyFunction( alph, energy, cAlphaDistScale, *in) != 0) {
				cerr <<"\n   Error: the given energy file '"<<energyFile <<"' is not valid !\n";
				return -1;
			}
			
		} else {
		/////// CONTACT BASED ENERGY FUNCTION ////////////////////////////////////
			  // do parsing
			if (initContactEnergyFunction( alph, energyMatrix, *in) != 0) {
				cerr <<"\n   Error: the given energy file '"<<energyFile <<"' is not valid !\n";
				return -1;
			}
			  // create energy function
			energy = new biu::ContactEnergyFunction(alph, energyMatrix, &lattice);
		}
		  // close stream if necessary
		if (energyFile.compare("STDIN") != 0) {
			in = &std::cin;
			inFile->close();
			delete inFile;
		}
	}
	// init string representation of the used alphabet for parameter output
	{
		biu::Alphabet::Sequence tmp;
		for (size_t i = 0; i < alph->getAlphabetSize(); i++) {
			tmp.push_back(alph->getElement(i));
		}
		alphString = alph->getString(tmp);
	}
		
	
	  // read sequence and check if conform to alphabet 
	{
		std::string seqStr = opts.getStrVal("seq");
		if (seqStr.size() < MIN_SEQ_LENGTH) {
			cerr	<<"\n   ERROR : minimal sequence length is "<<MIN_SEQ_LENGTH <<" \n\n";
			return -1;
		}
		if (!alph->isAlphabetString(seqStr)) {
			cerr	<<"\n   ERROR : Given sequence '" <<seqStr
					<<"' is not conform to alphabet from energy file '" <<energyFile <<"'\n\n";
			return -1;
		}
		seq = alph->getSequence(seqStr);
	}
	
	if (opts.argExist("deltaE")) {
		double tmp = opts.getDoubleVal("deltaE");
		if (tmp < 0.0) {
			cerr	<<"\n   ERROR : Given deltaE '" <<tmp
					<<"' is smaller than ZERO\n\n";
			return -1;
		}
		deltaE = tmp; 
	}
	
	if (opts.argExist("print")) {
		int tmp = opts.getIntVal("print");
		if (tmp < 0) {
			cerr	<<"\n   ERROR : Given number of structures to print '" <<tmp
					<<"' is smaller than ZERO\n\n";
			return -1;
		}
		toPrint = (size_t)tmp;
	}
	
	if (opts.argExist("extMax")) {
		int tmp = opts.getIntVal("extMax");
		if (tmp < 0) {
			cerr	<<"\n   ERROR : Given maximal number of structures to extend '" <<tmp
					<<"' is smaller than ZERO\n\n";
			return -1;
		}
		maxStoreSize = (size_t)tmp;
	}
	
	  // print only structures with best energy ?
	printBestOnly = opts.getBoolVal("printBest");
	  // extensibility check of last monomer ?
	if (opts.getBoolVal("noExtCheck"))
		minExtCheckLength = 0;
	  // determine output mode
	if (opts.getBoolVal("s") && opts.getBoolVal("v")) {
		cerr	<<"\n   ERROR : Either SILENT or VERBOSE output possible ;)\n\n";
		return -1;
	}
	  // verbose output ?
	if (opts.getBoolVal("v"))
		outMode = OUT_VERBOSE;
	  // silent output ?
	if (opts.getBoolVal("s"))
		outMode = OUT_SILENT;
	  // side chain model ?
	sideChain = opts.getBoolVal("sideChain");
	
	  // PDB file ?
	if (opts.argExist("best2PDB")) {
		if (opts.getStrVal("best2PDB").size() == 0) {
			cerr	<<"\n   ERROR : no PDB file name given\n\n";
			return -1;
		}
		pdbFile = opts.getStrVal("best2PDB");
	}
	
	//////////////////////////////////////////////////////////////
	// input parameter output 
	//////////////////////////////////////////////////////////////
	
	  // string stream for parameter output in results stream
	std::ostringstream infostream(std::ostringstream::out);
	
	if (outMode != OUT_SILENT) std::cout <<"\n  energy-file     = " <<energyFile;
	infostream <<" energy-file=" <<energyFile;
	if (outMode != OUT_SILENT) std::cout <<"\n  contact-energy  = " <<(opts.argExist("energyForDist")?"false":"true");
	infostream <<" contact-energy=" <<(opts.argExist("energyForDist")?"false":"true");
	if (outMode != OUT_SILENT) std::cout <<"\n  alphabet        = " <<alphString;
	infostream <<" alphabet=" <<alphString;
	if (outMode != OUT_SILENT) std::cout <<"\n  sequence        = " <<alph->getString(seq);
	infostream <<" sequence=" <<alph->getString(seq);
	if (outMode != OUT_SILENT) std::cout <<"\n  lattice         = " <<latDescr->getName();
	infostream <<" lattice=" <<latDescr->getName();
	if (outMode != OUT_SILENT) std::cout <<"\n  delta-E         = " <<deltaE;
	infostream <<" deltaE=" <<deltaE;
	if (outMode != OUT_SILENT) std::cout <<"\n  extend max.     = " <<maxStoreSize;
	infostream <<" max-to-extend=" <<maxStoreSize;
	if (outMode != OUT_SILENT) std::cout <<"\n  extens. check   = " <<(minExtCheckLength>0?"true":"false");
	infostream <<" ext-check=" <<(minExtCheckLength>0?"true":"false");
	if (outMode != OUT_SILENT) std::cout <<"\n  print max.      = " <<toPrint;
	infostream <<" print-max=" <<toPrint;
	if (outMode != OUT_SILENT) std::cout <<"\n  print best only = " <<(printBestOnly?"true":"false");
	infostream <<" print-best-only=" <<(printBestOnly?"true":"false");
	
	if (outMode != OUT_SILENT) std::cout <<"\n" <<std::endl;

	//////////////////////////////////////////////////////////////
	// initialise sequential folding
	//////////////////////////////////////////////////////////////
	
	  // stores sub structures with corresponding energy
	std::vector<StructStore*> structStore = std::vector<StructStore*>(2);
	structStore[0] = new StructStore();
	structStore[1] = new StructStore();
	
	  // length of the current sub structures in the StructStore
	size_t curIdx = sideChain ? 0 : 1;
	const biu::Move FIRST_MOVE =  latDescr->getAlphabet()->getElement(0);
	  // initialise structure store with substructure of length 2 with first move
	structStore[curIdx%2]->insert( ES_pair(0.0, latDescr->getAlphabet()->compress(
												 latDescr->getAlphabet()->getSequence(
												  latDescr->getAlphabet()->getString(
												   FIRST_MOVE )))) );
	  // increase index of current monomer index to append
	curIdx++;

	//////////////////////////////////////////////////////////////
	// perform sequential folding
	//////////////////////////////////////////////////////////////
	
	  // do chain elongation
	bool wasSuccessful = false;
	if (sideChain)
		wasSuccessful = chainGrowthSC( seq, structStore, curIdx, lattice, energy, deltaE, minExtCheckLength, outMode, maxStoreSize);
	else
		wasSuccessful = chainGrowth( seq, structStore, curIdx, lattice, energy, deltaE, minExtCheckLength, outMode, maxStoreSize);


	//////////////////////////////////////////////////////////////
	// result output
	//////////////////////////////////////////////////////////////
	
	  // than normal exit of the chain growth algorithm and no abortion 
	if (wasSuccessful) {
		
		const double minE = ( structStore[(curIdx-1)%2]->empty() ? 0.0 : structStore[(curIdx-1)%2]->begin()->first);
		
		  // get number of best structures only
		size_t numOfBest = 0;
		for (StructStore::const_iterator actStr=structStore[(curIdx-1)%2]->begin();
				actStr != structStore[(curIdx-1)%2]->end() && minE == actStr->first;
				actStr++)
		{
			numOfBest++;
		}

		// result summary
		std::cout 
			<<"\n number of structures      = " <<std::setw(12) <<structStore[(curIdx-1)%2]->size()
			<<"\n minimal energy structures = " <<std::setw(12) <<numOfBest
			<<"\n minimal energy found      = " <<std::fixed <<std::setprecision(OUT_PRECISION) <<std::setw(13+OUT_PRECISION) <<structStore[(curIdx-1)%2]->begin()->first
			<<std::endl;
		
		  // print structures
		size_t printed = toPrint;
		for (StructStore::const_iterator actStr=structStore[(curIdx-1)%2]->begin();
				actStr != structStore[(curIdx-1)%2]->end() && printed-- > 0;
				actStr++) 
		{
			  // check if all structures with minimal energy have been printed
			if (printBestOnly && minE < actStr->first)
				break;
			if (sideChain)
				std::cout
					<<"\n" <<abs2absSC(lattice.getString(latDescr->getAlphabet()->decompress(actStr->second, (2*seq.size())-1)), lattice)
					<<"  " <<std::fixed <<std::setprecision(OUT_PRECISION) <<actStr->first;
			else
				std::cout
					<<"\n" <<lattice.getString(latDescr->getAlphabet()->decompress(actStr->second, seq.size()-1))
					<<"  " <<std::fixed <<std::setprecision(OUT_PRECISION) <<actStr->first;
		}
		std::cout <<std::endl;
		
		  // write PDB file if necessary
		if (pdbFile.size() > 0) {
			std::ofstream fout;
			fout.open(pdbFile.c_str());
			if (!fout.is_open()) {
				std::cerr <<"\n  ERROR : can not open PDB output file '" <<pdbFile <<"' !\n";
			} else {
				  // scale factor to allow a reasonable viewing
				const double LAT_SCALE = 3.8;
				  // get move string of best structure
				biu::MoveSequence absMoves;
				if (sideChain)
					absMoves = latDescr->getAlphabet()->decompress(structStore[(curIdx-1)%2]->begin()->second, (2*seq.size())-1);
				else
					absMoves = latDescr->getAlphabet()->decompress(structStore[(curIdx-1)%2]->begin()->second, seq.size()-1);
				  // get coordinates of that structure
				biu::IPointVec pointsI;
				if (sideChain)
					abs2pointsSC( absMoves, pointsI, lattice );
				else
					pointsI = lattice.absMovesToPoints(absMoves);
				biu::DPointVec pointsD;
				for (size_t i=0; i<pointsI.size(); i++) {
					pointsD.push_back(biu::DblPoint(pointsI[i]).operator *=(LAT_SCALE));
				}
				std::cout <<std::endl;
				  // write to File
				std::vector<std::string> aaSeq = std::vector<std::string>(seq.size(),"HIS");
				writePDB(	std::string("created using LatSeF"),
							aaSeq,
							pointsD, 
							lattice,
							fout,
							sideChain,
							lattice);  
				fout.close();
			}
		}
	}
	
	  // final output flush
	std::cout <<std::endl;
	
	//////////////////////////////////////////////////////////////
	// clear memory
	//////////////////////////////////////////////////////////////
	
	delete structStore[0];
	delete structStore[1];
	
	if (latDescr != NULL) {
		delete latDescr; latDescr = NULL;
	}
	if (energy != NULL) {
		delete energy; energy = NULL;
	}
	if (energyMatrix != NULL) { 
		delete energyMatrix; energyMatrix = NULL;
	}
	if (alph != NULL) { 
		delete alph; alph = NULL; 
	}


	return 0;	
}

/*

===================================
   HPstructCG - CPSP-tool-library
===================================


 HP-sequence      : PHHPPHHPPPPPHPPPHPPH
 lattice type     : sqr
 delta E          : 0

    3. monomer appending to            1 structures with minimal E =        0
    4. monomer appending to            2 structures with minimal E =        0
    5. monomer appending to            5 structures with minimal E =        0
    6. monomer appending to           13 structures with minimal E =        0
    7. monomer appending to            6 structures with minimal E =       -1
    8. monomer appending to            2 structures with minimal E =       -2
    9. monomer appending to            4 structures with minimal E =       -2
   10. monomer appending to           11 structures with minimal E =       -2
   11. monomer appending to           30 structures with minimal E =       -2
   12. monomer appending to           81 structures with minimal E =       -2
   13. monomer appending to          223 structures with minimal E =       -2
   14. monomer appending to           22 structures with minimal E =       -3
   15. monomer appending to           29 structures with minimal E =       -3
   16. monomer appending to           56 structures with minimal E =       -3
   17. monomer appending to          151 structures with minimal E =       -3
   18. monomer appending to          401 structures with minimal E =       -3
   19. monomer appending to         1092 structures with minimal E =       -3
   20. monomer appending to         2927 structures with minimal E =       -3


 Folded structures and energies :

FFFLBBBBRRFFRRFLFLB       -6

 number of structures =            1
 minimal energy found =           -6
 min. E structures    =            1

====================================
*/
