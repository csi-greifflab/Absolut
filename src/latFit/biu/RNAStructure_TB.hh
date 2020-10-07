// $Id: RNAStructure_TB.hh,v 1.2 2016/08/08 12:42:00 mmann Exp $
#ifndef BIU_RNASTRUCTURE_TB_HH_
#define BIU_RNASTRUCTURE_TB_HH_


#include "biu/BioMolecule.hh"
#include "biu/AllowedBasePairs.hh"

namespace biu
{
	
		/**
		 * An object of the class RNAStructure_TB represents a RNA as a
		 * BioMolecule.
		 * 
		 * An object of this class is allowed to have only 
		 * 1- nested bonds
		 * 2- loop size greater than MIN_LOOP_LENGTH
		 * 3- certain matching base pairs
		 * -otherwise, isValid() returns false, and 
		 * 	methods other than getSequenceString, getStructureString should not be called
		 * 
		 * both MS1 (insertion & deletion) and MS2 (shifting) are supported
		 * 
		 * the validity of a move should be ensured seperatly (using the corresponding method)
		 * befor performing the move.
		 * 
		 * @author Salem Dekelbab
		 * @author Martin Mann
		 */
	class RNAStructure_TB : public BioMolecule
	{
	protected:
		
			// Before:
			// 	-A-B-C-D-E-F-G-H-I-J-K-L-M-N-O-P-Q-R-S-T-U-V-W-X-Y-Z-
			
			// After binding S with K
			//	-A-B-C-D-E-F-G-H-I-J-K-T-U-V-W-X-Y-Z-
			//	                     |
			//                      -S-L-M-N-O-P-Q-R-
			
			// 		- 	<=> 	next
			//		|	<=>		pair
		class TreeItem {
		public:
			
			size_t next;	//  next in sequence if there are no bonds, 
							//	next in the same level, where a new bond does not cause a non-nested structure
			size_t pair;	// (lower/upper) level when paired
							 
			
			bool 
			operator == (const TreeItem& op2) const {
				return ((this->next==op2.next)&&(this->pair==op2.pair));
			}
		};
		
		
		
		
	public:	//----------- abstract functions (BioMolecule)
			
		
			Sequence getSequence() const ;

			
			//! Returns the RNA structure in bracket notation.
			Structure getStructure() const ;
		
			/*! DUMMY IMPLEMENTATION THAT ALWAYS RETURNS 0.0 !!!
			 * @return 0.0
			 * */
			virtual double	getEnergy() const ;
	
			
			//! Returns the length of the Biomolecule, i.e. the number of 
			//! monomers.
			size_t	getLength() const ;
		
			
			/*! Returns true if an RNAStructure_TB object has a
				 *  valid structure with valid base pairs and
				 *  loops of a length of at least MIN_LOOP_LENGTH.
				 */
			bool isValid() const;
			
			
				/*! Returns a combination of bracket notation 
				 *  and sequence information of this RNA.
				 * 
				 *  BRACKETNOTATION(SEQUENCE) */
			std::string getStringRepresentation() const;
			
		
			//-----------(END) abstract functions (BioMolecule)
		
		
		
		
	protected:


			//! Minimal length of an RNA loop
		 static const size_t	MIN_LOOP_LENGTH;

		
			/*! The alphabet of the structure representation which contains
			 *  the left and right parenthesis and a dot.
			 */
		 static const Alphabet STRUCT_ALPH;
		
		
			//! whether or not the rnaSeq pointer is shared
		const bool		seqShared;

		
			/*! A pointer to a BasePair object containing the
			 *  valid RNA base pairs
			 */
		const AllowedBasePairs* 	bPair;

		
			// - Detailed status about the structure, i.e. why exactly the structure "is" invalid.
			//   Structure problem, Basepair problem or Loopsize problem 
			// - 0 if valid
			// - Status is used internally to determine whether an invalidity cause is determined or not yet
		mutable size_t strucStatus;
		
		
			//! RNA sequence
		Sequence*		rnaSeq;
		

			//! - RNA structure in bracket notation string
			//	- Storing it as string was prefered because the energy calculation requires a string,
			//		otherwise the Structure datatype could be used instead.  
		std::string bracketStructStr;
		

			//! RNA structure, tree representation, NULL if structure is invalid
		std::vector<TreeItem>* validTreeStruc; 
		
		
	public:
		
		
			//! structure element encoding a bond opening position
		static const Alphabet::AlphElem	STRUCT_BND_OPEN;
			//! structure element encoding a bond closing position
		static const Alphabet::AlphElem	STRUCT_BND_CLOSE;
			//! structure element encoding an unbound position
		static const Alphabet::AlphElem	STRUCT_UNBOUND;
		
			//! Constant that represents a invalid structure index position.
			//	i.e. means that it is not pair to any base
		static const size_t INVALID_INDEX;
		
		
			// Structure Status Bits 
		static const size_t BASEPAIRS_VALIDITY ; 			//000001
		static const size_t BASEPAIRS_VALIDITY_CALCULATED ;	//000010
		static const size_t LOOPSIZE_VALIDITY ;				//000100
		static const size_t LOOPSIZE_VALIDITY_CALCULATED ;	//001000
		static const size_t STRUCTURE_VALIDITY ;			//010000
		static const size_t STRUCTURE_VALIDITY_CALCULATED ;	//100000
		
		//	Do not change the following values without checking 
		//size_t moveType(size_t i, size_t j) const;
		static const size_t INSERT_MOVE;
		static const size_t SHIFT_MOVE;
		static const size_t REVERSE_SHIFT_MOVE;
		static const size_t INVALID_DELETE_MOVE;
		static const size_t DELETE_MOVE;
		
		
		
		// All global properties
		static const size_t getMinLoopLength() ;
		
		static const Alphabet* getStructureAlphabet() ;


		
			/*! Construction of an RNA molecule
			 * 
			 * Using this constructor no sequence sharing is possible.
			 * 
			 * @param rnaSeqStr is the RNA sequence
			 * @param rnaStructBracketDotStr	is the RNA structure in bracket
			 *  							notation.
			 * @param bPair is the handler for the allowed base pairs.
			 */
		RNAStructure_TB(	const std::string& rnaSeqStr, 
						const std::string& rnaStructBracketDotStr, 
						const AllowedBasePairs* const bPair);
			
		
			/*! Construction of an RNA molecule
			 * 
			 * @param rnaSeq the RNA sequence representation
			 * @param rnaStructBracketDot the RNA structure representation
			 *          in bracket notation.
			 * @param bPair is the handler for the allowed base pairs.
			 * @param seqIsShared whether or not the rnaSeqStr representation 
			 *          should be shared among all derived objects of this class
			 *          and no copy should be done
			 */
		RNAStructure_TB(	Sequence* rnaSeq, 
						const Structure* const rnaStructBracketDot, 
						const AllowedBasePairs* const bPair,
						const bool seqIsShared);
			     
				
		RNAStructure_TB(const RNAStructure_TB& rnaStruct);
				
		
		virtual ~RNAStructure_TB();
	
		
		RNAStructure_TB&	operator= (const RNAStructure_TB& rnaStruct2);
		bool	 			operator== (const RNAStructure_TB& rnaStruct2) const;
		bool	 			operator!= (const RNAStructure_TB& rnaStruct2) const;

		
		
			//! Sets the RNA structure in bracket notation.
			//! @param str the structure to set in bracket notation of correct 
			//!         length
		virtual void setStructure(const Structure& str);
		
		
		// String representation of the sequence
		std::string getSequenceString() const ;
		
		// String representation of the structure
		std::string getStructureString() const ;

		// Reference to string representation of the structure
		const std::string& getStructureStringRef() const ;

		
		//The following three methods calculate the return value only if not yet calculated (by initialization)
		// and only for once, unless structure is changed
		
		
		/*! Returns whether the RNA structure is in valid
			 *  bracket notation (each opening has a closing bracket).
			 */
		bool hasValidStructure() const ;
		
		
			/*! Returns whether all RNA base pairs of an
			 *  RNAStructure_TB object are valid.
			 */
		bool hasValidBasePairs() const ;
		
		
			/*! Returns whether all loops of an RNAStructure_TB object
			 *  have at least a length of MIN_LOOP_LENGTH.
			 */
		bool hasValidLoopSize() const ;
		
		
		// Single move operations, (insert, delete, shift)

		// Move validity check: i,j order does not matter, k is always the free base
		bool isValidInsertMove(const size_t i, const size_t j) const;
		bool isValidDeleteMove(const size_t i, const size_t j) const;
		
		/*!
		 * Checks if the given base pair (i,j) represents a valid insertion or
		 * deletion of a base pair.
		 * @return 0 if not valid, -1 if valid deletion, +1 if valid insertion
		 */
		int isValidSingleMove( const size_t i, const size_t j) const;
		
		bool isValidLeftShiftMove(const size_t i, const size_t j, const size_t k) const;
		bool isValidRightShiftMove(const size_t i, const size_t j, const size_t k) const;
		
	
		
		
		
		
		
		
		
		// Move Set 1: i,j order does not matter
		void insertBond(const size_t i, const size_t j);
		void deleteBond(const size_t i, const size_t j);
				
		//	Move Set 2: i,j are the already existing bond, k is the open base
		void leftShift(const size_t i, const size_t j, const size_t k);
		void rightShift(const size_t i, const size_t j, const size_t k);
		
		void shiftUnordered(const size_t i, const size_t j);
		
		void shiftToBond(const size_t i, const size_t k);
		
		void shiftToBond(const size_t i, const size_t j, const size_t k);
		
		
		
		size_t moveType(size_t i, size_t j) const ;
				
		void executeOrderedMove(size_t i, size_t j);

		bool executeOrderedMoveTry(size_t i, size_t j);
		
		
		/*!
		 * Calculates the next single move according to the internal order 
		 * 'after' the currently given single move (i,j).
		 * 
		 * To get the FIRST allowed single move call this function with i set
		 * to INVALID_INDEX. 
		 * 
		 * @param i IN/OUT the first position of the move to update
		 * @param j IN/OUT the second position of the move to update (>=i)
		 * @return whether or not another move was found and the given new (i,j)
		 *         is a valid move
		 */
		bool getNextSingleMove( size_t& i, size_t& j ) const ;
		
		
		
		//! Returns  the corresponding closing/opening bond index to the given
		//! opening/closing index
		//! 
		//! Returns INVALID_INDEX in error case, if there doesnt exist 
		//! a bond.
		size_t getCorrespondingBase(size_t basePos) const;
	
		
		/*! Returns whether a base pair between positions
		 *  first and second is valid.
		 */
		bool isAllowedBasePair(size_t pos1, size_t pos2) const;
		
		

		

	
	protected:

		
		  //! Initialises the bonds data structure. (called by setStructure)
		bool initTree(const Structure& str);
		
		
		// auxiliary function used by the move validity checking functions
		bool areOnSameLevel(const size_t i, const size_t j) const ;
		
		
		// auxiliary function used to set the brackets in the string representation 
		// for the not necessarily ordered i,j pair 
		void setStringBrackets(const size_t i, const size_t j);

		
		
		
		
		
	private:

		std::string decodedStatus() const;
		
		
		
		
	public:
		bool isValidUnorderedShift(const size_t i, const size_t j) const;
			
			
		bool isValidShift(const size_t i, const size_t j) const;
		
		
		void decodeMoveIndex(size_t const index,  size_t&  i,  size_t&  j) const;
						
		bool executeOrderedMoveTry(size_t const index);
		
		size_t getMoveIndexCount() const;

		
	};


} // namespace biu

#include "biu/RNAStructure_TB.icc"

#endif /*RNASTRUCTURE_TB_HH_*/
