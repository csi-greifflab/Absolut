// $Id: RNAStructure.hh,v 1.2 2016/08/08 12:41:58 mmann Exp $
#ifndef BIU_RNASTRUCTURE_HH_
#define BIU_RNASTRUCTURE_HH_


#include "biu/BioMolecule.hh"
#include "biu/AllowedBasePairs.hh"

namespace biu
{
		/**
		 * An object of the class RNAStructure represents a RNA as a
		 * BioMolecule.
		 * 
		 * An object of this class is allowed to have only nested bonds.
		 * 
		 * @author Martin Mann, Sebastian Will, Andreas Richter
		 */
	class RNAStructure : public BioMolecule
	{
	protected:
		
			//! RNA sequence
		Sequence*		rnaSeq;
		
			//! whether or not the rnaSeq pointer is shared
		const bool		seqShared;
		
			//! Minimal length of an RNA loop
		static const size_t	MIN_LOOP_LENGTH;

			/*! The alphabet of the structure representation which contains
			 *  the left and right parenthesis and a dot.
			 */
		static const Alphabet			STRUCT_ALPH;

			//! RNA structure in bracket notation
		Structure*		rnaStructBracketDot;
		
			/*! Mapping from index of an opening bond to the
			 *  index of the corresponding closing bond
			 *  Entry at index i is index of the closing bond
			 *  and INVALID_INDEX if index i is not an opening bond
			 */
		std::vector<size_t> rnaBonds;
		
			/*! A pointer to an BasePair object containing the
			 *  valid RNA base pairs
			 */
		const AllowedBasePairs* 	bPair;
		
	public:
	
			//! structure element encoding a bond opening position
		static const Alphabet::AlphElem	STRUCT_BND_OPEN;
			//! structure element encoding a bond closing position
		static const Alphabet::AlphElem	STRUCT_BND_CLOSE;
			//! structure element encoding an unbound position
		static const Alphabet::AlphElem	STRUCT_UNBOUND;
		
			//! Constant that represents a invalid structure index position.
		static const size_t INVALID_INDEX;

			/*! Construction of an RNA molecule
			 * 
			 *  It is explicitly allowed to construct invalid
			 *  RNAStructures, but each bond has to be
			 *  complete (each opening has a closing bracket).
			 * 
			 * Using this constructor no sequence sharing is possible.
			 * 
			 * @param rnaSeqStr is the RNA sequence
			 * @param rnaStructBracketDotStr	is the RNA structure in bracket
			 *  							notation.
			 * @param bPair is the handler for the allowed base pairs.
			 */
		RNAStructure(	const std::string& rnaSeqStr, 
						const std::string& rnaStructBracketDotStr, 
						const AllowedBasePairs* const bPair);
			     
			/*! Construction of an RNA molecule
			 * 
			 *  It is explicitly allowed to construct invalid
			 *  RNAStructures, but each bond has to be
			 *  complete (each opening has a closing bracket).
			 * 
			 * @param rnaSeq the RNA sequence representation
			 * @param rnaStructBracketDot the RNA structure representation
			 *          in bracket notation.
			 * @param bPair is the handler for the allowed base pairs.
			 * @param seqIsShared whether or not the rnaSeqStr representation 
			 *          should be shared among all derived objects of this class
			 *          and no copy should be done
			 */
		RNAStructure(	Sequence* rnaSeq, 
						const Structure* const rnaStructBracketDot, 
						const AllowedBasePairs* const bPair,
						const bool seqIsShared);
			     
		RNAStructure(const RNAStructure& rnaStruct);
		
		virtual ~RNAStructure();
		
		RNAStructure&	operator= (const RNAStructure& rnaStruct2);
		bool	 		operator== (const RNAStructure& rnaStruct2) const;
		bool	 		operator!= (const RNAStructure& rnaStruct2) const;

			/*! Returns whether the RNA structure is in valid
			 *  bracket notation (each opening has a closing bracket).
			 */
		bool			hasValidStructure() const;
		
			/*! Returns whether all RNA base pairs of an
			 *  RNAStructure object are valid.
			 */
		bool			hasValidBasePairs() const;
		
			/*! Returns whether all loops of an RNAStructure object
			 *  have at least a length of MIN_LOOP_LENGTH.
			 */
		bool			hasMinLoopLength() const;
		
			/*! Returns whether a base pair between positions
			 *  first and second is valid.
			 */
		bool isAllowedBasePair(size_t first, size_t second) const;

			//! Returns  the corresponding closing bond index to the given
			//! opening index
			//! 
			//! Returns INVALID_INDEX in error case, if there doesnt exist 
			//! a closing bond.
		size_t getClosingBond(size_t openingBond) const { 
			return rnaBonds[openingBond]; 
		}
		
		size_t getMinLoopLength() const {
			return MIN_LOOP_LENGTH; 
		}

		std::string getSequenceString() const {
			return bPair->getAlphabet()->getString(*rnaSeq);
		}

		std::string getStructureString() const {
			return STRUCT_ALPH.getString(*rnaStructBracketDot);
		}

		static const Alphabet* getStructureAlphabet() {
			return &STRUCT_ALPH;
		}
		
	// abstract functions (BioMolecule)
	
		Sequence getSequence() const {
			assertbiu(rnaSeq != NULL, "no sequence available");
			return *rnaSeq;
		}

			//! Returns the length of the Biomolecule, i.e. the number of 
			//! monomers.
		size_t	getLength() const {
			assertbiu(rnaSeq != NULL, "no sequence available");
			return rnaSeq->size();
		}
	
			//! Returns the RNA structure in bracket notation.
		Structure getStructure() const {
			assertbiu(rnaStructBracketDot != NULL, "no structure available");
			return *rnaStructBracketDot;
		}
		
			//! Sets the RNA structure in bracket notation.
			//! @param str the structure to set in bracket notation of correct 
			//!         length
		virtual void setStructure(const Structure& str) {
			assertbiu(rnaSeq->size() == str.size(), "structure has wrong length");
			delete rnaStructBracketDot;
			rnaStructBracketDot = new Structure(str);
		}
		
		const Structure & getStructureRef() const {
			assertbiu(rnaStructBracketDot != NULL, "no structure available");
			return *rnaStructBracketDot;
		}
		
		
			/*! Returns true if an RNAStructure object has a
			 *  valid structure with valid base pairs and
			 *  loops of a length of at least MIN_LOOP_LENGTH.
			 */
		virtual bool 		isValid() const;
		
			/*! Returns a combination of bracket notation 
			 *  and sequence information of this RNA.
			 * 
			 *  BRACKETNOTATION(SEQUENCE) */
		virtual std::string	getStringRepresentation() const;
		
	protected:
		
		  //! Initialises the bonds data structure. (called from constructors)
		void initBonds();
		
	};
	
} // namespace biu

#endif /*RNASTRUCTURE_HH_*/
