// $Id: RNAStructure.cc,v 1.2 2016/08/08 12:41:57 mmann Exp $

#include <biu/RNAStructure.hh>
#include <stack>
#include <biu/assertbiu.hh>
#include <limits>

namespace biu
{
	const size_t RNAStructure::MIN_LOOP_LENGTH = 3;
	const size_t RNAStructure::INVALID_INDEX =
		std::numeric_limits<size_t>::max();
	const Alphabet RNAStructure::STRUCT_ALPH("().", 1);
	const Alphabet::AlphElem	RNAStructure::STRUCT_BND_OPEN = RNAStructure::STRUCT_ALPH.getElement("(");
	const Alphabet::AlphElem	RNAStructure::STRUCT_BND_CLOSE = RNAStructure::STRUCT_ALPH.getElement(")");
	const Alphabet::AlphElem	RNAStructure::STRUCT_UNBOUND = RNAStructure::STRUCT_ALPH.getElement(".");

//////////////////////////////////////////////////////////
// RNAStructure
//////////////////////////////////////////////////////////

	RNAStructure::RNAStructure(	const std::string& rnaSeqStr,
	 							const std::string& rnaStructBracketDotStr,
	 							const AllowedBasePairs* const _bPair)
	 :	rnaSeq(NULL),
	 	seqShared(false),
	 	rnaStructBracketDot(NULL),
	 	rnaBonds(),
		bPair(_bPair) 
	{
		assertbiu(bPair != NULL, "bPair is not allowed to be NULL.");
		assertbiu(bPair->getAlphabet()->isAlphabetString(rnaSeqStr),
		 "RNA sequence contains characters which are not in the alphabet.");
		assertbiu(STRUCT_ALPH.isAlphabetString(rnaStructBracketDotStr),
		 "RNA structure contains characters which are not in the alphabet.");

		
		rnaSeq = new Sequence(bPair->getAlphabet()->getSequence(rnaSeqStr));
		rnaStructBracketDot = new Structure(
		 STRUCT_ALPH.getSequence(rnaStructBracketDotStr));

		assertbiu(rnaSeq->size() == rnaStructBracketDot->size(),
		 "RNA sequence and structure have to have the same length.");
	
		  // init bond data structure
		initBonds();
	}

	RNAStructure::RNAStructure(	Sequence* rnaSeq_,
	 							const Structure* const rnaStructBracketDot_,
	 							const AllowedBasePairs* const bPair_,
	 							const bool seqIsShared)
	 :	rnaSeq(NULL),
	 	seqShared(seqIsShared),
	 	rnaStructBracketDot(NULL),
		bPair(bPair_) 
	{
		assertbiu(rnaSeq_ != NULL, "no RNA sequence given.");
		assertbiu(rnaStructBracketDot_ != NULL, "no RNA structure given.");
		assertbiu(bPair != NULL, "bPair is not allowed to be NULL.");
		assertbiu(bPair->getAlphabet()->isAlphabetSequence(*rnaSeq_),
		 "RNA sequence contains characters which are not in the alphabet.");
		assertbiu(STRUCT_ALPH.isAlphabetSequence(*rnaStructBracketDot_),
		 "RNA structure contains characters which are not in the alphabet.");
		assertbiu(rnaSeq_->size() == rnaStructBracketDot_->size(),
		 "RNA sequence and structure have to have the same length.");

		
		  // init sequence
		if (seqShared) {
			rnaSeq = rnaSeq_;
		} else {
			rnaSeq = new Sequence(*rnaSeq_);
		}
		  // init structure
		rnaStructBracketDot = new Structure(*rnaStructBracketDot_);

		  // init bond data structure
		initBonds();
	}

	RNAStructure::RNAStructure(const RNAStructure& rnaStruct)
	 :	rnaSeq(rnaStruct.rnaSeq),
		seqShared(rnaStruct.seqShared),
		rnaStructBracketDot(new Structure(*(rnaStruct.rnaStructBracketDot))),
		rnaBonds(rnaStruct.rnaBonds), bPair(rnaStruct.bPair)
	{
		assertbiu(bPair != NULL, "bPair is not allowed to be NULL.");
		assertbiu(rnaSeq->size() == rnaStructBracketDot->size(),
		 "RNA sequence and structure have to have the same length.");
		
		if (!seqShared) {
			rnaSeq = new Sequence(*(rnaStruct.rnaSeq));
		}
	}

	RNAStructure::~RNAStructure() {
		if(!seqShared && rnaSeq != NULL) {
			delete rnaSeq; rnaSeq = NULL;
		}
		if (rnaStructBracketDot != NULL) {
			delete rnaStructBracketDot; rnaStructBracketDot = NULL;
		}
	}

	bool 		
	RNAStructure::isAllowedBasePair(	size_t first,
			 							size_t second)  const 
	{
		if (first == RNAStructure::INVALID_INDEX || second == RNAStructure::INVALID_INDEX)
			return false;
		assertbiu(first < rnaSeq->size() && second < rnaSeq->size(),
		 "First or second is not a base position in the RNA sequence.");
		return bPair->allowedBasePair((*rnaSeq)[first], (*rnaSeq)[second]);
	}

	RNAStructure& 	
	RNAStructure::operator= (const RNAStructure& rnaStruct2) {
		assertbiu(	seqShared == rnaStruct2.seqShared, 
					"sequence of both has to be shared or not");
		if (this != &rnaStruct2) {
			  // copy seq
			if (seqShared) { // copy only pointer
				rnaSeq = rnaStruct2.rnaSeq;
			} else { // copy
				rnaSeq->resize(rnaStruct2.rnaSeq->size());
				std::copy(	rnaStruct2.rnaSeq->begin(),
							rnaStruct2.rnaSeq->end(),
							rnaSeq->begin());
			}
			  // copy structure
			rnaStructBracketDot->resize(rnaStruct2.rnaStructBracketDot->size());
			std::copy(	rnaStruct2.rnaStructBracketDot->begin(),
						rnaStruct2.rnaStructBracketDot->end(),
						rnaStructBracketDot->begin());
			  // copy remaining data
			bPair = rnaStruct2.bPair;
			rnaBonds = rnaStruct2.rnaBonds;
		}
		return *this;
	}

	bool	 	
	RNAStructure::operator== (const RNAStructure& rnaStruct2)
			 const {
		return (rnaSeq == rnaStruct2.rnaSeq || *rnaSeq == *(rnaStruct2.rnaSeq))
			&& *rnaStructBracketDot == *(rnaStruct2.rnaStructBracketDot)
			&& rnaBonds == rnaStruct2.rnaBonds
			&& *bPair == *(rnaStruct2.bPair);
	}

	bool	 	
	RNAStructure::operator!= (const RNAStructure& rnaStruct2)
			 const {
		return (rnaSeq != rnaStruct2.rnaSeq && *rnaSeq != *(rnaStruct2.rnaSeq))
			|| *rnaStructBracketDot != *(rnaStruct2.rnaStructBracketDot)
			|| rnaBonds != rnaStruct2.rnaBonds
			|| *bPair != *(rnaStruct2.bPair);
	}

	bool		
	RNAStructure::hasValidStructure() const {
		int count = 0;	// counter for open bonds
		// each opening bracket has to have a closing bracket
		for (size_t i=0; i < rnaStructBracketDot->size(); i++) {
			if ((*rnaStructBracketDot)[i] == STRUCT_BND_OPEN ) {
				count++;
			}
			else if ((*rnaStructBracketDot)[i] == STRUCT_BND_CLOSE) {
				count--;
			}
			if (count < 0) {
				return false;
			}
		}
		return (count == 0);
	}

	bool		
	RNAStructure::hasValidBasePairs() const{
		for (size_t i=0; i < rnaBonds.size(); i++) {
			// existing bond has to be an allowed base pair
			if (rnaBonds[i] != RNAStructure::INVALID_INDEX
			    && !isAllowedBasePair(i, rnaBonds[i])) {
				return false;
			}
		}
		return true;
	}	

	bool		
	RNAStructure::hasMinLoopLength() const {
		for (size_t i=0; i < rnaBonds.size(); i++) {
			// existing bond has to have at least MIN_LENGTH_LOOP
			if (rnaBonds[i] != RNAStructure::INVALID_INDEX
			    && (rnaBonds[i]-i) <= RNAStructure::MIN_LOOP_LENGTH) {
				return false;
			}
		}
		return true;
	}


	bool		
	RNAStructure::isValid() const{
		assertbiu(rnaSeq->size() == rnaStructBracketDot->size(),
                 "RNA sequence and structure have to have the same length.");
		return (	hasValidStructure() && hasMinLoopLength()
					&& hasValidBasePairs() );
	}

	std::string	
	RNAStructure::getStringRepresentation() const {
		return ( getStructureString() + "(" + getSequenceString() + ")" );
	}
	
	void
	RNAStructure::initBonds() {
		  // resize and initialise bonds vector
		rnaBonds = std::vector<size_t>
			   (rnaStructBracketDot->size(), RNAStructure::INVALID_INDEX);

		  // do initialization only if structure is valid 
		if (!hasValidStructure()) {
			return;
		}
		  // std::stack openingBonds contains opening bonds to create
		  // vector rnaBonds during intitialisation
		std::stack<size_t>	openingBonds;
		
		  // do update for real initialisation
		for (size_t i=0; i <  rnaStructBracketDot->size(); i++) {
			if ( (*rnaStructBracketDot)[i] == STRUCT_BND_OPEN ) {
				openingBonds.push(i);
			}
			else if ( (*rnaStructBracketDot)[i] == STRUCT_BND_CLOSE ) {
				assertbiu(!openingBonds.empty(),
				 "Closing without opening bond.");
				// check whether pairing of bases allowed	
				rnaBonds[openingBonds.top()] = i;
				openingBonds.pop();
			}
		}
		assertbiu(openingBonds.empty(), "Opening without closing bond.");
	}

} // namespace biu
