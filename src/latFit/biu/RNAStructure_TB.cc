// $Id: RNAStructure_TB.cc,v 1.2 2016/08/08 12:41:57 mmann Exp $


#include "biu/RNAStructure_TB.hh"

#include <biu/assertbiu.hh>

#include <stack>
#include <cmath>
#include <iostream>
#include <limits>

#include "biu/Matrix.hh"


#include <inttypes.h>

#include "biu/RandomNumberFactory.hh"


namespace biu
{



	const size_t RNAStructure_TB::MIN_LOOP_LENGTH = 3;
	const Alphabet RNAStructure_TB::STRUCT_ALPH("().", 1);
	const Alphabet::AlphElem	RNAStructure_TB::STRUCT_BND_OPEN = RNAStructure_TB::STRUCT_ALPH.getElement("(");
	const Alphabet::AlphElem	RNAStructure_TB::STRUCT_BND_CLOSE = RNAStructure_TB::STRUCT_ALPH.getElement(")");
	const Alphabet::AlphElem	RNAStructure_TB::STRUCT_UNBOUND = RNAStructure_TB::STRUCT_ALPH.getElement(".");
	const size_t RNAStructure_TB::INVALID_INDEX = std::numeric_limits<size_t>::max();

	const size_t RNAStructure_TB::BASEPAIRS_VALIDITY = 1;
	const size_t RNAStructure_TB::BASEPAIRS_VALIDITY_CALCULATED = 2;
	const size_t RNAStructure_TB::LOOPSIZE_VALIDITY = 4;
	const size_t RNAStructure_TB::LOOPSIZE_VALIDITY_CALCULATED = 8;
	const size_t RNAStructure_TB::STRUCTURE_VALIDITY = 16;
	const size_t RNAStructure_TB::STRUCTURE_VALIDITY_CALCULATED = 32;


	//	Do not change the following values without checking
	//size_t RNAStructure_TB::moveType(size_t i, size_t j) const;
	const size_t RNAStructure_TB::INSERT_MOVE=0;
	const size_t RNAStructure_TB::SHIFT_MOVE=2;
	const size_t RNAStructure_TB::REVERSE_SHIFT_MOVE=4;
	const size_t RNAStructure_TB::INVALID_DELETE_MOVE=6;
	const size_t RNAStructure_TB::DELETE_MOVE=7;





	RNAStructure_TB::RNAStructure_TB(	const std::string& rnaSeqStr,
	 							const std::string& rnaStructBracketDotStr,
	 							const AllowedBasePairs* const _bPair)
	 :
	 	seqShared(false),
		bPair(_bPair),
		strucStatus(0),
		rnaSeq(NULL),
	 	bracketStructStr(rnaStructBracketDotStr),
	 	validTreeStruc(NULL)
		{

		assertbiu(bPair != NULL, "bPair is not allowed to be NULL.");
		assertbiu(bPair->getAlphabet()->isAlphabetString(rnaSeqStr),
		 "RNA sequence contains characters which are not in the alphabet.");
		assertbiu(STRUCT_ALPH.isAlphabetString(rnaStructBracketDotStr),
		 "RNA structure contains characters which are not in the alphabet.");


		rnaSeq = new Sequence(bPair->getAlphabet()->getSequence(rnaSeqStr));


		assertbiu(rnaSeq->size() == rnaStructBracketDotStr.size(),
		 "RNA sequence and structure have to have the same length.");


		Structure* tempStruct = new Structure(
									STRUCT_ALPH.getSequence(rnaStructBracketDotStr));


		assertbiu(bPair->getAlphabet()->isAlphabetSequence(*rnaSeq),
		 "RNA sequence contains characters which are not in the alphabet.");
		assertbiu(STRUCT_ALPH.isAlphabetSequence(*tempStruct),
		 "RNA structure contains characters which are not in the alphabet.");


		setStructure(*tempStruct);

		delete tempStruct;
		tempStruct = NULL;

	}

	RNAStructure_TB::RNAStructure_TB(	Sequence* rnaSeq_,
	 							const Structure* const rnaStructBracketDot_,
	 							const AllowedBasePairs* const bPair_,
	 							const bool seqIsShared)
	 :
	 	seqShared(seqIsShared),
		bPair(bPair_),
		strucStatus(0),
		rnaSeq(NULL),
	 	bracketStructStr(STRUCT_ALPH.getString(*rnaStructBracketDot_)),
	 	validTreeStruc(NULL)
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
		};

		setStructure(*rnaStructBracketDot_);
	}

	RNAStructure_TB::RNAStructure_TB(const RNAStructure_TB& rnaStruct)
	 :
		seqShared(rnaStruct.seqShared),
 		bPair(rnaStruct.bPair),
 		strucStatus(rnaStruct.strucStatus),
		rnaSeq(rnaStruct.rnaSeq),
		bracketStructStr(rnaStruct.bracketStructStr),
		validTreeStruc(NULL)
	{
		assertbiu(bPair != NULL, "bPair is not allowed to be NULL.");
		assertbiu(bracketStructStr.size() != 0, "bracketStructStr should not be an empty string");
		assertbiu(rnaSeq->size() == bracketStructStr.size(),
		 "RNA sequence and structure have to have the same length.");
		assertbiu(rnaSeq != NULL, "no RNA sequence given.");
		assertbiu(bPair->getAlphabet()->isAlphabetSequence(*rnaSeq),
		 "RNA sequence contains characters which are not in the alphabet.");

		if (rnaStruct.validTreeStruc != NULL) {
			validTreeStruc = new std::vector<TreeItem>(*(rnaStruct.validTreeStruc));
		}

		if (!seqShared ) {
			rnaSeq = new Sequence(*(rnaStruct.rnaSeq));
		}


	}

	RNAStructure_TB::~RNAStructure_TB() {
		if(!seqShared && rnaSeq != NULL) {
			delete rnaSeq;
			rnaSeq = NULL;
		}
		if(validTreeStruc != NULL){
			delete validTreeStruc;
			validTreeStruc = NULL;
		}

	}



	RNAStructure_TB& RNAStructure_TB::operator= (const RNAStructure_TB& rnaStruct2) {
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
			  // copy remaining data
			bPair = rnaStruct2.bPair;
			strucStatus = rnaStruct2.strucStatus;

			if (rnaStruct2.validTreeStruc != NULL){
				if (validTreeStruc == NULL)
					validTreeStruc = new std::vector<TreeItem>(rnaStruct2.validTreeStruc->size());
				else
					validTreeStruc->resize(rnaStruct2.validTreeStruc->size());

				std::copy(	rnaStruct2.validTreeStruc->begin(),
								rnaStruct2.validTreeStruc->end(),
								validTreeStruc->begin());

			} else if (validTreeStruc != NULL) {
				delete validTreeStruc;
				validTreeStruc = NULL;
			}

			bracketStructStr = rnaStruct2.bracketStructStr;
		}
		return *this;
	}




	bool RNAStructure_TB::operator== (const RNAStructure_TB& rnaStruct2)
			 const {
		assertbiu(bPair != NULL, "bPair should not be NULL");
		assertbiu(rnaSeq != NULL, "rnaSeq should not be NULL");
		assertbiu(bracketStructStr.size() != 0, "bracketStructStr should not be an empty string");
		bool res =
			(
			(rnaSeq == rnaStruct2.rnaSeq )
			&& *bPair == *(rnaStruct2.bPair)
			&& bracketStructStr == rnaStruct2.bracketStructStr
			)||(
			!(rnaSeq == rnaStruct2.rnaSeq )
			&&*rnaSeq == *(rnaStruct2.rnaSeq)
			&& *bPair == *(rnaStruct2.bPair)
			&& bracketStructStr == rnaStruct2.bracketStructStr
			);

//This is "the best" optimized form of:
//		res= (a) rnaSeq == rnaStruct2.rnaSeq || (b) *rnaSeq == *(rnaStruct2.rnaSeq)
//			(c)	&& *bPair == *(rnaStruct2.bPair)
//				&& bracketStructStr == rnaStruct2.bracketStructStr
//				;
//		(a|b)&c <==>(a&c)|(!a&b&c) a is cheap, b & c are expensive


		return res && (
				validTreeStruc==rnaStruct2.validTreeStruc 		//If equal, they would be NULL
				||	*validTreeStruc==*(rnaStruct2.validTreeStruc)
				);
	}


//	bool RNAStructure_TB::operator!= (const RNAStructure_TB& rnaStruct2) const {
//		return !operator==(rnaStruct2);
//	}


	void RNAStructure_TB::setStructure(const Structure& str) {
		assertbiu(rnaSeq->size() == str.size(), "structure has wrong length");

		strucStatus = 0;// i.e. All status-valvulated bits are 0
		if (validTreeStruc==NULL) {
			validTreeStruc = new std::vector<TreeItem>(str.size());
		} else {
			validTreeStruc->resize(str.size());
		}
		if (! initTree(str)) {
			delete validTreeStruc;
			validTreeStruc = NULL;
			bracketStructStr.resize(0);
			strucStatus = 0;
		} else {
			bracketStructStr = STRUCT_ALPH.getString(str);
		}
		return;

	}


	bool RNAStructure_TB::hasValidStructure() const {
		if (strucStatus & STRUCTURE_VALIDITY_CALCULATED) // If calculated, no need to recalculate
			return (strucStatus & STRUCTURE_VALIDITY);

		//Assume validity, then change if not
		strucStatus |= 	STRUCTURE_VALIDITY_CALCULATED | STRUCTURE_VALIDITY;

		size_t bondOpenings = 0;

		const char open = STRUCT_ALPH.getString(STRUCT_BND_OPEN).at(0);
		const char close = STRUCT_ALPH.getString(STRUCT_BND_CLOSE).at(0);

		for (size_t j=0; j < bracketStructStr.size(); j++) {

			if ( bracketStructStr[j] == open )  {
				bondOpenings++;
			}
			else if ( bracketStructStr[j] == close ) {
				if ( bondOpenings == 0 ) {
					strucStatus &= 	~STRUCTURE_VALIDITY;
					return false;
				}
				bondOpenings--;
			}
		}


		// check if no open brackets without closing
		if ( bondOpenings == 0) {
			return true;
		} else {
			strucStatus &= 	~STRUCTURE_VALIDITY;
			return false;
		}
	}

	bool RNAStructure_TB::hasValidBasePairs() const {

		if (strucStatus & BASEPAIRS_VALIDITY_CALCULATED) // If calclated, no need to recalculate
			return (strucStatus & BASEPAIRS_VALIDITY);

		//Assume validity, then change if not
		strucStatus |= 	BASEPAIRS_VALIDITY_CALCULATED | BASEPAIRS_VALIDITY;

		std::stack<size_t>	openingBonds;

		const char open = STRUCT_ALPH.getString(STRUCT_BND_OPEN).at(0);
		const char close = STRUCT_ALPH.getString(STRUCT_BND_CLOSE).at(0);

		for (size_t j=0; j < bracketStructStr.size(); j++) {

			if ( bracketStructStr[j] == open ) {
				openingBonds.push(j);
			} else if ( bracketStructStr[j] == close ) {
				if (!openingBonds.empty()) {
					if ( !isAllowedBasePair( j, openingBonds.top()) ) {
						strucStatus &= 	~BASEPAIRS_VALIDITY;
						return false;
					}
					openingBonds.pop();
				} else {	// A chance to detect structure invalidity
					strucStatus |= 	STRUCTURE_VALIDITY_CALCULATED;
					strucStatus &= ~STRUCTURE_VALIDITY;
				}
			}
		}

		//we are able to determine stuctural validity
		if (!(strucStatus & STRUCTURE_VALIDITY_CALCULATED)) {	// if not calculated
			strucStatus |= STRUCTURE_VALIDITY_CALCULATED;
			if (!openingBonds.empty())
				strucStatus &= ~STRUCTURE_VALIDITY;
			else
				strucStatus |= STRUCTURE_VALIDITY;
		}

		return true;
	}



	bool RNAStructure_TB::hasValidLoopSize() const {
		if (strucStatus & LOOPSIZE_VALIDITY_CALCULATED) // If calclated, no need to recalculate
			return (strucStatus & LOOPSIZE_VALIDITY);
		//Assume validity, then change if not
		strucStatus |= 	LOOPSIZE_VALIDITY_CALCULATED | LOOPSIZE_VALIDITY;

		size_t lastOpen = std::string::npos;

		const char open = STRUCT_ALPH.getString(STRUCT_BND_OPEN).at(0);
		const char close = STRUCT_ALPH.getString(STRUCT_BND_CLOSE).at(0);

		for (size_t j=0; j < bracketStructStr.size(); j++) {

			if ( bracketStructStr[j] == open ) {
				lastOpen = j;
				continue;
			} else if (lastOpen != std::string::npos && bracketStructStr[j] == close) {
				if ( (j - lastOpen) <= MIN_LOOP_LENGTH) {
					strucStatus &= ~LOOPSIZE_VALIDITY;
					return false;
				}
				lastOpen = std::string::npos;
			}
		}

		return true;
	}



	bool RNAStructure_TB::isValidInsertMove(const size_t i, const size_t j) const{
		assertbiu(validTreeStruc != NULL, "no valid structure present");
		assertbiu(i<validTreeStruc->size(),"i index should be smaller than the structure size");
		assertbiu(j<validTreeStruc->size(),"j index should be smaller than the structure size");
		assertbiu(i!=j,"i and j should not be equal");

		return (
				isAllowedBasePair(i, j) 						//Match (More probable to fail)
			&& (*validTreeStruc)[i].pair == INVALID_INDEX  		//Free
			&& (*validTreeStruc)[j].pair == INVALID_INDEX		//Free
			&& (size_t) abs(j-i) > RNAStructure_TB::MIN_LOOP_LENGTH 		//Loop size
			&& areOnSameLevel(i,j)
			);													//Nested  (On the same level) (High cost)
	}

	int
	RNAStructure_TB
	::isValidSingleMove( const size_t i, const size_t j) const {
		assertbiu(validTreeStruc != NULL, "no valid structure present");
		assertbiu(i<validTreeStruc->size(),"i index should be smaller than the structure size");
		assertbiu(j<validTreeStruc->size(),"j index should be smaller than the structure size");
		assertbiu(i!=j,"i and j should not be equal");

		  // check for valid deletion
		if (validTreeStruc->at(i).pair == j) {
			return -1;
		}
		  // check for valid insertion
		if (isValidInsertMove(i,j)) {
			return +1;
		}

		  // no valid indel
		return 0;
	}


	bool RNAStructure_TB::isValidLeftShiftMove(const size_t i, const size_t j, const size_t k) const{
		assertbiu(validTreeStruc != NULL, "no valid structure present");
		assertbiu(i<validTreeStruc->size(),"i index should be smaller than the structure size");
		assertbiu(j<validTreeStruc->size(),"j index should be smaller than the structure size");
		assertbiu(k<validTreeStruc->size(),"k index should be smaller than the structure size");
		assertbiu(i!=j,"i and j should not be equal");
		assertbiu(i!=k,"i and k should not be equal");
		assertbiu(k!=j,"k and j should not be equal");

		return(
					isAllowedBasePair(k, j) 						//Match (More probable to fail)
				&& (*validTreeStruc)[j].pair == i					//Paired (More probable to fail than free)
				&& (*validTreeStruc)[k].pair == INVALID_INDEX  		//Free
				&& (size_t) abs(j-k) > RNAStructure_TB::MIN_LOOP_LENGTH 		//Loop size (Less probable to fail)
				&& (areOnSameLevel(i,k) || areOnSameLevel(k,j)) 	//k on one of the levels (High cost)
				);
	}

	bool RNAStructure_TB::isValidRightShiftMove(const size_t i, const size_t j, const size_t k) const{
		assertbiu(i<validTreeStruc->size(),"i index should be smaller than the structure size");
		assertbiu(j<validTreeStruc->size(),"j index should be smaller than the structure size");
		assertbiu(k<validTreeStruc->size(),"k index should be smaller than the structure size");
		assertbiu(i!=j,"i and j should not be equal");
		assertbiu(i!=k,"i and k should not be equal");
		assertbiu(k!=j,"k and j should not be equal");

		return(
				 	isAllowedBasePair(k, i) 						//Match (More probable to fail)
				&& (*validTreeStruc)[i].pair == j					//Paired (More probable to fail than free)
				&& (*validTreeStruc)[k].pair == INVALID_INDEX  		//Free
				&& (size_t) abs(i-k) > RNAStructure_TB::MIN_LOOP_LENGTH 		//Loop size (Less probable to fail)
				&& (areOnSameLevel(i,k) || areOnSameLevel(k,j)) 	//k on one of the levels (High cost)
				);
	}



	void RNAStructure_TB::insertBond(const size_t i, const size_t j){

		assertbiu(isValidInsertMove(i,j),"i,j must be a valid insert move to be inserted!");

		setStringBrackets(i,j);			// Adjust string representation

		//Pair them
		(*validTreeStruc)[i].pair = j;
		(*validTreeStruc)[j].pair = i;

		//Push the structure between them down
		int in = (*validTreeStruc)[i].next;
		(*validTreeStruc)[i].next = (*validTreeStruc)[j].next;
		(*validTreeStruc)[j].next = in;
	}
	void RNAStructure_TB::deleteBond(const size_t i, const size_t j){
		assertbiu(isValidDeleteMove(i,j),"i,j must be a valid delete move to be deleted!");

		// Adjust string representation
		bracketStructStr[i]=bracketStructStr[j]='.';

		// Unpair them
		(*validTreeStruc)[i].pair = (*validTreeStruc)[j].pair = RNAStructure_TB::INVALID_INDEX;

		// Pull up the structure between them to the same level back
		int in = (*validTreeStruc)[i].next;
		(*validTreeStruc)[i].next = (*validTreeStruc)[j].next;
		(*validTreeStruc)[j].next = in;


	}




	void RNAStructure_TB::leftShift(const size_t i, const size_t j, const size_t k){

		assertbiu(isValidLeftShiftMove(i,j,k),"i,j,k must be a valid leftshift move to be shifted!");
		// Assertionsss

		shiftToBond(j,i,k);

	}

	void RNAStructure_TB::rightShift(const size_t i, const size_t j, const size_t k){

		assertbiu(isValidRightShiftMove(i,j,k),"i,j,k must be a valid rightshift move to be shifted!");

		shiftToBond(i,j,k);
	}



	bool RNAStructure_TB::initTree(const Structure& str) {
		assertbiu(str.size()>0,"str shouldn't be an empty structure");

		  // std::stack openingBonds contains opening bonds to create
		  // vector rnaBonds during intitialisation
		std::stack<size_t>	openingBonds;
		size_t i,j;
		  // do update for real initialisation

		//-1 the step will be done seperatly after that for efficiency
		for (j=0; j < str.size()-1; j++) {
			(*validTreeStruc)[j].next = j+1;	//Building the default unfolded next structure

			if ( str[j] == STRUCT_BND_OPEN )
				openingBonds.push(j);
			else if ( str[j] == STRUCT_BND_CLOSE ) {
				if (	openingBonds.empty()
					||	j-(i = openingBonds.top()) <= RNAStructure_TB::MIN_LOOP_LENGTH
					||	!isAllowedBasePair(i, j)
						) return false;			// Break but still uncertain why
				openingBonds.pop();

				//Pair them
				(*validTreeStruc)[j].pair = i;
				(*validTreeStruc)[i].pair = j;
				//Level down
				size_t in = (*validTreeStruc)[i].next;
				(*validTreeStruc)[i].next = (*validTreeStruc)[j].next;
				(*validTreeStruc)[j].next = in;
			}else (*validTreeStruc)[j].pair = INVALID_INDEX;
		}

		(*validTreeStruc)[j].next = 0; // Block was repeated to avoid an extra if
		(*validTreeStruc)[j].pair = INVALID_INDEX;
		if ( str[j] == STRUCT_BND_OPEN ) // Loop block specialized for last step
			return false;
		else if ( str[j] == STRUCT_BND_CLOSE ) {
			if (	openingBonds.empty()
				||	j-(i = openingBonds.top()) <= RNAStructure_TB::MIN_LOOP_LENGTH
				||	!isAllowedBasePair(i, j)
					) return false;
			openingBonds.pop();


			(*validTreeStruc)[j].pair = i;
			(*validTreeStruc)[i].pair = j;
			size_t in = (*validTreeStruc)[i].next;
			(*validTreeStruc)[i].next = (*validTreeStruc)[j].next;
			(*validTreeStruc)[j].next = in;
		}else (*validTreeStruc)[j].pair = INVALID_INDEX;

		if (openingBonds.empty()) {
			strucStatus =		// All calculated & valid
								BASEPAIRS_VALIDITY
							| 	BASEPAIRS_VALIDITY_CALCULATED
							|	LOOPSIZE_VALIDITY
							|	LOOPSIZE_VALIDITY_CALCULATED
							|	STRUCTURE_VALIDITY
							|	STRUCTURE_VALIDITY_CALCULATED;
			return true;
		}else {					// All calculated, all valid except for structure
			strucStatus =		BASEPAIRS_VALIDITY
							| 	BASEPAIRS_VALIDITY_CALCULATED
							|	LOOPSIZE_VALIDITY
							|	LOOPSIZE_VALIDITY_CALCULATED
							|	!STRUCTURE_VALIDITY
							|	STRUCTURE_VALIDITY_CALCULATED;
			return false;
		}

	}



	std::string RNAStructure_TB::decodedStatus() const{
		std::string x =   (strucStatus&BASEPAIRS_VALIDITY)!=0?"BP| ":"";
		x+= (strucStatus&BASEPAIRS_VALIDITY_CALCULATED)!=0?"BPC ":"";
		x+="\n";
		x+= (strucStatus&LOOPSIZE_VALIDITY)!=0?"L| ":"";
		x+= (strucStatus&LOOPSIZE_VALIDITY_CALCULATED)!=0?"LC ":"";
		x+="\n";
		x+= (strucStatus&STRUCTURE_VALIDITY)!=0?"S| ":"";
		x+= (strucStatus&STRUCTURE_VALIDITY_CALCULATED)!=0?"SC ":"";

		return x;

	}



	bool RNAStructure_TB::isValidUnorderedShift(const size_t i, const size_t j) const{
		assertbiu(i!=j,"i and j should not be equal");
		assertbiu(	(*validTreeStruc)[i].pair== INVALID_INDEX ||
					(*validTreeStruc)[j].pair== INVALID_INDEX ||
					(*validTreeStruc)[i].pair!=(*validTreeStruc)[j].pair
					,"i and j chouldn't be paired to the same base, invalid structure");

		if (!isAllowedBasePair(i,j)) return false;
		if ((*validTreeStruc)[i].pair!=INVALID_INDEX){ //i is the paired one
			return
					(*validTreeStruc)[j].pair==INVALID_INDEX
				&& 	(size_t) abs(i-j) > RNAStructure_TB::MIN_LOOP_LENGTH 		//Loop size (Less probable to fail)
				&& 	(areOnSameLevel(i,j) || areOnSameLevel(j,(*validTreeStruc)[i].pair));

		}else if ((*validTreeStruc)[j].pair!=INVALID_INDEX){//j is the paired one
			return
					(*validTreeStruc)[i].pair==INVALID_INDEX
				&& 	(size_t) abs(i-j) > RNAStructure_TB::MIN_LOOP_LENGTH 		//Loop size (Less probable to fail)
				&& 	(areOnSameLevel(i,j) || areOnSameLevel(i,(*validTreeStruc)[j].pair));

		}else return false;//both unpaired
	}


	bool RNAStructure_TB::isValidShift(const size_t i, const size_t j) const{
		assertbiu(i!=j,"i and j should not be equal");
		assertbiu(	(*validTreeStruc)[i].pair== INVALID_INDEX ||
					(*validTreeStruc)[j].pair== INVALID_INDEX ||
					(*validTreeStruc)[i].pair!=(*validTreeStruc)[j].pair
					,"i and j chouldn't be paired to the same base, Invalid structure!!");


			return
				isAllowedBasePair(i, j) 						//Match (More probable to fail)
				&& (*validTreeStruc)[j].pair!=INVALID_INDEX
				&& (*validTreeStruc)[i].pair==INVALID_INDEX
				&& (size_t) abs(i-j) > RNAStructure_TB::MIN_LOOP_LENGTH 		//Loop size (Less probable to fail)
				&& (areOnSameLevel(i,j) || areOnSameLevel(j,(*validTreeStruc)[i].pair));

	}



	void RNAStructure_TB::shiftUnordered(const size_t i, const size_t j){
		assertbiu(isValidShift(i,j),"i,j must be a valid shift move to be shifted!");

		if ((*validTreeStruc)[i].pair!=INVALID_INDEX)//i is paired
			shiftToBond(i,j);
		else// j is paired
			shiftToBond(j,i);
	}

	void RNAStructure_TB::shiftToBond(const size_t i, const size_t k){
		assertbiu( (*validTreeStruc)[i].pair!=k
				,"Warning: rearrangeWith(i,k) should be called before changing the .pair of i into k");
		// Assert every thing

		size_t j = (*validTreeStruc)[i].pair;
		size_t jn = (*validTreeStruc)[j].next;
		(*validTreeStruc)[j].next = (*validTreeStruc)[i].next;
		(*validTreeStruc)[i].next = (*validTreeStruc)[k].next;
		(*validTreeStruc)[k].next = jn;


		(*validTreeStruc)[i].pair = k;
		(*validTreeStruc)[k].pair = i;
		(*validTreeStruc)[j].pair =  RNAStructure_TB::INVALID_INDEX;

		setStringBrackets(i,k);
		bracketStructStr[j]='.';

	}

	void RNAStructure_TB::shiftToBond(const size_t i, const size_t j, const size_t k){
		// Assert i,j are pairs
		// Assert every thing


		size_t jn = (*validTreeStruc)[j].next;
		(*validTreeStruc)[j].next = (*validTreeStruc)[i].next;
		(*validTreeStruc)[i].next = (*validTreeStruc)[k].next;
		(*validTreeStruc)[k].next = jn;


		(*validTreeStruc)[i].pair = k;
		(*validTreeStruc)[k].pair = i;
		(*validTreeStruc)[j].pair =  RNAStructure_TB::INVALID_INDEX;

		setStringBrackets(i,k);
		bracketStructStr[j]='.';
	}



	size_t RNAStructure_TB::moveType(size_t i, size_t j) const {
		//000	<=>	both unpaired				<=>	insert;	//Order not important
		//010	<=> i is paired					<=>	shiftToBond	(i,j)
		//100	<=>	j is paired					<=>	shiftToBond	(j,i)
		//110	<=>	both paired					<=>	invalid delete(i,j)
		//111	<=>	both paired to each other	<=>	valid delete(i,j)	//Order not important
		size_t res = 0;

		  // direct return if paired with each other == valid delete
		if ((*validTreeStruc)[i].pair==j) return 7;

		  // check which is bound
		if ((*validTreeStruc)[i].pair!=INVALID_INDEX) res|=2;
		if ((*validTreeStruc)[j].pair!=INVALID_INDEX) res|=4;

		return res;

	}

	void RNAStructure_TB::executeOrderedMove(size_t i, size_t j){

		if ((*validTreeStruc)[i].pair==INVALID_INDEX){
			if ((*validTreeStruc)[j].pair==INVALID_INDEX){
				insertBond(i,j);
			}
			else{
				shiftToBond(j,i);
			}
		}
		else if ((*validTreeStruc)[j].pair==INVALID_INDEX){
			shiftToBond(i,j);
		}else{
			if ((*validTreeStruc)[j].pair==i) deleteBond(i,j);
			assertbiu((*validTreeStruc)[j].pair==i,"executing an move that is decoded into an invalid delete move!!");
		}
	}

	bool RNAStructure_TB::executeOrderedMoveTry(size_t i, size_t j){
		if ((*validTreeStruc)[i].pair==INVALID_INDEX){
			if ((*validTreeStruc)[j].pair==INVALID_INDEX){
				if (isValidInsertMove(i,j)){
					insertBond(i,j);
					return true;
				}
			}
			else{
				if (isValidShift(j,i)){
					shiftToBond(j,i);
					return true;
				}
			}
		}
		else if ((*validTreeStruc)[j].pair==INVALID_INDEX){
			if (isValidShift(i,j)){
				shiftToBond(i,j);
				return true;
			}
		}else{
			if (isValidDeleteMove(i,j)){
				deleteBond(i,j);
				return true;
			}
		}

		return false;
	}


	bool
	RNAStructure_TB
	::getNextSingleMove( size_t& i, size_t& j ) const
	{
		assertbiu(validTreeStruc != NULL, "no valid structure present");
		assertbiu(i<=j, "index i is greater than j");
		assertbiu(&i!=&j, "referenced variable i and j are the same !!!");
		  // check if first single move is to find
		if ( i == INVALID_INDEX ) {
			i = 0;
			j = 0;
		} else {
			  // last move was a base pair
			if (validTreeStruc->at(i).pair == j ) {
				i++;
				j = size_t(i);
			}
			  // was no base pair deletion .. maybe an insertion
			else {
				j = size_t(validTreeStruc->at(j).next);
				  // check if j out of bound
				if (j <= i || j >= validTreeStruc->size()) {
					i++;
					j = size_t(i);
				}
			}
		}

		size_t minJpos = i+MIN_LOOP_LENGTH+1;

		while ( minJpos < validTreeStruc->size() ) {

			  // check if j was not set yet
			if (i==j) {
				if (validTreeStruc->at(i).pair != INVALID_INDEX
						&& validTreeStruc->at(i).pair > i )
				{
					j = size_t(validTreeStruc->at(i).pair);
					return true;
				}
				j = size_t(validTreeStruc->at(i).next);
				  // set j to minimal distance allowed for a base pair
				while (i<j && j < validTreeStruc->size() && j < minJpos) {
					j = size_t(validTreeStruc->at(j).next);
				}
			}
			  // otherwise assure that j has min-loop-distance from i
			else if (j < minJpos) {
				j = size_t(validTreeStruc->at(i).next);
				  // set j to minimal distance allowed for a base pair
				while (i<j && j < validTreeStruc->size() && j < minJpos) {
					j = size_t(validTreeStruc->at(j).next);
				}
			}

			  // check if j is a valid index
			if ( i < j && j < validTreeStruc->size()) {
				  // check if i and j are paired or can form a base pair
				if (validTreeStruc->at(i).pair == j
					|| (validTreeStruc->at(j).pair
							== INVALID_INDEX && isAllowedBasePair(i,j)))
				{
					return true;
				}

				  // update j to next possible pairing partner of i
				j = size_t(validTreeStruc->at(j).next);
				if (i>=j || j > validTreeStruc->size()) {
					i++;
					j = size_t(i);
				}
			} else {
				i++;
				j = size_t(i);
			}
			  // minimal position of j to allow for a base pair
			minJpos = i+MIN_LOOP_LENGTH+1;
		}
		return false;
	}


	void RNAStructure_TB::decodeMoveIndex(size_t const index,  size_t&  i,  size_t&  j) const
	{

		assertbiu(index<getMoveIndexCount(),"Index should be in [0 .. n.(n-1)/2[ ");

		size_t n =	rnaSeq->size();

		if ((n & 1)==0){	//even

			j = index % (n-1);
			i = index / (n-1);
			if (i>j) {
				i=n-i-1;
				j=n-j-1;

			} else {
				//	i = i;
				j=j+1;
			}

		} else {			//odd
			j = index % n;
			i = index / n;
			if (i>=j) {
				i=n-i-2;
				j=n-j-1;

			} else {
				//Same! i and j
			}
		}
	}

					//The previous code is a compact form of the following sample:
//					----------------------------------------------------------
//					int i,j,n;
//					n=4;
//						for (i=0; i < n/2; ++i){
//					for (j = 0; j < n-1; ++j) {
//							if (i>j){
//								std::cout <<n-i-1<<"-"<<n-j-1<<"\n";
//							}else{
//								std::cout <<i<<" "<<j+1<<"\n";
//							}
//						}
//					}
//					std::cout << "------\n";
//					n=5;
//					//
//						for (i=0; i < (n-1)/2; ++i){
//					for (j = 0; j < n; ++j) {
//							if (i>=j){
//								std::cout <<n-i-2<<"-"<<n-j-1<<"\n";
//							}else{
//								std::cout <<i<<" "<<j<<"\n";
//							}
//						}
//					}

//					Output:
//					0 1
//					2-3
//					0 2
//					1 2
//					0 3
//					1 3
//					-----
//					3-4
//					2-4
//					0 1
//					2-3
//					0 2
//					1 2
//					0 3
//					1 3
//					0 4
//					1 4





	bool RNAStructure_TB::executeOrderedMoveTry(size_t const index){
		size_t i,j;
		decodeMoveIndex(index,i,j);
		return  executeOrderedMoveTry(i,j);
	}

	inline size_t RNAStructure_TB::getMoveIndexCount() const{
		return rnaSeq->size()*(rnaSeq->size()-1)/2;
	}


} // namespace biu

