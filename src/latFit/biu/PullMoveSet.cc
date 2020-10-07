// $Id: PullMoveSet.cc,v 1.2 2016/08/08 12:41:59 mmann Exp $

#include "biu/assertbiu.hh"
#include "biu/PullMoveSet.hh"
#include "biu/NeighborVector.hh"
#include "biu/LatticeProtein_Ipnt.hh"

#include <iostream>
namespace biu
{
	PullMoveSet::PullMoveDecoder::PullMoveDecoder(const LatticeModel* const lattice_)
	 :	lattice(lattice_)
	 	, triangularLattice(lattice->getDescriptor()->isPossibleRing(3))
	 	, neighDirections()
	 	, neighDirectionsAccess()
	 	, stdPulls()
	 	, endPulls()
	 	, directionNumber(0)
	 	, stdPullNumber(0)
	 	, endPullNumber(0)
	{
		assertbiu(lattice != NULL, "no lattice model given");
		const LatticeNeighborhood& neighbors = lattice->getNeighborhood();
// std::cerr <<"\n lattice = " <<lattice->getDescriptor()->getName() <<" : \n";

		  // ensure that a pull of 2 monomers is possible if lattice has no triangles
		#ifndef NDEBUG
			if (!isTriangular()) {
				assertbiu(lattice->getDescriptor()->isPossibleRing(4), "lattice does not support selfavoiding rings of size 4");
			}
		#endif

		  // initialize neighboring vectors with first vector
		neighDirections.resize( neighbors.size(),*(neighbors.begin()) );
		directionNumber = neighDirections.size();
		  // copy reference of all neighbor vectors for direct lookup
		size_t i=0;
		for ( LatticeNeighborhood::const_iterator it = neighbors.begin();
				it != neighbors.end(); it++)
		{
			  // get reference for this neighboring vector
			neighDirections[i] = (*it);
			neighDirectionsAccess[*it] = i;
// std::cerr <<" neighDirections["<<i<<"] = " <<neighDirections[i]<<"\n";
			i++;
		}
// std::cerr <<" directionNumber = " <<directionNumber<<"\n";
		// set up stdPulls
		// structure: direction i->(x) : pull number : (c, l)

		  // resize container for filling
		stdPulls.resize(neighDirections.size());
		// for every direction between two neighbors
		for (size_t dir = 0; dir < neighDirections.size(); dir++)
		{
// std::cerr <<" stdPulls["<<dir<<"] = " ;
			std::vector< IPointPair >& stdMoves = stdPulls[dir];
			// for each neighbor of i
			for (size_t newC = 0; newC != neighDirections.size(); newC++)
			{
				  // for each possible C position
				if ( newC != dir) {
					  // calculate triangular or cubic pulls
					if(isTriangular()) {
						  // if i and (i+1)=~newL are neighbored save c
						if ( lattice->areNeighbored( neighDirections.at(dir), neighDirections.at(newC) )) {
							IPointPair clpair;
							clpair.first = neighDirections.at(newC);
							stdMoves.push_back( clpair );
						}
					} else {
					  // generate possible L positions
						for (size_t newLdir = 0; newLdir < neighDirections.size(); newLdir++)
						{
							  // calculate new L position
							IntPoint newL = neighDirections.at(dir) + neighDirections.at(newLdir);
								  // if neighbor of i and neighbor of (x) are neighbored
								  // save new c,l pair
							if ( newL != IntPoint(0,0,0)
								&& lattice->areNeighbored( neighDirections.at(newC), newL ))
							{
								IPointPair clpair;
								clpair.first = neighDirections.at(newC);
								clpair.second = newL;
								stdMoves.push_back( clpair );
	// std::cerr <<clpair.first <<" " <<clpair.second <<" # " ;
							}
						}
					}
				}
			}
// std::cerr <<"\n" ;
		}
		stdPullNumber = stdPulls.at(0).size();
// std::cerr <<" stdPullNumber = " <<stdPullNumber<<"\n" ;

#ifndef NDEBUG
	for (size_t i=0; i<stdPulls.size(); i++) {
		assertbiu( stdPulls.at(i).size() == stdPullNumber, "there is no equal number of pull moves for all neighboring vectors");
	}
#endif

		endPullNumber = 0;
// std::cerr <<" endPulls = " ;
		// set up endPulls
		for (LatticeNeighborhood::const_iterator nb1 = neighbors.begin();
			nb1 != neighbors.end(); nb1++)
		{
			if (isTriangular()) {
				IPointPair clpair;
				clpair.first = (*nb1);
				endPulls.push_back(clpair);
			} else {
				for (LatticeNeighborhood::const_iterator nb2 = neighbors.begin();
					nb2 != neighbors.end(); nb2++)
				{
					if ((*nb1 + *nb2) != IntPoint(0,0,0))
					{
						IPointPair clpair;
						clpair.first = (*nb1);
						clpair.second = (*nb1 + *nb2);
						endPulls.push_back(clpair);
	// std::cerr <<clpair.first <<" " <<clpair.second <<" # " ;
					}
				}
			}
		}
// std::cerr <<"\n" ;
		endPullNumber = endPulls.size();
// std::cerr <<" endPullNumber = "<<endPullNumber<<"\n" ;
//		endPullNumber = (directionNumber-1)*(directionNumber-1);
	}

	void
	PullMoveSet::PullMoveDecoder::lookupMove(const size_t& todoSize,
			size_t& moveIndex, bool& pullFront,
			bool& stdPull, size_t& movePosition) const
	{
// std::cerr <<" lookupMove call ( "<<todoSize <<", " <<moveIndex <<", " <<pullFront <<", " <<stdPull <<", " <<movePosition <<")\n";
		// *********************************
		// code table:
		// *********************************
		// legal moves:
		// 		[0, endPullNumber*2 + (proteinLength-1)*stdPullNumber*2 - 1]
		// end pulls @ position 0:
		//		[endPullNumber, 2*endPullNumber - 1]
		// end pulls @ last position:
		//		[0, endPullNumber - 1]
		// standart pulls:
		//		[2*endPullNumber,
		//			endPullNumber*2 + (proteinLength-1)*stdPullNumber*2 - 1]

		if (moveIndex < endPullNumber)
		{
			// end pull at end position
			pullFront=true;
			stdPull=false;
			movePosition=todoSize-1;
		}
		else if (moveIndex < endPullNumber*2)
		{
			// end pull at position 0
			pullFront=false;
			stdPull=false;
			movePosition=0;
			moveIndex-=endPullNumber;
		}
		else
		{
			// do std pulls
			stdPull = true;
			moveIndex -= endPullNumber*2;
			if (moveIndex < stdPullNumber)
			{
				// pull first element & pull front
				pullFront = true;
				movePosition = 0;
			}
			else if (moveIndex < stdPullNumber*2)
			{
				// pull last element & pull back
				moveIndex -= stdPullNumber;
				pullFront = false;
				movePosition = todoSize-1;
			}
			else
			{
				// do standard pull @ element in [1, todoSize-1]
				moveIndex -= stdPullNumber*2;
				movePosition = 1 + moveIndex/(2*stdPullNumber);
				pullFront = (moveIndex % (2*stdPullNumber)) < stdPullNumber ?
						true : false;
				moveIndex = moveIndex % stdPullNumber;
			}
		}
// std::cerr <<" lookupMove done ( "<<todoSize <<", " <<moveIndex <<", " <<pullFront <<", " <<stdPull <<", " <<movePosition <<")\n";
	}

	PullMoveSet::IPointPair
	PullMoveSet::PullMoveDecoder::lookupEndMove(const IntPoint& rootpoint,
			const size_t& moveIndex) const
	{
// std::cerr <<" lookupEndMove done ( "<<rootpoint <<", " <<moveIndex <<") ... ";
		  // copy pull coordinates
		IPointPair movePositions = endPulls.at(moveIndex);
		  // shift to root position
		movePositions.first  += rootpoint;
		movePositions.second += rootpoint;
		  // return
// std::cerr <<"done\n";
		return movePositions;
	}

	PullMoveSet::IPointPair
	PullMoveSet::PullMoveDecoder::lookupStdMove(const IntPoint& pullPoint,
			const IntPoint& fixedPoint,
			const size_t& moveIndex) const
	{
// std::cerr <<" lookupEndMove done ( "<<pullPoint <<", "<<fixedPoint <<", " <<moveIndex <<") ... ";
		const IntPoint moveDirection = fixedPoint - pullPoint;

		assertbiu( neighDirectionsAccess.find( moveDirection ) != neighDirectionsAccess.end(), "(fixedPoint - pullPoint) is no known neighboring vector");

// std::cerr <<"done\n";
		return stdPulls.at( neighDirectionsAccess.find( moveDirection )->second ).at(moveIndex);

//		for (size_t i = 0;	i < directionNumber; i++)
//		{
//			// check where to look in stdPulls
//			if (neighDirections.at(i) == moveDirection)
//			{
//				return stdPulls.at(i).at(moveIndex);
//			}
//		}
//
//		return IPointVec();
	}

	size_t
	PullMoveSet::PullMoveDecoder::getEndMoveNumber() const
	{
		return endPullNumber;
	}

	size_t
	PullMoveSet::PullMoveDecoder::getPullMoveNumber() const
	{
		return stdPullNumber;
	}

	size_t
	PullMoveSet::PullMoveDecoder::getMoveNumber(const size_t proteinLength) const
	{
		return	endPullNumber*2 +							// # of end pulls
				(proteinLength-1)*stdPullNumber*2;			// # of std pulls
	}

	const bool
	PullMoveSet::PullMoveDecoder::isTriangular() const {
		return triangularLattice;
	}

	const LatticeModel* const
	PullMoveSet::PullMoveDecoder::getLattice(void) const {
		return lattice;
	}

	PullMoveSet::PullMoveSet(const LatticeModel* lattice)
	 :	LatticeMoveSet(lattice),
	 	decoder(new PullMoveSet::PullMoveDecoder(lattice)),
	 	decoderIsShared(false)
	{
		assertbiu( lattice != NULL, "no lattice model given");
	 	undoRec.lastChangedObject = NULL;
	}

	PullMoveSet::PullMoveSet(PullMoveSet::PullMoveDecoder* decoder_,
					const bool decoderIsShared_ )
	 :	LatticeMoveSet(decoder_->getLattice()),
	 	decoder(decoder_),
	 	decoderIsShared(decoderIsShared_)
	{
		assertbiu(decoder_ != NULL, "no PullMoveDecoder given");
		if (!decoderIsShared) {
			decoder = new PullMoveSet::PullMoveDecoder(*decoder_);
		}
	 	undoRec.lastChangedObject = NULL;
	}

	PullMoveSet*
	PullMoveSet::clone() {
		return new PullMoveSet(*this);
	}

	PullMoveSet::PullMoveSet(const PullMoveSet& moveSet)
		: 	LatticeMoveSet(moveSet.lattice),
			undoRec(moveSet.undoRec),
			decoder(moveSet.decoder),
			decoderIsShared(moveSet.decoderIsShared)
	{
		if (!decoderIsShared) {
			decoder = new PullMoveSet::PullMoveDecoder(*(moveSet.decoder));
		}
	}

	PullMoveSet&
	PullMoveSet::operator=(const PullMoveSet& moveSet2) {
		if (this != &moveSet2) {
			assertbiu(decoderIsShared == moveSet2.decoderIsShared, "both move sets have either to share or not");
			lattice = moveSet2.lattice;
			if (decoderIsShared) {
				decoder = moveSet2.decoder;
			} else {
				delete decoder;
				decoder = new PullMoveSet::PullMoveDecoder(*(moveSet2.decoder));
			}
			undoRec = moveSet2.undoRec;
		}
		return *this;
	}

	PullMoveSet::~PullMoveSet()
	{
		if (!decoderIsShared)
			delete decoder;
	}

	LatticeProtein*
	PullMoveSet::applyMove(const LatticeProtein* const todo,
			const size_t moveIndex)  {

		// apply the pull move inplace to a copy of todo
		LatticeProtein_Ipnt* l = new LatticeProtein_Ipnt(*todo);
		LatticeProtein* l2 = applyMoveInPlace(l, moveIndex);

		// delete l2 if we don't need it!
		if (l2 == NULL) delete l;

		return l2;
	}

	LatticeProtein*
	PullMoveSet::applyMoveInPlace(LatticeProtein* _todo,
			const size_t moveIndex_)
	{
		size_t moveIndex = moveIndex_;

		LatticeProtein_Ipnt* todo = dynamic_cast<LatticeProtein_Ipnt*>(_todo);
		assertbiu(todo!=NULL,
				"Downcasting to LatticeProtein_Ipnt failed miserably.");

		const size_t todoSize = todo->getSequence().size();

		assertbiu(moveIndex >=0 && moveIndex < getMoveNumber(todo),
			"The pull move moveIndex has to be in [0, number of available moves)."
			<< getMoveNumber(todo));

		undoRec.hasChanged = false;

		// lookup PullFront, stdPull and the pull position
		bool pullFront;
		bool stdPull;
		size_t position;
		decoder->lookupMove(todoSize, moveIndex, pullFront, stdPull, position);

		// RelativeInt objects are used to be able to compute both pull
		// directions the same code
		const RelativeInt relPos(position, pullFront);
		const RelativeInt relLastElement((pullFront ? 0 : todoSize-1), pullFront);

		// const pointer to const reference to lattice protein IPoints
		const IPointVec* const todoConstPoints = todo->getPointsRef();
		// pointer to mutable LatticeProtein IPoints
		IPointVec* const todoPoints = todo->getPointsRef();

			// do triangular end pull
		if (decoder->isTriangular() && !stdPull) {
			IPointPair pullPoints =
							decoder->lookupEndMove(todoConstPoints->at(position),
													moveIndex);

				// move legal? -> pull positions not occupied
			for (IPointVec::const_iterator it = todoConstPoints->begin();
					it != todoConstPoints->end(); it++)
			{
					// if position occupied, leave
				if ( (*it == pullPoints.first) )
					return NULL;
			}

				// save undo information
			undoRec.hasChanged = true;
			undoRec.lastChangedObject = todo;
			undoRec.position = position;
			undoRec.pullFront = pullFront;

			// check how far to pull
			if (lattice->areNeighbored(	pullPoints.first,
										todoPoints->at(relPos-1)) ) {
				undoRec.lastPullPosition = relPos+0;
			}
			else {
				undoRec.lastPullPosition = relLastElement.getValue(); // default
				for (RelativeInt i(position, pullFront); i>= relLastElement+2; i--)
				{
					if (lattice->areNeighbored( todoPoints->at(i.getValue()),
												todoPoints->at(i-2))) {
						undoRec.lastPullPosition = i-1;
						break;
					}
				}
			}
			undoRec.lostPos0 = todoPoints->at(undoRec.lastPullPosition);

			// pull tail
			for (RelativeInt i(undoRec.lastPullPosition, pullFront);
					i <= relPos -1; i++)
			{
				todoPoints->at(i.getValue()) = todoPoints->at(i+1);
			}
			  // overwrite pulled positions
			todoPoints->at(relPos-0) = pullPoints.first;

			return todo;
		}

			// do triangular std pull
		if (decoder->isTriangular() && stdPull) {
			// look up pullPoint & fixedPoint
			const IntPoint pullPoint = todoConstPoints->at(position);
			const IntPoint fixedPoint = todoConstPoints->at(relPos + 1);

			// look up l
			const IPointPair clPair =
				decoder->lookupStdMove(pullPoint, fixedPoint, moveIndex);
			const IntPoint l = pullPoint + clPair.first;

			for (IPointVec::const_iterator it = todoConstPoints->begin();
			it != todoConstPoints->end(); it++)
			{
				// if L occupied, leave
				if (*it == l) return NULL;
			}

			// compute how much needs to be moved
			int lastPullPosition=-1;					// last element to be moved
			bool lastPullPositionFound = false;

			// if moving first element while pulling front
			// or if moving last element while pulling back
			// or if moving second element
			if (relPos <= relLastElement + 1)
			{
				lastPullPosition = relLastElement.getValue();
				lastPullPositionFound = true;
			}
			// if moving anything above second element
			else
			{
				// if l and position-1 are neighbored
				if (lattice->areNeighbored(l, todoConstPoints->at(relPos-1)))
				{
					// just change relPos
					lastPullPosition = relPos+0;
					lastPullPositionFound = true;
				}
			}
			// if moving anything above second element
			if (relPos >= relLastElement + 2
					&& !lastPullPositionFound)
			{
				// check beyond position-1
				for (RelativeInt i(position, pullFront);
					(relLastElement <= i-2) && !lastPullPositionFound; i--)
				{
					if (lattice->areNeighbored(
							todoPoints->at(i.getValue()),
							todoPoints->at(i-2)))
					{
						// pull until position i-1
						lastPullPosition = i-1;
						lastPullPositionFound = true;
					}
				}
			}
			if (!lastPullPositionFound) lastPullPosition = relLastElement.getValue();
				assertbiu(lastPullPosition!=-1,
						"Unhandled case in PullMoveSet move calculation.");

			// save values for undo:
			undoRec.lastChangedObject = todo;
			undoRec.hasChanged = true;
			undoRec.pullFront = pullFront;
			undoRec.position = position;
			undoRec.lastPullPosition = lastPullPosition;

			// move
			const RelativeInt relLastPullPosition(lastPullPosition, pullFront);
			if (relLastPullPosition <= relPos-1)
			{
				// save for undo
				undoRec.lostPos0 = todoConstPoints->at(relLastPullPosition.getValue());
				// pull tail if existent
				for (RelativeInt i(lastPullPosition, pullFront); i<=relPos-1; i++)
				{
					todoPoints->at(i.getValue()) = todoConstPoints->at(i+1);
				}
				todoPoints->at(relPos.getValue()) = l;
			}
			else if (lastPullPosition == relPos.getValue())
			{
				// save for undo
				undoRec.lostPos0 = todoConstPoints->at(relPos.getValue());
				todoPoints->at(relPos.getValue()) = l;
			}

			return todo;
		}

		// do end pull
		if (!stdPull)
		{
			IPointPair pullPoints =
				decoder->lookupEndMove(todoConstPoints->at(position),
										moveIndex);

			// move legal? -> pull positions not occupied
			for (IPointVec::const_iterator it = todoConstPoints->begin();
					it != todoConstPoints->end(); it++)
			{
				// if position occupied, leave
				if ( (*it == pullPoints.first) || (*it == pullPoints.second) )
					return NULL;
			}

			// save undo information
			undoRec.hasChanged = true;
			undoRec.lastChangedObject = todo;
			undoRec.position = position;
			undoRec.pullFront = pullFront;

			// check how far to pull
			if (lattice->areNeighbored(	pullPoints.first,
										todoPoints->at(relPos-2)) ) {
				undoRec.lastPullPosition = relPos-1;
			}
			else {
				undoRec.lastPullPosition = relLastElement.getValue(); // default
				for (RelativeInt i(position, pullFront); i>= relLastElement+3; i--)
				{
					if (lattice->areNeighbored( todoPoints->at(i.getValue()),
												todoPoints->at(i-3))) {
						undoRec.lastPullPosition = i-2;
						break;
					}
				}
			}
			undoRec.lostPos0 = todoPoints->at(undoRec.lastPullPosition);
			undoRec.lostPos1 = todoPoints->at(RelativeInt(
													undoRec.lastPullPosition,
													pullFront)+1);

			// pull tail
			for (RelativeInt i(undoRec.lastPullPosition, pullFront);
					i <= relPos -2; i++)
			{
				todoPoints->at(i.getValue()) = todoPoints->at(i+2);
			}
			  // overwrite pulled positions
			todoPoints->at(relPos-1) = pullPoints.first;
			todoPoints->at(relPos.getValue()) = pullPoints.second;

			return todo;
		}

		// look up pullPoint & fixedPoint
		const IntPoint pullPoint = todoConstPoints->at(position);
		const IntPoint fixedPoint = todoConstPoints->at(relPos + 1);

		// look up c & l
		const IPointPair clPair =
			decoder->lookupStdMove(pullPoint, fixedPoint, moveIndex);
		const IntPoint c = pullPoint + clPair.first;
		const IntPoint l = pullPoint + clPair.second;

		// check if C and L are valid positions
		// if C position occupied by element in pull direction
		// or if moving first element while pulling front
		// or if moving last element while pulling back
//		if ( (pullFront && position == 0)
//				|| (!pullFront && position == todoSize-1)
//				|| (c==todoConstPoints->at(relPos - 1)) )
//		// just check for L position
//		{
//			for (IPointVec::const_iterator it = todoConstPoints->begin();
//				it != todoConstPoints->end(); it++)
//			{
//	 			// if L position occupied, leave
//				if ((*it) == l) return NULL;
//			}
//		}
		// check both C and L positions
//		else
//		{
			for (IPointVec::const_iterator it = todoConstPoints->begin();
			it != todoConstPoints->end(); it++)
			{
				// if C or L position occupied, leave
				if (((*it) == l) || ((*it) == c)) return NULL;
			}
//		}

		// compute how much needs to be moved
		int lastPullPosition=-1;					// last element to be moved
		bool lastPullPositionFound = false;

		// if moving first element while pulling front
		// or if moving last element while pulling back
		// or if moving second element
		if (relPos <= relLastElement + 1)
		{
			lastPullPosition = relLastElement.getValue();
			lastPullPositionFound = true;
		}
		// if moving anything above second element
		else
		{
			// if c and position-2 are neighbored
			if (lattice->areNeighbored(c, todoConstPoints->at(relPos-2)))
			{
				// just change position and position-1
				lastPullPosition = relPos-1;
				lastPullPositionFound = true;
			}
		}
		// if moving anything above third element
		if (relPos >= relLastElement + 3
				&& !lastPullPositionFound)
		{
			// check beyond position-2
			for (RelativeInt i(position, pullFront);
					(relLastElement <= i-3) && !lastPullPositionFound; i--)
			{
				if (lattice->areNeighbored(
						todoPoints->at(i.getValue()),
						todoPoints->at(i-3)))
				{
					// pull until position i-2
					lastPullPosition = i-2;
					lastPullPositionFound = true;
				}
			}
		}
		if (!lastPullPositionFound) lastPullPosition = relLastElement.getValue();
		assertbiu(lastPullPosition!=-1,
				"Unhandled case in PullMoveSet move calculation.");

		// save values for undo:
		undoRec.lastChangedObject = todo;
		undoRec.hasChanged = true;
		undoRec.pullFront = pullFront;
		undoRec.position = position;
		undoRec.lastPullPosition = lastPullPosition;

		// move
		const RelativeInt relLastPullPosition(lastPullPosition, pullFront);
		if (relLastPullPosition <= relPos-2)
		{
			// save for undo
			undoRec.lostPos0 = todoConstPoints->at(relLastPullPosition.getValue());
			undoRec.lostPos1 = todoConstPoints->at(relLastPullPosition+1);
			// pull tail if existent
			for (RelativeInt i(lastPullPosition, pullFront); i<=relPos-2; i++)
			{
				todoPoints->at(i.getValue()) = todoConstPoints->at(i+2);
			}
			todoPoints->at(relPos-1) = c;
			todoPoints->at(relPos.getValue()) = l;
		}
		else if (lastPullPosition == relPos-1)
		{
			// save for undo
			undoRec.lostPos0 = todoConstPoints->at(relPos-1);
			todoPoints->at(relPos-1) = c;
			undoRec.lostPos1 = todoConstPoints->at(relPos.getValue());
			todoPoints->at(relPos.getValue()) = l;
		}
		else if (lastPullPosition == relPos.getValue())
		{
			// save for undo
			undoRec.lostPos1 = todoConstPoints->at(relPos.getValue());
			todoPoints->at(relPos.getValue()) = l;
		}

		return todo;
	}

	size_t
	PullMoveSet::getMoveNumber(const LatticeProtein * const lp) const
	{
		assertbiu(lp!=NULL,
					"PullMoveSet::getMoveNumber: " <<
					"lp is not allowed to be NULL.");
		return decoder->getMoveNumber(lp->getSequence().size());
	}

	LatticeProtein*
	PullMoveSet::undoLastMove(LatticeProtein* _toUndo)
	{
		LatticeProtein_Ipnt* toUndo =
			dynamic_cast<LatticeProtein_Ipnt*>(_toUndo);
		assertbiu(toUndo!=NULL,
					"PullMoveSet::undoLastMove: Downcasting to " <<
					"LatticeProtein_Ipnt failed.");

		IPointVec* const toUndoPoints = toUndo->getPointsRef();

		if (toUndo == undoRec.lastChangedObject && undoRec.hasChanged)
		{
			undoRec.hasChanged = false;

			const RelativeInt relPos(undoRec.position, undoRec.pullFront);
			const RelativeInt relLastPullPos
								(undoRec.lastPullPosition, undoRec.pullFront);

				// undo triangular pull moves
			if (decoder->isTriangular()) {
				if ((int)undoRec.lastPullPosition == relPos.getValue()) {
					toUndoPoints->at(relPos.getValue()) = undoRec.lostPos0;
				}
				else if (relLastPullPos <= relPos-1)
				{
					// undo pulling whole tail
					for (RelativeInt i(undoRec.position, undoRec.pullFront);
							i>=relLastPullPos+1; i--)
					{
						toUndoPoints->at(i.getValue()) = toUndoPoints->at(i-1);
					}
					toUndoPoints->at(relLastPullPos.getValue()) = undoRec.lostPos0;
				}

				return toUndo;
			}

				// undo regular pull moves
			if ((int)undoRec.lastPullPosition == relPos.getValue())
			{
				toUndoPoints->at(relPos.getValue()) = undoRec.lostPos1;
			}
			if ((int)undoRec.lastPullPosition == relPos-1)
			{
				toUndoPoints->at(relPos.getValue()) = undoRec.lostPos1;
				toUndoPoints->at(relPos-1) = undoRec.lostPos0;
			}
			if (relLastPullPos <= relPos-2)
			{
				// undo pulling whole tail
				for (RelativeInt i(undoRec.position, undoRec.pullFront);
						i>=relLastPullPos+2; i--)
				{
					toUndoPoints->at(i.getValue()) = toUndoPoints->at(i-2);
				}
				toUndoPoints->at(relLastPullPos.getValue()) = undoRec.lostPos0;
				toUndoPoints->at(relLastPullPos + 1) = undoRec.lostPos1;
			}
		}

		return toUndo;
	}

	const PullMoveSet::PullMoveDecoder* const
	PullMoveSet::getDecoder(void) const {
		return decoder;
	}
}
