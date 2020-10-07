// $Id: PivotMoveSet.cc,v 1.2 2016/08/08 12:42:00 mmann Exp $

#include <biu/PivotMoveSet.hh>

namespace biu
{
	PivotMoveSet::PivotMoveSet(const LatticeModel* lattice)
	 :	LatticeMoveSet(lattice) , 
	 	lastOverwrittenMove(lattice->getNeighborhood().size()-1),
	 	lastOverwritePos(0),
	 	lastChangedObject(NULL)
	{} 

	PivotMoveSet::~PivotMoveSet()
	{}
	
	PivotMoveSet* 
	PivotMoveSet::clone() {
		return new PivotMoveSet(*this);
	}

// abstract functions (LatticeMoveSet)
	LatticeProtein*
	PivotMoveSet::applyMove(const LatticeProtein* const todo, 
							size_t moveIndex)  {
			// apply the pivot move inplace to a copy of todo
		return applyMoveInPlace(new LatticeProtein_Ipnt(*todo), moveIndex);
	}

	LatticeProtein* 
	PivotMoveSet::applyMoveInPlace(LatticeProtein* todo, 
	 size_t _moveIndex)  {

		assertbiu(_moveIndex >=0 && _moveIndex < getMoveNumber(todo),
	"The pivot move moveIndex has to be in [0, number of available moves).");

			//compute position & moveIndex from _moveIndex
		size_t position = (_moveIndex / (lattice->getNeighborhood().size()-1))+1;
		size_t moveIndex = _moveIndex % (lattice->getNeighborhood().size()-1);
		
			// get the mutated relative move
		Move mutatedMove = lattice->getNeighborhood().getElementByIndex(moveIndex).getMove();

			// if the current relative move at 'position' is equal to the chosen
			// mutated move, choose the move with the highest move index
			// this one is not chosen randomly
			// by that it is ensured that always a random neighbor is returned
		if (todo->getMoveSeqRel()[position] == mutatedMove) {
			mutatedMove = lattice->getNeighborhood().getElementByIndex(
					lattice->getNeighborhood().size() - 1).getMove();
		}

		lastChangedObject = NULL;

		if ( todo->getMoveSeqRel().at(position) != mutatedMove) {
				// save for next undo
			lastOverwrittenMove	= todo->getMoveSeqRel().at(position);
			lastOverwritePos 	= position;
			lastChangedObject	= todo;
			
			MoveSequence relMoves = todo->getMoveSeqRel();
				// change the relative move string at position to the mutated move
			relMoves.at(position) = mutatedMove;
			
			LatticeProtein_Ipnt* lpi = dynamic_cast<LatticeProtein_Ipnt*>(todo);
			assertbiu(lpi!=NULL, "Downcasting to LatticeProtein_Ipnt failed.");
			
				// update points accordingly
			IPointVec* points = lpi->getPointsRef();
			*points = lattice->relMovesToPoints(relMoves);
			
				// check selfavoidingness
			if (!todo->isSelfAvoiding())
				undoLastMove(todo);
		}

		return todo;
	}

	size_t
	PivotMoveSet::getMoveNumber( const LatticeProtein* const lp ) const {
			// number of pivot moves is the number of relative moves minus one
			// for the current relative move multiplied by the number of movable
			// elements
		return (lattice->getNeighborhood().size()-1)*(lp->getLength()-2);
	}


	LatticeProtein* 
	PivotMoveSet::undoLastMove(LatticeProtein* toUndo) {
		if (toUndo == lastChangedObject 
			&& toUndo->getMoveSeqRel()[lastOverwritePos] != lastOverwrittenMove) 
		{
			MoveSequence relMoves = toUndo->getMoveSeqRel();
			relMoves[lastOverwritePos] = lastOverwrittenMove;
			LatticeProtein_Ipnt* lpi = dynamic_cast<LatticeProtein_Ipnt*>(toUndo);
			assertbiu(lpi!=NULL, "Downcasting to LatticeProtein_Ipnt failed.");
			IPointVec* points = lpi->getPointsRef();
			*points = lattice->relMovesToPoints(relMoves);
		}
		
		return toUndo;
	}
}
