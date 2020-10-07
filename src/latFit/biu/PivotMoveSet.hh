// $Id: PivotMoveSet.hh,v 1.2 2016/08/08 12:41:58 mmann Exp $
#ifndef BIU_PIVOTMOVESET_HH_
#define BIU_PIVOTMOVESET_HH_


#include "biu/LatticeMoveSet.hh"
#include "biu/LatticeProtein_Ipnt.hh"

namespace biu
{
		/**
		 * This class handles pivot moves on LatticeProtein objects.
		 * The pivot moves are implemented as point mutations on the
		 * relative move string.
		 * 
		 * @author Martin Mann, Sebastian Will, Andreas Richter
		 */
	class PivotMoveSet : public LatticeMoveSet
	{
	protected:
	
			//! the last move overwritten by this moveset
		Move			lastOverwrittenMove;
			//! last position that was overwritten
		size_t			lastOverwritePos;
			//! pointer to the last object that was changed
		LatticeProtein*	lastChangedObject;
		
	public:
		PivotMoveSet(const LatticeModel* lattice); 

		virtual ~PivotMoveSet();
		
		PivotMoveSet* clone();

	// abstract functions (LatticeMoveSet)
			/*! Applies a pivot move to the given LatticeProtein on the 
			 * specified position and returns a new LatticeProtein object 
			 * as result.
			 * 
			 * @param toChange the initial state
			 * @param moveIndex the move index (has to be less than the return 
			 * 					value of getMoveNumber())
			 * @return a new LatticeProtein object neighbored to the todo object
			 *
			 *  In case of invalid moves, the behaviour
			 *  of the program is not specified.
			 */
		virtual LatticeProtein* applyMove(const LatticeProtein* const toChange, 
									 size_t moveIndex) ;

			/*! Applies a pivot move to the given LatticeProtein on the 
			 * specified position INPLACE and returns the modified object
			 * as result.
			 *
			 * @param toChange the initial state to modify
			 * @param moveIndex the move index (has to be less than the return 
			 * 					value of getMoveNumber())
			 * @return the modified LatticeProtein 
			 *
			 *  In case of invalid moves, the behaviour
			 *  of the program is not specified.
			 */
		virtual LatticeProtein* 
		applyMoveInPlace(LatticeProtein* toChange, size_t moveIndex) ;
		
		virtual size_t 
		getMoveNumber( const LatticeProtein* const lp ) const;
		
		virtual LatticeProtein* undoLastMove(LatticeProtein* toUndo);
	};

} // namespace biu

#endif /*PIVOTMOVESET_HH_*/
