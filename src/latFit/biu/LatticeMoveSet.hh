// $Id: LatticeMoveSet.hh,v 1.2 2016/08/08 12:41:56 mmann Exp $
#ifndef BIU_LATTICEMOVESET_H_
#define BIU_LATTICEMOVESET_H_


#include "biu/LatticeModel.hh"
#include "biu/LatticeProtein.hh"

namespace biu
{

		/*! This class provides an interface for move set implementation
		 * to define the neighborhood in the energy landscape of LatticeProtein
		 * objects.
		 * 
		 * An object of this class is able to generate a neighbored LatticeProtein
		 * by applying one of its moves.
		 * 
		 * @author Martin Mann, Sebastian Will, Andreas Richter
		 */
	class LatticeMoveSet
	{
	protected:
			/*! the lattice model this move set bases on */
		const LatticeModel* lattice;
			
	public:
			/*! Constructs a LatticeMoveSet on a given lattice model.
			 * @param lattice has to be != NULL
			 */
		LatticeMoveSet(const LatticeModel* lattice);
		
		virtual LatticeMoveSet* clone() = 0;
		
		virtual ~LatticeMoveSet();
		
			/*! Applies a special move to the given LatticeProtein on the 
			 * specified position and returns a new LatticeProtein object 
			 * as result.
			 * 
			 * @param todo the initial state
			 * @param moveIndex the move index (has to be less than the return 
			 * 					value of getMoveNumber())
			 * @return a new LatticeProtein object neighbored to the todo object
			 *
			 *  In case of invalid moves, the behaviour
			 *  of the program is not specified.
			 */
		virtual LatticeProtein* applyMove(const LatticeProtein* const todo, 
									size_t moveIndex) = 0;
									
			/*! Applies a special move to the given LatticeProtein on the 
			 * specified position INPLACE and returns the modified object
			 * as result.
			 * 
			 * @param todo the initial state to modify
			 * @param moveIndex the move index (has to be less than the return 
			 * 					value of getMoveNumber())
			 * @return the modified LatticeProtein or NULL if the move is not
			 *         valid.
			 */
		virtual LatticeProtein* applyMoveInPlace(LatticeProtein* const todo,
									size_t moveIndex) = 0;
		
			/*! Returns the number of available moves by this MoveSet. */
		virtual size_t getMoveNumber( const LatticeProtein* const lp) const = 0;
		
		/*!
		 * Undo of last move performed on LatticeProtein toUndo. Undo
		 * is performed in place!
		 */
		virtual LatticeProtein* undoLastMove(LatticeProtein* toUndo) = 0;
	};

} // namespace biu

#include "LatticeMoveSet.icc"

#endif /*LATTICEMOVESET_H_*/
