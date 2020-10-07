// $Id: LatticeModel.hh,v 1.2 2016/08/08 12:41:59 mmann Exp $
#ifndef BIU_LATTICEMODEL_H_
#define BIU_LATTICEMODEL_H_


#include "biu/Point.hh"
#include "biu/LatticeDescriptor.hh"
#include <set>

namespace biu
{
	

		/*! A lattice model handles the lattice type speficic calculations
		 *  based on the specified biu::LatticeDescriptor object.
		 * 
		 * This class is limited to an integer lattice point representation.
		 * 
		 * @author Martin Mann, Sebastian Will, Andreas Richter
		 */
	class LatticeModel
	{
	protected:
			/*! The lattice property handler (!= NULL). */
		LatticeDescriptor const* const latDescriptor;
		
			//! short reference to neighborhood of the lattice 
		const LatticeNeighborhood & latNeighborhood;
		
	public:
			//! construction
			//! @param _latDescriptor This lattice property handler has to
		    //! be != NULL.
		LatticeModel(const LatticeDescriptor* const _latDescriptor); 
		
		LatticeModel(const LatticeModel& toCopy);
		virtual ~LatticeModel();
		
			//! Returns a constant reference to the current lattice descriptor.
		LatticeDescriptor const* const getDescriptor() const {
			return latDescriptor;
		}
		
			//! Direct access to lattice neighborhood.
		const LatticeNeighborhood & getNeighborhood() const {
			return latNeighborhood;
		}
		
		
		bool operator == (const LatticeModel& lm2) const;
		bool operator != (const LatticeModel& lm2) const;
		
		
		//////////////////////////////////////////////////////
		// conversion of move sequences
		//////////////////////////////////////////////////////
		
			//! Convert a sequence of relative moves 
			//! to a sequence of absolute moves
		MoveSequence relMovesToAbsMoves( const MoveSequence &relMoves ) const;
		
			//! Convert a sequence of absolute moves
			//! to a sequence of relative moves
		MoveSequence absMovesToRelMoves( const MoveSequence &absMoves ) const;
		
		//////////////////////////////////////////////////////
		// conversion moves <=> coordinates
		//////////////////////////////////////////////////////
			
			//! Converts two points into an absolute move. 
			//!
			//! They have to be neighbored in the lattice. 
		virtual Move getAbsMove(	const IntPoint& lastPoint, 
									const IntPoint& actPoint ) const;

			//! Converts a vector of points to an absolute move string.
			//!
			//! Consecutive positions have to be neighbored in the lattice!
		virtual MoveSequence pointsToAbsMoves( const IPointVec& points ) const;

			//! Converts an absolute/relative move string representation into
			//! an absolute/relative move sequence.
		virtual MoveSequence parseMoveString( const std::string& moveStr) const;

			//! Converts a vector of points to a relative move string.
			//!
			//! Consecutive positions have to be neighbored in the lattice! 
		virtual MoveSequence pointsToRelMoves( const IPointVec& points) const;

			//! Applys an absolute move to a point, returns result.
		virtual IntPoint applyAbsMove(	const IntPoint& actPoint, 
										const Move& absMove ) const;

			//! Converts an absolute move string to coordinates. 
		virtual IPointVec absMovesToPoints( const MoveSequence& absMoves) const;

			//! Converts a relative move string to coordinates. 
		virtual IPointVec relMovesToPoints( const MoveSequence& relMoves) const;

			//! Converts a relative move string to coordinates written to
			//! the given vector toFill. 
		virtual IPointVec& relMovesToPoints( const MoveSequence& relMoves, 
			IPointVec& toFill) const;
			 
			//! Returns the string representation of the move sequence.
		virtual std::string	getString( const MoveSequence& moveSeq ) const;
		
		//////////////////////////////////////////////////////
		// nachbarschaft im gitter
		//////////////////////////////////////////////////////
		
			//! Calculates all neighbored points to the given center based
			//!  on the current lattice descriptor.
		virtual IPointSet	getAllNeighPoints( const IntPoint& center ) const;

			//! Return whether or not two points in the lattice are neighbored,
			//! based on the current lattice descriptor.
			//! @param first the first point
			//! @param second the second point
			//! @return true if (second-first) is element of the neighborhood,
			//!         false otherwise
		virtual bool areNeighbored( const IntPoint &first, 
									const IntPoint &second ) const;
		
	
	};

} // namespace biu


#include "LatticeModel.icc"

#endif /*LATTICEMODEL_H_*/
