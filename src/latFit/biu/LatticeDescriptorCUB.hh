// $Id: LatticeDescriptorCUB.hh,v 1.2 2016/08/08 12:42:01 mmann Exp $
#ifndef BIU_LATTICEDESCRIPTORCUB_HH_
#define BIU_LATTICEDESCRIPTORCUB_HH_


#include "biu/LatticeDescriptor.hh"
#include <cstdlib>

namespace biu
{

		/**
		 * Lattice property handler for the integer cubic lattice.
		 * 
		 * The neighboring vectors are defined by all permutations of (1,0,0)
		 * and its negations. Thus 6 neighboring vectors are possible.
		 * 
		 * @author Martin Mann, Sebastian Will, Andreas Richter
		 */
	class LatticeDescriptorCUB : public LatticeDescriptor
	{
	protected:

			// abstract function implementation
		virtual unsigned int getNeighborDataSize() const;
		virtual const NeighborData *getNeighborData() const;
		virtual unsigned int getAutomorphismDataSize() const;
		virtual const AutomorphismData *getAutomorphismData() const;
	public:
	
		LatticeDescriptorCUB();
		LatticeDescriptorCUB(const LatticeDescriptorCUB& toCopy);
		virtual ~LatticeDescriptorCUB();
			//! Returns whether or not two points in the lattice are neighbored
			//! using a specialised implementation for the cubic lattice
			//! @param first the first point
			//! @param second the second point
			//! @return true if (second-first) is element of the neighborhood,
			//!         false otherwise
		virtual bool areNeighbored( const IntPoint &first, 
									const IntPoint &second ) const;
		
		  /*! Calculates the base vector scaling 
		   * that scales all neighboring vectors to the given length. 
		   * 
		   * @param neighVecLength the length that the neighbor vectors should
		   *         be scaled to
		   * 
		   * @return the multiplicator for the base vectors to scale the lattice
		   */
		virtual
		double
		getBaseScale( const double neighVecLength ) const;


			//! Checks if the given point is a valid node of the lattice or not
			//! @param p the point to check
			//! @return true (all IntPoints are reachable)
		virtual bool isLatticeNode( const IntPoint & p ) const;
		
		
		  /*! Checks whether or not a given ring size is possible as a 
		   * selfavoiding walk with equal start/end position.
		   * @param ringSize the size of the ring to check for (at least 3)
		   * @return true if the ring size is even and greater than 3; 
		   *         false otherwise
		   */
		virtual 
		bool
		isPossibleRing( const size_t ringSize ) const;

	};

}


namespace biu
{
		/*! 
		 * A specialized LatticeNeighborhood to handle the
		 * the NeighborVector objects of the cubic lattice and the 
		 * access to them.
		 */
	class LatticeNeighborhoodCUB : public LatticeNeighborhood
	{
	public:
		LatticeNeighborhoodCUB(	const MoveAlphabet* moveAlph, 
								const NeighSet& neighbors)
		 :	LatticeNeighborhood(moveAlph, neighbors) 
		{
		}
		
		LatticeNeighborhoodCUB(	const LatticeNeighborhood& nh ) 
		:	LatticeNeighborhood( nh )
		{
		}
		
								
		virtual ~LatticeNeighborhoodCUB()
		{}
		
		
			//! Returns whether or not a vector v belongs to the neighborhood.
		virtual bool isElement( const IntPoint& v) const {
			int sum = abs(v.getX());
			return	sum < 2 
					&& (sum+=abs(v.getY())) < 2 
					&& (sum+=abs(v.getZ())) == 1;
/*
			return (v.getX()* v.getX() + v.getY()*v.getY() + v.getZ()*v.getZ()) 
					== 1;
			return (abs( v.getX() + v.getY() + v.getZ() ) == 1) && 
			((abs(v.getX()) == 1 && v.getY() == 0 && v.getZ() == 0) ||
			(v.getX() == 0 && 
				((abs(v.getY()) == 1 && v.getZ() == 0) ||
				 (v.getY() == 0 && abs(v.getZ()) == 1))) 
			);
*/
		}
	};

} // namespace biu


#endif /*LATTICEDESCRIPTORCUB_HH_*/
