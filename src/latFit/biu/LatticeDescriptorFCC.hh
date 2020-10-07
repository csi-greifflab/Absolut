// $Id: LatticeDescriptorFCC.hh,v 1.2 2016/08/08 12:41:58 mmann Exp $
#ifndef BIU_LATTICEDESCRIPTORFCC_HH_
#define BIU_LATTICEDESCRIPTORFCC_HH_


#include "biu/LatticeDescriptor.hh"

namespace biu
{
	
		/**
		 * Lattice property handler for the face centered cubic lattice (FCC).
		 * 
		 * @author Martin Mann, Sebastian Will, Andreas Richter
		 */
	class LatticeDescriptorFCC : public LatticeDescriptor
	{
	protected:
		virtual unsigned int getNeighborDataSize() const;
		virtual const NeighborData *getNeighborData() const;
		virtual unsigned int getAutomorphismDataSize() const;
		virtual const AutomorphismData *getAutomorphismData() const;
	public:
		LatticeDescriptorFCC();
		virtual ~LatticeDescriptorFCC();
			//! Returns whether or not two points in the lattice are neighbored
			//! using a specialised implementation for the face centered cubic 
			//! lattice
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
			//! @return true if the point has even coordinate sum;
			//!         false otherwise
		virtual bool isLatticeNode( const IntPoint & p ) const; 

		
		  /*! Checks whether or not a given ring size is possible as a 
		   * selfavoiding walk with equal start/end position.
		   * @param ringSize the size of the ring to check for (at least 3)
		   * @return true if the ring size is at least 3; 
		   *         false otherwise
		   */
		virtual 
		bool
		isPossibleRing( const size_t ringSize ) const;

	};

} // namespace biu

#endif /*LATTICEDESCRIPTORFCC_HH_*/
