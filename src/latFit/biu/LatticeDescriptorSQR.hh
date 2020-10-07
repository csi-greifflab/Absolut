// $Id: LatticeDescriptorSQR.hh,v 1.2 2016/08/08 12:41:56 mmann Exp $
#ifndef BIU_LATTICEDESCRIPTORSQR_HH_
#define BIU_LATTICEDESCRIPTORSQR_HH_


#include "biu/LatticeDescriptor.hh"

namespace biu
{

		/**
		 * Lattice property handler for the integer square lattice.
		 * 
		 * @author Martin Mann, Sebastian Will, Andreas Richter
		 */
	
	class LatticeDescriptorSQR : public LatticeDescriptor
	{
	protected:
		virtual unsigned int getNeighborDataSize() const;
		virtual const NeighborData *getNeighborData() const;
		virtual unsigned int getAutomorphismDataSize() const;
		virtual const AutomorphismData *getAutomorphismData() const;
	public:
		LatticeDescriptorSQR();
		virtual ~LatticeDescriptorSQR();
			//! Returns whether or not two points in the lattice are neighbored
			//! using a specialised implementation for the square lattice
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

			//! Checks if the given point is a valid node of the 2D-lattice or not
			//! @param p the point to check
			//! @return true if the point is element of the XY-plane (Z==0);
			//!    false otherwise
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

} // namespace biu

#endif /*LATTICEDESCRIPTORSQR_HH_*/
