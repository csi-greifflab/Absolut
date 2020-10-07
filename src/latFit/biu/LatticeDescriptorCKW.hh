#ifndef BIU_LATTICEDESCRIPTORCKW_HH_
#define BIU_LATTICEDESCRIPTORCKW_HH_

#include "biu/LatticeDescriptor.hh"

namespace biu
{

	/*!
	 * Lattice descriptor of the chess knight walk lattice also known as the
	 * 210 lattice.
	 * 
	 * The neighborhood vertices are all permutations of (2,1,0) and its 
	 * negations. Thus 24 neighboring vectors are possible.
	 * 
	 * @author Martin Mann http://www.bioinf.uni-freiburg.de/~mmann/
	 * 
	 */
	class LatticeDescriptorCKW : public LatticeDescriptor
	{
	protected:

			// abstract function implementation
		virtual unsigned int getNeighborDataSize() const;
		virtual const NeighborData *getNeighborData() const;
		virtual unsigned int getAutomorphismDataSize() const;
		virtual const AutomorphismData *getAutomorphismData() const;
	public:
	
		 /*! Construction
		  */
		LatticeDescriptorCKW();
		
		 /*! Copy construction.
		  * @param toCopy the CKW lattice descriptor to copy
		  */
		LatticeDescriptorCKW( const LatticeDescriptorCKW& toCopy);
		
		 /*! Destruction
		  */
		virtual ~LatticeDescriptorCKW();
		
			//! Returns whether or not two points in the lattice are neighbored
			//! using a specialised implementation for the chess knight walk
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
			//! @return true (all IntPoints are reachable)
		virtual bool isLatticeNode( const IntPoint & p ) const; 
		
		  /*! Checks whether or not a given ring size is possible as a 
		   * selfavoiding walk with equal start/end position.
		   * @param ringSize the size of the ring to check for (at least 3)
		   * @return true if a SAW ring of that size is possible; 
		   *         false otherwise
		   */
		virtual 
		bool
		isPossibleRing( const size_t ringSize ) const {
			return false;
		}

	};

}

#endif /*LATTICEDESCRIPTORCKW_HH_*/
