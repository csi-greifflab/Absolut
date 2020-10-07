// $Id: LatticeNeighborhood.hh,v 1.2 2016/08/08 12:41:56 mmann Exp $
#ifndef BIU_LATTICE_NEIGHBORHOOD_HH_
#define BIU_LATTICE_NEIGHBORHOOD_HH_


#include "biu/NeighborVector.hh"

// get information on available hashes or maps
#include <biu/HashMap.hh>

// include best available hash or map
#if HAVE_UNORDERED_MAP == 1
	#include <unordered_map>
#else
#if HAVE_TR1_UNORDERED_MAP == 1
	#include <tr1/unordered_map>
#else
#if HAVE_GNU_HASH_MAP == 1
	#include <ext/hash_map>
#else
	#include <map>
#endif
#endif
#endif

namespace biu
{
		/*! 
		 * A LatticeNeighborhood handles the
		 * the NeighborVector objects of a lattice and the 
		 * access to them.
		 */
	class LatticeNeighborhood
	{
	protected:
		

		/*! An adaption of Daniel J. Bernstein's string hash function 
		 *  taken from http://www.cs.yorku.ca/~oz/hash.html.
		 * 
		 * Maybe not the best choice! Suggestions are welcome!
		 */
		struct hash_IntPoint {
			size_t operator()(const IntPoint& p) const
			{
				size_t hash = 5381;
				
				hash = ((hash << 5) + hash) + (size_t)p.getX(); // hash * 33 + coordinate X
				hash = ((hash << 5) + hash) + (size_t)p.getY(); // hash * 33 + coordinate Y
				hash = ((hash << 5) + hash) + (size_t)p.getZ(); // hash * 33 + coordinate Z
	
				return hash;
			}
			     
		};

		// set typedef for best available hash or map
		#if HAVE_UNORDERED_MAP == 1
			typedef std::unordered_map< IntPoint, const NeighborVector*, hash_IntPoint > P2N_MAP;
		#else
		#if HAVE_TR1_UNORDERED_MAP == 1
			typedef std::tr1::unordered_map< IntPoint, const NeighborVector*, hash_IntPoint > P2N_MAP;
		#else
		#if HAVE_GNU_HASH_MAP == 1
			typedef __gnu_cxx::hash_map< IntPoint, const NeighborVector*, hash_IntPoint > P2N_MAP;
		#else
			typedef std::map< IntPoint, const NeighborVector* > P2N_MAP;
		#endif
		#endif
		#endif
		
	protected:

		const MoveAlphabet* moveAlph;	//!< the underlying move string alphabet
		NeighSet			neighSet;	//!< the NeighborVector data
		
			//! a crossreference from move index to NeighborVectors
		std::vector<const NeighborVector*>	neighVec;	

			//! a crossreference from vectors to NeighborVectors
		P2N_MAP vec2neigh;
		
	public:
		LatticeNeighborhood(	const MoveAlphabet* moveAlph_, 
								const NeighSet& neighbors);
								
		LatticeNeighborhood(	const LatticeNeighborhood& nh );
								
		virtual ~LatticeNeighborhood();
		
			//! Returns the number of elements of the neighborhood.
		virtual 
		unsigned int size() const;
		
			//! Returns whether or not a vector belongs to the neighborhood.
			//! @param vector the vector to check
			//! @return true if the vector is part of the neighborhood,
			//!         false otherwise
		virtual 
		bool isElement( const IntPoint& vector) const;
		
			//! Returns the corresponding NeighborVector.
			//! @param index has to be the index of an element of this
			//!				neighborhood (has to be less than the return value
			//!				of size())
			//! @return the neighbor vector of the given index
		virtual 
		const NeighborVector& getElementByIndex(unsigned int index) const;

			//! Returns the corresponding NeighborVector.
			//! @param vector has to be element of this neighborhood
			//! @return the neighbor vector that corresponds to the given vector
		virtual 
		const NeighborVector& getElement(const IntPoint& vector) const;
		
			//! Returns the corresponding NeighborVector.
			//! @param move has to be element of the underlying move string 
			//!				alphabet
			//! @return the neighbor vector that corresponds to the given move
		virtual 
		const NeighborVector& getElement(const Move& move) const;
		
			//! a constant iterator to access the elements of the neighborhood
		typedef NeighSet::const_iterator const_iterator;
		
			//! Returns a constant iterator to the first element of the 
			//! neighborhood
		virtual 
		const_iterator begin() const;
		
			//! Returns a constant iterator behind the last element of the 
			//! neighborhood
		virtual 
		const_iterator end() const;
	};

} // namespace biu

#include "LatticeNeighborhood.icc"

#endif /*LATTICE_NEIGHBORHOOD_HH_*/
