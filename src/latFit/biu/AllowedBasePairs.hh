// $Id: AllowedBasePairs.hh,v 1.2 2016/08/08 12:41:56 mmann Exp $
#ifndef BIU_ALLOWEDBASEPAIRS_HH_
#define BIU_ALLOWEDBASEPAIRS_HH_

#include <set>

#include "biu/Alphabet.hh"
#include "biu/Matrix.hh"

namespace biu
{
	/*! This class handels possible base pairs of an RNA structure.
	 *
	 * @author Salem Dekelbab
	 * @author Martin Mann
	 * @author Sebastian Will
	 * @author Andreas Richter
	 */	
	class AllowedBasePairs {
	private:
			//! Allowed base pairs of an RNA structure
		typedef std::pair<Alphabet::AlphElem, Alphabet::AlphElem> BPair;
		
			//! the Alphabet the allowed base pairs are basing on
		const Alphabet* alph;
		
			//! Matrix of the allowed combinations of base pairs
		biu::Matrix< bool > allowedPairs ;
		
	public:
			/*! Construction 
			 * @param alphabet has to be a structure alphabet for a RNA
			 * @param bps the allowed combinations of alphabet elements as
			 * comma separated list of alphabet element, e.g. "AU,GU".
			 * @param symmetric when true, it builds a symmetric matrix, 
			 * 			so supplying AU is enough for UA to be valid as well
			 */
		AllowedBasePairs(	const Alphabet* alphabet
							, const std::string& bps
							, const bool symmetric = true );
		
		~AllowedBasePairs();

		bool operator== (const AllowedBasePairs& abp2) const;
		bool operator!= (const AllowedBasePairs& abp2) const;
		
			/*! Returns if two characters form a valid base pair, e.g.
			 *  a Watson-Crick (AU, CG) or a non-standard (GU) pair.
			 */
		bool allowedBasePair(	const Alphabet::AlphElem& first,
			 	 				const Alphabet::AlphElem& second) const;
			 	 				
			 //! Returns the alphabet the base pairs are basing on
		const Alphabet* getAlphabet() const;
		
		
	};	


} // namespace biu

#include "biu/AllowedBasePairs.icc"

#endif /*ALLOWEDBASEPAIRS_HH_*/
