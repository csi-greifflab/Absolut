// $Id: RandomNumberFactory.hh,v 1.2 2016/08/08 12:41:58 mmann Exp $
#ifndef BIU_RANDOMNUMBERFACTORY_HH_
#define BIU_RANDOMNUMBERFACTORY_HH_

#include "biu/RandomNumberGenerator.hh"

namespace biu
{

	/*!
	 * A wrapper class for a central but variable random number generator 
	 * access.
	 * 
	 * @author Martin Mann
	 */
	class RandomNumberFactory 
	{
	private:
		//! The current random number generator in use. (default = RNG_ISO)
		static RandomNumberGenerator * rng;
		
	public:
		RandomNumberFactory();
		virtual ~RandomNumberFactory();
		
		
		//////////////  RNG SETUP  /////////////////////
		
		//! Access to the random number generator.
		//! @return the current RNG
		static RandomNumberGenerator& getRNG(void);
		
		//! Setting a new random number generator.
		//! @param rng the RNG to copy
		static void setRNG(RandomNumberGenerator &rng);
		
		//! Setting a new random number generator.
		//! @param rng the RNG to copy
		static void setRNG(RandomNumberGenerator *rng);
		
		
		
		
		//////////////  DIRECT ACCESS FUNCTION WRAPPER /////////////////////

		/*!
		 * Returns the next random number in the series. It's value will be in 
		 * [0, getMaxRN()].
		 * @return rng->getRN()
		 */
		static unsigned int getRN();
		
		/*!
		 * Returns the next random number in the series. It's value will be in 
		 * [0, max).
		 * @return (rng->getRN()%max)
		 */
		static unsigned int getRN(unsigned int max);
		
		/*!
		 * Returns the largest value the rand function will return.
		 * @return rng->getMaxRN()
		 */
		static unsigned int getMaxRN();
		
		//! pointer object to random generator function for RNF usage with STL
		//! algorithms like stl::random_shuffle
		static unsigned int (*pt_getRN)(unsigned int);
		
	}; // class RNF
	
	//! shortcut typedef
	typedef RandomNumberFactory RNF;
	
	
	
} // namespace biu



#endif /*RANDOMNUMBERFACTORY_HH_*/
