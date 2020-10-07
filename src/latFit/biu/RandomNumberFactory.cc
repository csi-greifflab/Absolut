// $Id: RandomNumberFactory.cc,v 1.2 2016/08/08 12:41:57 mmann Exp $


#include "biu/RandomNumberFactory.hh"
#include <biu/assertbiu.hh>

namespace biu
{

	// default random number generator
	RandomNumberGenerator* RandomNumberFactory::rng = new RNG_ISO();
	
	// pointer object to random generator function for RNF usage with STL
	// algorithms like stl::random_shuffle
	unsigned int (*RandomNumberFactory::pt_getRN)(unsigned int) = &RandomNumberFactory::getRN;


	RandomNumberFactory::RandomNumberFactory()
	{
	}
	
	RandomNumberFactory::~RandomNumberFactory()
	{
	}
	
	RandomNumberGenerator&
	RandomNumberFactory::getRNG() {
		return *rng;
	}
	
	void
	RandomNumberFactory::setRNG(RandomNumberGenerator& rng_) {

		delete rng; // delete old RNG
		
		rng = rng_.copy(); // set new one
	}
	
	void
	RandomNumberFactory::setRNG(RandomNumberGenerator* rng_) {

		assertbiu( rng_ != NULL, "given RNG is not available (NULL)");
		delete rng; // delete old RNG
		
		rng = rng_->copy(); // set new one
	}
	
	unsigned int
	RandomNumberFactory::getRN(void) {
		return rng->getRN();
	}
	
	unsigned int 
	RandomNumberFactory::getRN(unsigned int max) {
		assertbiu(max != 0, "maximal value == 0 would cause division by zero");
		return rng->getRN()%max;
	}
	
	unsigned int
	RandomNumberFactory::getMaxRN(void) {
		return rng->getMaxRN();
	}
	
} // namespace 
