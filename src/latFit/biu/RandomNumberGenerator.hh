// $Id: RandomNumberGenerator.hh,v 1.2 2016/08/08 12:42:00 mmann Exp $
#ifndef BIU_RANDOMNUMBERGENERATOR_HH_
#define BIU_RANDOMNUMBERGENERATOR_HH_


#include "Random123/ars.h"



namespace biu {

/*!
 * The abstract class works as a random number generator for unsigned integers.
 * 
 * @author Daniel Maticzka
 * @author Martin Mann
 */
class RandomNumberGenerator
{
protected:
	unsigned int seed;
	
public:
	/*!
	 * Creates a RandomNumberGenerator object and initialises the seed. 
	 * @param seed_ the initial seed value defaults to 1.
	 */
	RandomNumberGenerator(unsigned int seed_ = 1) : seed(seed_) {};
	RandomNumberGenerator(RandomNumberGenerator& c) : seed(c.seed) {};
	virtual ~RandomNumberGenerator() {};

	/*!
	 * Specifies the seed value used for random number generation. 
	 * @param seed the new seed value to use
	 */
	virtual void setSeed(unsigned int seed) = 0;
	
	/*!
	 * Returns the next random number in the series. It's value will be in 
	 * [0, getMaxRN()].
	 * @return the next random number
	 */
	virtual unsigned int getRN() = 0;
	
	/*!
	 * Returns the largest value the rand function will return.
	 * @return larges possible random number
	 */
	virtual unsigned int getMaxRN() = 0;
	
	/*!
	 * Creates a new copy of this object.
	 * @return new copy
	 */
	virtual RandomNumberGenerator* copy(void) = 0;
};

/*!
 * An object of this class works as a random number generator. This  
 * implementation maps to the rand & srand ISO C random number functions.
 * The random numbers derived from a particular seed value will differ
 * on different CPU types and C librarys.
 * For a defined set of random numbers, one of it's subclasses should be used.
 * 
 * @author Daniel Maticzka
 */
class RNG_ISO : public RandomNumberGenerator
{
private:
	using RandomNumberGenerator::seed;
	
public:
	/*!
	 * Creates a RandomNumberGenerator object. 
	 * @param seed the initial seed value defaults to 1.
	 */
	RNG_ISO(unsigned int seed = 1);
	virtual ~RNG_ISO();

	/*!
	 * Specifies the seed value used for random number generation. Because this
	 * implementation maps to ISO C rand & srand this will not be 
	 * specific to an object, it will affect the "global" seed value.
	 */
	virtual void setSeed(unsigned int _seed);
	
	/*!
	 * Returns the next random number in the series. It's value will be in 
	 * [0, getMaxRN()].
	 */
	virtual unsigned int getRN();
	
	/*!
	 * Returns the largest value the rand function will return.
	 */
	virtual unsigned int getMaxRN();

	/*!
	 * Creates a new Copy of this object.
	 */
	virtual RandomNumberGenerator* copy(void);
};

/**
 * This subclass of RandomNumberGenerator uses a linear congruent generator
 * to generate pseudo random numbers based on a seed value. The series of
 * random numbers generated does not depend on the CPU used.
 * 
 * @author Daniel Maticzka
 */
class RNG_LCG : public RandomNumberGenerator
{
private:
	using RandomNumberGenerator::seed;

public:
	/*!
	 * Creates a RNG_LCG object.
	 * @param seed the initial seed value defaults to 1.
	 */
	RNG_LCG(unsigned int seed = 1);
	virtual ~RNG_LCG();
	
	/*!
	 * Specifies the seed value used for random number generation.
	 */
	virtual void setSeed(unsigned int _seed);
	
	/*!
	 * Returns the next random number in the series. It's value will be in 
	 * [0, getMaxRN()].
	 */
	virtual unsigned int getRN();
	
	/*!
	 * Returns the largest value the rand function will return.
	 */
	virtual unsigned int getMaxRN();

	/*!
	 * Creates a new Copy of this object.
	 */
	virtual RandomNumberGenerator* copy(void);
};

/**
 * This subclass of RandomNumberGenerator uses a counter-based RNG
 * from the Random123 library (local version).
 * A CBRNG is stateless, generating a
 * number given an input key and a input counter value. Incrementing
 * the counter and giving the CBRNG the new (key, counter) pair
 * generates a different random number. The library is
 * CPU-architecture independent. See
 * https://www.deshawresearch.com/resources_random123.html
 *
 * The implemented RNG, ARS4x32, performs a simplified version
 * of AES, and requires the CPU to support AES-NI instructions.
 *
 * @author Victor Zhao
 */
class RNG_ARS4x32 : public RandomNumberGenerator
{
private:

	//!The particular RNG algorithm from r123 is ARS.
	//!   ARS4x32 generates 4 32-bit uints at a time.
    // Philippe 2019-09-25 adapt to newer version GSL typedef r123::ARS4x32 ars4x32;
    //typedef r123::r123array4x32 ars4x32;
    typedef r123::Array4x32 ars4x32;

	//! the random number generator
	ars4x32 generator;
	//! for ARS4x32, both key_type and ctr_type are
	// !  4x32-bit r123 arrays.
    //ars4x32::key_type
    ars4x32::value_type key;
    //ars4x32::ctr_type ctr;
    ars4x32::value_type  ctr;

    //! the last random number drawn (holds 4 random numbers per update)
    // ars4x32::ctr_type randomNumbers;
    r123::Array4x32 randomNumbers;
	//! next index within randomNumbers to return
	uint32_t randomNumberIndex;


public:
	/*!
	 * Creates a RNG_ARS4x32 object.
	 * @param seed Integer seed value, defaults to 1. Forms part of the CBRNG key.
	 */
	RNG_ARS4x32(unsigned int seed = 1);

	virtual ~RNG_ARS4x32();
	
	/*!
	 * Implements setSeed function
	 */
	virtual void setSeed(unsigned int _seed);
	
	/*!
	 * Returns the next random number in the series. It's value will be in 
	 * [0, getMaxRN()].
	 */
	virtual unsigned int getRN();
	
	/*!
	 * Returns the largest possible unsigned int32.
	 */
	virtual unsigned int getMaxRN();

	/*!
	 * Creates a new Copy of this object.
	 */
	virtual RandomNumberGenerator* copy(void);
};

} // end namespace biu

#endif /*RANDOMNUMBERGENERATOR_HH_*/
