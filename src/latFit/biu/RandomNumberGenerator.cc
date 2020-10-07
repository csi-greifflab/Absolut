// $Id: RandomNumberGenerator.cc,v 1.2 2016/08/08 12:41:56 mmann Exp $

#include "biu/RandomNumberGenerator.hh"
#include <stdlib.h>
#include <biu/assertbiu.hh>

namespace biu {


	
	
	RNG_ISO::RNG_ISO(unsigned int seed_) : RandomNumberGenerator(seed_)
	{
		setSeed(seed_);
	}
	
	RNG_ISO::~RNG_ISO()
	{}
	
	void 
	RNG_ISO::setSeed(unsigned int _seed)
	{
		seed = _seed;
		srand(_seed);
	}
	
	unsigned int 
	RNG_ISO::getRN()
	{
		return rand();
	}
	
	unsigned int 
	RNG_ISO::getMaxRN() {
		return RAND_MAX;
	}
	
	RandomNumberGenerator* 
	RNG_ISO::copy(void) {
		return new RNG_ISO(*this);
	}
	
	
	
	
	
	
	
	
	
	
	RNG_LCG::RNG_LCG(unsigned int seed_) : RandomNumberGenerator(seed_)
	{
		setSeed(seed_);
	}
	
	RNG_LCG::~RNG_LCG()
	{}
	
	
	void 
	RNG_LCG::setSeed(unsigned int _seed)
	{
		assertbiu(_seed > 0, "seed value has to be > zero!");
		seed = _seed;
	}
	
	unsigned int 
	RNG_LCG::getRN()
	{
		
		// Wikipedia: using 64 bit integer arithmetics to compute 32 bit integer
		// V_{j+1} = (279470273 * V_j) % 4294967291
		//	
		//        uint64_t z;
		//        z = a;
		//        z *= 279470273; 
		//        z %= 4294967291U;
		//        a = z;
		//        return a;
		
		// same algorithm using shift operations
		// taken from rand_r.c (glibc 2.5)
		unsigned int next = seed;
		int result;
	
		next *= 1103515245;
		next += 12345;
		result = (unsigned int) (next / 65536) % 2048;
	
		next *= 1103515245;
		next += 12345;
		result <<= 10;
		result ^= (unsigned int) (next / 65536) % 1024;
	
		next *= 1103515245;
		next += 12345;
		result <<= 10;
		result ^= (unsigned int) (next / 65536) % 1024;
	
		seed = next;
		
		return result;
	}

	unsigned int 
	RNG_LCG::getMaxRN() {
		return 2147483647;
	}
	
	RandomNumberGenerator* 
	RNG_LCG::copy(void) {
		return new RNG_LCG(*this);
	}







	RNG_ARS4x32::RNG_ARS4x32(unsigned int seed_) :
		RandomNumberGenerator(seed_)
	{
		setSeed(seed_);
	}
	
	RNG_ARS4x32::~RNG_ARS4x32()
	{}
	
	
	void
	RNG_ARS4x32::setSeed(unsigned int _seed)
	{
//		key = { _seed, _seed, _seed, _seed}};
//		key[0] = key[1] = key[2] = key[3] = _seed;
		key[0] = _seed;
		key[1] = _seed;
		key[2] = _seed;
		key[3] = _seed;

		// Might as well start from 0.
		// The periodicity of the counter is 2^128.
		// Since generator(ctr, key) generates 4 32-bit ints,
		//   the periodicity of RNG_ARS4x32 is 2^130.
//		ctr = = {{0,0,0,0}};
//		ctr[0] = ctr[1] = ctr[2] = ctr[3] = 0;
		ctr[0] = 0;
		ctr[1] = 0;
		ctr[2] = 0;
		ctr[3] = 0;

		// reset random number index for recalculation
		randomNumberIndex=4;
	}

	
	unsigned int
	RNG_ARS4x32::getRN()
	{
		// Return random uint32. generator(ctr, key)
		//  produces four ints at a time, so it is called
		//  every four calls to getRN()
		if (!(randomNumberIndex % 4)) {
            randomNumbers = generator(ctr, key);
			randomNumberIndex -= 4;
			ctr.incr();
		}
		return (unsigned int)randomNumbers[randomNumberIndex++ % 4];
	}

	unsigned int 
	RNG_ARS4x32::getMaxRN() {
		return 4294967295; // 2^32 - 1
	}
	
	RandomNumberGenerator* 
	RNG_ARS4x32::copy(void) {
		return new RNG_ARS4x32(*this);
	}

	

} // namespace biu
