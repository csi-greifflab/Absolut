// $Id: QuasiRandomNumberGenerator.hh,v 1.2 2016/08/08 12:42:00 mmann Exp $
#ifndef BIU_QUASIRANDOMNUMBERGENERATOR_HH_
#define BIU_QUASIRANDOMNUMBERGENERATOR_HH_

#include "biu/qrng/gsl_qrng.h"

namespace biu {

/*!
 * An object of this class works as a quasi random number generator. It is used
 * to generate a low-discrepancy sequence with the property that any 
 * subsequence is almost uniformly distributed. This class is mainly for use
 * by quasi-Monte Carlo algorithms. This class is a wrapper for the 
 * sobol sequence generator from the GNU Scientific Library.
 * 
 * @author Daniel Maticzka
 */
class QuasiRandomNumberGenerator
{
private:
	//! a pointer to the generated gsl_qrng
	gsl_qrng * q;
	
public:
	QuasiRandomNumberGenerator();
	virtual ~QuasiRandomNumberGenerator();
	
	/*!
	 * Returns a new quasi random number.
	 * @return 	The next item of the quasirandom sequence. It's value is 
	 * 			in the range [0;1].
	 */
	double getQuasiRN();
};

} // end namespace biu

#endif /*QUASIRANDOMNUMBERGENERATOR_HH_*/
