// $Id: QuasiRandomNumberGenerator.cc,v 1.2 2016/08/08 12:42:01 mmann Exp $

#include "biu/QuasiRandomNumberGenerator.hh"
#include "biu/qrng/gsl_qrng.h"

using namespace biu;

QuasiRandomNumberGenerator::QuasiRandomNumberGenerator()
{
	q = gsl_qrng_alloc (gsl_qrng_sobol, 1);
}

QuasiRandomNumberGenerator::~QuasiRandomNumberGenerator()
{
	gsl_qrng_free(q);
}

double 
QuasiRandomNumberGenerator::getQuasiRN()
{
	double v[1];
	gsl_qrng_get(q, v);
	return v[0];
}
