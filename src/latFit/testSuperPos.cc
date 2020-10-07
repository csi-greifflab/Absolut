
#include <stdio.h>

#include "biu/SuperPos_Kabsch.hh"

#include <cmath>
#include <iostream>
#include <iomanip>

#include "kabsch.hh"

using namespace biu;

void kabschTest() 
{
	
	unsigned int size = 4;
	
	gsl_matrix *X = gsl_matrix_alloc(size, 3);
	gsl_matrix *Y = gsl_matrix_alloc(size, 3);
	
	gsl_matrix_set( X, 0, 0, 0);
	gsl_matrix_set( X, 0, 1, 0);
	gsl_matrix_set( X, 0, 2, 0);
	
	gsl_matrix_set( X, 1, 0, 1);
	gsl_matrix_set( X, 1, 1, 0);
	gsl_matrix_set( X, 1, 2, 0);
	
	gsl_matrix_set( X, 2, 0, 1);
	gsl_matrix_set( X, 2, 1, 1);
	gsl_matrix_set( X, 2, 2, 0);
	
	gsl_matrix_set( X, 3, 0, 2);
	gsl_matrix_set( X, 3, 1, 1);
	gsl_matrix_set( X, 3, 2, 0);

	
	
	
	gsl_matrix_set( Y, 0, 0, 3);
	gsl_matrix_set( Y, 0, 1, 3);
	gsl_matrix_set( Y, 0, 2, 3);
	
	gsl_matrix_set( Y, 1, 0, 3);
	gsl_matrix_set( Y, 1, 1, 4);
	gsl_matrix_set( Y, 1, 2, 3);
	
	gsl_matrix_set( Y, 2, 0, 2);
	gsl_matrix_set( Y, 2, 1, 4);
	gsl_matrix_set( Y, 2, 2, 3);
	
	gsl_matrix_set( Y, 3, 0, 2);
	gsl_matrix_set( Y, 3, 1, 5);
	gsl_matrix_set( Y, 3, 2, 3);

//	gsl_matrix_memcpy(Y,X);
	
	{
		double diffSum = 0.0;
		std::cout <<"\n";
		for (size_t i=0; i<size; i++) {
			

			double x = gsl_matrix_get( X, i, 0) - gsl_matrix_get( Y, i, 0);
			double y = gsl_matrix_get( X, i, 1) - gsl_matrix_get( Y, i, 1);
			double z = gsl_matrix_get( X, i, 2) - gsl_matrix_get( Y, i, 2);
			
			
			std::cout <<" " <<i <<" :   " <<std::fixed <<std::setprecision(3)
				<<gsl_matrix_get( X, i, 0) <<"|" <<gsl_matrix_get( X, i, 1) <<"|" <<gsl_matrix_get( X, i, 2)
				<<"\t  #  "
				<<gsl_matrix_get( Y, i, 0) <<"|" <<gsl_matrix_get( Y, i, 1) <<"|" <<gsl_matrix_get( Y, i, 2)
				<<std::endl;

			double dist = sqrt( x*x + y*y + z*z );
			
			diffSum += pow( dist, 2 );
		}
		double cRMSD = sqrt (diffSum / double(size));
		
		std::cout <<"\n cRMSD before : " <<cRMSD <<std::endl;
	}
	
	 // rotation matrix
	gsl_matrix *U = gsl_matrix_alloc(3, 3);
	 // translation vector
	gsl_vector *t = gsl_vector_alloc(3);
	 // scaling
	double *s = NULL;
	
	std::cout <<"\n running kabsch : " <<kabsch( size, X, Y, U, t, s ) <<std::endl;
	
	
	{
		double diffSum = 0.0;
		std::cout <<"\n";
		for (size_t i=0; i<size; i++) {
			double xRot = 	  gsl_matrix_get( U, 0, 0) * gsl_matrix_get( X, i, 0)
							+ gsl_matrix_get( U, 0, 1) * gsl_matrix_get( X, i, 1)
							+ gsl_matrix_get( U, 0, 2) * gsl_matrix_get( X, i, 2);
			double x = xRot - gsl_matrix_get( Y, i, 0);
			double yRot = 	  gsl_matrix_get( U, 1, 0) * gsl_matrix_get( X, i, 0)
							+ gsl_matrix_get( U, 1, 1) * gsl_matrix_get( X, i, 1)
							+ gsl_matrix_get( U, 1, 2) * gsl_matrix_get( X, i, 2);
			double y = yRot - gsl_matrix_get( Y, i, 1);
			double zRot = 	  gsl_matrix_get( U, 2, 0) * gsl_matrix_get( X, i, 0)
							+ gsl_matrix_get( U, 2, 1) * gsl_matrix_get( X, i, 1)
							+ gsl_matrix_get( U, 2, 2) * gsl_matrix_get( X, i, 2);
			double z = zRot- gsl_matrix_get( Y, i, 2);
			
			std::cout <<" " <<i <<" :   " <<std::fixed <<std::setprecision(3)
				<<gsl_matrix_get( X, i, 0) <<"|" <<gsl_matrix_get( X, i, 1) <<"|" <<gsl_matrix_get( X, i, 2)
				<<"\t  #  "
				<<xRot <<"|" <<yRot <<"|" <<zRot
				<<"\t  #  "
				<<gsl_matrix_get( Y, i, 0) <<"|" <<gsl_matrix_get( Y, i, 1) <<"|" <<gsl_matrix_get( Y, i, 2)
				<<std::endl;

			
			double dist = sqrt( x*x + y*y + z*z );
			
			diffSum += pow( dist, 2 );
		}
		double cRMSD = sqrt (diffSum / double(size));
		
		std::cout <<"\n cRMSD after : " <<cRMSD <<std::endl;
	}
	
	std::cout <<"\n U:\n";
	for (size_t i=0; i<3; i++) {
		std::cout <<std::fixed <<std::setprecision(3)<<gsl_matrix_get(U,i,0) <<" " <<gsl_matrix_get(U,i,1) <<" " <<gsl_matrix_get(U,i,2) <<std::endl;
	}
	std::cout <<"\n t:\n";
	for (size_t i=0; i<3; i++) {
		std::cout <<gsl_vector_get(t,i) <<" ";
	}
	std::cout <<std::endl;
	
	gsl_matrix_free( X );
	gsl_matrix_free( Y );
	gsl_matrix_free( U );
	gsl_vector_free( t );
}


int main (int argc, char** argv ) {
	
	
//	kabschTest();
	
	std::cout <<" running biu::SuperPos_Kabsch::test();" <<std::endl;
	biu::SuperPos_Kabsch::test();
	
	return 0;
}
