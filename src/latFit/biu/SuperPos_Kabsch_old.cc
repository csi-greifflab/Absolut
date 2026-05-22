
/*
 * Kabsch's algorithm: compute least-square-best-fit transformation
 * Copyright (C) 2003, Arno Formella (formella@ei.uvigo.es)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * This is a gsl-based implementation of Kabsch's algorithm presented in:
 *
 *   W. Kabsch, A solution for the best rotation to relate two sets of vectors,
 *   Acta Cryst. (1976), A32, 922-923
 *
 *   W. Kabsch, A discussion of the solution for the best rotation to relate
 *   two sets of vectors, Acta Cryst. (1978), A34, 827-828
 *
 * The code is C and C++ compilable.
 *
 * More information about GSL:
 *
 * The project homepage is http://www.gnu.org/software/gsl/
 * The development site is http://sources.redhat.com/gsl/
 *
 * bugfix from: Adrien Brilhault
 */

#include <stdio.h>

#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

/* gsl does not provide it */
static inline void kabsch_gsl_vector_cross(
		const gsl_vector *a,
		const gsl_vector *b,
		gsl_vector *c ) 
{
	double a0=gsl_vector_get(a,0);
	double a1=gsl_vector_get(a,1);
	double a2=gsl_vector_get(a,2);
	double b0=gsl_vector_get(b,0);
	double b1=gsl_vector_get(b,1);
	double b2=gsl_vector_get(b,2);
	gsl_vector_set(c,0,a1*b2-b1*a2);
	gsl_vector_set(c,1,a2*b0-b2*a0);
	gsl_vector_set(c,2,a0*b1-b0*a1);
}

#define NORM_EPS 0.00000001

int kabsch_superpositioning(
		unsigned int size, /* the number of points */
		gsl_matrix *X,     /* the points to be moved */
		gsl_matrix *Y,     /* the points to move to */
		gsl_matrix *U,     /* the rotation matrix */
		gsl_vector *t,     /* the translation vector */
		double *s          /* the optimal scaling, if != 0 */
) {
	unsigned int i,j,k;
	int U_ok=1;
	double n=1.0/size;
	gsl_vector *cx=gsl_vector_alloc(3);     /* centroid of X */
	gsl_vector *cy=gsl_vector_alloc(3);     /* centroid of Y */
	gsl_matrix *R=gsl_matrix_alloc(3,3);    /* Kabsch's R */
	gsl_matrix *RTR=gsl_matrix_alloc(3,3);  /* R_trans * R (and Kabsch's bk) */
	gsl_eigen_symmv_workspace *espace=gsl_eigen_symmv_alloc(3);
	gsl_matrix *evec=gsl_matrix_alloc(3,3); /* eigenvectors (and Kabsch's ak) */
	gsl_vector *eval=gsl_vector_alloc(3);   /* vector of eigenvalues */

	/* compute centroid of X */
	gsl_vector_set_zero(cx);
	for(i=size;i>0;) {
		gsl_vector_const_view row=gsl_matrix_const_row(X,--i);
		gsl_vector_add(cx,&row.vector);
	} 
	gsl_vector_scale(cx,n);

	/* compute centroid of Y */
	gsl_vector_set_zero(cy);
	for(i=size;i>0;) {
		gsl_vector_const_view row=gsl_matrix_const_row(Y,--i);
		gsl_vector_add(cy,&row.vector);
	} 
	gsl_vector_scale(cy,n);

	/* move X to origin */
	for(i=size;i>0;) {
		gsl_vector_view row=gsl_matrix_row(X,--i);
		gsl_vector_sub(&row.vector,cx);
	}
	/* move Y to origin */
	for(i=size;i>0;) {
		gsl_vector_view row=gsl_matrix_row(Y,--i);
		gsl_vector_sub(&row.vector,cy);
	}

	if(size==1) {
		/* just one point, so U is trival */
		gsl_matrix_set_identity(U);
	}
	else {
		/* compute R */
		gsl_matrix_set_zero(R);
		for(k=size;k>0;) {
			--k;
			for(i=3;i>0;) {
				--i;
				for(j=3;j>0;) {
					--j;
					gsl_matrix_set(R,i,j,
							gsl_matrix_get(R,i,j)+
							gsl_matrix_get(Y,k,i)*gsl_matrix_get(X,k,j)
					);
				}
			}
		}

		/* compute RTR = R_trans * R */
		gsl_matrix_set_zero(RTR);
		gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,R,R,0.0,RTR);

		/* compute orthonormal eigenvectors */
		gsl_eigen_symmv(RTR,eval,evec,espace);  /* RTR will be modified! */
		gsl_eigen_symmv_sort(eval,evec,GSL_EIGEN_SORT_VAL_DESC);
		if(gsl_vector_get(eval,1)>NORM_EPS) {
			/* compute ak's (as columns of evec) and bk's (as columns of RTR) */
			double norm_b0,norm_b1,norm_b2;
			gsl_vector_const_view a0=gsl_matrix_const_column(evec,0);
			gsl_vector_const_view a1=gsl_matrix_const_column(evec,1);
			gsl_vector_view a2=gsl_matrix_column(evec,2);
			gsl_vector_view b0=gsl_matrix_column(RTR,0);
			gsl_vector_view b1=gsl_matrix_column(RTR,1);
			gsl_vector_view b2=gsl_matrix_column(RTR,2);
			kabsch_gsl_vector_cross(&a0.vector,&a1.vector,&a2.vector); /* a2 = a0 x a1 */
			gsl_blas_dgemv(CblasNoTrans,1.0,R,&a0.vector,0.0,&b0.vector);
			norm_b0=gsl_blas_dnrm2(&b0.vector);
			gsl_blas_dgemv(CblasNoTrans,1.0,R,&a1.vector,0.0,&b1.vector);
			norm_b1=gsl_blas_dnrm2(&b1.vector);
			if(norm_b0>NORM_EPS&&norm_b1>NORM_EPS) {
				gsl_vector_scale(&b0.vector,1.0/norm_b0);         /* b0 = ||R * a0|| */
				gsl_vector_scale(&b1.vector,1.0/norm_b1);         /* b1 = ||R * a1|| */
				kabsch_gsl_vector_cross(&b0.vector,&b1.vector,&b2.vector);  /* b2 = b0 x b1 */

				norm_b2=gsl_blas_dnrm2(&b2.vector);
				if(norm_b2>NORM_EPS) {
					/* we reach this point only if all bk different from 0 */
					/* compute U = B * A_trans (use RTR as B and evec as A) */
					gsl_matrix_set_zero(U); /* to avoid nan */
					gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,RTR,evec,0.0,U);
				}
				else {
					U_ok=0;
					gsl_matrix_set_identity(U);
				}
			}
			else {
				U_ok=0;
				gsl_matrix_set_identity(U);
			}
		}
		else {
			U_ok=0;
			gsl_matrix_set_identity(U);
		}
	}

	if(s) {
		/* let us compute the optimal scaling as well */
		/* s = <Y,UX> / <UX,UX> */
		*s=1.0;
		if(U_ok&&size>1) {
			double dom=0.0;
			double nom=0.0;
			double dom_i,nom_i;
			gsl_vector *Uxi=gsl_vector_alloc(3);
			for(i=size;i>0;) {
				gsl_vector_const_view row_x=gsl_matrix_const_row(X,--i);
				gsl_vector_const_view row_y=gsl_matrix_const_row(Y,i);
				gsl_vector_set_zero(Uxi);
				gsl_blas_dgemv(CblasNoTrans,1.0,U,&row_x.vector,1.0,Uxi);
				gsl_blas_ddot(&row_y.vector,Uxi,&nom_i);
				nom+=nom_i;
				gsl_blas_ddot(Uxi,Uxi,&dom_i);
				dom+=dom_i;
			}
			*s=nom/dom;
			gsl_vector_free(Uxi);
		}
		gsl_vector_scale(cx,*s);
	}
	/* compute t = cy - s * U * cx  */
	gsl_vector_memcpy(t,cy);
	gsl_blas_dgemv(CblasNoTrans,-1.0,U,cx,1.0,t);

	gsl_vector_free(eval);
	gsl_matrix_free(evec);
	gsl_eigen_symmv_free(espace);
	gsl_matrix_free(RTR);
	gsl_matrix_free(R);
	gsl_vector_free(cy);
	gsl_vector_free(cx);

	return U_ok;
}

#include "SuperPos_Kabsch.hh"


#include <biu/assertbiu.hh>
#include <biu/LatticeProteinUtil.hh>

#include <cmath>

std::vector<double>
p2v( const biu::DblPoint& p) {
	std::vector<double> v(3);
	v[0] = p.getX();
	v[1] = p.getY();
	v[2] = p.getZ();
	return v;
}

biu::DblPoint
v2p( const std::vector<double> & v) {
	return biu::DblPoint(v.at(0),v.at(1),v.at(2));
}


namespace biu {

	
	biu::DblPoint
	SuperPos_Kabsch::
	rotate(	const biu::SquareMatrix<double,3>& U
			, const biu::DblPoint & p )
	{
		return v2p(U * p2v(p));
	}
	
	biu::DblPoint
	SuperPos_Kabsch::
	applyAutomorphism(	const biu::Automorphism& am
						, const biu::DblPoint & p )
	{
		return biu::DblPoint(
					am[0][0]*p.getX() + am[0][1]*p.getY() + am[0][2]*p.getZ(),
					am[1][0]*p.getX() + am[1][1]*p.getY() + am[1][2]*p.getZ(),
					am[2][0]*p.getX() + am[2][1]*p.getY() + am[2][2]*p.getZ()
				);
	}

	
	SuperPos_Kabsch::
	Transformation
	SuperPos_Kabsch::
	superposition(	biu::DPointVec & pos1
					, const biu::DPointVec & pos2
					, const bool scale )
	{
		assertbiu(pos1.size() == pos2.size(), "pos1 and pos2 are of different size");
		
		  // simple checks
		if (pos1.size() == 0) {
			return Transformation();
		}
		if (pos1.size() == 1) {
			pos1[0] = pos2[0];
			return Transformation();
		}
		
		gsl_matrix *p1 = gsl_matrix_alloc(pos1.size(), 3);
		gsl_matrix *p2 = gsl_matrix_alloc(pos1.size(), 3);
		
		for (size_t i = 0; i<pos1.size(); i++) {
			  // copy pos1
			gsl_matrix_set( p1, i, 0, pos1.at(i).getX());
			gsl_matrix_set( p1, i, 1, pos1.at(i).getY());
			gsl_matrix_set( p1, i, 2, pos1.at(i).getZ());
			  // copy pos2
			gsl_matrix_set( p2, i, 0, pos2.at(i).getX());
			gsl_matrix_set( p2, i, 1, pos2.at(i).getY());
			gsl_matrix_set( p2, i, 2, pos2.at(i).getZ());
		}
		
		 // rotation matrix to fill
		gsl_matrix *U = gsl_matrix_alloc(3, 3);
		 // translation vector to fill
		gsl_vector *t = gsl_vector_alloc(3);
		 // scaling 
		double s = 1.0;
		
		  // run the Kabsch algorithm
		kabsch_superpositioning( pos1.size(), p1, p2, U, t, (scale ? &s : NULL) );
		
		  // prepare the return container for the transformation information
		Transformation trans;
		  // set translation information
		trans.translation.setX( gsl_vector_get( t, 0 ) );
		trans.translation.setY( gsl_vector_get( t, 1 ) );
		trans.translation.setZ( gsl_vector_get( t, 2 ) );
		  // set rotation matrix
		for (size_t i=0; i<3; i++) {
			trans.rotation[i][0] = gsl_matrix_get( U, i, 0);
			trans.rotation[i][1] = gsl_matrix_get( U, i, 1);
			trans.rotation[i][2] = gsl_matrix_get( U, i, 2);
		}
		using namespace biu;
		  // set scaling
		trans.scale = s;
		
		  // copy data back
		for (size_t i = 0; i<pos1.size(); i++) {
			  // copy pos1 (U*s*p1)
			pos1[i] = trans.translation + rotate( trans.rotation, pos1[i])*trans.scale;
		}

		  // deallocation
		gsl_matrix_free( p1 );
		gsl_matrix_free( p2 );
		gsl_matrix_free( U );
		gsl_vector_free( t );

		return trans;
	}
	
	SuperPos_Kabsch::
	Transformation
	SuperPos_Kabsch::
	superposition(	biu::IPointVec & pos1
					, const biu::IPointVec & pos2 )
	{
		assertbiu(pos1.size() == pos2.size(), "pos1 and pos2 are of different size");

		  // copy data for superpositioning
		biu::DPointVec p1(pos1.size());
		biu::DPointVec p2(pos2.size());
		for (size_t i=0; i<pos1.size(); i++) {
			p1[i] = biu::DblPoint(pos1[i]);
			p2[i] = biu::DblPoint(pos2[i]);
		}
		
		  // call superpositioning and store transformation information
		Transformation trans = SuperPos_Kabsch::superposition( p1, p2 );
		
		  // transform results into integer coordinates
		for (size_t i=0; i<pos1.size(); i++) {
			  // pos1 update
			pos1[i].setX( int(round(p1.at(i).getX())) );
			pos1[i].setY( int(round(p1.at(i).getY())) );
			pos1[i].setZ( int(round(p1.at(i).getZ())) );
		}
	
		  // return transformation information
		return trans;
	}


	double
	SuperPos_Kabsch::
	bestsuperposition(	biu::DPointVec & pos1
					, const biu::DPointVec & pos2
					, const biu::AutomorphismVec & automorphisms )
	{
		assertbiu(pos1.size() == pos2.size(), "pos1 and pos2 are of different size");
		
		  // simple checks
		if (pos1.size() == 0) {
			return 0.0;
		}
		if (pos1.size() == 1) {
			pos1[0] = pos2[0];
			return 0.0;
		}
		
		  // holds the best superpositioning found
		biu::DPointVec bestPos1(pos1);		
		double bestCRMSD = LatticeProteinUtil::cRMSD(pos1,pos2);

		  // holds the superpositioning of the current automorphism applied
		biu::DPointVec curPos1(pos1.size());
		biu::DPointVec curPos2(pos2.size());
		double curCRMSD = bestCRMSD;
		
		for (size_t a=0; a<automorphisms.size(); a++) {
			  // apply automorphism and generate points for superpositioning
			for (size_t i=0; i<pos1.size(); i++) {
				curPos1[i] = SuperPos_Kabsch::applyAutomorphism(automorphisms.at(a),pos1.at(i));
				curPos2[i] = pos2.at(i);
			}
			  // find superpositioning
			SuperPos_Kabsch::superposition(curPos1,curPos2,false);
			  // evaluate current superpositioning
			curCRMSD = LatticeProteinUtil::cRMSD(curPos1,curPos2);
			  // update best if necessary
			if (curCRMSD < bestCRMSD) {
				bestPos1 = curPos1;
				bestCRMSD = curCRMSD;
			}
		}
		
		  // copy best superpositioning back to output parameters
		pos1 = bestPos1;
		return bestCRMSD;
	}		

}

#include <iostream> // for test routine
#include <iomanip> // for test routine

void
print (const biu::SquareMatrix<double,3>& m)
{
	for (size_t i = 0; i < m.numRows(); i++)
	{
		for (size_t j = 0; j < m.numColumns(); j++)
		{
			std::cout << m[i][j] << " ";
		}
		std::cout << "\n";
	}
}

namespace biu {

	void
	SuperPos_Kabsch::
	test(void)
	{
		SuperPos_Kabsch::Transformation trans;
		
		biu::IPointVec ti1;
		ti1.push_back( biu::IntPoint( 0, 0, 0) );
		ti1.push_back( biu::IntPoint( 1, 0, 0) );
		ti1.push_back( biu::IntPoint( 1, 1, 0) );
		ti1.push_back( biu::IntPoint( 2, 1, 0) );
		
		biu::IPointVec ti2;
		ti2.push_back( biu::IntPoint( 3, 3, 3) );
		ti2.push_back( biu::IntPoint( 3, 4, 3) );
		ti2.push_back( biu::IntPoint( 2, 4, 3) );
		ti2.push_back( biu::IntPoint( 2, 5, 3) );
		
		biu::DPointVec td1(ti1.size());
		biu::DPointVec td2(ti2.size());
		for (size_t i=0; i<ti1.size(); i++) {
			td1[i] = biu::DblPoint(ti1[i]);
			td2[i] = biu::DblPoint(ti2[i]);
		}
		
		biu::IPointVec ti3;
		ti3.push_back( biu::IntPoint( 0, 1, 0) );
		ti3.push_back( biu::IntPoint( 0, 0, 0) );
		ti3.push_back( biu::IntPoint( 1, 0, 0) );
		ti3.push_back( biu::IntPoint( 1, 1, 0) );
		ti3.push_back( biu::IntPoint( 2, 1, 0) );
		ti3.push_back( biu::IntPoint( 2, 2, 0) );
		
		biu::DPointVec td3(ti3.size());
		for (size_t i=0; i<ti3.size(); i++) {
			td3[i] = biu::DblPoint(ti3[i]);
		}
		
		biu::DPointVec d1, d2;
		biu::IPointVec i1, i2;
		
		  // run test
		std::cout <<std::fixed <<std::setprecision(3);

		std::cout <<" d1 = d2 " <<std::endl;
		d1 = td1;
		d2 = td1;
		for (size_t i=0; i<d1.size(); i++) {
			std::cout <<" " <<i <<" : (" <<d1.at(i) <<")  /  (" <<d2.at(i) <<")\n";
		}
		std::cout <<" cRMSD - initial    : " << LatticeProteinUtil::cRMSD( d1, d2 ) <<std::endl;
		std::cout <<" dRMSD - initial    : " << LatticeProteinUtil::dRMSD( d1, d2 ) <<std::endl;
		trans = SuperPos_Kabsch::superposition( d1, d2 );
		std::cout <<"  rotation =\n"; print(trans.rotation); 
		std::cout	<<"  translation = " <<trans.translation
				<<"\n  scaling     = " <<trans.scale
				<<std::endl;
		std::cout <<" cRMSD - afterwards : " << LatticeProteinUtil::cRMSD( d1, d2 ) <<std::endl;
		std::cout <<" dRMSD - afterwards : " << LatticeProteinUtil::dRMSD( d1, d2 ) <<std::endl;
		for (size_t i=0; i<d1.size(); i++) {
			std::cout <<" " <<i <<" : (" <<d1.at(i) <<")  /  (" <<d2.at(i) <<")\n";
		}
		std::cout <<std::endl;
		std::cout <<" i1 = i2 " <<std::endl;
		i1 = ti1;
		i2 = ti1;
		for (size_t i=0; i<i1.size(); i++) {
			std::cout <<" " <<i <<" : (" <<i1.at(i) <<")  /  (" <<i2.at(i) <<")\n";
		}
		std::cout <<" cRMSD - initial    : " << LatticeProteinUtil::cRMSD( i1, i2 ) <<std::endl;
		std::cout <<" dRMSD - initial    : " << LatticeProteinUtil::dRMSD( i1, i2 ) <<std::endl;
		SuperPos_Kabsch::superposition( i1, i2 );
		std::cout <<" cRMSD - afterwards : " << LatticeProteinUtil::cRMSD( i1, i2 ) <<std::endl;
		std::cout <<" dRMSD - afterwards : " << LatticeProteinUtil::dRMSD( i1, i2 ) <<std::endl;
		for (size_t i=0; i<i1.size(); i++) {
			std::cout <<" " <<i <<" : (" <<i1.at(i) <<")  /  (" <<i2.at(i) <<")\n";
		}
		std::cout <<std::endl;

		  // run test
		std::cout <<" d1 , d2 " <<std::endl;
		d1 = td1;
		d2 = td2;
		for (size_t i=0; i<d1.size(); i++) {
			std::cout <<" " <<i <<" : (" <<d1.at(i) <<")  /  (" <<d2.at(i) <<")\n";
		}
		std::cout <<" cRMSD - initial    : " << LatticeProteinUtil::cRMSD( d1, d2 ) <<std::endl;
		trans = SuperPos_Kabsch::superposition( d1, d2 );
		std::cout <<"  rotation =\n"; print(trans.rotation); 
		std::cout	<<"  translation = " <<trans.translation
				<<"\n  scaling     = " <<trans.scale
				<<std::endl;
		std::cout <<" cRMSD - afterwards : " << LatticeProteinUtil::cRMSD( d1, d2 ) <<std::endl;
		for (size_t i=0; i<d1.size(); i++) {
			std::cout <<" " <<i <<" : (" <<d1.at(i) <<")  /  (" <<d2.at(i) <<")\n";
		}
		std::cout <<std::endl;
		std::cout <<" i1 , i2 " <<std::endl;
		i1 = ti1;
		i2 = ti2;
		for (size_t i=0; i<i1.size(); i++) {
			std::cout <<" " <<i <<" : (" <<i1.at(i) <<")  /  (" <<i2.at(i) <<")\n";
		}
		std::cout <<" cRMSD - initial    : " << LatticeProteinUtil::cRMSD( i1, i2 ) <<std::endl;
		SuperPos_Kabsch::superposition( i1, i2 );
		std::cout <<" cRMSD - afterwards : " << LatticeProteinUtil::cRMSD( i1, i2 ) <<std::endl;
		for (size_t i=0; i<i1.size(); i++) {
			std::cout <<" " <<i <<" : (" <<i1.at(i) <<")  /  (" <<i2.at(i) <<")\n";
		}
		std::cout <<std::endl;
		
		  // run test
		std::cout <<" d2 , d1 " <<std::endl;
		d1 = td2;
		d2 = td1;
		for (size_t i=0; i<d1.size(); i++) {
			std::cout <<" " <<i <<" : (" <<d1.at(i) <<")  /  (" <<d2.at(i) <<")\n";
		}
		std::cout <<" cRMSD - initial    : " << LatticeProteinUtil::cRMSD( d1, d2 ) <<std::endl;
		trans = SuperPos_Kabsch::superposition( d1, d2 );
		std::cout <<"  rotation =\n"; print(trans.rotation); 
		std::cout	<<"  translation = " <<trans.translation
				<<"\n  scaling     = " <<trans.scale
				<<std::endl;
		std::cout <<" cRMSD - afterwards : " << LatticeProteinUtil::cRMSD( d1, d2 ) <<std::endl;
		for (size_t i=0; i<d1.size(); i++) {
			std::cout <<" " <<i <<" : (" <<d1.at(i) <<")  /  (" <<d2.at(i) <<")\n";
		}
		std::cout <<std::endl;
		std::cout <<" i2 , i1 " <<std::endl;
		i1 = ti2;
		i2 = ti1;
		for (size_t i=0; i<i1.size(); i++) {
			std::cout <<" " <<i <<" : (" <<i1.at(i) <<")  /  (" <<i2.at(i) <<")\n";
		}
		std::cout <<" cRMSD - initial    : " << LatticeProteinUtil::cRMSD( i1, i2 ) <<std::endl;
		SuperPos_Kabsch::superposition( i1, i2 );
		std::cout <<" cRMSD - afterwards : " << LatticeProteinUtil::cRMSD( i1, i2 ) <<std::endl;
		for (size_t i=0; i<i1.size(); i++) {
			std::cout <<" " <<i <<" : (" <<i1.at(i) <<")  /  (" <<i2.at(i) <<")\n";
		}
		std::cout <<std::endl;
		
		  // run test
		std::cout <<" d1 , d1+(1,1,1) " <<std::endl;
		d1 = td1;
		d2 = td1;
		for (size_t i=0; i<d2.size(); i++) {
			d2[i] += biu::DblPoint(1,1,1);
		}
		for (size_t i=0; i<d1.size(); i++) {
			std::cout <<" " <<i <<" : (" <<d1.at(i) <<")  /  (" <<d2.at(i) <<")\n";
		}
		std::cout <<" cRMSD - initial     : " << LatticeProteinUtil::cRMSD( d1, d2 ) <<std::endl;
		std::cout <<" GDT_TS - initial    : " << LatticeProteinUtil::GDT_TS( d1, d2 ) <<std::endl;
		std::cout <<" GDT_HA - initial    : " << LatticeProteinUtil::GDT_HA( d1, d2 ) <<std::endl;
		trans = SuperPos_Kabsch::superposition( d1, d2 );
		std::cout <<"  rotation =\n"; print(trans.rotation); 
		std::cout	<<"  translation = " <<trans.translation
				<<"\n  scaling     = " <<trans.scale
				<<std::endl;
		std::cout <<" cRMSD - afterwards  : " << LatticeProteinUtil::cRMSD( d1, d2 ) <<std::endl;
		std::cout <<" GDT_TS - afterwards : " << LatticeProteinUtil::GDT_TS( d1, d2 ) <<std::endl;
		std::cout <<" GDT_HA - afterwards : " << LatticeProteinUtil::GDT_HA( d1, d2 ) <<std::endl;
		for (size_t i=0; i<d1.size(); i++) {
			std::cout <<" " <<i <<" : (" <<d1.at(i) <<")  /  (" <<d2.at(i) <<")\n";
		}
		std::cout <<std::endl;
		std::cout <<" i1 , i1+(1,1,1) " <<std::endl;
		i1 = ti1;
		i2 = ti1;
		for (size_t i=0; i<i2.size(); i++) {
			i2[i] += biu::IntPoint(1,1,1);
		}
		for (size_t i=0; i<i1.size(); i++) {
			std::cout <<" " <<i <<" : (" <<i1.at(i) <<")  /  (" <<i2.at(i) <<")\n";
		}
		std::cout <<" cRMSD - initial    : " << LatticeProteinUtil::cRMSD( i1, i2 ) <<std::endl;
		SuperPos_Kabsch::superposition( i1, i2 );
		std::cout <<" cRMSD - afterwards : " << LatticeProteinUtil::cRMSD( i1, i2 ) <<std::endl;
		for (size_t i=0; i<i1.size(); i++) {
			std::cout <<" " <<i <<" : (" <<i1.at(i) <<")  /  (" <<i2.at(i) <<")\n";
		}
		std::cout <<std::endl;
		
		  // run test
		std::cout <<" d1 , d2+(1,1,1) " <<std::endl;
		d1 = td1;
		d2 = td2;
		for (size_t i=0; i<d2.size(); i++) {
			d2[i] += biu::DblPoint(1,1,1);
		}
		for (size_t i=0; i<d1.size(); i++) {
			std::cout <<" " <<i <<" : (" <<d1.at(i) <<")  /  (" <<d2.at(i) <<")\n";
		}
		std::cout <<" cRMSD - initial    : " << LatticeProteinUtil::cRMSD( d1, d2 ) <<std::endl;
		std::cout <<" dRMSD - initial    : " << LatticeProteinUtil::dRMSD( d1, d2 ) <<std::endl;
		std::cout <<" GDT_TS - initial    : " << LatticeProteinUtil::GDT_TS( d1, d2 ) <<std::endl;
		std::cout <<" GDT_HA - initial    : " << LatticeProteinUtil::GDT_HA( d1, d2 ) <<std::endl;
		SuperPos_Kabsch::superposition( d1, d2 );
		std::cout <<" cRMSD - afterwards : " << LatticeProteinUtil::cRMSD( d1, d2 ) <<std::endl;
		std::cout <<" dRMSD - afterwards : " << LatticeProteinUtil::dRMSD( d1, d2 ) <<std::endl;
		std::cout <<" GDT_TS - afterwards : " << LatticeProteinUtil::GDT_TS( d1, d2 ) <<std::endl;
		std::cout <<" GDT_HA - afterwards : " << LatticeProteinUtil::GDT_HA( d1, d2 ) <<std::endl;
		for (size_t i=0; i<d1.size(); i++) {
			std::cout <<" " <<i <<" : (" <<d1.at(i) <<")  /  (" <<d2.at(i) <<")\n";
		}
		std::cout <<std::endl;
		std::cout <<" i1 , i2+(1,1,1) " <<std::endl;
		i1 = ti1;
		i2 = ti2;
		for (size_t i=0; i<i2.size(); i++) {
			i2[i] += biu::IntPoint(1,1,1);
		}
		for (size_t i=0; i<i1.size(); i++) {
			std::cout <<" " <<i <<" : (" <<i1.at(i) <<")  /  (" <<i2.at(i) <<")\n";
		}
		std::cout <<" cRMSD - initial    : " << LatticeProteinUtil::cRMSD( i1, i2 ) <<std::endl;
		SuperPos_Kabsch::superposition( i1, i2 );
		std::cout <<" cRMSD - afterwards : " << LatticeProteinUtil::cRMSD( i1, i2 ) <<std::endl;
		for (size_t i=0; i<i1.size(); i++) {
			std::cout <<" " <<i <<" : (" <<i1.at(i) <<")  /  (" <<i2.at(i) <<")\n";
		}
		std::cout <<std::endl;
		
		  // run test
		std::cout <<" d3 , d3(-X) " <<std::endl;
		d1 = td3;
		d2 = td3;
		for (size_t i=0; i<d2.size(); i++) {
			d2[i].setX( d2[i].getX() * -1.0);
		}
		for (size_t i=0; i<d1.size(); i++) {
			std::cout <<" " <<i <<" : (" <<d1.at(i) <<")  /  (" <<d2.at(i) <<")\n";
		}
		std::cout <<" cRMSD - initial    : " << LatticeProteinUtil::cRMSD( d1, d2 ) <<std::endl;
		std::cout <<" dRMSD - initial    : " << LatticeProteinUtil::dRMSD( d1, d2 ) <<std::endl;
		SuperPos_Kabsch::superposition( d1, d2 );
		std::cout <<" cRMSD - afterwards : " << LatticeProteinUtil::cRMSD( d1, d2 ) <<std::endl;
		std::cout <<" dRMSD - afterwards : " << LatticeProteinUtil::dRMSD( d1, d2 ) <<std::endl;
		for (size_t i=0; i<d1.size(); i++) {
			std::cout <<" " <<i <<" : (" <<d1.at(i) <<")  /  (" <<d2.at(i) <<")\n";
		}
		std::cout <<std::endl;
		
		
		
	}

} // namespace biu

