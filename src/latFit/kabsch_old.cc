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

#include "kabsch.hh"

#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

/* gsl does not provide it */
static inline void gsl_vector_cross(
  const gsl_vector *a,
  const gsl_vector *b,
  gsl_vector *c
) {
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

int kabsch(
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
      gsl_vector_cross(&a0.vector,&a1.vector,&a2.vector); /* a2 = a0 x a1 */
      gsl_blas_dgemv(CblasNoTrans,1.0,R,&a0.vector,0.0,&b0.vector);
      norm_b0=gsl_blas_dnrm2(&b0.vector);
      gsl_blas_dgemv(CblasNoTrans,1.0,R,&a1.vector,0.0,&b1.vector);
      norm_b1=gsl_blas_dnrm2(&b1.vector);
      if(norm_b0>NORM_EPS&&norm_b1>NORM_EPS) {
        gsl_vector_scale(&b0.vector,1.0/norm_b0);         /* b0 = ||R * a0|| */
        gsl_vector_scale(&b1.vector,1.0/norm_b1);         /* b1 = ||R * a1|| */
        gsl_vector_cross(&b0.vector,&b1.vector,&b2.vector);  /* b2 = b0 x b1 */

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
