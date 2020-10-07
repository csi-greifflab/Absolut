/*
 * Kabsch's algorithm: compute least-square-best transformation
 * Copyright (C) 2003, Arno Formella (formella@ei.uvigo.es)
 *
 * Source code is take from 
 * 
 * http://trevinca.ei.uvigo.es/~formella/inv/aaa/ca_en.html
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
 * The code is C++-robust.
 *
 * More information about GSL
 * ==========================
 *
 * The project homepage is http://www.gnu.org/software/gsl/
 * The development site is http://sources.redhat.com/gsl/
 * 
 */

#ifndef KABSCH_H
#define KABSCH_H

#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_matrix_double.h>

#if defined(__cplusplus)
extern "C" {
#endif

/*
   preconditions:
     size > 0
     X and Y point to (size x 3)-matrices specifying 3D points x_i and y_i
     U points to a (3 x 3)-matrix
     t points to a (3)-vector
     s points to double
   postconditions if return 1:
     X and Y will be centralized at their respective origins
     U will hold the rotation
     t will hold the translation
     s will hold the scaling value
     such that:
       sum_i (U * s * x_i + t - y_i)^2  is minimal
   postconditions if return 0:
     X and Y will be centralized at their respective origins
     U will hold the identity
     t will hold the difference vector of the centroids t = cy -cx
     s will hold 1.0
   note:
     if s == NULL, no scaling is computed
*/
int kabsch(
  unsigned int size,
  gsl_matrix *X,
  gsl_matrix *Y,
  gsl_matrix *U,
  gsl_vector *t,
  double *s
);

#if defined(__cplusplus)
}
#endif

#endif
