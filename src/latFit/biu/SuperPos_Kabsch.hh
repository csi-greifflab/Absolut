#ifndef SUPERPOS_KABSCH_HH_
#define SUPERPOS_KABSCH_HH_

#include <utility>

#include <biu/Point.hh>
#include <biu/NeighborVector.hh>

namespace biu {
	
	/*! Allows for the superpositioning via the algorithm introduced by
	 * Kabsch. The class is a BIU wrapper around the implementation by Arno
	 * Formella (see below).
	 * 
	 * @author Martin Mann - http://www.bioinf.uni-freiburg.de
	 * 
	 *  
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
	class SuperPos_Kabsch {
		
	public:
		
		  /*! Defines the neccessary transformation done during the 
		   * superpositioning. It describes the translation
		   * vector, the rotation matrix, and the scaling
		   * needed.
		   */
		struct Transformation {
			  //! the translation vector applied
			biu::DblPoint translation;
			  //! the rotation matrix applied
			biu::SquareMatrix<double,3> rotation;
			  //! the vector length scaling applied 
			double scale;
		};
		
	protected:
	public:
		
		
		  /*! Applies the an automorphism (rotation/reflection) on a given point
		   * 
		   * @param am the automorphism to be applied
		   * @param p the point to rotate/reflect
		   * @return the resulting coordinate of p
		   */
		static
		biu::DblPoint
		applyAutomorphism(	const biu::Automorphism& am
							, const biu::DblPoint & p );
		
		
		
		  /*! Rotates a given point according to a rotation matrix.
		   * 
		   * @param U the rotation matrix to be applied
		   * @param p the point to rotate/reflect
		   * @return the resulting coordinate of p
		   */
		static
		biu::DblPoint
		rotate(	const biu::SquareMatrix<double,3>& U
							, const biu::DblPoint & p );
		
		  /*!
		   * Performs a superpositioning of the two point vectors utilizing the
		   * Kabsch algorithm. The provided point vectors are overwritten with 
		   * new coordinates.
		   * 
		   * Result: for all i in [0,pos1.size()]:
		   * 
		   * pos1[i] = T + (s * U * pos1[i])
		   * 
		   * with the returned Transformation values
		   * - T the translation vector
		   * - U the rotation matrix
		   * - s the scaling factor
		   * 
		   * @param pos1 INOUT first position data that is changed accordingly
		   * @param pos2 IN second position data (equal size to pos1)
		   * @param scale IN whether or not an appropriate scaling of the data 
		   *              should be done
		   * 
		   * @return the translation vector, rotation matrix and scaling factor
		   *   used to superposition pos1 to pos2
		   */
		static
		Transformation
		superposition(	biu::DPointVec & pos1
						, const biu::DPointVec & pos2
						, const bool scale = false );
		
		  /*!
		   * Performs a superpositioning of the two point vectors utilizing the
		   * Kabsch algorithm. The provided point vectors are overwritten with 
		   * new coordinates.
		   * 
		   * Result: for all i in [0,pos1.size()]:
		   * 
		   * pos1[i] = T + (s * U * pos1[i])
		   * 
		   * with the returned Transformation values
		   * - T the translation vector
		   * - U the rotation matrix
		   * - s the scaling factor
		   * 
		   * @param pos1 INOUT first position data that is changed accordingly
		   * @param pos2 IN second position data (equal size to pos1)
		   * 
		   * @return the translation vector, rotation matrix and scaling factor
		   *   used to superposition pos1 to pos2
		   */
		static
		Transformation
		superposition(	biu::IPointVec & pos1
						, const biu::IPointVec & pos2 );
		

		  /*!
		   * Calculates the best superposition of all automorphisms of pos1
		   * to pos2 according to a coordinate root mean square deviation 
		   * (cRMSD) evaluation.
		   * 
		   * @param pos1 first position data that is changed accordingly
		   * @param pos2 second position data (equal size to pos1)
		   * @param automorphisms the rotation and reflection matrices needed 
		   *         to calculate the automorphisms of pos1
		   * 
		   * @return the lowest cRMSD found that correspondes to the cRMSD of
		   *          the final pos1 and pos2 vectors.
		   */
		static
		double
		bestsuperposition(	biu::DPointVec & pos1
							, const biu::DPointVec & pos2
							, const biu::AutomorphismVec & automorphisms );
		
		
		  /*!
		   * Runs several test scenarios on the interface.
		   * 
		   * For internal use only ...
		   * 
		   */
		static
		void
		test( void );
	};


} // namespace biu


#endif /*SUPERPOS_KABSCH_HH_*/
