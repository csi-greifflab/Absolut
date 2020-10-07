#ifndef BIU_ROTATOR3D_HH_
#define BIU_ROTATOR3D_HH_

#include "biu/Point.hh"
#include "biu/Matrix.hh"

namespace biu
{


	/*!
	 * Handles the rotation of 3D points in 3 dimensions around given angles.
	 * 
	 * @author Martin Mann http://www.bioinf.uni-freiburg.de/~mmann/
	 * 
	 */
	class Rotator3D
	{
	public:
		
		 //! the approximate value of the constant PI used
		static const double PI;
		
		 //! typedef for matrix of double values
		typedef biu::Matrix<double> DblMatrix;

	protected:
		
		 //! value of (2*PI)
		static const double TWO_PI;
		
		 //! rotation angle around X axis in radian
		double rotX;
		 //! rotation angle around Y axis in radian
		double rotY;
		 //! rotation angle around Z axis in radian
		double rotZ;
		
		 //! the rotation matrix used to rotate around X axis for rotX in radian 
		DblMatrix mX;
		 //! the rotation matrix used to rotate around Y axis for rotX in radian 
		DblMatrix mY;
		 //! the rotation matrix used to rotate around Z axis for rotX in radian 
		DblMatrix mZ;
		
	public:

		/*! Construction of a new rotation handler that does no rotation at all.
		 */
		Rotator3D();
		
		/*! Construction of a new rotation handler using the given rotation 
		 * angles.
		 * @param rotX the rotation angle around X axis in radian 
		 * @param rotY the rotation angle around Y axis in radian 
		 * @param rotZ the rotation angle around Z axis in radian 
		 */
		Rotator3D(const double rotX, const double rotY, const double rotZ);
		
		/*! Destruction of the rotator
		 */
		virtual ~Rotator3D();

		
		//////////// ROTATION HANDLING ////////////////////////////////////
		
		
		/*! Performs the rotation of a given point around X, Y, and Z axis for 
		 * the current angles set.
		 * @param toRotate the point to rotate
		 * @return the rotated point
		 */
		biu::DblPoint
		rotate( const biu::DblPoint& toRotate ) const;
		
		
		//////////// SETTER ////////////////////////////////////
		
		/*!
		 * Sets the rotation angles to be applied for coming rotations.
		 * @param rotX the rotation angle around X axis in radian 
		 * @param rotY the rotation angle around Y axis in radian 
		 * @param rotZ the rotation angle around Z axis in radian 
		 */
		virtual void
		setRotation( const double rotX, const double rotY, const double rotZ );
		
		/*! Sets the rotation angle around X axis
		 * @param val the new rotation angle in radian
		 */
		void 
		setRotationX( const double val );
		
		/*! Sets the rotation angle around Y axis
		 * @param val the new rotation angle in radian
		 */
		void 
		setRotationY( const double val );

		/*! Sets the rotation angle around Z axis
		 * @param val the new rotation angle in radian
		 */
		void 
		setRotationZ( const double val );
		
		
		//////////// GETTER ////////////////////////////////////
		
		
		/*! Access to the currently used rotation angle around X axis
		 * @return the current rotation angle in radian
		 */
		const double
		getRotationX( void ) const;

		/*! Access to the currently used rotation angle around Y axis
		 * @return the current rotation angle in radian
		 */
		const double
		getRotationY( void ) const;
		
		/*! Access to the currently used rotation angle around Z axis
		 * @return the current rotation angle in radian
		 */
		const double
		getRotationZ( void ) const;

		/*! Access to the currently used rotation angle around X axis
		 * @return the current rotation angle in degree
		 */
		const double
		getRotationDegreeX( void ) const;
		
		/*! Access to the currently used rotation angle around Y axis
		 * @return the current rotation angle in degree
		 */
		const double
		getRotationDegreeY( void ) const;
		
		/*! Access to the currently used rotation angle around Z axis
		 * @return the current rotation angle in degree
		 */
		const double
		getRotationDegreeZ( void ) const;
		
	};

}

#endif /*ROTATOR3D_HH_*/
