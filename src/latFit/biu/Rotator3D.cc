#include "Rotator3D.hh"

#include "biu/assertbiu.hh"

#include <cmath>

namespace biu
{

	////////////////////////////////////////////////////////////////////////////
	

	const double Rotator3D::PI = 3.14159265;

	const double Rotator3D::TWO_PI = 2*Rotator3D::PI;

	
	////////////////////////////////////////////////////////////////////////////
	
	Rotator3D::
	Rotator3D()
	 :	rotX(0.0)
	 	, rotY(0.0)
	 	, rotZ(0.0)
	 	, mX(3,3,0.0)
	 	, mY(3,3,0.0)
	 	, mZ(3,3,0.0)
	{
		 // initialize rotation matrices
		setRotationX(rotX);
		setRotationY(rotY);
		setRotationZ(rotZ);
	}

	
	////////////////////////////////////////////////////////////////////////////
	
	Rotator3D::
	Rotator3D( const double rotX_, const double rotY_, const double rotZ_ )
	 :	rotX(0.0)
	 	, rotY(0.0)
	 	, rotZ(0.0)
	 	, mX(3,3,0.0)
	 	, mY(3,3,0.0)
	 	, mZ(3,3,0.0)
	{
		 // initialize rotation matrices
		setRotationX(rotX_);
		setRotationY(rotY_);
		setRotationZ(rotZ_);
	}
	
	
	////////////////////////////////////////////////////////////////////////////
	
	Rotator3D::~Rotator3D()
	{
	}

	
	////////////////////////////////////////////////////////////////////////////
	
	void
	Rotator3D::
	setRotation( const double rotX_, const double rotY_, const double rotZ_ ) {
		 // update data structures
		setRotationX(rotX_);
		setRotationY(rotY_);
		setRotationZ(rotZ_);
	}
	
	////////////////////////////////////////////////////////////////////////////
	

	// setter for rotX
	void
	Rotator3D::
	setRotationX( const double val ) 
	{ 
		assertbiu( std::abs(rotX) <= TWO_PI, "angle for X rotation should be within [-2*PI,2*PI]");
		  // update rotation angle in radian
		rotX = val;
		  // update rotation matrix 
		mX[0][0] = 1;
		mX[1][1] = cos(rotX);
		mX[1][2] = -sin(rotX);
		mX[2][1] = sin(rotX);
		mX[2][2] = cos(rotX);

	}
	
	////////////////////////////////////////////////////////////////////////////
	
	
	// setter for rotY
	void 
	Rotator3D::
	setRotationY( const double val ) 
	{ 
		assertbiu( std::abs(rotY) <= TWO_PI, "angle for Y rotation should be within [-2*PI,2*PI]");
		  // update rotation angle in radian
		rotY = val;
		  // update rotation matrix 
		mY[0][0] = cos(rotY);
		mY[0][2] = sin(rotY);
		mY[1][1] = 1;
		mY[2][0] = -sin(rotY);
		mY[2][2] = cos(rotY);
	}
	
	
	////////////////////////////////////////////////////////////////////////////
	
	// setter for rotZ
	void 
	Rotator3D::
	setRotationZ( const double val ) 
	{ 
		assertbiu( std::abs(rotZ) <= TWO_PI, "angle for Z rotation should be within [-2*PI,2*PI]");
		  // update rotation angle in radian
		rotZ = val;
		  // update rotation matrix 
		mZ[0][0] = cos(rotZ);
		mZ[0][1] = -sin(rotZ);
		mZ[1][0] = sin(rotZ);
		mZ[1][1] = cos(rotZ);
		mZ[2][2] = 1;
	}
	
	
	////////////////////////////////////////////////////////////////////////////
	
		// getter for rotX
	const double
	Rotator3D::
	getRotationX( void ) const { return rotX; }
	
	////////////////////////////////////////////////////////////////////////////
	
		// getter for rotY
	const double
	Rotator3D::
	getRotationY( void ) const { return rotY; }
	
	////////////////////////////////////////////////////////////////////////////
	
	// getter for rotZ
	const double
	Rotator3D::
	getRotationZ( void ) const { return rotZ; }
	
	////////////////////////////////////////////////////////////////////////////

	biu::DblPoint
	Rotator3D::
	rotate( const biu::DblPoint& toRotate ) const
	{
		 // perform successive rotation around the three axes
		return (mZ * ( mY * ( mX * toRotate ) ) );
	}
	
	////////////////////////////////////////////////////////////////////////////
	
		// getter for rotX in degree
	const double
	Rotator3D::
	getRotationDegreeX( void ) const 
	{ 
		 // convert from radian to degree
		return (rotX / TWO_PI) * 360.0 ; 
	}
	
	////////////////////////////////////////////////////////////////////////////
	
		// getter for rotY in degree
	const double
	Rotator3D::
	getRotationDegreeY( void ) const 
	{ 
		 // convert from radian to degree
		return (rotY / TWO_PI) * 360.0 ; 
	}
	
	////////////////////////////////////////////////////////////////////////////
	
		// getter for rotZ in degree
	const double
	Rotator3D::
	getRotationDegreeZ( void ) const 
	{ 
		 // convert from radian to degree
		return (rotZ / TWO_PI) * 360.0 ; 
	}
	
	////////////////////////////////////////////////////////////////////////////

		
} // namespace biu
