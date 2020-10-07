// $Id: Point.hh,v 1.2 2016/08/08 12:42:00 mmann Exp $
#ifndef BIU_POINT_HH_
#define BIU_POINT_HH_


#include <cmath>
#include <set>
#include <vector>
#include <iostream>
#include "biu/assertbiu.hh"

namespace biu
{

	template<class T>	
	class Point3D;

	///////////////////
	// TYPEDEF
	///////////////////

		//!An IntPoint represents an integer position in 3-dimensional space.
	typedef Point3D<int> IntPoint;
	typedef std::set<IntPoint> IPointSet;	//!< a std::set of IntPoint objects
	typedef std::vector<IntPoint> IPointVec;	/*!< a std::vector of IntPoint 
												 *   objects */
	
		 //! A DblPoint represents a positions in 3-dimensional space.
	typedef Point3D<double> DblPoint;
	typedef std::set<DblPoint> DPointSet;	//!< a std::set of DblPoint objects
	typedef std::vector<DblPoint> DPointVec;	/*!< a std::vector if DblPoint 
												 *   objects */


		/**
		 * The template Point3D represents a point in 3-dimensional space 
		 * for different base types.
		 * 
		 * @author Martin Mann, Sebastian Will, Andreas Richter
		 */
	template<class T>
	class Point3D
	{
	private:
			/*! the coordinates */	
		T x, y, z;		
		
			// for comparison of 2 floating point coordinates
		bool equal(const T &a, const T &b) const {
			return fabs(a-b)<1.0e-10;
		}
			
			// specialisation for 2 int coordinates
		template<int> bool 
		equal(const int &a, const int &b) const {
			return a==b;
		}

	public:
			/*! Constructs a point placed in the lattice center (0,0,0). */
		Point3D() : x(T(0)), y(T(0)), z(T(0)) 
		{}
		
		Point3D(const T _x, const T _y, const T _z) : x(_x), y(_y), z(_z) 
		{}
		
		Point3D(const Point3D &p2): x(p2.x), y(p2.y), z(p2.z) 
		{}

		Point3D(const std::vector<T> &v): x(v[0]),y(v[1]),z(v[2]) {
			assertbiu( (v.size()==3),				"size of the initializing vector has to be 3");
		}

		~Point3D() {}
		
			// conversion
		operator std::vector<T> () const {
			std::vector<T> v(3);
			v[0]=x; v[1]=y; v[2]=z;
			return v;
		}

			// operators
		Point3D& operator+= (const Point3D &p2) {
			x+=p2.x; y+=p2.y; z+=p2.z;
    		return *this;
		}
		
		Point3D& operator+= (const T &d) {
			x+=d; y+=d; z+=d;
			return *this;
		}
		
		Point3D& operator-= (const Point3D &p2) {
        	x-=p2.x; y-=p2.y; z-=p2.z;
    		return *this;
		}
		
		Point3D& operator-= (const T &d) {
			x-=d; y-=d; z-=d;
			return *this;
		}
		
		Point3D& operator*= (const T &d) {
			x*=d; y*=d; z*=d;
			return *this;
		}
		
		Point3D& operator/= (const T &d) {
	        x/=d; y/=d; z/=d;
        	return *this;
		}
		
		Point3D operator+ (const Point3D &p2) const {
			return Point3D (x+p2.x, y+p2.y, z+p2.z);
		}
		
		Point3D operator- (const Point3D &p2) const {
			return Point3D (x-p2.x, y-p2.y, z-p2.z);
		}
		
		Point3D operator* (const T &d) const {
			return Point3D (d*x, d*y, d*z);
		}
		
		Point3D operator/ (const T &d) const {
			return Point3D (x/d, y/d, z/d);
		}
		
		Point3D& operator -() {
			x=-x; y=-y; z=-z;
			return *this;
		}
		
		Point3D& operator= (const Point3D &p2) {
			if (this != &p2) {
				x=p2.x; y=p2.y; z=p2.z;
			}
    		return *this;
		}
		
		bool	operator== (const Point3D &p2) const {
			return equal(x,p2.x) && equal(y,p2.y) && equal(z,p2.z);
		}
		
		bool	operator!= (const Point3D &p2) const {
			return !(equal(x,p2.x) && equal(y,p2.y) && equal(z,p2.z));
		}
		
		operator Point3D<double>() const {
			return Point3D<double>((double) x, (double) y, (double) z);
		}
		
			/*! Calculates the less operator in the order: x-component, 
			 * y-component, z-component */
		bool	operator< (const Point3D &p2) const {
			return (x<p2.x) ||	(equal(x,p2.x) 
								&& (y<p2.y || ( equal(y,p2.y) && z<p2.z)));
		}
		
			/*! Calculates the greater operator in the order: z-component, 
			 * y-component, x-component */
		bool	operator> (const Point3D &p2) const {
			return (z>p2.z) ||	(equal(z,p2.z) 
								&& (y>p2.y || ( equal(y,p2.y) && x>p2.x)));
		}
		
	    		/*! Calculates the euclidian distance of two Point3D objects*/
		double	distance (const Point3D &p2) const {
			double dist;
        		dist =	pow( (double)(x - p2.x) , 2) 
        				+ pow( (double)(y - p2.y) , 2) 
        				+ pow( (double) (z - p2.z) , 2);
        		return sqrt(dist);
		}
		
			/*! Returns if the sum of x-, y- and z-component is even or not. */
		bool	isEven() const {return (x+y+z)%2 == 0;}
		
			/*! Returns the length of the vector 
			 * (between lattice center and this point). */
		double	vectorLength() const {return sqrt(double(x*x + y*y + z*z));}
		
		friend std::ostream& operator<< (std::ostream &out, const Point3D& p) {
	                out << p.x << "  " << p.y << "  " << p.z;
                	return out;
        }

		T	getX() const { return x; }
		T	getY() const { return y; }
		T	getZ() const { return z; }
		
		void setX(const T& x_) { x = x_; }
		void setY(const T& y_) { y = y_; }
		void setZ(const T& z_) { z = z_; }

	};


} // namespace biu


#endif /*POINT_HH_*/
