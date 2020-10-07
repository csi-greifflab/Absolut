
#ifndef UTIL_POINT_HH_
#define UTIL_POINT_HH_

#include "biu/Point.hh"

	std::ostream& 
	operator<< (std::ostream &out, const biu::IPointSet& s);

	std::ostream& 
	operator<< (std::ostream &out, const biu::IPointVec& s);

	std::ostream& 
	operator<< (std::ostream &out, const biu::DPointSet& s);

	std::ostream& 
	operator<< (std::ostream &out, const biu::DPointVec& s);


namespace biu {
	
	
	
	/**
	 * Checks whether or not the P set sub is a subset of s.
	 */
	template<class T>
	bool
	isSubset(const std::set< Point3D<T> > &s, const std::set< Point3D<T> > &sub);
	
	
} // namespace biu


// implementations
#include "point.icc"

#endif /*UTIL_POINT_HH_*/
