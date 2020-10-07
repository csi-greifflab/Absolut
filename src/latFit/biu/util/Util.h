// Util.h: 
/*

	<author>	Martin Mann					</author>
	<created>	9.10.2005					</created>

	<info>									</info>
*/
//
//////////////////////////////////////////////////////////////////////

#if !defined(IT_UTIL_H__INCLUDED)
#define IT_UTIL_H__INCLUDED

#include <iostream>
#include <sstream>
#include <ostream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <functional>

#ifdef WIN32
	#pragma warning( disable : 4172 )
#endif

namespace biu {
  namespace util {

	//////////////////////////////////
	// TYPEDEFs
	//////////////////////////////////


	typedef std::vector<int>			vecInt;
	typedef std::set<int>				setInt;
	typedef std::map<int, int>			mapInt;
	typedef std::vector<vecInt>			vecVecInt;
	typedef std::vector<setInt>			vecSetInt;
	typedef std::map<int, vecInt>		mapIntVecInt;
	typedef std::vector<double>			vecDbl;
	typedef std::vector<std::string>	vecStr;
	typedef std::pair<int,int>			pairII;
	
	/////////////////////////////////
	// globale funktionen
	/////////////////////////////////
	
	std::ostream& operator<< (std::ostream& os ,const  mapInt& c);
	std::ostream& operator<< (std::ostream& os ,const  vecInt& c);
	std::ostream& operator<< (std::ostream& os ,const  setInt& c);
	
	std::ostream& operator<< (std::ostream& os ,const  vecSetInt& c);
	
	
	/////////////////////////////////
	//! print class to add simple output support to classes by implementing
	//! the toString() method.
	/////////////////////////////////
	class CPrintable {
	public:
		virtual ~CPrintable();
			//! string conversion method to implement by inherited subclasses.
		virtual std::string toString() const =0;
	};
	
	/////////////////////////////////
	// ausgabe fuer CPrintable und abgeleitete klassen
	/////////////////////////////////
	std::ostream& operator<< (std::ostream& os, const CPrintable &p);
	
	
	///////////////////////////////////////////
	//! exception class
	///////////////////////////////////////////
	class CException : public CPrintable {
	private:	
			//! error message
		std::string outString;		
	public:
		// construction
		CException() : outString("\nException raised\n"), retValue(-1)  { errorCount++; };
		CException(const std::string &outString) : outString(outString), retValue(-1) { errorCount++; }
		CException(const std::string &outString, const int retValue) : outString(outString), retValue(retValue) { errorCount++; }
		
			//! destruction
		virtual ~CException() {}
	
			//! returns error message
		std::string toString() const { return outString; }
	
			//! available return value
		int retValue;
			//! number of raised exceptions
		static int errorCount;
	};

  } // namespace util	
}  // namespace biu

#endif // !defined(IT_UTIL_H__INCLUDED)
