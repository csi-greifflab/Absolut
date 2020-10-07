// Util.cpp: Implementierung der Klasse CUtil.
/*

	<author>	Martin Mann					</author>
	<created>	9.10.2005					</created>

	<info>									</info>
*/
//
//////////////////////////////////////////////////////////////////////

#include "biu/util/Util.h"



namespace biu {
  namespace util {
	
	int CException::errorCount = 0;
	
	
	///////////////////////////////////////////////////////////////////
	// Ausgabeoperatoren der typen
	///////////////////////////////////////////////////////////////////
	
	std::ostream& operator<< (std::ostream& os ,const  mapInt& c) {
		os <<"{ ";
		for (mapInt::const_iterator i=c.begin(); i!=c.end(); i++)
			os <<(i!=c.begin()?", (":"(") <<i->first <<","<<i->second <<")";
		os <<"}";
		return os ;
	}
	
	
	std::ostream& operator<< (std::ostream& os ,const  setInt& c) {
		os <<"{ ";
		for (setInt::const_iterator i=c.begin(); i!=c.end(); i++)
			os  <<(i!=c.begin()?", ":"") <<(*i);
		os <<"}";
		return os ;
	}
	
	std::ostream& operator<< (std::ostream& os ,const  vecInt& c) {
		os <<"[ ";
		for (vecInt::const_iterator i=c.begin(); i!=c.end(); i++)
			os  <<(i!=c.begin()?", ":"") <<(*i);
		os <<"]";
		return os ;
	}
	
	std::ostream& operator<< (std::ostream& os ,const  vecSetInt& c) {
		os <<"[ ";
		for (vecSetInt::const_iterator i=c.begin(); i!=c.end(); i++)
			os  <<"\t"<<(*i) <<"\n";
		os <<"]";
		return os ;
	}
	
	
	std::ostream& operator<< (std::ostream& os, const CPrintable& p) {
		if (&p == NULL)
			os <<"NULL";
		else
			os <<p.toString();
		return os;
	}
		
	CPrintable::~CPrintable() {}

  } // namespace util
} // namespace biu
