
#include "biu/util/biu-string.h"
#include <iostream>

namespace biu {
  namespace util {
	
	//-------------------------------------------------------------------------------
	// search in a target all occurrences of query
	// @author Michael Hiller, Sebastian Will (leichte Änderungen)
	//-------------------------------------------------------------------------------
	bool findQueryString(std::string target, std::string query,std::ostream & out) {
	   std::string::size_type loc=0;
	   bool found=false;
	   
	   while ( (loc=target.find( query, loc )) != std::string::npos ) {
	       out <<"query: " <<query <<" found at " <<loc <<std::endl;
	       found=true;
	       loc++;
	   }
	   
	   if (!found) {
	       out <<"query: " <<query <<" not found" <<std::endl;
	   }
	   
	   return found;
	}

  } // namespace util	
}  // namespace biu
