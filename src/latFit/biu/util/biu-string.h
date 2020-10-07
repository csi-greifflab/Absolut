#ifndef BIUSTRING_HH_
#define BIUSTRING_HH_

/* **********************************
 *  BioInformatics Utilities
 * **********************************
 *  String tool library
 * **********************************
 * 
 * 
 * *********************************/


#include <string>
#include <iostream>

namespace biu {
  namespace util {

	/************************************************************
	 * search in a target all occurrences of query and writes 
	 * them to the given out stream
	 ***********************************************************/
	bool findQueryString(std::string target, std::string query, 
		std::ostream & out = std::cout);	

  } // namespace util	
}  // namespace biu

#endif /*BIUSTRING_HH_*/
