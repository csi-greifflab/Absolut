/*
 *  Main authors:
 *     Martin Mann <mmann@informatik.uni-freiburg.de>
 *
 *  Copyright:
 *     Martin Mann, 2008
 *
 *
 *  See the file "LICENSE" for information on usage and
 *  redistribution of this file, and for a
 *     DISCLAIMER OF ALL WARRANTIES.
 *
 */

#include "Exception.hh"

	unsigned int Exception::errorCount = 0;


	Exception::Exception(): 
		outString("\nException raised\n"), retValue(-1)  
	{ 
		errorCount++; 
	}
	
	Exception::Exception(const std::string &outString) : 
		outString(outString), retValue(-1) 
	{ 
		errorCount++; 
	}
	
	Exception::Exception(const std::string &outString, const int retValue) : 
		outString(outString), retValue(retValue) 
	{ 
		errorCount++; 
	}
	
	Exception::~Exception() {
	}

	std::ostream& operator<< (std::ostream& os, const Exception& p) {
		if (&p == NULL)
			os <<"NULL";
		else
			os <<p.toString();
		return os;
	}
