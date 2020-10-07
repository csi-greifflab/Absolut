/*
 *  Main authors:
 *     Martin Mann <mmann@informatik.uni-freiburg.de>
 *
 *  Copyright:
 *     Martin Mann, 2008
 *
 *  See the file "LICENSE" for information on usage and
 *  redistribution of this file, and for a
 *     DISCLAIMER OF ALL WARRANTIES.
 *
 */

#ifndef EXCEPTION_HH_
#define EXCEPTION_HH_

#include <string>
#include <iostream>

	class Exception
	{
	private:
		std::string outString;
	public:
		Exception();
		Exception(const std::string &outString);
		Exception(const std::string &outString, const int retValue);
		virtual ~Exception();

		std::string toString() const { return outString; }
	
		int retValue;
		static unsigned int errorCount;
	};
	
	std::ostream& operator<< (std::ostream& os, const Exception &p);
	

#endif /*EXCEPTION_HH_*/
