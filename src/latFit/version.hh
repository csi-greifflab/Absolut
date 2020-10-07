#ifndef VERSION_HH_IN_
#define VERSION_HH_IN_

#define BIN_PACKAGE_NAME "@PACKAGE_NAME@"
#define BIN_PACKAGE_VERSION "@PACKAGE_VERSION@"

#include <iostream>

void
giveVersion()
{
	std::cout	<<"\n " <<BIN_PACKAGE_NAME <<" package version " <<BIN_PACKAGE_VERSION 
				<<"\n"
				<<std::endl;
}

#endif /*VERSION_HH_IN_*/
