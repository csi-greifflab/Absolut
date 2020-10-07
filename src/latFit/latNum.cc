
#include <iostream>
#include <iomanip>
#include <algorithm>

#include <biu/OptionParser.hh>
#include <biu/LatticeDescriptorSQR.hh>
#include <biu/LatticeDescriptorCUB.hh>
#include <biu/LatticeDescriptorFCC.hh>
#include <biu/LatticeDescriptorCKW.hh>
#include <biu/LatticeModel.hh>

#include <biu/LatticeProteinUtil.hh>

#include "version.hh"

/*!
 * Allows for the enumeration of all lattice protein structures.
 * 
 * Done via an exhaustive enumeration of all possible absolute move strings.
 * Each is checked if self-avoiding. Symmetric structures are identified by
 * the comparison of the move string to its normalized version: only if equal
 * the move string is non-symmetric and thus reported.
 *
 * @author (c) 2010 Martin Mann http://www.bioinf.uni-freiburg.de/~mmann/
 * 
 */

enum OUTMODE { OUT_SILENT, OUT_NORMAL, OUT_VERBOSE };

size_t extendAndCheck_BM( const size_t curPos
		, biu::IPointVec & P
		, biu::MoveSequence & M
		, const biu::LatticeDescriptor& latDescr
		, OUTMODE outMode )
{
	  // recursion abortion
	if (curPos == P.size()) {
		  // check if M equals its normalization, i.e. is a non-symmetric structure
		if ( std::equal( M.begin(), M.end(), latDescr.normalizeSequence( M ).begin()) ) {
			if (outMode != OUT_SILENT) {
				std::cout <<latDescr.getString( M ) <<"\n";
			}
			return 1;
		} else {
			return 0;
		}
	}

	  // neighborhood access
	const biu::LatticeNeighborhood & neigh = latDescr.getNeighborhood();

	  // special handling of first
	if (curPos == 1) {
		  // use first neighboring vector
		P[curPos] = P[curPos-1] + (neigh.getElementByIndex(0));
		M[curPos-1] = neigh.getElementByIndex(0).getMove();
		return extendAndCheck_BM( curPos + 1, P, M, latDescr, outMode );
	} else {

		size_t structureNumber = 0;
		for (biu::LatticeNeighborhood::const_iterator n = neigh.begin();
				n != neigh.end(); n++)
		{
			  // update position information
			P[curPos] = P[curPos-1] + (*n);
			  // check if selfavoiding
			bool isSelfAvoiding = true;
			for (size_t i = curPos-1; isSelfAvoiding && i > 0; i--) {
				isSelfAvoiding = (P[curPos] != P[i-1]);
			}
			if (isSelfAvoiding) {
				  // recursive call if selfavoiding
				M[curPos-1] = n->getMove();
				structureNumber += extendAndCheck_BM( curPos + 1, P, M, latDescr, outMode );
			}
		}
		return structureNumber;
	}
}

size_t extendAndCheck_SC( const size_t curPos
		, biu::IPointVec & P
		, biu::MoveSequence & M
		, const biu::LatticeDescriptor& latDescr
		, OUTMODE outMode  )
{
	  // recursion abortion
	if (curPos == P.size()) {
		  // check if M equals its normalization, i.e. is a non-symmetric structure
		if ( M == latDescr.normalizeSequence( M ) ) {
			if (outMode != OUT_SILENT) {
				std::cout <<latDescr.getString( M ) <<"\n";
			}
			return 1;
		} else {
			return 0;
		}
	}

	  // neighborhood access
	const biu::LatticeNeighborhood & neigh = latDescr.getNeighborhood();

	  // special handling of first
	if (curPos == 1) {
		  // use first neighboring vector
		P[curPos] = P[curPos-1] + (neigh.getElementByIndex(0));
		M[curPos-1] = neigh.getElementByIndex(0).getMove();
		return extendAndCheck_SC( curPos + 1, P, M, latDescr, outMode );
	} else {

		size_t structureNumber = 0;
		const size_t lastPos = (curPos%2 == 0) ? curPos-2 : curPos-1;
		for (biu::LatticeNeighborhood::const_iterator n = neigh.begin();
				n != neigh.end(); n++)
		{
			  // update position information
			P[curPos] = P[lastPos] + (*n);
			  // check if selfavoiding
			bool isSelfAvoiding = true;
			for (size_t i = curPos; isSelfAvoiding && i > 0; i--) {
				isSelfAvoiding = (P[curPos] != P[i-1]);
			}
			if (isSelfAvoiding) {
				  // recursive call if selfavoiding
				M[curPos-1] = n->getMove();
				structureNumber += extendAndCheck_SC( curPos + 1, P, M, latDescr, outMode );
			}
		}
		return structureNumber;
	}
}



//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// DEFINITIONS
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////

 /*! The parameter setup for this tools.
  * 
  * @param allowedArgs the provided parameters to fill
  * @param infoText the additional information text to setup 
  */
void
initAllowedArguments(biu::OptionMap & allowedArgs, std::string &infoText );

//////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// IMPLEMENTATIONS
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

int
main( int argc, char** argv ) 
{

	//////////////////////////////////////////////////////////////
	// parameter parsing and checking
	//////////////////////////////////////////////////////////////

	biu::OptionMap allowedArgs;
	std::string infoText;
	initAllowedArguments(allowedArgs,infoText);	// init

		// parse programm arguments
	biu::COptionParser opts = biu::COptionParser(	allowedArgs, argc,
													argv, infoText);
		// check if arguments parseable and all mandatory arguments given
	if (!opts.noErrors()) {
		return -1;
	}
		// help output
	if (opts.getBoolVal("help")) {
		opts.coutUsage();
		return 0;
	}
	  // version information requested
	if (opts.getBoolVal("version")) {
		giveVersion();
		return 0;
	}
	
	//////////////////////////////////////////////////////////////////////
	// setup parameters etc.
	//////////////////////////////////////////////////////////////////////
	
	
	OUTMODE outMode = OUT_NORMAL;
//	if (opts.argExist("s") && opts.argExist("v")) {
//		std::cerr <<"\n   ERROR : cannot be silent and verbose at the same time!\n";
//		return -1;
//	}
	if (opts.argExist("s")) {
		outMode = OUT_SILENT;
	}
//	if (opts.argExist("v")){
//		outMode = OUT_VERBOSE;
//	}
	
	biu::LatticeDescriptor* latDescr = NULL;
	std::string latStr = opts.getStrVal("lat");
	if (latStr.compare("SQR") == 0)
		latDescr = new biu::LatticeDescriptorSQR();
	else if (latStr.compare("CUB") == 0)
		latDescr = new biu::LatticeDescriptorCUB();
	else if (latStr.compare("FCC") == 0)
		latDescr = new biu::LatticeDescriptorFCC();
	else if (latStr.compare("210") == 0)
		latDescr = new biu::LatticeDescriptorCKW();
	else {
		std::cerr <<"\n   ERROR : Unknown lattice type '"+latStr+"'\n\n";
		return -1;
	}
	  // create lattice
	biu::LatticeModel lattice(latDescr);
	
	  // check if side chain structure
	const bool sidechain = opts.argExist("sideChain");

	  // get and convert structures
	if (opts.getIntVal("len") < 2) {
		std::cerr <<"\n   ERROR : structure length has to be >= 2\n\n";
		return -1;
	}
	const size_t length = (size_t)opts.getIntVal("len");


	//////////////////////////////////////////////////////////////////////
	// setup data structures
	//////////////////////////////////////////////////////////////////////


		// generate points vector to fill
	biu::IPointVec P( sidechain ? (2*length) : length );
	biu::MoveSequence M( sidechain ? ((2*length)-1) : (length-1) );
	
	
	//////////////////////////////////////////////////////////////////////
	// run enumeration
	//////////////////////////////////////////////////////////////////////
	
	
	P[0] = biu::IntPoint(0,0,0);

	size_t structureNumber = 0;
	
	if (sidechain) {
		// sidechain models
		structureNumber = extendAndCheck_SC( 1, P, M, *latDescr, outMode );
	} else {
		// backbone-only models
		structureNumber = extendAndCheck_BM( 1, P, M, *latDescr, outMode );
	}

	if (outMode != OUT_SILENT) {
		std::cout <<"\n number of structures of length "
					<<length
					<<" = ";
	}
	std::cout	<<structureNumber
				<<std::endl;
	

	
	//////////////////////////////////////////////////////////////////////
	// clear data structures
	//////////////////////////////////////////////////////////////////////
	
	delete latDescr;

	return 0;
}


void
initAllowedArguments(biu::OptionMap & allowedArgs, std::string &infoText )
{
	allowedArgs.push_back(biu::COption(
							"len", false, biu::COption::INT,
							"the length of the structures to enumerate"));
	allowedArgs.push_back(biu::COption(
							"lat", false, biu::COption::STRING,
							"lattice of the structures : SQR - 2D square, CUB - 3D cubic, FCC - 3D face-centered-cubic, 210 - 3D chess knights walk)"));
	allowedArgs.push_back(biu::COption(
							"sideChain", true, biu::COption::BOOL,
							"use when given structures are side chain structures"));
	allowedArgs.push_back(biu::COption(
							"s", true, biu::COption::BOOL,
							"silent : no output except results"));
	allowedArgs.push_back(biu::COption(
							"help", true, biu::COption::BOOL,
							"program parameters and help"));
	allowedArgs.push_back(biu::COption(
							"version", true, biu::COption::BOOL,
							"version information of this program"));

	infoText =	std::string(
"LatNum allows for the enumeration of all lattice protein structures of a given length.");

} // initArguments


//////////////////////////////////////////////////////////////////////////




