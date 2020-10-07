
#include <iostream>
#include <iomanip>

#include <biu/OptionParser.hh>
#include <biu/LatticeDescriptorSQR.hh>
#include <biu/LatticeDescriptorCUB.hh>
#include <biu/LatticeDescriptorFCC.hh>
#include <biu/LatticeDescriptorCKW.hh>
#include <biu/LatticeModel.hh>

#include <biu/LatticeProteinUtil.hh>
#include <biu/SuperPos_Kabsch.hh>

#include "version.hh"

/*!
 * Allows for the structural comparison of two lattice protein structures.
 * 
 * It calculates the best superpositioning among all symmetric structures and
 * allows for the evaluation of the different structural measures:
 * 
 * - dRMSD
 * - cRMSD
 * - GDT_TS
 * - GDT_HA
 * 
 * For the superpositioning the Kabsch algorithm is used.
 * 
 * @author (c) 2009 Martin Mann http://www.bioinf.uni-freiburg.de/~mmann/
 * 
 */


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


void
printPoints( const biu::DPointVec & p, std::ostream& out = std::cout ) {
	for (size_t i=0; i<p.size(); i++) {
		out <<" (" <<p.at(i) <<"),";
	}
	out <<std::endl;
}



//////////////////////////////////////////////////////////////////////////

void
printEvaluation( const biu::DPointVec & pos1
				, const biu::DPointVec & pos2
				, std::ostream& out = std::cout
				, const size_t outPrec = 3)
{
	out
		<<std::fixed <<std::setprecision(outPrec)
		<<" cRMSD  = " <<biu::LatticeProteinUtil::cRMSD( pos1, pos2 )
					<<" Angstroms" <<std::endl
		<<" dRMSD  = " <<biu::LatticeProteinUtil::dRMSD( pos1, pos2 )
					<<" Angstroms" <<std::endl
		<<" GDT_TS = " <<biu::LatticeProteinUtil::GDT_TS( pos1, pos2 )
					<<std::endl
		<<" GDT_HA = " <<biu::LatticeProteinUtil::GDT_HA( pos1, pos2 )
					<<std::endl
		;

}





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
	
	const double cAlphaDist = opts.getDoubleVal("cAdist");
	if (cAlphaDist <= 0.0) {
		std::cerr <<"\n   ERROR : C_alpha distance has to be greater than zero!\n";
		return -1;
	}
	const int outPrec = opts.getIntVal("outPrec");
	if (outPrec < 0) {
		std::cerr <<"\n   ERROR : output precision has to be at least zero!\n";
		return -1;
	}
	
	enum OUTMODE { OUT_SILENT, OUT_NORMAL, OUT_VERBOSE };
	OUTMODE outMode = OUT_NORMAL;
	if (opts.argExist("s") && opts.argExist("v")) {
		std::cerr <<"\n   ERROR : cannot be silent and verbose at the same time!\n";
		return -1;
	}
	if (opts.argExist("s")) {
		outMode = OUT_SILENT;
	}
	if (opts.argExist("v")){
		outMode = OUT_VERBOSE;	
	}
	
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
	const bool sideChain = opts.argExist("sideChain");

	  // get and convert structures
	const std::string abs1 = opts.getStrVal("abs1");
	const std::string abs2 = opts.getStrVal("abs2");
	  // the according move sequences
	biu::MoveSequence moves1;
	biu::MoveSequence moves2;
	
	const std::string whiteChars = " \t\n\r";
	{
		size_t i = abs1.find_first_not_of( whiteChars , 0);
		size_t l = abs1.find_last_not_of( whiteChars );
		if (i<l) {
			moves1 = biu::LatticeProteinUtil::toMoveSequence(
						abs1.substr(i, l-i+1)
						, *latDescr
						, sideChain
					);
		} else {
			std::cerr <<"\n   ERROR : no first move string given!\n";
			return -1;
		}
	}
	{
		size_t i = abs2.find_first_not_of( whiteChars , 0);
		size_t l = abs2.find_last_not_of( whiteChars );
		if (i<l) {
			moves2 = biu::LatticeProteinUtil::toMoveSequence(
						abs2.substr(i, l-i+1)
						, *latDescr
						, sideChain
					);
		} else {
			std::cerr <<"\n   ERROR : no second move string given!\n";
			return -1;
		}
	}
	if (abs1.size() != abs2.size()) {
		std::cerr <<"\n   ERROR : given move strings are of different lengths!\n";
		return -1;
	}
	
	if (outMode >= OUT_VERBOSE) {
		std::cout <<" abs1 : " 
				<<biu::LatticeProteinUtil::toString( moves1, *latDescr, sideChain )
				<<std::endl;
		std::cout <<" abs2 : " 
				<<biu::LatticeProteinUtil::toString( moves2, *latDescr, sideChain )
				<<std::endl;
	}

		// get points
	biu::DPointVec pos1 = biu::LatticeProteinUtil::toDblPoints(
						moves1
						, *latDescr
						, sideChain
						, cAlphaDist
					);
	biu::DPointVec pos2 = biu::LatticeProteinUtil::toDblPoints(
						moves2
						, *latDescr
						, sideChain
						, cAlphaDist
					);

	if (outMode >= OUT_VERBOSE) {
		std::cout <<" pos1 : ";
		printPoints(pos1, std::cout);
		std::cout <<" pos2 : ";
		printPoints(pos2, std::cout);
	}
	
	//////////////////////////////////////////////////////////////////////
	// run superpositioning
	//////////////////////////////////////////////////////////////////////
	
	if (outMode >= OUT_VERBOSE) {
		std::cout <<"\n ==> superpositioning :\n";
	}
	
	// TODO: try reflection ! maybe trigger via parameter !
	
	biu::SuperPos_Kabsch::bestsuperposition( 
									pos1
									, pos2
									, latDescr->getAutomorphisms() );
	
	if (outMode >= OUT_VERBOSE) {
		std::cout <<" sup1 : ";
		printPoints(pos1, std::cout);
		std::cout <<" sup2 : ";
		printPoints(pos2, std::cout);
	}
	
	//////////////////////////////////////////////////////////////////////
	// calculate distances
	//////////////////////////////////////////////////////////////////////
	
	if (outMode >= OUT_VERBOSE) {
		std::cout <<"\n ==> distance :\n";
	}
	  // print distance evaluation of whole proteins
	printEvaluation( pos1, pos2, std::cout, outPrec );
	
	if (sideChain) {
		std::cout <<"\n ==> backbone data only :\n";
		  // get backbone data only
		biu::DPointVec pos1bb(pos1.size());
		biu::DPointVec pos2bb(pos2.size());
		  // copy backbone data
		for (size_t i=0; i<pos1.size(); i+=2) {
			pos1bb[i/2] = pos1[i];
			pos2bb[i/2] = pos2[i];
		}
		  // superposition backbone data
		biu::SuperPos_Kabsch::bestsuperposition( 
										pos1bb
										, pos2bb
										, latDescr->getAutomorphisms() );
		  // print distance evaluation of protein backbones
		printEvaluation( pos1bb, pos2bb, std::cout, outPrec );
	}
	
	if (sideChain) {
		std::cout <<"\n ==> sidechain data only :\n";
		  // get sidechain data only
		biu::DPointVec pos1sc(pos1.size());
		biu::DPointVec pos2sc(pos2.size());
		  // copy sidechain data
		for (size_t i=1; i<pos1.size(); i+=2) {
			pos1sc[i/2] = pos1[i];
			pos2sc[i/2] = pos2[i];
		}
		  // superposition sidechain data
		biu::SuperPos_Kabsch::bestsuperposition( 
										pos1sc
										, pos2sc
										, latDescr->getAutomorphisms() );
		  // print distance evaluation of protein sidechains
		printEvaluation( pos1sc, pos2sc, std::cout, outPrec );
	}
	
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
							"abs1", false, biu::COption::STRING,
							"the first structure in absolute move string representation"));
	allowedArgs.push_back(biu::COption(
							"abs2", false, biu::COption::STRING,
							"the second structure in absolute move string representation"));
	allowedArgs.push_back(biu::COption(
							"lat", false, biu::COption::STRING,
							"lattice of the structures : SQR - 2D square, CUB - 3D cubic, FCC - 3D face-centered-cubic, 210 - 3D chess knights walk)"));
	allowedArgs.push_back(biu::COption(
							"sideChain", true, biu::COption::BOOL,
							"use when given structures are side chain structures"));
	allowedArgs.push_back(biu::COption(
							"cAdist", true, biu::COption::DOUBLE,
							"the C_alpha atom distance to scale the lattice protein to", "3.8"));
	allowedArgs.push_back(biu::COption(
							"outPrec", true, biu::COption::INT,
							"output precision, i.e. number of decimal places given", "3"));
	allowedArgs.push_back(biu::COption(
							"v", true, biu::COption::BOOL,
							"do verbose output"));
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
"LatMap allows for the structural comparison of two lattice protein structures."
"\n\n"
"It calculates the best superpositioning among all symmetric structures and"
" allows for the evaluation of the different structural measures:"
"\n\n"
" - dRMSD  : distance root mean square deviation\n"
" - cRMSD  : coordinate root mean square deviation\n"
" - GDT_TS : global distance test - total score\n"
" - GDT_HA : global distance test - high accuracy score\n"
"\n"
"For superpositioning the Kabsch algorithm is used.\n") ;

} // initArguments


//////////////////////////////////////////////////////////////////////////




