// $Id: OptionParser.hh,v 1.2 2016/08/08 12:41:56 mmann Exp $
#ifndef BIU_OPTIONPARSER_HH_
#define BIU_OPTIONPARSER_HH_

// OptionParser.h: Schnittstelle für die Klasse COptionParser.
/*

	<author>	Martin Mann					</author>
	<created>	9.10.2005					</created>

	<info>									</info>
*/
//
//////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>
#include <sstream>

// get information on available hashes or maps
#include "biu/HashMap.hh"

// include best available hash or map
#if HAVE_UNORDERED_MAP == 1
	#include <unordered_map>
#else
#if HAVE_TR1_UNORDERED_MAP == 1
	#include <tr1/unordered_map>
#else
#if HAVE_GNU_HASH_MAP == 1
	#include <ext/hash_map>
#else
	#include <map>
#endif
#endif
#endif


namespace biu {

/** 
 * Parameter description for COptionParser
 * 
 * @author Martin Mann <mmann@@informatik.uni-freiburg.de>
 */
class COption {

public:

	static std::string DEF_INIT;

	std::string option;			//!< the parameter
	bool optional;				//!< if optional in the list
	int retType;				//!< the type of the parameter
	std::string description;	//!< the help description
	std::string strValue;		//!< the value if given
	std::string defValue;		//!< the default value if optional
	bool exist;					//!< is true if the parameter was given

		//! supported parameter types
	enum TYPES { STRING, CHAR, INT, FLOAT, DOUBLE, BOOL , TYPES_SIZE};
	
		//! inits the parameter names for help output
	static std::vector<std::string> initTypeNames() {
		std::vector<std::string> ret = std::vector<std::string>(TYPES_SIZE, "TYPE");

		ret[STRING]	= "str ";
		ret[CHAR]	= "char";
		ret[INT]	= "int ";
		ret[FLOAT]	= "flt ";
		ret[DOUBLE]	= "dbl ";
		ret[BOOL]	= "bool";

		return ret;
	}
		//! a list of the type names
	static std::vector<std::string> TYPE_NAME;

		//! construction
	COption (std::string _option, bool _optional, int _retType, std::string _description, std::string _defaultValue=DEF_INIT) :
		option(_option),
		optional(_optional),
		retType(_retType),
		description(_description),
		strValue(_defaultValue),
		defValue(_defaultValue),
		exist(false)
		{}

	~COption() {}
};


///////////////////////////////////////////////////////////////
// TYPEDEFS
///////////////////////////////////////////////////////////////

	//! The list of the possible options for the parameter list.
typedef std::vector<COption> OptionMap;
	//! Constant iterator for an OptionMap.
typedef OptionMap::const_iterator OptionMapIt;

/**
 *  Class for type assured parameter parsing and help output
 *  generation.
 * 
 *  @author Martin Mann <mmann@@informatik.uni-freiburg.de>
 */
class COptionParser {
	
private:


	// set typedef for best available hash or map
	#if HAVE_UNORDERED_MAP == 1
		typedef std::unordered_map< std::string, int > STR2INT_MAP;
	#else
	#if HAVE_TR1_UNORDERED_MAP == 1
		typedef std::tr1::unordered_map< std::string, int > STR2INT_MAP;
	#else
	#if HAVE_GNU_HASH_MAP == 1
		typedef __gnu_cxx::hash_map< std::string, int, hash_string > STR2INT_MAP;
	#else
		typedef std::map< std::string, int > STR2INT_MAP;
	#endif
	#endif
	#endif


	OptionMap opt;				//!< the list of possible parameter options
	std::string programName;	//!< the name of the program
	std::string infoText;		//!< additional informations for help output

	bool errorOccured;			//!< == true if the parameters are not parseable
	unsigned int maxOName;
	STR2INT_MAP mopt;	//! mapping of parameter names to the option indices in opt

		//! possible error codes
	enum ERRORS { ERR_NO_OPT, ERR_WR_USE, ERR_WR_VAL, ERR_NO_ARG };

		//! prints a formatted error output
		//! @param error the error code for the main return value
		//! @param optionName the name of the option not found
		//! @param errormsg the error message to display
	void coutError(int error, std::string optionName, std::string errormsg);

		//! tests whether or not the parameters are parseable
		//! @param val the string to cast
		//! @param type the variable to cast to
		//! @return true if a cast is possible, false otherwise
	bool isCastable(std::string val, int type) const;

		//! checks the parameters
	void parseOpt(int argc, char** argv);

		//! help function for help output generation
	void coutLineBreaking(std::string text, std::ostream &os,const int emptyHeadSize,const int lineLength) const;

public:

		//! the maximal line length for the help output
	static int OUTPUT_LINE_LENGTH;

		//! construction
	COptionParser(OptionMap _options, int argc, char** argv, std::string infoText);

		//! == true if an error occured during parameter parsing
	bool noErrors();

		//! prints the program help output to std::cout
	void coutUsage() const;

		//! == true if the corresponding parameter was given or has a default
		//! value
	bool argExist(std::string option);

	// Parse-Funktionen

	std::string	getStrVal(std::string arg);
	char	getCharVal(std::string arg);
	int		getIntVal(std::string arg);
	float	getFloatVal(std::string arg);
	double	getDoubleVal(std::string arg);
	bool	getBoolVal(std::string arg);

};

} // namespace biu


/*  Beispiel

int main(int argc, char ** argv) {

	OptionMap allowedArgs;
	allowedArgs.push_back(COption("text", false, COption::STRING, "input sequence file with given degeneration"));
	allowedArgs.push_back(COption("max", false, COption::INT, "maximal degeneration that is used"));
	...

	COptionParser opts = COptionParser(allowedArgs, argc, argv, "infotext");

	if (opts.noErrors()) {
	    std::string text = opts.getStrVal("text");

		...

		return 0;

	} else {

		...

		return -1;
	}
}

*/

#endif // BIU_OPTIONPARSER_HH_
