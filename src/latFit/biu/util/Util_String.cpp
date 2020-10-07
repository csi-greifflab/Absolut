// Util_String.cpp: Implementierung der Klasse Util_String.
/*

	<author>	Martin Mann					</author>
	<created>	9.10.2005					</created>

	<info>									</info>
*/
//
//////////////////////////////////////////////////////////////////////

#include "biu/util/Util_String.h"


namespace biu {
  namespace util {

///////////////////////////////////////////////////////////////////
// wandelt int in string um
///////////////////////////////////////////////////////////////////
std::string Util_String::int2str(const int& number) {
	std::ostringstream oss;
	oss << number;
	return oss.str();
}

///////////////////////////////////////////////////////////////////
// wandelt integer string in int um
///////////////////////////////////////////////////////////////////
int Util_String::str2int(const std::string& numString) {
	if (numString.empty())
		return 0;
	std::istringstream iss(numString);
	int retInt=0;
	iss >> retInt;
	return retInt;
}


///////////////////////////////////////////////////////////////////
// liefert anzahl der vorkommen von <c> in <str>
///////////////////////////////////////////////////////////////////
int Util_String::countChar( const std::string &str, const char c) {
	std::string::size_type  i=0;
	int ret = 0;
	for ( ;i<str.size();i++) 
		if (str[i] == c) 
			ret++;
	return ret;
}


///////////////////////////////////////////////////////////////////
// liefert laenge der laengsten wiederholung von <c> in <str>
///////////////////////////////////////////////////////////////////
int Util_String::maxSubseq( const std::string &str, const char c) {
	int m=0, act=0;
	std::string::size_type i=str.find_first_not_of(c);
	for(;i<str.find_last_not_of(c);i++) {
//		std::cerr <<i <<" = " <<str[i] ;
		if (str[i]==c)
			act++;
		else 
			act=0;
		if (act > m) m = act;
//		std::cerr <<" act = " <<act <<" max = " << m <<"\n";
	}
	return m;
}


///////////////////////////////////////////////////////////////////
// konvertiert den string in grossbuchstaben
///////////////////////////////////////////////////////////////////
std::string Util_String::str2upperCase(const std::string &str) {
	std::string ret(str);
	char diff = 'A'-'a';
	for (std::string::size_type i=0; i<ret.size(); i++)
		if (ret[i]>96 && ret[i]<123)
			ret[i] = ret[i]+diff;
	return ret;
}



///////////////////////////////////////////////////////////////////
// prueft ob string nur aus elementen des alphabets besteht
///////////////////////////////////////////////////////////////////
bool Util_String::isAlphStr(const std::string str, const std::string alph) {
	return str.find_first_not_of(alph,0) < str.size();
}

///////////////////////////////////////////////////////////////////
// entfernt fuehrende und endende "\n" und " " und "\t"
///////////////////////////////////////////////////////////////////
std::string Util_String::chompStr(const std::string &str) {
	std::string ret(str);
		// remove leading whitespaces
	ret.erase(0,ret.find_first_not_of("\n \t\r"));
		// remove tailing whitespaces
	ret.erase(ret.find_last_not_of("\n \t\r")+1);
	return ret;
}


///////////////////////////////////////////////////////////////////
// inverts the string
///////////////////////////////////////////////////////////////////
std::string Util_String::invert(const std::string& str) {
	std::string ret(str);
	std::string::size_type i = 0;
	for (std::string::const_reverse_iterator it = str.rbegin(); it!= str.rend(); it++)
		ret[i++] = *it;
	return ret;
}

  } //namespace util
} // namespace biu
