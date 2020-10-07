// $Id: Alphabet.hh,v 1.2 2016/08/08 12:42:01 mmann Exp $
#ifndef BIU_ALPHABET_HH_
#define BIU_ALPHABET_HH_


#include <vector>
#include <string>

// get information on available hashes or maps
#include "biu/HashMap.hh"

// include best available hash or map
#if HAVE_UNORDERED_MAP == 1
	#include <unordered_map>

	template< class K, class T >
	bool
	operator==( const std::unordered_map<K,T> h1, const std::unordered_map<K,T> h2) 
	{
		  // check if equal size
		if (h1.size() != h2.size()) {
			return false;
		}
		  // check if all elements of h1 are present in h2
		typename std::unordered_map<K,T>::const_iterator i2;
		for (typename std::unordered_map<K,T>::const_iterator i=h1.begin(); i!=h1.end(); i++){
			if ((i2 = h2.find(i->first))==h2.end() || (i2->second != i->second) )
				return false;  // not present --> not equal
		}
		  // they are equal
		return true;
	}
	template< class K, class T >
	bool
	operator!=( const std::unordered_map<K,T> h1, const std::unordered_map<K,T> h2) 
	{
		return !(h1 == h2);
	}
#elif HAVE_TR1_UNORDERED_MAP == 1
	#include <tr1/unordered_map>

	template< class K, class T >
	bool
	operator==( const std::tr1::unordered_map<K,T> h1, const std::tr1::unordered_map<K,T> h2) 
	{
		  // check if equal size
		if (h1.size() != h2.size()) {
			return false;
		}
		  // check if all elements of h1 are present in h2
		typename std::tr1::unordered_map<K,T>::const_iterator i2;
		for (typename std::tr1::unordered_map<K,T>::const_iterator i=h1.begin(); i!=h1.end(); i++){
			if ((i2 = h2.find(i->first))==h2.end() || (i2->second != i->second) )
				return false;  // not present --> not equal
		}
		  // they are equal
		return true;
	}
	template< class K, class T >
	bool
	operator!=( const std::tr1::unordered_map<K,T> h1, const std::tr1::unordered_map<K,T> h2) 
	{
		return !(h1 == h2);
	}
#elif HAVE_GNU_HASH_MAP == 1
	#include <ext/hash_map>
#else
	#include <map>
#endif

namespace biu
{
		/*! This class handels an alphabet for the sequence.
		 *
		 *  Each member of the alphabet maps to an unique int index. 
		 *
		 * @author Martin Mann, Sebastian Will, Andreas Richter
		 */
	class Alphabet
	{
	public:
			//! The internal representation of an alphabet element.
		typedef size_t	AlphElem;
			//! The internal representation of a compressed alphabet element.
		typedef unsigned char	CAlphElem;
		
			//! The internal representation of a sequence of alphabet elements.
		typedef std::vector<AlphElem> Sequence;
		
			//! The internal representation of a compressed sequence of alphabet elements.
		typedef std::vector<CAlphElem> CSequence;
		
	private:

		// set typedef for best available hash or map
		#if HAVE_UNORDERED_MAP == 1
			typedef std::unordered_map< std::string, AlphElem > STR2ALPH_MAP;
		#elif HAVE_TR1_UNORDERED_MAP == 1
			typedef std::tr1::unordered_map< std::string, AlphElem > STR2ALPH_MAP;
		#elif HAVE_GNU_HASH_MAP == 1
			typedef __gnu_cxx::hash_map< std::string, AlphElem, hash_string > STR2ALPH_MAP;
		#else
			typedef std::map< std::string, AlphElem > STR2ALPH_MAP;
		#endif


			//! a mapping from a string representation to the corresponding 
			//! alphabet element
		STR2ALPH_MAP string2alph;	
		
			//! a mapping from alphabet elements to their string representation
		std::vector<std::string> alph2string;	
		
			//! the length of a string representation of an alphabet element.
		size_t elementLength;
		
			//! Holds the multiplyer for the ith position of a sequence
			//! when compressed into a CSequence. Each CAlphElem can encode
			//! compressBase.size() many AlphElem elements.
		std::vector<int> compressBase;
		
	public:
	
			/*! Initialises an Alphabet object with the string representations 
			 * given by the alphabetString. The string is members given in 
			 * alphabetString.
			 * 
			 * The string it sliced into substrings of the length elementLength_
			 * , therefor all alphabet elements have to have the same length.
			 */
		Alphabet(	const std::string& alphabetString, 
					const size_t elementLength);
		
			/*! Construct an alphabet from a vector of strings
			*
			* The indexing of elements in the vector has to be preserved by
			* the alphabet
			*/
		Alphabet(const std::vector<std::string> & alphabetStrings);
		
		virtual ~Alphabet();

		bool operator== (const Alphabet& alph2) const;
		bool operator!= (const Alphabet& alph2) const;

			//! Returns the size of the alphabet. 
			//! @return the number of elements in the alphabet
		size_t getAlphabetSize() const;
		
			//! Returns the length of a alphabet element.
			//! @return the string length of one alphabet element
		size_t getElementLength() const;
	
		//////////////////////////////////////////////////	
		// string to sequence to string
		//////////////////////////////////////////////////
	
			//! Returns an internal sequence representation of the string.
			//! @param seqString the string to encode
			//! @return the encoded sequence
		Sequence getSequence(const std::string& seqString) const;
		
			//! Returns the index the alphabet member maps to.
			//! @param alphElemStr the string to encode
			//! @return internal representation of the alphabet element string
		AlphElem getElement(const std::string& alphElemStr) const;
		
			//! Converts an internal sequence representation into a string.
			//! @param sequence the sequence to decode
			//! @return the string encoded by the sequence
		std::string	getString(const Alphabet::Sequence& sequence) const;
		
			//! Returns the string representation of the alphabet member.
			//! @param elem the internal alphabet element representation
			//! @return the string encoded by the element
		std::string getString(const Alphabet::AlphElem& elem) const;
		
			//! Compresses a sequence.
			//! @param sequence the sequence to compress
			//! @return the compressed sequence representation
		CSequence compress(const Alphabet::Sequence& sequence) const;
		
			//! Compresses a string.
			//! @param sequence the string to compress
			//! @return the compressed sequence representation
		CSequence compressS(const std::string& sequence) const;
		
			//! Decompresses a sequence.
			//! @param sequence		the sequence to decompress
			//! @param seqLength 	the length of the uncompressed sequence encoded
			//! @return the sequence that was compressed
		Sequence decompress(const CSequence& sequence, const size_t seqLength) const;
		
			//! Decompresses a sequence and generates directely a string.
			//! @param sequence		the sequence to decompress
			//! @param seqLength 	the length of the uncompressed sequence encoded 
			//! @return the string that was compressed
		std::string decompressS(const CSequence& sequence, const size_t seqLength) const;
		
	
		//////////////////////////////////////////////////
		// miscellaneous
		//////////////////////////////////////////////////
	
			//! Returns whether or not a string contains only elements
			//! of the alphabet or not
			//! @param str the string to check
			//! @return true if the string is valid for this alphabet, false
			//! otherwise
		bool isAlphabetString(const std::string& str) const;
		
			//! Returns whether or not a sequencec contains only elements
			//! compatible with the alphabet or not
			//! @param seq the sequence to check
			//! @return true if the sequence is valid for this alphabet, false
			//! otherwise
		bool isAlphabetSequence(const Sequence& seq) const;
		
			//! Returns the internal index of the AlphElem. The index
			//! starts with 0 and is < getAlphabetSize()
			//! @param elem the alphabet element to get the index for
			//! @return the index of this element
		size_t getIndex(const AlphElem& elem) const;
		
			//! Returns the internal index of the AlphElem corresponding
			//! to the string. The index
			//! starts with 0 and is < getAlphabetSize()
			//! @param elemStr the alphabet element string to get the index for
			//! @return the index of this element
		size_t getIndex(const std::string& elemStr) const;
		
			//! Returns the corresponding AlphElem to the index.
			//! @param index alphabet element index (has to be < getAlphabetSize())
			//! @return the indexed element
		AlphElem getElement(const size_t index) const;
		
	};

} // namespace biu

  // include definitions
#include "biu/Alphabet.icc"

#endif /*ALPHABET_HH_*/
