#ifndef SC_OUTSTREAMABS_HH_
#define SC_OUTSTREAMABS_HH_

#include <iostream>
#include <string>

#include "SC_MinE.hh"

namespace ell
{

	/*! A StateCollector that prints the energy of each reported state 
	 * and the move string representation to the given outstream.
	 * 
	 * @author Daniel Maticzka
	 */
	class SC_OutAbs : public SC_MinE
	{
	protected:
		  //! the stream to write the states to
		std::ostream& out;
		const size_t cutoff;
		
	public:
	
		SC_OutAbs(	std::ostream& out, size_t cutoff);
		
		virtual ~SC_OutAbs();
	
		  //! This function is used to track all accepted intermediate States.
		  //! @param s the accepted State
		virtual void add(const State& s);
			
	};

}

#endif /*SC_OUTSTREAMABS_HH_*/
