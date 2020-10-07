#include "SC_OutEnergy.hh"

namespace ell
{

	SC_OutEnergy::SC_OutEnergy(	std::ostream& out_)
	 :	SC_MinE(), out(out_)
	{
	}
	
	SC_OutEnergy::~SC_OutEnergy()
	{
	}

	  // This function is used to track all added intermediate States.
	  // @param s the added State
	void 
	SC_OutEnergy::add(const State& s) {
		  //  call handler of superclass
		SC_MinE::add(s);
		  // print to stream
		
		out << s.getEnergy() << std::endl;
	}

}
