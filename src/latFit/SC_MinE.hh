#ifndef SC_MINE_HH_
#define SC_MINE_HH_


#include "ell/StateCollector.hh"

namespace ell
{

	/*! A StateCollector that stores the minimal energy seen so far.
	 * 
	 * @author Martin Mann
	 */
	class SC_MinE : public SC_Counting
	{
	protected:
		  //! the minimal energy seen so far
		double minE;
		
	public:
	
		SC_MinE();
		
		virtual ~SC_MinE();
		
		virtual void add(const State& s);
	
		  //! access to the minimal energy seen so far
		virtual double getMinE() const;
		
	};

}

#endif /*SC_MINE_HH_*/
