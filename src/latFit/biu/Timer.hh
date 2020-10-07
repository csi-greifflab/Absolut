// $Id: Timer.hh,v 1.2 2016/08/08 12:41:58 mmann Exp $
#ifndef BIU_TIMER_HH_
#define BIU_TIMER_HH_

#include <ctime>

namespace biu
{

		/**
		 * Timer class to measure runtime in miliseconds.
		 *
		 * @author Martin Mann <mmann@@informatik.uni-freiburg.de>
		 */
	class Timer {
		private:
				//! starting time
			clock_t t0;
		public:
				//! Sets starting time.
			void start(void){
				t0 = clock();
			}
				//! Returns time consumption in miliseconds until now from last
				//! start() call on.
			double stop(void) {
				return (static_cast<double>(clock()-t0) / CLOCKS_PER_SEC) * 1000.0;
			}
	};
	
} // namespace biu
#endif /*TIMER_HH_*/
