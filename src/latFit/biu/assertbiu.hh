// $Id: assertbiu.hh,v 1.2 2016/08/08 12:42:00 mmann Exp $
#ifndef BIU_ASSERT_BIU_HH
#define BIU_ASSERT_BIU_HH

#ifdef NDEBUG
	#define assertbiu(EXPR,MSG) ((void)0)
#else
	#include <cassert>
	#include <iostream>
	
	/** assertbiu(EXPR,MSG) calls the macro "assert(EXPR)"
	 * if EXPR is true, with the additional error message MSG
	 */
	#define assertbiu(EXPR,MSG) \
            ( (EXPR) || \
              (	(std::cerr<<"\n\tASSERT FAILED : " <<MSG <<"\n\n") \
              	&& (assert(EXPR), 0)) \
            )

#endif // NDEBUG

#endif // ASSERT_BIU_HH
