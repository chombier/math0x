#ifndef MATH0X_FUNC_COVECTOR_H
#define MATH0X_FUNC_COVECTOR_H

#include <math0x/func/array.h>
#include <math0x/covector.h>

namespace math0x {
	namespace func {
		
		template<class F, int M>
		using covector = array<F, 
		                       math0x::covector< domain<F>, M >,
		                       math0x::covector< range<F>, M >,
		                       M >;

		template<class F, int M>
		using covector_tie = array<F, math0x::covector< range<F>, M >, M >;
		
	}
}


#endif
