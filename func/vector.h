#ifndef MATH0X_FUNC_VECTOR_H
#define MATH0X_FUNC_VECTOR_H

#include <math0x/func/array.h>
#include <math0x/vector.h>

namespace math0x {
	namespace func {
		
		template<class F, int M>
		using vector = array<F, 
		                     math0x::vector< domain<F>, M >,
		                     math0x::vector< range<F>, M >,
		                     M >;

	template<class F, int M>
	using vector_tie = array<F, math0x::vector< range<F>, M >, M >;
	
	}
}


#endif
