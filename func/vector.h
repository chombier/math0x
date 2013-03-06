#ifndef MATH0X_FUNC_VECTOR_H
#define MATH0X_FUNC_VECTOR_H

#include <math0x/func/array.h>
#include <math0x/vector.h>

namespace math0x {
	namespace func {
		
		template<class F, int M = -1 >
		using vector = array<F, 
		                     math0x::vector< domain<F>, M >,
		                     math0x::vector< range<F>, M >,
		                     M >;

		template<class F, int M = -1>
		using vector_tie = array_tie<F, math0x::vector< range<F>, M >, M >;
		

		// f: NN -> (Domain -> Range)
		template<class Fun>
		auto make_vector_tie(NN size, Fun&& fun) -> func::vector_tie< decltype(fun(0)) > {
			return {size, std::forward<Fun>(fun)};
		}
		
	}
}


#endif
