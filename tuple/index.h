#ifndef MATH0X_TUPLE_INDEX_H
#define MATH0X_TUPLE_INDEX_H

#include <tuple>

#include <math0x/meta.h>

// TODO wtf ?
namespace math0x { 
	namespace tuple {

		template<class F, int ... I>
		struct michel {
			typedef std::tuple< decltype( std::declval<F>().template operator()<I>() )... > type;
		};
		

		// indices for static access
		template<int ... I>
		struct index {
    
			static constexpr unsigned N = sizeof...(I);
			
			// maps fun.operator<I>() for all values of I to a tuple. called
			// in reverse order (at least on gcc)
			template<class F>
			static typename michel<F, I...>::type  map(const F& fun) {
				return std::make_tuple( fun.template operator()<I>()... );
			}
    
			// calls fun successively for all values of I
			template<class F>
			static void apply(const F& fun) {
				meta::noop( (fun.template operator()<N - 1 - I>(), 0)... );
			}


			typedef index< I..., N > next_type;
		};

	}

}

#endif

