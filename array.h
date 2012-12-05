#ifndef MATH0X_ARRAY_H
#define MATH0X_ARRAY_H

#include <math0x/types.h>
#include <math0x/tuple/range.h>

#include <array>
#include <vector>

namespace math0x {

	// array container, possibly fixed size. initialize with a function
	// f: NN -> U

	// array container, fixed size:
	template<class U, int M = -1> 
	class array {
		typedef	std::array<U, M> data_type;
		data_type data;

		template<class F, int... I>
		static data_type make_data(F&& f, tuple::index<I...> ) {
			return data_type{{ f(I)... }};
		}
		 
		template<class F>
		static data_type make_data(F&& f) {
			return make_data(std::forward<F>(f), typename tuple::impl::range<M>::type{} );
		}
		
	public:
		
		U& operator()(NN i) { 
			assert( i < size() );
			return data[i];
		}

		const U& operator()(NN i) const { 
			assert( i < size() );
			return data[i];
		}

		NN size() const { return data.size(); }
		
		template<class F>
		array(NN size, F&& f) 
			: data(make_data(std::forward<F>(f) ) ) {
			assert( int(size) == M );
			meta::noop( size );
		}

		
		
	};
	

	// array container, dynamic size
	template<class U>
	class array<U, -1> {
		std::vector<U> data;
	public:

		U& operator()(NN i) { 
			assert( i < size() );
			return data[i];
		}

		const U& operator()(NN i) const { 
			assert( i < size() );
			return data[i];
		}

		NN size() const { return data.size(); }
		
		array( NN size = 0 ) : data(size) { }

		template<class F>
		array(NN size, F&& f) {
			data.reserve(size);
			for(NN i = 0; i < size; ++i) {
				data.push_back(f(i));
			}
		}
		
	};
	
	

}


#endif
