#ifndef MATH0X_NLLS_H
#define MATH0X_NLLS_H

// a convenient way of building levmar function/desired value for
// non-linear least-squares

#include <math0x/func/any.h>
#include <math0x/func/vector.h>

namespace math0x {

	template<class Domain, class Range>
	class nlls {
		NN samples;
	public:
		typedef vector<Range> target_type;
		typedef func::any<Domain, target_type > fun_type;

	private:
		fun_type fun;
		target_type target;
	public:
		
		
		nlls(NN samples) 
		: samples(samples) { 	
			
		}
		
		// f: NN -> ( Domain -> Range )
		template<class F>
		void models(F&& f) { 
			typedef decltype(f(0)) fun_type;
			
			static_assert( std::is_same< func::domain<fun_type>, Domain >::value,
			               "domain mismatch" );

			static_assert( std::is_same< func::range<fun_type>, Range >::value,
			               "range mismatch" );
			
			fun = func::make_vector_tie(samples, std::forward<F>(f)); 
		}
		
		// f:: NN -> Range
		template<class F>
		void observations(F&& f) {
			target.resize( samples );

			typedef decltype( f(0) ) fun_type;

			static_assert( std::is_same< meta::decay<fun_type>, Range >::value,
			               "range mismatch" );
			
			each(target, [&](NN i) {
					target(i) = f(i);
				});
		}
		
		const fun_type& models() const { return fun; }
		const target_type& observations() const { return target; }
		
	};
	
}

#endif
