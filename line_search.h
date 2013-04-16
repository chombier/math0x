#ifndef MATH0X_LINE_SEARCH_H
#define MATH0X_LINE_SEARCH_H

#include <math0x/iter.h>
#include <math0x/func.h>


namespace math0x {

	// secant method on alpha -> f(x + alpha * dir ) 
	struct line_search {
		math0x::iter iter;
		
		template<class F>
		void solve(func::domain<F>& x, const F& f, 
		           const lie::algebra< func::domain<F> >& dir) const {
			
			lie::group< func::domain<F> > dmn(x);
			auto alg = dmn.alg();
			auto exp = dmn.exp();
			
			real alpha_prev = 0;
			real f_alpha_prev = func::push<F>(f, x)(dir);

			real alpha = 1.0;
			real f_alpha = func::push<F>(f, dmn.prod(x, exp(dir) ) )(dir);
			
			iter([&] {
					real num = (alpha - alpha_prev) * f_alpha;
					real den = (f_alpha - f_alpha_prev);
					
					alpha_prev = alpha;
					f_alpha_prev = f_alpha;
					
					alpha -= num / den;
					f_alpha = func::push<F>(f, dmn.prod(x, exp( alg.scal(alpha, dir) ) ) )( dir );
					
					return std::abs(f_alpha);
				});
			
			x = dmn.prod(x, exp( alg.scal(alpha, dir) ) );
		}

	};

}


#endif
