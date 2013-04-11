#ifndef MATH0X_NLCG_H
#define MATH0X_NLCG_H

#include <math0x/euclid.h>
#include <math0x/func.h>

#include <math0x/iter.h>


namespace math0x {

	// non-linear conjugate gradient
	struct nlcg {
		
		math0x::iter iter;

		real omega{1};

		template<class F>
		void solve( func::domain<F>& x, 
		            const F& f,
		            const func::range<F>& b) const {
			
			// domain structure
			lie::group< func::domain<F> > dmn(x);
			auto dmn_alg = dmn.alg();
			auto dmn_coalg = *dmn_alg;
			
			// range structure
			lie::group< func::range<F> > rng(b);
			auto rng_alg = rng.alg();
			auto rng_coalg = *rng_alg;
			
			// exp/log
			lie::exp< func::domain<F> > exp = dmn.exp();
			lie::log< func::range<F> > log = rng.log();
			
			// descent/conjugate direction
			lie::algebra< func::domain<F> > d, s, delta;

			d = dmn_alg.zero();
			s = d;

			// inv(b)
			func::range<F> b_inv = rng.inv(b);
			
			// temporaries
			lie::algebra< func::range<F> > df_s, e;

			func::range<F> fx;

			real old = 1.0;
			
			fx = f(x);
			e = log( rng.prod(b_inv, fx) );
			
			iter( [&] {
					
					// backup
					delta = d;
					
					// descent direction
					d = dmn_alg.minus( dmn_coalg.transpose( func::dT(f)(x)( rng_alg.transpose(e)) ) );
					
					real d_norm2 = dmn_alg.norm2(d);
					
					// fletcher-reeves
					real beta = d_norm2 / old;
					
					// polak-ribiere
					beta -= dmn_alg.dot(d, delta) / old;
					
					// direction reset
					beta = std::max(0.0, beta);

					s = dmn_alg.sum(d, dmn_alg.scal(beta, s));
					
					df_s = func::d(f)(x)(s);
					
					// line-search based on linearization
					real alpha = - rng_alg.dot(e, df_s) / rng_alg.norm2( df_s );

					delta = dmn_alg.scal(omega * alpha, s);
					x = dmn.prod(x, exp( delta  ) );
					
					old = d_norm2;
					
					fx = f(x);
					e = log( rng.prod(b_inv, fx) );
					
					return std::min( dmn_alg.norm(delta),
					                 rng_alg.norm(e) );
				} );

		}


	};


}


#endif
