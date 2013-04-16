#ifndef MATH0X_NLCG_H
#define MATH0X_NLCG_H

#include <math0x/euclid.h>
#include <math0x/func.h>

#include <math0x/iter.h>


namespace math0x {

	// non-linear conjugate gradient
	struct nlcg {
		
		math0x::iter iter;

		template<class F>
		void solve( func::domain<F>& x, const F& f) const {
			
			// domain structure
			lie::group< func::domain<F> > dmn(x);
			
			auto dmn_alg = dmn.alg();
			auto dmn_coalg = dmn.coalg();
			
			auto exp = dmn.exp();
			
			// dimension
			NN n = dmn_alg.dim();

			// descent/conjugate direction
			typedef euclid::coords< lie::algebra< func::domain<F> > > vec;
			vec d, s, d_prev;
			
			d = vec::Zero( n );
			d_prev = vec::Zero( n );
			s = vec::Zero( n );

			// position delta
			lie::algebra< func::domain<F> > delta;
			delta = dmn_alg.zero();

			real old = 1.0;
			real theta_prev = f(x);
			
			iter( [&] {
					
					// descent direction 
					dmn_coalg.get(d, func::pull<F>(f, x)(1.0) ); 
					
					real d_norm2 = d.squaredNorm();
					
					// fletcher-reeves
					real beta = d_norm2 / old;
					
					// polak-ribiere
					beta -= d.dot( d_prev ) / old;
					
					// direction reset TODO wtf is this ?
					beta = std::max(0.0, beta);

					// conjugation
					s = d + beta * s;

					// line-search along s
					dmn_alg.set(delta, s);
					
					line_search(x, f, delta);
					
					d_prev.swap( d );
					old = d_norm2;
						
					// stop criteria
					real theta = f(x);
					real res = std::abs( theta - theta_prev ) / theta;
					theta_prev = theta;

					return res;
				});
			
		}


		
		// nlls: minimizes || log(inv(v).f(x)) ||^2
		// lambda is a levenberg-marquart-like damping parameter
		template<class F>
		void solve( func::domain<F>& x, 
		            const F& f,
		            const func::range<F>& b,
		            real lambda = 0) const {
			
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
			
			// dimension
			NN n = dmn_alg.dim();
			
			// descent/conjugate direction
			typedef euclid::coords< lie::algebra< func::domain<F> > > vec;
			vec d, s, d_prev;

			d = vec::Zero( n );
			d_prev = vec::Zero( n );
			s = vec::Zero( n );
			
			lie::algebra< func::domain<F> > delta;
			delta = dmn_alg.zero();
			
			// inv(b)
			func::range<F> b_inv = rng.inv(b);
			
			// temporaries
			lie::algebra< func::range<F> > df_s, e;

			func::range<F> fx;

			real old = 1.0;
			
			fx = f(x);
			e = log( rng.prod(b_inv, fx) );
			real theta_prev = rng_alg.norm2( e );
			
			iter( [&] {
					
					// descent direction 
					dmn_coalg.get(d, func::pull<F>(f, x)( rng_alg.transpose(e)) ); // TODO alloc
					
					real d_norm2 = d.squaredNorm();
					
					// fletcher-reeves
					real beta = d_norm2 / old;
					
					// polak-ribiere
					beta -= d.dot( d_prev ) / old;
					
					// direction reset TODO wtf is this ?
					beta = std::max(0.0, beta);

					// conjugation
					s = d + beta * s;
					
					dmn_alg.set(delta, s);
					df_s = func::push<F>(f, x)(delta);
					
					// line-search based on linearization
					real alpha = -rng_alg.dot(e, df_s) / rng_alg.norm2( df_s );
					dmn_alg.set( delta, (alpha / (1 + lambda)) * s );

					// step estimate
					x = dmn.prod(x, exp( delta ) );
					
					// update stuff
					fx = f(x);
					
					e = log( rng.prod(b_inv, fx) );
					
					d_prev.swap( d );
					old = d_norm2;
						
					// stop criteria
					real theta = rng_alg.norm2( e );
					real res = std::abs( theta - theta_prev ) / theta;
					theta_prev = theta;
				
					return res;
				} );

		}


	};


}


#endif
