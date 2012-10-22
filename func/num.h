#ifndef MATH0X_FUNC_NUM_H
#define MATH0X_FUNC_NUM_H

#include <math0x/func.h>
#include <math0x/lie.h>
#include <math0x/func/form.h>

namespace math0x {

	namespace func {

		// numerical differentiation
		template<class F>
		struct num {
			typedef num base;

			F of;
			RR step;
		
			// TODO perfect forwarding ?
			range<F> operator()(const domain<F>& x) const {
				return of(x);
			}		
			
			struct push  {
				
				F of;
				domain<F> at;
				range<F> value;
				
				lie::group< domain<F> > dmn;
				lie::group< range<F> > rng;
				
				range<F> inv_value; 
				
				lie::exp< domain<F> > exp;
				lie::log< range<F> > log;
				
				RR step;

				push(const num& of, const domain<F>& at)
					: of(of.of),
					  at(at),
					  value( of(at) ),
					  dmn(at),
					  rng(value),
					  inv_value( rng.inv( value ) ),
					  exp( dmn.exp() ),
					  log( rng.log()),
					  step( of.step )
				{
					
				}
				     
				 
				lie::algebra< range<F> > operator()(const lie::algebra< domain<F> >& v) const {
					range<F> forward = rng.prod(inv_value,  of( dmn.prod(at, exp( dmn.alg().scal(step, v)) ) ) );
					range<F> backward = rng.prod(inv_value,  of( dmn.prod(at, exp( dmn.alg().scal(-step, v)) ) ) );
				
					auto rng_alg = rng.alg();
					
					return rng_alg.scal(0.5 / step, rng_alg.sum(log(forward), rng_alg.minus(log(backward))));
				}				
				
			};

			
			// TODO find more efficient ?
			struct pull {
				
				push df;

				euclid::space< lie::algebra< domain<F> > > dmn_alg;
				euclid::space< lie::coalgebra< domain<F> > > dmn_coalg;
				
				pull(const num& of, const domain<F>& at) 
					: df(of, at),
					  dmn_alg( df.dmn.alg() ),
					  dmn_coalg( *dmn_alg )
				{
					
				}
				
				lie::coalgebra< domain<F> > operator()( const lie::coalgebra<range<F> >& f) const {
					
					lie::algebra< domain<F> > v = dmn_alg.zero();
					
					lie::coalgebra< domain<F> > res = dmn_coalg.zero();
					
					auto coord = make_form(f) << df;
					
					for(NN i = 0, n = dmn_alg.dim(); i < n; ++i ) {
						dmn_alg.coord(i, v) = 1;

						dmn_coalg.coord(i, res) = coord(v);
						
						dmn_alg.coord(i, v) = 0;
					}

					return res;
				};


			};

			
			
		};

		template<class F>
		num< meta::decay<F> > make_num(F&& of, RR step) { return { std::forward<F>(of), step}; }
		
	}
}
#endif
