#ifndef MATH0X_LEVMAR_H
#define MATH0X_LEVMAR_H

#include <math0x/func.h>
#include <math0x/iter.h>
#include <math0x/coords.h>

namespace math0x {

	// levenberg-marquardt
	struct levmar {
    
		iter outer;			// iterations for gauss-newton (outer loop)
		iter inner;			// iterations for minres (inner loop)
    
		RR lambda;			// TODO template ?
    
		// defaults to vanilla gauss-newton
		levmar(RR lambda = 0) 
			: lambda(lambda) {
      
		}
    
		// don't assemble jacobians
		template<class F>
		iter dense(func::domain< meta::decay<F> >& x,
		           const F& f,
		           const func::range< meta::decay<F> >& y) const {
			typedef F f_type;
      
			typedef func::domain< f_type > domain;
			typedef func::range< f_type > range;
      
			typedef lie::algebra< domain > domain_algebra;
			typedef lie::algebra< range > range_algebra;
      
			typedef euclid::coords<domain_algebra> domain_coords;
			typedef euclid::coords<range_algebra> range_coords;
			
			// lie groups
			lie::group< domain > dmn(x);
			lie::group< range > rng(y);
			
			// lie algebras
			auto dmn_alg = dmn.alg();
			auto rng_alg = rng.alg();
			
			// functions
			auto residual = rng.log() << func::L( rng.inv(y) ) << f;
			auto exp = dmn.exp();
			auto Jf = func::J(f, dmn_alg, rng_alg);
			
			// vars
			range_algebra r = residual( x );
			domain_algebra delta = dmn_alg.zero();
			
			// jacobian matrix and JTJ's diagonal
			typename func::jacobian<F>::range J;
			domain_coords diag;
			
			// temporaries
			domain_coords dmn_tmp, rhs, dx;
			range_coords rng_tmp;
			
			// tangent normal equations
			minres< euclid::space<domain_algebra>::static_dim > normal;
			normal.iter = inner;
			
			// FIXME g++-4.7 chokes when this is defined inside the lambda
			// below
			auto JTJ = [&](const domain_coords& u) -> const domain_coords& {
				rng_tmp.noalias() = J * u;
				dmn_tmp.noalias() = J.transpose() * rng_tmp;

				if( this->lambda ) dmn_tmp.array() += lambda * (diag.array() * u.array());
				
				return dmn_tmp;
			};

			return outer( [&] ()-> RR  {
	  
					J = Jf(x);
					if( this->lambda ) diag = J.colwise().squaredNorm().transpose();					
					
					rng_alg.get(rng_tmp, r);
					rhs.noalias() = -J.transpose() * rng_tmp;
	  
					normal.solve(dx, JTJ, rhs);
					dmn_alg.set(delta, dx);
					
					x = dmn.prod(x, exp( delta ) );
					r = residual(x);
	  
					return 1.0; //rng_alg.norm(r);
				});
      
		}
		
	};


}


#endif
