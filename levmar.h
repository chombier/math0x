#ifndef MATH0X_LEVMAR_H
#define MATH0X_LEVMAR_H

#include <math0x/func.h>
#include <math0x/iter.h>


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
      
      lie::group< domain > dmn(x);
      lie::group< range > rng(y);
      
      auto dmn_alg = dmn.alg();
      auto rng_alg = rng.alg();
      
      auto residual = rng.log() << func::L( rng.inv(y) ) << f;
      
      range_algebra r = residual( x );
      domain_algebra delta = dmn_alg.zero();
      
      auto Jf = func::J(f, dmn_alg, rng_alg);
      
      typename func::jacobian<F>::range J;

      typedef vector<RR, euclid::space<domain_algebra>::static_dim> domain_vec;
      typedef vector<RR, euclid::space<range_algebra>::static_dim> range_vec;
      
      domain_vec dmn_tmp, rhs, dx;
      range_vec rng_tmp;
      
      auto exp = dmn.exp();

      // tangent normal equations
      minres<euclid::space<domain_algebra>::static_dim> normal;
      normal.iter = inner;
      
      return outer( [&]()  {
	  
	  J = Jf(x);
	  
	  // TODO actually use lambda
	  auto JTJ = [&](const domain_vec& x) {
	    return x;
	    // rng_tmp.noalias() = J * x;
	    // dmn_tmp.noalias() = J.transpose() * rng_tmp;
	    // return dmn_tmp;
	  };
	  
	  rng_alg.get(rng_tmp, r);
	  rhs.noalias() = -J.transpose() * rng_tmp;
	  
	  normal.solve(dx, JTJ, rhs);
	  dmn_alg.set(delta, dx);
	  
	  x = dmn.prod(x, exp(delta) );
	  r = residual(x);
	  
	  return 1.0; //rng_alg.norm(r);
	});
      
    }
		
  };


}


#endif
