#ifndef MATH0X_FUNC_JACOBIAN_H
#define MATH0X_FUNC_JACOBIAN_H

#include <math0x/func.h>
#include <math0x/lie.h>
#include <math0x/matrix.h>

namespace math0x {

  namespace func {

    template<class F>
    struct jacobian {
      
      F of;
      
      typedef func::domain<F> of_domain;
      typedef func::range<F> of_range;

      typedef lie::algebra<of_range> range_algebra;
      typedef lie::algebra<of_domain> domain_algebra;

      static_assert( std::is_same<euclid::field<range_algebra>,
				  euclid::field<domain_algebra> >::value, 
		     "field mismatch" );

      typedef euclid::field<domain_algebra> field;

      euclid::space< domain_algebra > dmn_alg;
      euclid::space< range_algebra > rng_alg;
      
      // jacobian(const F& of,
      // 	       const euclid::space< domain_algebra >& dmn_alg,
      // 	       const euclid::space< range_algebra >& rng_alg)
      // 	: of(of),
      // 	  dmn_alg(dmn_alg),
      // 	  rng_alg(rng_alg) {

      // }
      
      typedef of_domain domain;
      typedef matrix< field, 
		      euclid::space< domain_algebra >::static_dim, 
		      euclid::space< range_algebra >::static_dim > range;
      
      range operator()(const domain& x) const {

	range res; res.resize( rng_alg.dim(), dmn_alg.dim() );

	domain_algebra delta = dmn_alg.zero();
	
	func::push<F> df_x(of, x);
	
	// TODO parallelize ?
	for(NN i = 0, n = dmn_alg.dim(); i < n; ++i) {
	  dmn_alg.coord(i, delta) = 1;
	  
	  rng_alg.get( res.col(i), df_x(delta) );
	  
	  dmn_alg.coord(i, delta) = 0;
	};
	
	return res;
      }
      
    };


    template<class F>
    jacobian< meta::decay<F> > J(F&& f, 
				 const euclid::space< lie::algebra< func::domain< meta::decay<F> > > >& dmn = {},
				 const euclid::space< lie::algebra< func::range< meta::decay<F> > > > & rng = {} ) {
      return {std::forward<F>(f), dmn, rng};
    }
    
    // convenience
    template<class F>
    jacobian< meta::decay<F> > J(F&& f, 
				 const lie::group< func::domain< meta::decay<F> > >& dmn = {},
				 const lie::group< func::range< meta::decay<F> > >& rng = {} ) {
      return {std::forward<F>(f), dmn.alg(), rng.alg()};
    }


  }

}


#endif
