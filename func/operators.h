#ifndef MATH0X_FUNC_OPERATORS_H
#define MATH0X_FUNC_OPERATORS_H

#include <math0x/func.h>
#include <math0x/func/tie.h>
#include <math0x/func/comp.h>

#include <math0x/func/sum.h>
#include <math0x/func/minus.h>
#include <math0x/func/scal.h>
#include <math0x/func/poly.h>

#include <math0x/macro.h>

namespace math0x { 
  namespace func {
  

	  // TODO requirement checks !
	  
	  namespace impl {


		  template<class Arg>
		  struct unary_traits {
			  
			  typedef func::range<Arg> range_type;

			  typedef func::comp< func::minus<range_type>, Arg > minus_type;
			  typedef func::comp< func::scal<range_type>, Arg > scal_type;
			  typedef func::comp< func::poly< range_type >, Arg > pow_type;
		  };

		  template<class LHS, class RHS>
		  struct binary_traits {
			  typedef func::range<LHS> range_type;
			  // TODO assert range type is consistent ?

			  typedef func::comp< func::sum< range_type >, func::tie<LHS, RHS> > sum_type;
			  typedef func::comp< func::sum< range_type >, func::tie<LHS, typename unary_traits<RHS>::minus_type > > diff_type;
			  
		  };
		  
		  
	  }
	  
	  // unary
	  template<class Arg>
	  meta::second< decltype( func::requires< meta::decay<Arg> >() ),
	                typename impl::unary_traits< meta::decay<Arg> >::scal_type >
	  operator*( euclid::field< func::range< meta::decay<Arg> > > lambda, Arg&& arg ) {
		  return  scal< func::range< meta::decay<Arg> > >{lambda} << std::forward<Arg>(arg);
	  }

	  template<class Arg>
	  meta::second< decltype( func::requires< meta::decay<Arg> >() ),
	                typename impl::unary_traits< meta::decay<Arg> >::scal_type >
	  operator*( Arg&& arg, euclid::field< func::range< meta::decay<Arg> > > lambda ) {
		  return lambda * std::forward<Arg>(arg);
	  }

	  template<class Arg>
	  meta::second< decltype( func::requires< meta::decay<Arg> >() ),
	                typename impl::unary_traits< meta::decay<Arg> >::scal_type >
	  operator/( Arg&& arg, euclid::field< func::range< meta::decay<Arg> > > lambda ) {
		  assert( lambda );
		  return (1.0 / lambda) * std::forward<Arg>(arg);
	  }


	  template<class Arg>
	  meta::second< decltype( func::requires< meta::decay<Arg> >() ),
	                typename impl::unary_traits< meta::decay<Arg> >::minus_type >
	  operator-(Arg&& arg ) {
		  return minus< func::range< meta::decay<Arg> > >{} << std::forward<Arg>(arg);
	  }
	
	  template<class Arg>
	  meta::second< decltype( func::requires< meta::decay<Arg> >() ),
	                typename impl::unary_traits< meta::decay<Arg> >::pow_type >
	  operator^( Arg&& arg, NN n) {
		  return  poly< func::range< meta::decay<Arg> > >{n} <<  std::forward<Arg>(arg);
	  }



	  // binary
    template<class LHS, class RHS>
    meta::second< decltype( func::requires< meta::decay<LHS>, meta::decay<RHS> >() ),
                  typename impl::binary_traits< meta::decay<LHS>, meta::decay<RHS> >::sum_type > 
    operator+(LHS&& lhs, RHS&& rhs) {
	    return sum< func::range< meta::decay<RHS> > >{} << make_tie( std::forward<LHS>(lhs), 
	                                                                 std::forward<RHS>(rhs));
    }

	  template<class LHS, class RHS>
    meta::second< decltype( func::requires< meta::decay<LHS>, meta::decay<RHS> >() ),
                  typename impl::binary_traits< meta::decay<LHS>, meta::decay<RHS> >::diff_type > 
    operator-(LHS&& lhs, RHS&& rhs) {
		  return std::forward<LHS>(lhs) + (-std::forward<RHS>(rhs));
	  }

	


    // template<class Arg>
    // auto operator*( euclid::field< func::range< meta::decay<Arg> > > lambda, Arg&& arg) ->
    //   macro_auto( scal< func::range< meta::decay<Arg> > >{lambda} << std::forward<Arg>(arg) );

    // template<class Arg>
    // auto operator*(  Arg&& arg, euclid::field< func::range< meta::decay<Arg> > > lambda) ->
    //   macro_auto( lambda * std::forward<Arg>(arg) );
  
    // template<class Arg>
    // auto operator/( Arg&& arg, euclid::field< func::range< meta::decay<Arg> > > lambda) ->
    //   macro_auto( std::forward<Arg>(arg) * (1.0 / lambda) );
  
	
  }

}
#endif
