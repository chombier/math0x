#ifndef MATH0X_FUNC_OPS_H
#define MATH0X_FUNC_OPS_H

#include <math0x/func/tie.h>
#include <math0x/func/comp.h>

#include <math0x/func/sum.h>
#include <math0x/func/minus.h>
#include <math0x/func/scal.h>
#include <math0x/func/poly.h>

#include <math0x/macro.h>

namespace math0x { 
  namespace func {
  
    template<class LHS, class RHS>
    auto operator+(LHS&& lhs, RHS&& rhs) -> 
      macro_auto( (sum< func::range< meta::decay<RHS> > >{}) << make_tie( std::forward<LHS>(lhs), 
									  std::forward<RHS>(rhs)) );
    template<class Arg>
    auto operator-(Arg&& arg) ->
      macro_auto( minus< func::range< meta::decay<Arg> > >{} << std::forward<Arg>(arg) );
  
    template<class LHS, class RHS>
    auto operator-(LHS&& lhs, RHS&& rhs) -> 
      macro_auto( std::forward<LHS>(lhs) + (- std::forward<RHS>(rhs) ) );


    template<class Arg>
    auto operator*( euclid::field< func::range< meta::decay<Arg> > > lambda, Arg&& arg) ->
      macro_auto( scal< func::range< meta::decay<Arg> > >{lambda} << std::forward<Arg>(arg) );

    template<class Arg>
    auto operator*(  Arg&& arg, euclid::field< func::range< meta::decay<Arg> > > lambda) ->
      macro_auto( lambda * std::forward<Arg>(arg) );
  
    template<class Arg>
    auto operator/( Arg&& arg, euclid::field< func::range< meta::decay<Arg> > > lambda) ->
      macro_auto( std::forward<Arg>(arg) * (1.0 / lambda) );
  
    template<class Arg>
    auto operator^(  Arg&& arg, NN n) ->
      macro_auto( poly< func::range< meta::decay<Arg> > >{n} <<  std::forward<Arg>(arg) );
  
  }

}
#endif
