#ifndef MATH0X_FUNC_GET_H
#define MATH0X_FUNC_GET_H

#include <math0x/lie.h>
#include <math0x/tuple.h>
#include <math0x/func/part.h>

namespace math0x {

  namespace func {
    
    template<int I, class> struct get;
    template<int I, class> struct partial;

    // get the ith element from a tuple
    template<int I, class ... Args>
    struct get<I, std::tuple<Args...> > {
      typedef get base;
      
      typedef std::tuple<Args...> domain;
      typedef tuple::element<I, domain> range;

      range operator()(const domain& x) const {
	return std::get<I>(x);
      }


      struct push : get<I, lie::algebra<domain> > {
      	push(const get&, const domain& ) { }
      };


      struct pull : partial<I, lie::coalgebra<domain> > {
      	pull(const get&, const domain& at)
      	  : pull::base( lie::group<domain>(at).coalg().zero() ) { 

      	}
      };
      
    };

  }
}



#endif
