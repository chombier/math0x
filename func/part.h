#ifndef MATH0X_FUNC_PART_H
#define MATH0X_FUNC_PART_H

#include <math0x/tuple.h>
#include <math0x/func/get.h>

namespace math0x {

  namespace func {
    
    template<int I, class> struct partial;
    template<int I, class> struct get;

    // partial function on a tuple
    template<int I, class ... Args>
    struct partial<I, std::tuple<Args...> > {
      
      typedef partial base;

      typedef std::tuple<Args...> range;
      typedef tuple::element<I, range> domain;
      
      range at;
      
      range operator()(const domain& x) const {
	range res = at;
	std::get<I>(res) = x;
	return res;
      }


      struct push : partial<I, lie::algebra<range> > {

	push(const partial& of, const domain& at) 
	  : push::base{ lie::group<range>(of.at).alg().zero() }
	{
	  
	}

      };

       struct pull : get<I, lie::coalgebra<range> > {

	pull(const partial& of, const domain& at) 
	  : pull::base{ lie::group<range>(of.at).alg().zero() }
	{
	  
	}

      };

      
    };


  }

}


#endif
