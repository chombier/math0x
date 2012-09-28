#ifndef MATH0X_FUNC_VAL_H
#define MATH0X_FUNC_VAL_H

#include <math0x/lie.h>
#include <math0x/meta.h>

namespace math0x { 
  namespace func {

    // constant value
    template<class Domain, class Range>
    struct value {
      typedef value base;
    
      Range data;

      const Range& operator()(const Domain& ) const { return data; }
    
      struct push : value< lie::algebra<Domain>, lie::algebra<Range> > {
	push(const value& , const Domain& )
	  : push::base{ lie::group<Range>(data).alg().zero() } {

	}
      };

      struct pull : value< lie::coalgebra<Range>, lie::coalgebra<Domain> > {
	pull(const value& , const Domain& at)
	  : pull::base{ lie::group<Domain>(at).coalg().zero() } {
	
	}
      };
    
    };


    template<class Domain, class Range>
    value<Domain, meta::decay<Range> > val(Range&& data) { return { std::forward<Range>(data) }; }

  }

}
#endif
