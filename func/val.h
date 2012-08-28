#ifndef GROUP_FUNC_VAL_H
#define GROUP_FUNC_VAL_H

#include <group/lie.h>
#include <group/meta.h>

namespace func {

  // constant value
  template<class Domain, class Range>
  struct value {
    typedef value self;
    
    Range data;

    const Range& operator()(const Domain& ) const { return data; }
    
    struct push : value< lie::alg<Domain>, lie::alg<Range> > {
      push(const value& of, const Domain& )
	: push::self{ lie::group<Range>(data).alg().zero() } {

      }
    };

    struct pull : value< lie::coalg<Range>, lie::coalg<Domain> > {
      pull(const value& , const Domain& at)
	: pull::self{ lie::group<Domain>(at).coalg().zero() } {
	
      }
    };
    
  };


  template<class Domain, class Range>
  value<Domain, meta::decay<Range> > val(Range&& data) { return { std::forward<Range>(data) }; }

}


#endif
