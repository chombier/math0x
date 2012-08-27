#ifndef GROUP_FUNC_SUM_H
#define GROUP_FUNC_SUM_H

#include <group/euclid.h>

#include <group/func/tie.h>
#include <group/func/id.h>

namespace func {

  template<class E>
  struct sum {
    typedef sum self;
    
    euclid::space<E> space;

    sum(const euclid::space<E>& space = euclid::space<E>() ) : space(space) { }

    E operator()(const std::tuple<E, E>& x) const {
      return space.sum( std::get<0>(x), std::get<1>(x) );
    }
    
    struct push : sum {
      push(const sum& of, const E& ) : push::self(of) { }
    };
    
    
    struct pull : tuple_tie< id< euclid::dual<E> > > {

      pull(const sum&, const E& )  { }

    };
    
  };


}

#endif
