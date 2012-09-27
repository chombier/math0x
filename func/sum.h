#ifndef GROUP_FUNC_SUM_H
#define GROUP_FUNC_SUM_H

#include <group/euclid.h>

#include <group/func/tie.h>
#include <group/func/id.h>

namespace func {

  template<class E>
  struct sum {
    typedef sum base;
    
    euclid::space<E> space;

    E operator()(const std::tuple<E, E>& x) const {
      return space.sum( std::get<0>(x), std::get<1>(x) );
    }
    
    struct push;
    struct pull;
    
  };



  template<class E>
  struct sum<E>::push : sum {
    push(const sum& of, const E& ) : push::base(of) { }
  };



  template<class E>
  struct sum<E>::pull : tie< id< euclid::dual<E> > > {

    pull(const sum&, const E& )  { }

  }; 

}

#endif
