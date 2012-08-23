#ifndef GROUP_FUNC_H
#define GROUP_FUNC_H

#include <group/types.h>
#include <group/meta.h>

namespace func {

  namespace impl {

    template<class ... F> void _void(F...);

    // F qualifies as a function type if it has one of the following
    // operators:
    template<class Domain, class Range, class F>
    bool requires(Range (F::*)(const Domain& ) const);
      
    template<class Domain, class Range, class F>
    bool requires(Range (F::*)(Domain&& ) const);

  }
  
  // this one is only defined when each F qualifies as a function
  template<class ... F>
  auto requires() -> decltype( impl::_void( impl::requires(&F::operator())...)  );
  
  // the domain over which the function is defined
  template<class F>
  using domain = typename traits<F>::domain;
  
  // the range of the function
  template<class F>
  using range = typename traits<F>::range;
  
  // pushforward (derivative) type
  template<class F>
  using push = typename traits<F>::push;
  
  // pullback (derivative transpose) type
  template<class F>
  using pull = typename traits<F>::pull;
  
  // default automatic function traits
  namespace impl {
          
    template<class Range, class G, class Domain>
    meta::decay<Range> range(Range (G::*)(const Domain&) const );
      
    template<class Range, class G, class Domain>
    meta::decay<Range> range(Range (G::*)(Domain&&) const );
      
    template<class Range, class G, class Domain>
    Domain domain(Range (G::*)(const Domain&) const );

    template<class Range, class G, class Domain>
    Domain domain(Range (G::*)(Domain&&) const );
      
    // we require member typename ::self to get successive push/pull
    // rights. other, pull::pull refers to the first pull due to
    // c++ name injection
    template<class F>
    using self = typename F::self;
      
    template<class F>
    typename self<F>::push push(F*);

    template<class F>
    typename self<F>::pull pull(F*);
      
    // TODO default push/pull ?
  }


  template<class F>
  struct traits<F, decltype( requires<F>() ) >  {
      

    typedef decltype( impl::range(&F::operator()) ) range;
    typedef decltype( impl::domain(&F::operator()) ) domain;
      
    typedef decltype( impl::push<F>( 0 ))  push;
    typedef decltype( impl::pull<F>( 0 ))  pull;
  };
  
  

}


#endif
