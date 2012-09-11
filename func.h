#ifndef GROUP_FUNC_H
#define GROUP_FUNC_H

#include <group/types.h>
#include <group/meta.h>

#include <stdexcept>

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
  
  // default automatic function traits
  namespace impl {
          
    template<class Range, class G, class Domain>
    meta::decay<Range> range(Range (G::*)(const Domain&) const, ...  );

    template<class Range, class G, class Domain>
    meta::decay<Range> range(Range (G::*)(Domain&&) const, G* );

    template<class Range, class G, class Domain>
    Domain domain(Range (G::*)(const Domain&) const, ... );

    template<class Range, class G, class Domain>
    Domain domain(Range (G::*)(Domain&&) const, G* );
 
    // we require member typename ::self to get successive push/pull
    // rights. other, pull::pull refers to the first pull due to
    // c++ name injection
    template<class F>
    using self = typename F::self;
      
    template<class F>
    typename self<F>::push push(F*);

    template<class F>
    typename self<F>::pull pull(F*);
      
    template<class F>
    struct default_push : func::error< lie::alg< func::domain<F> >,
				       lie::alg< func::range<F> >,
				       std::logic_error > {
      default_push(const F&, const func::domain<F>& )
	: default_push::self{ std::logic_error("no pushforward lol") }  {
	
      }      
      
    };

    template<class F>
    struct default_pull : func::error< lie::coalg< func::range<F> >,
				       lie::coalg< func::domain<F> >,
				       std::logic_error > {
      default_pull(const F&, const func::domain<F>& )
	: default_pull::self{ std::logic_error("no pullback lol") } {
	
      }      
      
    };

    
    // default pushforward
    template<class F>
    default_push<F> push( ... );

    // default pullback
    template<class F>
    default_pull<F> pull( ... );
    
  }

  template<class F, class >
  struct traits  {

    typedef decltype( impl::range( &F::operator(), 0) ) range;
    typedef decltype( impl::domain(&F::operator()) ) domain;
      
    typedef decltype( impl::push<F>( 0 )) push;
    typedef decltype( impl::pull<F>( 0 )) pull;
  };
  
  

}


#endif
