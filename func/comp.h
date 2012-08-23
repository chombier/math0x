#ifndef GROUP_FUNC_COMP_H
#define GROUP_FUNC_COMP_H

#include <group/func.h>

// function composition
namespace func {

  template<class Outer, class Inner>
  struct comp {
    typedef comp self;
    
    Outer outer;
    Inner inner;
  
    static_assert( std::is_same< domain<Outer>, range<Inner> >::value,
		   "function domains must agree" );

    range<Outer> operator()(const domain<Inner>& x) const {
      return outer( inner( x ) );
    }
    
    struct push : comp< func::push<Outer>, 
			func::push<Inner> > {

      push(const comp& of, const domain<Inner>& at) 
	: push::self{ func::push<Outer>(of.outer, of.inner(at)),
	  func::push<Inner>(of.inner, at ) } {
      }
      
    };

    struct pull : comp< func::pull<Inner>, 
			func::pull<Outer> > {

      pull( const comp& of, const domain<Inner>& at) 
	: pull::self{ func::pull<Inner>(of.inner, at),
	  func::pull<Outer>(of.outer, of.inner(at)) } {
	
      }
      
    };


  };

  namespace impl {

    // only defined when both types are function types
    template<class Outer, class Inner, class = void> struct comp;
    
    template<class Outer, class Inner>
    struct comp<Outer, Inner, decltype( func::requires<Outer, Inner>() ) > {
      typedef func::comp<Outer, Inner> type;
    };
    
  }

  template<class Outer, class Inner>
  typename impl::comp< meta::decay<Outer>,
		       meta::decay<Inner> >::type operator<<(Outer&& outer, Inner&& inner) {
    return { std::forward<Outer>(outer), std::forward<Inner>(inner) };
  }

  template<class Outer, class Inner>
  typename impl::comp< meta::decay<Outer>,
		       meta::decay<Inner> >::type operator>>(Inner&& inner, Outer&& outer) {
    return { std::forward<Outer>(outer), std::forward<Inner>(inner) };
  }
  
}


#endif