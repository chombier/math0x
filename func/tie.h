#ifndef GROUP_FUNC_TIE_H
#define GROUP_FUNC_TIE_H

#include <group/tuple.h>

namespace func {

  // given functions f1 : A -> B1, f2 : A -> B2, ..., fn: A -> Bn,
  // builds the function f: A -> B1 x B2 x ... x Bn such that:

  // f(a) = (f1(a), f2(a), ..., fn(a))
  template<class ... Args>
  struct tuple_tie {
    typedef tuple_tie self;

    typedef std::tuple< Args... > args_type;
    args_type args;

    typedef func::domain< tuple::element<0, args_type > > domain;
    typedef std::tuple< func::range<Args> ... > range;

    static_assert( std::is_same< tuple::repeat< domain, sizeof...(Args) >,
				 std::tuple< func::domain<Args>... > >::value,
		   "functions should share domain" );

    typedef tuple::range<Args...> each;
    
    struct call {
      
      const args_type& args;
      const domain& x;

      template<int I>
      func::range< tuple::element<I, args_type> > operator()() const {
	return std::get<I>(args)(x);
      }
      
    };
    
    range operator()(const domain& x) const {
      return each::map(call{args, x} );
    }
    
    
    struct push : tuple_tie< func::push<Args>... > {
      
      push( const tuple_tie& of,
	    const domain& at ) 
	: push::self{ func::push<Args>(of.args, at)... } {

      }

    };

    struct pull  {

      typedef euclid::space< lie::coalg<domain> > space_type;
      space_type space ;
      
      typedef std::tuple< func::pull<Args> ... > args_type;
      args_type args;
      
      struct call {
	
	lie::coalg<domain>& res;
	
	const space_type& space;
	const args_type& args;
	const lie::coalg<range>& fx;
	
	
	template<int I>
	void operator()() const {
	  res = space.sum(res, std::get<I>(args)(fx) );
	}
	
      };

      lie::coalg<domain> operator()(const lie::coalg<range>& fx) const {
	lie::coalg<domain> res = space.zero();
	
	each::apply( call{res, space, args, fx} );

	return res;
      };

      
      // pull(const tuple_tie& of, const domain& at)
      // 	: space( *euclid::space<domain>(at) ),
      // 	  args( func::pull<Args>(of.args, at)... ) {
	
      // }
      
    };

  };

  template<class ... Args>
  tuple_tie< meta::decay<Args>... > tie( Args&& ... args) {
    return { std::make_tuple( std::forward<Args>(args)... )  };
  }

}


#endif
