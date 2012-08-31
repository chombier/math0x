#ifndef GROUP_FUNC_TUPLE_H
#define GROUP_FUNC_TUPLE_H

#include <group/tuple.h>
#include <group/meta.h>

namespace func {
  
  // function tuple: given f1: A1 -> B1, ..., fn: An -> Bn,
  // constructs a function f: (A1, ..., An) -> (B1, ..., Bn) such that
  // f(x1, ... xn) = (f1(x1), ..., fn(xn) )
  template<class ... Args>
  struct tuple {
    typedef tuple self;
    
    typedef std::tuple<Args...> args_type;
    args_type args;

    typedef ::tuple::range<Args...> each;
    
    typedef std::tuple< func::domain<Args>... > domain;
    typedef std::tuple< func::range<Args>... > range;

    struct call_lvalue {
      const args_type& args;
      const domain& x;
      
      template<int I>
      auto operator()() const -> decltype( std::get<I>(args)(std::get<I>(x))) {
	return std::get<I>(args)(std::get<I>(x));
      }
      
    };
    
    struct call_rvalue {
      const args_type& args;
      domain&& x;
      
      template<int I>
      auto operator()() const -> decltype( std::get<I>(args)(std::get<I>(x))) {
	return std::get<I>(args)( std::move( std::get<I>(x) ) );
      }
      
    };
    

    range operator()(const domain& x) const {
      return each::map( call_lvalue{args, x} );
    }

    range operator()(domain&& x) const {
      return each::map( call_rvalue{args, std::move(x) } );
    }
    
    
    struct push : tuple< func::push<Args>... > {

      struct get {
	
	const args_type& args;
	const domain& at;

	template<int I>
	func::push< ::tuple::element<I, args_type> > operator()() const {
	  return { std::get<I>(args), std::get<I>(at) };
	}
	
      };

      push(const tuple& of, const domain& at) 
	: push::self{ each::map( get{of.args, at} )  } {
	
      }

    };

    struct pull : tuple< func::pull<Args>... > {

      struct get {
	
	const args_type& args;
	const domain& at;

	template<int I>
	func::pull< ::tuple::element<I, args_type> > operator()() const {
	  return { std::get<I>(args), std::get<I>(at) };
	}
	
      };

      pull(const tuple& of, const domain& at) 
	: pull::self{ each::map( get{of.args, at} )  } {
	
      }
      
    };
    
    
  };


  template<class ... Args>
  tuple< meta::decay<Args>... > make_tuple(Args&& ... args) {
    return { std::make_tuple(std::forward<Args>(args)...) };
  }
  
}


#endif
