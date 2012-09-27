#ifndef GROUP_FUNC_REF_H
#define GROUP_FUNC_REF_H

#include <group/func.h>

namespace func {

  // reference wrapper (use it to avoid costly copies)
  template<class F>
  struct reference {
    typedef reference base;
    
    const F& to;

    auto operator()(const domain<F>& x) const -> decltype( to(x) ) { return to(x); }
    auto operator()(domain<F>&& x) const -> decltype( to( std::move(x) ) ) { return to( std::move(x) ); }
    
    struct push : func::push<F> {

      push(const reference& of, const domain<F>& at) 
	: func::push<F>( of.to, at ) { }

    };

    struct pull : func::pull<F> {
      
      pull(const reference& of, const domain<F>& at) 
	: func::pull<F>( of.to, at ) { }
      
    };
    

  };

  template<class F>
  reference< F> ref( const F& to) { return {to}; }
  
}


#endif
