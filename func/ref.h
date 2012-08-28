#ifndef GROUP_FUNC_REF_H
#define GROUP_FUNC_REF_H

#include <group/func.h>

namespace func {

  // reference wrapper (use to avoid spurious copies)
  template<class F>
  struct reference {
    typedef reference self;
    
    const F& to;

    range<F> operator()(const domain<F>& x) const { return to(x); }
    range<F> operator()(domain<F>&& x) const { return to( std::move(x) ); }

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
