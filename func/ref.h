#ifndef MATH0X_FUNC_REF_H
#define MATH0X_FUNC_REF_H

#include <math0x/func.h>
#include <math0x/macro.h>

namespace math0x { 
	namespace func {

		// reference wrapper (use it to avoid costly copies)
		template<class F>
		struct reference {
			// typedef reference base;
    
			const F& to;

			// auto operator()(const domain<F>& x) const -> macro_returns( to(x) );
			// auto operator()(domain<F>&& x) const -> macro_returns( to( std::move(x) ) );
			
			range<F> operator()(const domain<F>& x) const { return to(x); }
			range<F> operator()(domain<F>&& x) const { return to( std::move(x)); }
			
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

}
#endif
