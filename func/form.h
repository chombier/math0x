#ifndef MATH0X_FUNC_FORM_H
#define MATH0X_FUNC_FORM_H

#include <math0x/euclid.h>
#include <math0x/func/line.h>

namespace math0x { 
	namespace func {

		template<class E> struct line;

		// linear form over a vector space: E -> field
		template<class E>
		struct form {
			typedef form base;
    
			euclid::dual<E> value;
			
			form(const euclid::dual<E>& value)
				: value( value )
			{
				
			}

			form(euclid::dual<E>&& value)
				: value( std::move(value) ) {
      
			}
    
			euclid::field<E> operator()(const E& x) const {
      
				euclid::range< euclid::dual<E> > rvalue(value);
				euclid::range< E > rx(x);

				euclid::field<E> res = 0;

				for( ; !rx.empty(); rx.pop(), rvalue.pop() ) {
					res += rvalue.front() * rx.front();
				}
				
				return res;      
			}
			
			struct push;
			struct pull;
			
		};

		template<class E>
		struct form<E>::push : form {
			push(const form& of, const E& ) : form(of) { }
			push(const form& of) : form(of) { }
		};


		template<class E>
		struct form<E>::pull : line< euclid::dual<E> > {
      
			pull( const form& of, const E& )
				: pull::base(of.value) { }
      
			pull(const form& of) : pull::base(of.value) { } 
			
		};
		
		template<class E>
		form< euclid::dual< meta::decay<E> > > make_form(E&& x) { return { std::forward<E>(x) }; }
		
	}

}
#endif
