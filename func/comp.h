#ifndef MATH0X_FUNC_COMP_H
#define MATH0X_FUNC_COMP_H

#include <math0x/func.h>


namespace math0x { 
	namespace func {

		namespace impl {

			// we use external check class to get more meaningful error
			// messages
			template<class A, class B>
			struct comp_check {

				static_assert( std::is_same< A, B >::value,
				               "function domain/range must agree" );
				
			};

		}


		// function composition
		template<class Outer, class Inner>
		struct comp {
			typedef comp base;
    
			Outer outer;
			Inner inner;
			
			static constexpr impl::comp_check< domain<Outer>, range<Inner> > trigger = {};
			
			auto operator()(const domain<Inner>& x) const -> decltype( outer(inner( x ) ) ) {
				return outer( inner( x ) );
			}
    
			auto operator()(domain<Inner>&& x) const -> decltype( outer(inner( std::move(x)))) {
				return outer( inner( std::move(x) ) );
			}
    
			struct push : comp< func::push<Outer>, 
			                    func::push<Inner> > {

				push(const comp& of, const domain<Inner>& at) 
					: push::base{ func::push<Outer>(of.outer, of.inner(at)),
						func::push<Inner>(of.inner, at ) } {
				}

				push(const comp& of, domain<Inner>&& at) 
					: push::base{ func::push<Outer>(of.outer, of.inner(at)),
						func::push<Inner>(of.inner, std::move(at) ) } {
					
				}
      
			};

			struct pull : comp< func::pull<Inner>, 
			                    func::pull<Outer> > {

				pull( const comp& of, const domain<Inner>& at) 
					: pull::base{ func::pull<Inner>(of.inner, at),
						func::pull<Outer>(of.outer, of.inner(at)) } {
	
				}

				pull( const comp& of, domain<Inner>&& at) 
					: pull::base{ func::pull<Inner>(of.inner, at),
						func::pull<Outer>(of.outer, of.inner( std::move(at) )) } {
					
				}
      
			};


		};

		// TODO move to ops.h ?
		namespace impl {
    
			// only defined when both types are function types
			template<class Outer, class Inner, class = void> struct comp;
    
			template<class Outer, class Inner>
			struct comp<Outer, Inner// , decltype( func::requires<Outer, Inner>() )
			            > {
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

}

#endif
