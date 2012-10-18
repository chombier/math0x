#ifndef MATH0X_FUNC_H
#define MATH0X_FUNC_H

#include <math0x/types.h>
#include <math0x/meta.h>

#include <stdexcept>
#include <typeinfo>

namespace math0x { 
	namespace func {

		namespace impl {

			// converts return type to void
			template<class ... F> void to_void(F...);

			// F qualifies as a function type if it has one of the following
			// operators:
			template<class Domain, class Range, class F>
			bool requires(Range (F::*)(const Domain& ) const);
      
			template<class Domain, class Range, class F>
			bool requires(Range (F::*)(Domain&& ) const);

		}
  
		// this one is only defined when each F qualifies as a function
		template<class ... F>
		auto requires() -> decltype( impl::to_void( impl::requires(&F::operator())...)  );
		
		// default automatic function traits
		namespace impl {
          
			// range type
			
			// lower priority
			template<class Range, class G, class Domain>
			meta::decay<Range> range(Range (G::*)(const Domain&) const, meta::priority<0>  );

			// higher priority
			template<class Range, class G, class Domain>
			meta::decay<Range> range(Range (G::*)(Domain&&) const, meta::priority<1> );
			
			// same here for domain type
			template<class Range, class G, class Domain>
			Domain domain(Range (G::*)(const Domain&) const, meta::priority<0> );

			template<class Range, class G, class Domain>
			Domain domain(Range (G::*)(Domain&&) const, meta::priority<1> );
			
			// push/pull types
			
			// due to c++ name injection rules, if your typename F::push
			// derives from class A, A needs to have a member type A::base
			// equal to A, otherwise push<push<F>> will refer to push<F>
			// instead of push<A>

			template<class F> struct default_push;
			template<class F> struct default_pull;

			template<class F>
			default_push<F> push( meta::priority<0> );
			
			// this fails when we try to get push::push
			template<class F>
			meta::enable_if< !std::is_same< typename F::push, F>::value, 
			                 typename F::push> push( meta::priority<1> );
		
			template<class F>
			typename F::base::push push( meta::priority<2> );

			template<class F>
			default_pull<F> pull( meta::priority<0> );

			// this fails when we try to get pull::pull
			template<class F>
			meta::enable_if< !std::is_same< typename F::pull, F>::value, 
			                 typename F::pull> pull( meta::priority<1> );
			
			template<class F>
			typename F::base::pull pull( meta::priority<2> );
			

			// default pushforward/pullback
			template<class F>
			struct default_push : func::error< lie::algebra< func::domain<F> >,
			                                   lie::algebra< func::range<F> >,
			                                   std::logic_error > {
				default_push(const F&, const func::domain<F>& )
					: default_push::base{ std::logic_error( meta::name<F>() + " has no pushforward") }  {
					
				}      
      
			};

			template<class F>
			struct default_pull : func::error< lie::coalgebra< func::range<F> >,
			                                   lie::coalgebra< func::domain<F> >,
			                                   std::logic_error > {
				default_pull(const F&, const func::domain<F>& )
					: default_pull::base{ std::logic_error( meta::name<F>() + " has no pullback") } {
					
				}      
      
			};
			
		}

		template<class F>
		struct traits<F> // <F, decltype( requires<F>() ) >
		{
			static meta::priority<8> priority(); 
			
			typedef decltype( impl::range( &F::operator(), priority() ) ) range;
			typedef decltype( impl::domain(&F::operator(), priority() ) ) domain;
			
			typedef decltype( impl::push<F>( priority() )) push;
			typedef decltype( impl::pull<F>( priority() )) pull;
		};

	}
}

#endif
