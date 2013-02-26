#ifndef MATH0X_FUNC_TIE_H
#define MATH0X_FUNC_TIE_H

#include <math0x/tuple.h>

namespace math0x { 
	namespace func {

		// TODO provide a tie with dynamic sized arrays ? i.e. using
		// math0x::array ?

		// given functions f1 : A -> B1, f2 : A -> B2, ..., fn: A -> Bn,
		// builds the function f: A -> B1 x B2 x ... x Bn such that:

		// f(a) = (f1(a), f2(a), ..., fn(a))
		template<class ... Args>
		struct tie {
			typedef tie base;

			typedef std::tuple< Args... > args_type;
			args_type args;

			tie(const args_type& args) : args(args) { }
			tie(args_type&& args = {}) : args(std::move(args)) { }
			
			template<int I>
			using type = math0x::tuple::element<I, args_type>;

			typedef func::domain< type<0> > domain;
			typedef std::tuple< func::range<Args> ... > range;
    
			static constexpr int size = sizeof...(Args);

			static_assert( std::is_same< math0x::tuple::repeat< domain, size  >,
			                             std::tuple< func::domain<Args>... > >::value,
			               "functions should share domain" );

			typedef math0x::tuple::range< sizeof...(Args) > each;
    
			struct call {
      
				const args_type& args;
				const domain& x;

				template<int I>
				func::range< type<I> > operator()() const {
					return std::get<I>(args)(x);
				}
      
			};
    
			range operator()(const domain& x) const {
				return each::map( call{args, x} );
			}
    
    
			struct push : tie< func::push<Args>... > {
      
				struct get_push {
					
					const tie::args_type& args;
					const domain& at;
					
					template<int I>
					func::push< type<I> > operator()() const {
						return {std::get<I>(args), at};
					}
					
				};
				
				

				push( const tie& of,
				      const domain& at ) 
					: push::base( each::map( get_push{of.args, at} ) )
				{
					
				}

			};

			struct pull  {

				typedef euclid::space< lie::coalgebra<domain> > space_type;
				space_type space;
      
				typedef std::tuple< func::pull<Args> ... > args_type;
				args_type args;
      
				struct call {
	
					lie::coalgebra<domain>& res;
	
					const space_type& space;
					const args_type& args;
					const lie::coalgebra<range>& fx;
	
					template<int I>
					void operator()() const {
						res = space.sum(res, std::get<I>(args)( std::get<I>(fx) ) );
					}
					
				};

				
				lie::coalgebra<domain> operator()(const lie::coalgebra<range>& fx) const {
					lie::coalgebra<domain> res = space.zero();
	
					each::apply( call{res, space, args, fx} );

					return res;
				};


				struct get_pull {
					
					const tie::args_type& args;
					const domain& at;
					
					template<int I>
					func::pull< type<I> > operator()() const {
						return {std::get<I>(args), at};
					}
					
				};

      
				pull(const tie& of, const domain& at)
					: space( *lie::group<domain>(at).alg() ),
					  args( each::map(get_pull{of.args, at}) ) {
					
				}
      
			};

		};

		template<class ... Args>
		tie< meta::decay<Args>... > make_tie( Args&& ... args) {
			return { std::make_tuple( std::forward<Args>(args)... )  };
		}

	}

}
#endif
