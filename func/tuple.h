#ifndef MATH0X_FUNC_TUPLE_H
#define MATH0X_FUNC_TUPLE_H

#include <math0x/tuple.h>
#include <math0x/meta.h>
#include <math0x/func.h>

namespace math0x { 
	namespace func {
  
		// function tuple: given f1: A1 -> B1, ..., fn: An -> Bn,
		// constructs a function f: (A1, ..., An) -> (B1, ..., Bn) such that
		// f(x1, ... xn) = (f1(x1), ..., fn(xn) )
		template<class ... Args>
		struct tuple {
			typedef tuple base;
    
			typedef std::tuple<Args...> args_type;
			args_type args;

			typedef math0x::tuple::range< sizeof...(Args) > each;
    
			typedef std::tuple< func::domain<Args>... > domain_type;
			typedef std::tuple< func::range<Args>... > range_type;
    
			template<int I>
			using type = math0x::tuple::element<I, args_type>;
    

			// TODO better return types ?
			struct call_lvalue {
				const args_type& args;
				const domain_type& x;
      
				template<int I>
				func::range< type<I> > operator()() const {
					return std::get<I>(args)( std::get<I>(x) );
				}
      
			};
    
			struct call_rvalue {
				const args_type& args;
				domain_type&& x;
      
				template<int I>
				func::range< type<I> > operator()() const {
					return std::get<I>(args)( std::move( std::get<I>(x) ) );
				}
      
			};
    

			range_type operator()(const domain_type& x) const {
				return each::map( call_lvalue{args, x} );
			}

			range_type operator()(domain_type&& x) const {
				return each::map( call_rvalue{args, std::move(x) } );
			}
    
    
			struct push : tuple< func::push<Args>... > {

				struct get {
	
					const args_type& args;
					const domain_type& at;

					template<int I>
					func::push< type<I> > operator()() const {
						return { std::get<I>(args), std::get<I>(at) };
					}
	
				};
      
				push(const tuple& of, const domain_type& at) 
					: push::base{ each::map( get{of.args, at} )  } {
	
				}

			};


			struct pull : tuple< func::pull<Args>... > {

				struct get {
	
					const args_type& args;
					const domain_type& at;

					template<int I>
					func::pull< type<I> > operator()() const {
						return { std::get<I>(args), std::get<I>(at) };
					}
	
				};

				pull(const tuple& of, const domain_type& at) 
					: pull::base{ each::map( get{of.args, at} )  } {
	
				}
      
			};
    
    
		};


  
		template<class ... Args>
		tuple< meta::decay<Args>... > make_tuple(Args&& ... args) {
			return { std::make_tuple(std::forward<Args>(args)...) };
		}
  
	}

}
#endif
