#ifndef MATH0X_FUNC_GET_H
#define MATH0X_FUNC_GET_H

#include <math0x/types.h>
#include <math0x/lie.h>
#include <math0x/tuple.h>
#include <math0x/vector.h>
#include <math0x/covector.h>
#include <math0x/func/part.h>

namespace math0x {

	namespace func {
    
		// get the ith element from a tuple
		template<class ... Args, int I>
		struct get<std::tuple<Args...>, I> {

			static_assert( I >= 0, "tuple index must be positive" );
			
			typedef get base;
      
			typedef std::tuple<Args...> domain;
			typedef math0x::tuple::element<I, domain> range;


			const range& operator()(const domain& x) const {
				return std::get<I>(x);
			}

			range operator()(domain&& x) const {
				return std::move( std::get<I>(x) );
				
				// TODO is it ok to move here ?
				// return std::get<I>(x);
			}
      
      
			struct push : get<lie::algebra<domain>, I > {
				push(const get&, const domain& ) { }
			};


			struct pull : part<lie::coalgebra<domain>, I> {
				pull(const get&, const domain& at)
					: pull::base( (*lie::group<domain>(at).alg()).zero() ) { 
					
				}
			};
      
		};


		// ith-element from a vector
		template<class U, int M >
		struct get< vector<U, M> > {
			typedef get base;

			NN i;
			
			get(NN i) : i(i) {
				assert( M == -1 || i < M );
			}
			
			typedef vector<U, M> domain;

			const U& operator()(const domain& x) const {
				return x(i);
			}

			U operator()(domain&& x) const {
				// TODO is it safe to move here ?
				return std::move( x(i) );
			}

			
			struct push : get<lie::algebra<domain> > {
				push(const get& of, const domain& ) : push::get(of.i){ }
			};
			

			struct pull : part<lie::coalgebra<domain> > {
				pull(const get& of, const domain& at)
					: pull::base( (*lie::group<domain>(at).alg()).zero(), of.i ) { 
					
				}
			};

		};

		// TODO remove copypasta with base class
		template<class U, int M >
		struct get< covector<U, M> > {
			typedef get base;
			
			NN i;
			
			get(NN i) : i(i) {
				assert( M == -1 || i < M );
			}
			
			typedef covector<U, M> domain;

			const U& operator()(const domain& x) const {
				return x(i);
			}

			U operator()(domain&& x) const {
				// TODO is it safe to move here ?
				return std::move( x(i) );
			}

			
			struct push : get<lie::algebra<domain> > {
				push(const get& of, const domain& ) : push::get(of.i){ }
			};
			

			struct pull : part<lie::coalgebra<domain> > {
				pull(const get& of, const domain& at)
					: pull::base( (*lie::group<domain>(at).alg()).zero(), i ) { 
					
				}
			};

		};



	}
}



#endif
