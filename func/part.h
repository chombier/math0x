#ifndef MATH0X_FUNC_PART_H
#define MATH0X_FUNC_PART_H

#include <math0x/types.h>

#include <math0x/tuple.h>
#include <math0x/vector.h>
#include <math0x/covector.h>

namespace math0x {

	namespace func {
		
		// partial function on a tuple
		template<class ... Args, int I>
		struct part<std::tuple<Args...>, I> {
			
			typedef part base;

			typedef std::tuple<Args...> range;
			typedef math0x::tuple::element<I, range> domain;
      
			part(const range& at) : at(at) {}
			part(range&& at) : at(std::move(at)) {}
			
			range at;
      
			range operator()(const domain& x) const {
				range res = at;
				std::get<I>(res) = x;
				return res;
			}

			range operator()(domain&& x) const {
				range res = at;
				std::get<I>(res) = std::move(x);
				return res;
			}


			struct push : part<lie::algebra<range>, I> {

				push(const part& of, const domain&) 
					: push::base{ lie::group<range>(of.at).alg().zero() }
				{
	  
				}

			};

			struct pull : get<lie::coalgebra<range>, I> {

				pull(const part&, const domain&) { }

			};

      
			};

		template<int I, class ... Args>
		part<std::tuple<Args...>, I> make_part(const std::tuple<Args...>& at) {
			return {at};
		}

		template<int I, class ... Args>
		part<std::tuple<Args...>, I > make_part(std::tuple<Args...>&& at) {
			return {std::move(at)};
		}

		// template<int I, class ... Args>
		// part<std::tuple<meta::decay<Args>...>, I> make_part(Args&&... args) {
		// 	return { std::make_tuple(std::forward<Args>(args)...) };
		// }
    


		

		// partial function on vectors
		template<class U, int M>
		struct part< vector<U, M> > {
			typedef part base;
			
			typedef vector<U, M> range;
			
			range at;
			NN i;
      
			part(const range& at, NN i) : at(at), i(i) {
				assert( i < at.size() );
			}
			
			range operator()(const U& x) const {
				range res = at;
				res(i) = x;
				return res;
			}

			range operator()(U&& x) const {
				range res = at;
				res(i) = std::move(x);
				return res;
			}
			

			struct push : part<lie::algebra<range> > {

				push(const part& of, const U&) 
					: push::base( lie::group<range>(of.at).alg().zero(), of.i)
				{
	  
				}

			};

			struct pull : get<lie::coalgebra<range> > {

				pull(const part& of, const U&) 
					: pull::base( of.i )
				{
	  
				}

			};

		};

		// TODO remove copypasta with base class
		template<class U, int M>
		struct part< covector<U, M> > {
			typedef part base;
			
			typedef covector<U, M> range;
			
			range at;
			NN i;
      
			part(const range& at, NN i) : at(at), i(i) {
				assert( i < at.size() );
			}
			
			range operator()(const U& x) const {
				range res = at;
				res(i) = x;
				return res;
			}

			range operator()(U&& x) const {
				range res = at;
				res(i) = std::move(x);
				return res;
			}
			

			struct push : part<lie::algebra<range> > {

				push(const part& of, const U&) 
					: push::base( lie::group<range>(of.at).alg().zero(), i)
				{
	  
				}

			};

			struct pull : get<lie::coalgebra<range> > {

				pull(const part& of, const U&) 
					: pull::base( of.i )
				{
	  
				}

			};

		};



	}

}


#endif
