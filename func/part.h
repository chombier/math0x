#ifndef MATH0X_FUNC_PART_H
#define MATH0X_FUNC_PART_H

#include <math0x/types.h>
#include <math0x/tuple.h>

namespace math0x {

	namespace func {
		
		// partial function on a tuple
		template<class ... Args, int I>
		struct partial<std::tuple<Args...>, I> {
			
			typedef partial base;

			typedef std::tuple<Args...> range;
			typedef math0x::tuple::element<I, range> domain;
      
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


			struct push : partial<lie::algebra<range>, I> {

				push(const partial& of, const domain&) 
					: push::base{ lie::group<range>(of.at).alg().zero() }
				{
	  
				}

			};

			struct pull : get<lie::coalgebra<range>, I> {

				pull(const partial& of, const domain&) 
					: pull::base{ lie::group<range>(of.at).alg().zero() }
				{
	  
				}

			};

      
		};

		template<int I, class ... Args>
		partial<std::tuple<Args...>, I> part(const std::tuple<Args...>& at) {
			return {at};
		}

		template<int I, class ... Args>
		partial<std::tuple<Args...>, I > part(std::tuple<Args...>&& at) {
			return {std::move(at)};
		}

		template<int I, class ... Args>
		partial<std::tuple<meta::decay<Args>...>, I> part(Args&&... args) {
			return { std::make_tuple(std::forward<Args>(args)...) };
		}
    
	}

}


#endif
