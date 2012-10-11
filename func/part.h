#ifndef MATH0X_FUNC_PART_H
#define MATH0X_FUNC_PART_H

#include <math0x/tuple.h>
#include <math0x/func/get.h>

namespace math0x {

	namespace func {
    
		template<int I, class> struct partial;
		template<int I, class> struct get;

		// partial function on a tuple
		template<int I, class ... Args>
		struct partial<I, std::tuple<Args...> > {
      
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


			struct push : partial<I, lie::algebra<range> > {

				push(const partial& of, const domain&) 
					: push::base{ lie::group<range>(of.at).alg().zero() }
				{
	  
				}

			};

			struct pull : get<I, lie::coalgebra<range> > {

				pull(const partial& of, const domain&) 
					: pull::base{ lie::group<range>(of.at).alg().zero() }
				{
	  
				}

			};

      
		};

		template<int I, class ... Args>
		partial<I, std::tuple<Args...> > part(const std::tuple<Args...>& at) {
			return {at};
		}

		template<int I, class ... Args>
		partial<I, std::tuple<Args...> > part(std::tuple<Args...>&& at) {
			return {std::move(at)};
		}

		template<int I, class ... Args>
		partial<I, std::tuple<meta::decay<Args>...> > part(Args&&... args) {
			return { std::make_tuple(std::forward<Args>(args)...) };
		}
    
	}

}


#endif
