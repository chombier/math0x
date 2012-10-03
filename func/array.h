#ifndef MATH0X_FUNC_ARRAY_H
#define MATH0X_FUNC_ARRAY_H

#include <math0x/func.h>
#include <math0x/array.h>
#include <math0x/each.h>

namespace math0x {
	namespace func {

		template<class F, class Domain, class Range = Domain, int M = -1>
		struct array {
			typedef array base;
			
			math0x::array<F, M> data;
			
			template<class Fun>
			array(NN size, const Fun& fun) 
				: data(size, fun) { }
			
			Range operator()(const Domain& x) const {
				Range res; res.resize( data.size() );
				
				each(data, [&](NN i) {
						res(i) = data(i)(x(i));
					});
				
				return res;
			}
			

			struct push : array< func::push<F>, 
			                     lie::algebra<Domain>,
			                     lie::algebra<Range>,
			                     M > {

				push(const array& of, const Domain& at) 
					: push::base( of.data.size(), 
					              [&](NN i) {
						              return func::push<F>(of.data(i), at(i));
					              }) {
				}
				
			};

			struct pull : array< func::pull<F>, 
			                     lie::coalgebra<Range>,
			                     lie::coalgebra<Domain>,
			                     M > {

				pull(const array& of, const Domain& at) 
					: pull::base( of.data.size(), 
					              [&](NN i) {
						              return func::pull<F>(of.data(i), at(i));
					              }) {
				}
				
			};

			
		};

	}
}


#endif
