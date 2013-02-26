#ifndef MATH0X_FUNC_ARRAY_H
#define MATH0X_FUNC_ARRAY_H

#include <math0x/func.h>
#include <math0x/array.h>
#include <math0x/each.h>

namespace math0x {
	namespace func {

		template<class F, class Domain, class Range, int M>
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

		
		template<class F, class Range, int M>
		struct array_tie {
			typedef array_tie base;
			
			math0x::array<F, M> data;
			
			template<class Fun>
			array_tie(NN size, const Fun& fun) 
				: data(size, fun) { }
			
			
			Range operator()(const domain<F>& x) const {
				Range res; res.resize( data.size() );
				
				each(data, [&](NN i) {
						res(i) = data(i)(x);
					});
				
				return res;
			}
			

			struct push : array_tie< func::push<F>, 
			                         lie::algebra<Range>,
			                         M > {
				
				push(const array_tie& of, const domain<F>& at) 
					: push::base( of.data.size(), 
					              [&](NN i) {
						              return func::push<F>(of.data(i), at(i));
					              }) {
				}
				
			};
			
			
			struct pull {
				
				math0x::array< func::pull<F>, M > data;
				euclid::space< lie::coalgebra<domain<F> > > coalg;
				
				pull(const array_tie& of, const domain<F>& at) 
					: data( of.data.size(), 
					        [&](NN i) {
						        return func::pull<F>(of.data(i), at(i));
					        }),
					  coalg( *lie::group_of(at).alg() )
				{
				}
				
				
				lie::coalgebra<domain<F> > operator()(const lie::coalgebra<Range>& x) const {
					lie::coalgebra<domain<F> > res = coalg.zero();

					for(unsigned i = 0, n = data.size(); i < n; ++i) {
						res = coalg.sum(res, data(i)(x(i)));
					}
					
					return res;
				}
				
			};
			
			
		};





	}
}


#endif
