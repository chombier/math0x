#ifndef MATH0X_FUNC_VAL_H
#define MATH0X_FUNC_VAL_H

#include <math0x/lie.h>
#include <math0x/meta.h>

namespace math0x { 
	namespace func {

		// constant value
		template<class Domain, class Range>
		struct value {
			typedef value base;
    
			Range data;
			value(const Range& data) : data(data) { }
			value(Range&& data) : data(std::move(data) ) { }

			const Range& operator()(const Domain& ) const { return data; }
    
			struct push : value< lie::algebra<Domain>, lie::algebra<Range> > {
				push(const value& of, const Domain& )
					: push::base{ lie::group_of(of.data).alg().zero() } {

				}
			};

			struct pull : value< lie::coalgebra<Range>, lie::coalgebra<Domain> > {
				pull(const value& , const Domain& at)
					: pull::base{ (*lie::group_of(at).alg()).zero() } {
	
				}
			};
    
		};


		template<class Domain, class Range>
		value<Domain, meta::decay<Range> > val(Range&& data) { 
			return { std::forward<Range>(data) }; 
		}

	}

}
#endif
