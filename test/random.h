#ifndef MATH0X_TEST_RANDOM_H
#define MATH0X_TEST_RANDOM_H

#include <math0x/lie.h>
#include <math0x/coords.h>
#include <math0x/vector.h>

namespace math0x {

	namespace test {

		template<class G>
		lie::algebra<G> random_tangent(const lie::group<G>& group = {} ) {
			euclid::coords< lie::algebra<G> > vec;
			
			vec.resize( group.alg().dim() );
			vec.setRandom();
			
			lie::algebra<G> res = group.alg().zero();
			group.alg().set(res, vec);

			return res;
		}
		
		
		template<class G>
		lie::coalgebra<G> random_cotangent(const lie::group<G>& group = {} ) {
			euclid::coords< lie::coalgebra<G> > vec;
			
			auto coalg = *group.alg();
			vec.resize( coalg.dim() );
			vec.setRandom();
			
			lie::coalgebra<G> res = coalg.zero();
			coalg.set(res, vec);
			
			return res;
		}


		template<class G>
		G random(const lie::group<G>& group = {} ) {
			return group.exp()( random_tangent( group ) );
		}
		
	}

}


#endif
