#ifndef MATH0X_TEST_LIE_H
#define MATH0X_TEST_LIE_H

#include <math0x/test/func.h>

// TODO includes !
#include <math0x/func/subgroup.h>

namespace math0x {
	namespace test {

		template<class G>
		void lie( const math0x::lie::group<G>& group) {
			
			debug("testing Lie group structure for", meta::name<G>());
			
			G h = test::random( group );		
			math0x::lie::algebra<G> v = test::random_tangent( group );
			
			math0x::lie::group< std::tuple<G, G> > pair = std::make_tuple( group, group );

			using namespace func;
			test::func( make_prod(group), pair, group );
			test::func( make_inv(group), group, group );
			
			test::func( L(h), group, group );
			test::func( R(h), group, group );

			test::func( subgroup<G>(v, group), {}, group );
			
			math0x::lie::group< lie::algebra<G> > alg( group.alg().zero() );
			
			test::func( group.exp(), alg, group );
			test::func( group.log(), group, alg );
			
		};
		
	}
}


#endif
