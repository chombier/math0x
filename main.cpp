
// #include <math0x/euclid.h>
// #include <math0x/real.h>

// #include <math0x/tuple.h>
// #include <math0x/tuple/stream.h>

// #include <math0x/debug.h>
// #include <math0x/lie.h>

// #include <math0x/vector.h>
// #include <math0x/covector.h>

// #include <math0x/func/push.h>
// #include <math0x/func/pull.h>
// #include <math0x/func/line.h>
// #include <math0x/func/form.h>
// #include <math0x/func/comp.h>

// #include <math0x/func/norm2.h>
// #include <math0x/func/dot.h>

// #include <math0x/func/tie.h>
// #include <math0x/func/tuple.h>
// #include <math0x/func/sum.h>
// #include <math0x/func/minus.h>

// #include <math0x/func/any.h>
// #include <math0x/func/val.h>
// #include <math0x/func/ref.h>

// #include <math0x/func/scal.h>
// #include <math0x/func/poly.h>

// #include <math0x/quaternion.h>
// #include <math0x/SO3.h>
// #include <math0x/SE3.h>
// #include <math0x/vector.h>

// #include <math0x/func/error.h>

// #include <math0x/func/trans.h>
// #include <math0x/func/inv.h>
// #include <math0x/func/prod.h>

// // #include <math0x/func/ops.h>
// #include <math0x/func/val.h>

// #include <math0x/iter.h>
// #include <math0x/minres.h>

// #include <math0x/func/part.h>
// #include <math0x/func/get.h>

// #include <math0x/func/jacobian.h>
// #include <math0x/levmar.h>
// #include <math0x/array.h>


#include <math0x/test/func.h>
#include <math0x/test/lie.h>
#include <math0x/test/euclid.h>

#include <math0x/func/apply.h>
#include <math0x/func/default.h>

#include <math0x/SO3.h>
#include <math0x/SE3.h>
#include <math0x/vector.h>
#include <math0x/tuple.h>
#include <math0x/real.h>



#include <math0x/func/get.h>
#include <math0x/func/part.h>
// #include <math0x/extern.h>

using namespace math0x;
using namespace func;
 
// math0x_lambda( test, [] { return id<RR>(); } );


int main(int, char** ) {
	srand ( time(NULL) );

	lie::group< SO<3> > SO3;
	lie::group< SE<3> > SE3;

	euclid::space<vec3> RR3e;
	lie::group<vec3> RR3g;
 	
	euclid::space< std::tuple<vec3, vec3> > RR3xRR3e;
	lie::group< std::tuple<SO<3>, SE<3>> > SO3xSE3;
	
	// TODO dynamic sized vectors
	
	// 3-vectors
	test::euclid( RR3e );
	
	// rotations
	test::lie( SO3 );
	test::func( make_apply( SO3 ));
	
	// rigids
	test::lie( SE3 );
	test::func( make_apply( SE3 ));
	
	// tuples
	test::euclid( RR3xRR3e );
	test::lie( SO3xSE3 );
	
	// some more tuple tests
	test::func( get< std::tuple<SO<3>, SE<3> >, 0 >{} );
	test::func( get< vec3 >{0} );

	test::func( part< std::tuple<SO<3>, SE<3> >, 0 >{SO3xSE3.id()} );
	test::func( part< vec3 >{vec3::Zero(), 0} );
	
	
  return 0;
}
