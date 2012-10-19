
#include <math0x/euclid.h>
#include <math0x/real.h>

#include <math0x/tuple.h>
#include <math0x/tuple/stream.h>

#include <math0x/debug.h>
#include <math0x/lie.h>

#include <math0x/vector.h>
#include <math0x/covector.h>

#include <math0x/func/push.h>
#include <math0x/func/pull.h>
#include <math0x/func/line.h>
#include <math0x/func/form.h>
#include <math0x/func/comp.h>

#include <math0x/func/norm2.h>
#include <math0x/func/dot.h>

#include <math0x/func/tie.h>
#include <math0x/func/tuple.h>
#include <math0x/func/sum.h>
#include <math0x/func/minus.h>

#include <math0x/func/any.h>
#include <math0x/func/val.h>
#include <math0x/func/ref.h>

#include <math0x/func/scal.h>
#include <math0x/func/poly.h>

#include <math0x/quaternion.h>
#include <math0x/SO3.h>
#include <math0x/SE3.h>
#include <math0x/vector.h>

#include <math0x/func/error.h>

#include <math0x/func/trans.h>
#include <math0x/func/inv.h>
#include <math0x/func/prod.h>

// #include <math0x/func/ops.h>
#include <math0x/func/val.h>

#include <math0x/iter.h>
#include <math0x/minres.h>

#include <math0x/func/part.h>
#include <math0x/func/get.h>

#include <math0x/func/jacobian.h>
#include <math0x/levmar.h>
#include <math0x/array.h>


#include <math0x/test/func.h>
#include <math0x/test/lie.h>

#include <math0x/func/num.h>
#include <math0x/func/apply.h>

namespace michel {
	
	template<bool b>
	using enable_if = typename std::enable_if<b>::type;

	template<class F, class = void> 
	struct get_push {
		typedef math0x::func::impl::default_push<F> type;
	};
	

	template<class F>
	struct get_push<F, enable_if< !std::is_same<typename F::push, F>::value> > {
		typedef typename F::push type;
	};
	

	template<class F>
	using push = typename get_push<F>::type;


}

int main(int, char** ) {
	srand ( time(NULL) );

	using namespace math0x;
	using namespace func;

	typedef SO<3> SO3;
	
	// vectors
	{
		typedef vec3 space;
		typedef euclid::dual<space> dual;
		typedef euclid::field<space> field;
		
		space x = test::random<space>();
		dual f = test::random<dual>();
		
		field lambda = test::random<field>();
		
		test::func< line<space> > (x);
		test::func< form<space> > (f);
		test::func< scal<space> > (lambda); 
		
		test::func< sum<space> > ();
		test::func< minus<space> > ();
		test::func< dot<space> > (); 
		test::func< norm2<space> > (); 

	}

	// rotations
	{
		typedef SO3 G;

		lie::group<G> group;
		
		test::lie( group );
		test::func( make_apply(group));
	}
	
	lie::group< SE<3> > se3;
  vec6 twist = 1.5 * test::random_tangent< SE<3> >() ;
  
  linear(twist) *= 10;
  
  debug( twist.transpose(), (twist - se3.log()( se3.exp()(twist) )).norm() );

  return 0;
}
