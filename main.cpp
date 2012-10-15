
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


int main(int, char** ) {
  using namespace math0x;
  
  euclid::space<RR> E;
  
  RR x = E.zero();
  RR z = E.sum(E.minus(x), x);

  typedef tuple::repeat<RR, 2> pair_type;
  
  euclid::space<pair_type> Pair;

  typedef std::tuple< pair_type, pair_type > double_pair_type;
  euclid::space<double_pair_type> DPair;

  double_pair_type dpair = DPair.zero();
  
  DPair.coord(2, dpair) = 1;
  dpair = DPair.sum(dpair, dpair);

  debug( dpair );
  
  vec3 u, v;

  u = vec3::Ones();
  
  func::line<vec3> lu(u);
  
  v = func::d(lu)(1.0)(1.0);
  
  func::id<vec3> id;

  using func::operator<<;
  auto g = id << id << lu;

  auto dg = func::d(g)(1.0);

  RR n = func::norm2<vec3>()( u );
  debug( "bob", n );

  // func::tuple_tie< func::id<RR>, func::id<RR> > h = { std::make_tuple(lu, lu) };

  auto bob = func::make_tie( lu, lu );
  auto marcel = func::make_tuple( lu, lu );
  
  debug( "bob", bob(2.0) );
  
  // auto dTbob = func::dT( bob )( 1.0 );

  func::any<RR, vec3> michel = lu;
  func::any< euclid::dual<vec3>, RR> dTmichel = func::dT( michel )(1.0);
  
  debug("michel", dTmichel( u.transpose() ));
  
  if( michel ) michel.reset();
  
  auto gg = func::val<vec3>(10.0);

  auto ff = func::ref( michel );
  

  SO<3> g1, g2;

  {
    func::id<RR> x;
    
    // auto p = x + x + x;

    // auto p = 0.5 * (x^2) +  2 * x + func::val<RR>(1.0);
    
    // debug( "domain:", typeid( func::domain< decltype(p) > ).name());
    // debug( "range:", typeid( func::range< decltype(p) > ).name());
  }

  // typedef vector< SO<3>, 3 > test_type;
  // test_type test;

  // lie::group< test_type > test_lie;
 
  // test_type zz = test_lie.inv( test );


  // auto log = test_lie.log();
  
  // func::any< lie::alg< test_type >, lie::alg<test_type> > ww = func::d(log)(test);
  
  // // ww( test_lie.alg().zero() );
  
  // // lie::group< std::tuple< SO<3> > > joe;
  // // auto tt = joe.alg().zero();
  // // auto ss = joe.ad( joe.prod( joe.id(), joe.id() ) );


  // lie::adT<SO<3> > adT{ SO<3>() };
  
  // auto mitch = func::make_tuple(adT);
  
  // func::domain< decltype(mitch) > henri;
  // func::range< decltype(mitch) > bidou;
  
  // mitch( henri );
  
  lie::group< std::tuple<SO<3> > > ryan;

  debug(ryan.Ad( ryan.id() )( ryan.alg().zero() ));
  debug(ryan.AdT( ryan.id() )( (*ryan.alg()).zero() ));

  lie::group<SO<3> > so3;

  {
	  SO<3> q = so3.exp()( vec3::UnitY() );
	  
	  debug( so3.log()( q ).transpose() );
  }

  so3.exp()(so3.alg().zero());
  so3.log()(so3.id());

  typedef std::tuple< SO<3> > prout_type;
  
  // lie::log< prout_type  > prout;
  
  // func::domain< lie::log< prout_type  > > p;

  // why do these fail !?
  ryan.log()( ryan.id() );
  ryan.exp()( ryan.alg().zero() );

  debug("iter");
  math0x::iter it(10, 1e-4);

  it([&] {
      return 0.1;
    });

 
  vec3 ex = vec3::UnitX();
  vec3 ey = vec3::UnitY();

  lie::Ad<vec3> zob( vec3::Zero()) ;
  
  zob( vec3::Zero() );

  SO<3> q;
  
  auto map = func::apply< SO<3> >() << func::part<0>(q, ex);
  
  debug("levmar");
  levmar opt;

  opt.outer.bound = 10;
  opt.outer.epsilon = 1e-7;
  
  opt.inner.bound = 10;
  opt.inner.epsilon = 0;
  
  opt.dense(q, map, ey);
  // debug(so3.log()(id3).transpose());
  
  // RR c = (*hermite<RR>::ptr)( 1.0 );
  

  array< lie::Ad<RR> > henri(5, [&](NN  ) {
		  return lie::Ad<RR>(0.0);
	  });
  
  
  array< lie::Ad<RR>, 5 > roger(5, [&](NN i )  {
	   return lie::Ad<RR>(0.0);
  });
  
  
  auto dexp = func::d(ryan.exp())(ryan.alg().zero()) (ryan.alg().zero() );
  
  return z == E.zero();
}
