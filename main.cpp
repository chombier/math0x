
#include <group/euclid.h>
#include <group/real.h>

#include <group/tuple.h>
#include <group/tuple/stream.h>

#include <group/debug.h>
#include <group/lie.h>

#include <group/vector.h>
#include <group/covector.h>

#include <group/func/push.h>
#include <group/func/pull.h>
#include <group/func/line.h>
#include <group/func/form.h>
#include <group/func/comp.h>

#include <group/func/norm2.h>
#include <group/func/dot.h>

#include <group/func/tie.h>
#include <group/func/tuple.h>
#include <group/func/sum.h>
#include <group/func/minus.h>

#include <group/func/any.h>
#include <group/func/val.h>
#include <group/func/ref.h>

#include <group/func/scal.h>
#include <group/func/poly.h>

#include <group/quaternion.h>
#include <group/SO.h>
#include <group/vector.h>

#include <group/func/error.h>

#include <group/func/trans.h>
#include <group/func/inv.h>
#include <group/func/prod.h>

#include <group/func/ops.h>
#include <group/func/val.h>

int main(int, char** ) {

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

    auto p = 0.5 * (x^2) +  2 * x + func::val<RR>(1.0);
    
    debug( "domain:", typeid( func::domain< decltype(p) > ).name());
    debug( "range:", typeid( func::range< decltype(p) > ).name());
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
  so3.exp()(so3.alg().zero());
  so3.log()(so3.id());

  typedef std::tuple< SO<3> > prout_type;
  
  // lie::log< prout_type  > prout;
  
  // func::domain< lie::log< prout_type  > > p;

  // why do these fail !?
  ryan.log()( ryan.id() );
  ryan.exp()( ryan.alg().zero() );

  // RR c = (*hermite<RR>::ptr)( 1.0 );
  
  return z == E.zero();
}
