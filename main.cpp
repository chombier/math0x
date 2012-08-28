
#include <group/euclid.h>
#include <group/real.h>

#include <group/tuple.h>
#include <group/tuple/stream.h>

#include <group/debug.h>
#include <group/lie.h>

#include <group/vector.h>
#include <group/eigen.h>

#include <group/func/push.h>
#include <group/func/pull.h>
#include <group/func/line.h>
#include <group/func/form.h>
#include <group/func/comp.h>

#include <group/func/norm2.h>
#include <group/func/dot.h>

#include <group/func/tie.h>
#include <group/func/sum.h>
#include <group/func/minus.h>

#include <group/func/any.h>
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

 
  
  typedef Eigen::Vector3d vec3;

  vec3 u, v;

  u = vec3::Ones();
  
  func::line<vec3> lu(u);
  
  v = func::d(lu)(1.0)(1.0);
  
  func::id<vec3> id;

  auto g = id << id << lu;

  auto dg = func::d(g)(1.0);

  RR n = func::norm2<vec3>()( u );
  debug( "bob", n );

  // func::tuple_tie< func::id<RR>, func::id<RR> > h = { std::make_tuple(lu, lu) };

  auto bob = func::tie( lu, lu );

  debug( "bob", bob(2.0) );
  
  // auto dTbob = func::dT( bob )( 1.0 );

  func::any<RR, vec3> michel = lu;
  func::any< euclid::dual<vec3>, RR> dTmichel = func::dT( michel )(1.0);
  
  debug("michel", dTmichel( u.transpose() ));
  
  if( michel ) michel.reset();
  
  auto gg = func::val<vec3, RR>(10.0);

  return z == E.zero();
}
