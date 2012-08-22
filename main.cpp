
#include <group/euclid.h>
#include <group/real.h>

int main(int, char** ) {

  euclid::space<RR> E;

  RR x = E.zero();
  RR z = E.sum(E.minus(x), x);

  return z == E.zero();
}
