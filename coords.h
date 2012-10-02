#ifndef MATH0X_COORDS_H
#define MATH0X_COORDS_H

#include <math0x/euclid.h>
#include <math0x/vector.h>

namespace math0x {
	namespace euclid {

		// coordinate type for an euclidean space, based on static_dim
		template<class E>
		using coords = vector< field<E>, space<E>::static_dim >;

	}
}

#endif
