#ifndef MATH0X_EACH_H
#define MATH0X_EACH_H

#include <math0x/types.h>

namespace math0x {

	template<class What, class F>
	void each(const What& v, F&& f) {
		for(NN i = 0, n = v.size(); i < n; ++i) {
			f(i);
		}
	};

}


#endif
