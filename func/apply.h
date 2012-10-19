#ifndef MATH0X_FUNC_APPLY_H
#define MATH0X_FUNC_APPLY_H

#include <math0x/types.h>
#include <math0x/lie.h>

namespace math0x {
	namespace func {

		template<class G>
		apply< G > make_apply(const lie::group<G>& group) { return {group}; }
		
	}
}


#endif
