#ifndef MATH0X_SO3_ALG_H
#define MATH0X_SO3_ALG_H

#include <math0x/types.h>
#include <math0x/vector.h>

namespace math0x {

	template<class U>
	struct so<3, U> {
		typedef so base;
		
		typedef vector<U, 3> twist_type;
		twist_type twist;
		
		so( const twist_type& twist ) : twist(twist) { }
		
		typedef math0x::matrix<U, 3, 3> matrix_type;
		
		matrix_type matrix() const {
			matrix_type res;

			res << 
				0, -twist.z(), twist.y(),
				twist.z(), 0, -twist.x(),
				-twist.y(), twist.x(), 0;
			
			return res;
		}
		
		typedef vector<U, 3> domain;
		typedef vector<U, 3> range;

		range operator()(const domain& x) const {
			return twist.cross(x);
		}

		struct push : so {
			push(const so& of, const domain& ) : push::base(of) { }
			push(const so& of) : push::base(of) { }
		};

		struct pull {
			twist_type twist;
			
			pull(const so& of, const domain& ) : twist(of.twist) { }
			pull(const so& of) : twist(of.twist) { }

			euclid::dual<range> operator()(const euclid::dual<domain>& f) const {
				return -twist.cross(f.transpose()).transpose();
			}
		};
	};

}


#endif
