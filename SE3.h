#ifndef MATH0X_SE3_H
#define MATH0X_SE3_H

#include <math0x/types.h>
#include <math0x/lie.h>

namespace math0x { 

	// 3-dimensional euclidean group
	template<class U>
	struct SE<3, U> {
		typedef vector<U, 3> translation_type;
		translation_type translation;

		typedef SO<3, U> rotation_type;
		rotation_type rotation;
		
		SE(const translation_type& translation = translation_type::Zero(),
		   const rotation_type& rotation = {} ) 
			: translation( translation ),
			  rotation( rotation ) {

		}
		
		SE operator*(const SE& other) const {
			return { translation + rotation( other.translation ),
					rotation * other.rotation };
		}
		
		SE inv() const {
			SO<3, U> rot_inv = rotation.inv();
			return { -rot_inv(translation), rot_inv};
		}
		
	};


	namespace lie {

		// TODO traits
		template<class U>
		struct traits< SE<3, U> > {

			typedef vector<6, U> algebra;
			typedef euclid::dual< algebra > coalgebra;
			
			


		};

	}




}





#endif
