#ifndef MATH0X_SE3_H
#define MATH0X_SE3_H

#include <math0x/types.h>
#include <math0x/lie.h>

#include <math0x/SO3.h>

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

	template<class Array6>
	auto angular(Array6&& x) -> decltype( std::forward<Array6>(x).template head<3>() ) {
		static_assert( meta::decay<Array6>::SizeAtCompileTime == 6, "only for 6 dimensional" );
		return std::forward<Array6>(x).template head<3>();
	}
	
	template<class Array6>
	auto linear(Array6&& x) -> decltype( std::forward<Array6>(x).template tail<3>() ) {
		static_assert( meta::decay<Array6>::SizeAtCompileTime == 6, "only for 6 dimensional" );
		return std::forward<Array6>(x).template tail<3>();
	}
	
	namespace lie {

		template<class U>
		struct traits< SE<3, U> > {

			typedef vector<U, 6> algebra;
			typedef euclid::dual< algebra > coalgebra;
			
			typedef SE<3, U> SE3;
			typedef SO<3, U> SO3;

			typedef SE3 G;
			
			traits(const SE3& )  { }
			traits() {  }

			G id() const { return {}; }
			G inv(const G& g) { return g.inv(); }
			G prod( const G& a, const G& b) const { return a * b; }
			
			euclid::space< algebra > alg() const { return {}; }

			struct Ad {
				Ad(const G& at) {

				}

				algebra operator()(const algebra& x) const {

				}
				
			};

			struct AdT {
				AdT(const G& at) {

				}

				coalgebra operator()(const coalgebra& x) const {
					
				}
				
			};
			

			struct exp {
				exp(const group<G>& = {}) { }
				
			
				
				G operator()( const algebra& x ) const {
					typedef lie::exp< SO3 > sub_type;
					
					typename sub_type::push dsub(sub_type{}, angular(x));
					
					SO3 rotation = dsub.value_inv().inv();
					lie::algebra< SO3 > translation = rotation( dsub( linear(x) ) );
					
					return {translation, rotation};
				}
				
				struct push {

					push(const exp& of, const algebra& at) { }

					algebra operator()(const algebra& h) const {

					}

				};

				struct pull {

					pull(const exp& of, const algebra& at) { }
					
					coalgebra operator()(const coalgebra& h) const {
						
					}
					
				};


			};

			struct log {
				log(const group<G>& = {} ) { }

				algebra operator()( const G& g ) const {
					typedef lie::log< SO3 > sub_type;
					
					typename sub_type::push dsub(sub_type{}, g.rotation);
					
					algebra res;
					
					angular(res) = dsub.value();
					linear(res) = dsub( g.rotation.inv()( g.translation ) );
					
					return res;
				}

				struct push {
				
					
					push(const log& of, const G& at) { }
					
					algebra operator()(const algebra& h) const {
					

						
					}
					
				};
				
				struct pull {

					pull(const log& of, const G& at) { }
					
					coalgebra operator()(const coalgebra& h) const {
						
					}
					
				};
				
			};


		};

	}




}





#endif
