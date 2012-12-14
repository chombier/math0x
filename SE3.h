#ifndef MATH0X_SE3_H
#define MATH0X_SE3_H

#include <math0x/types.h>
#include <math0x/lie.h>

#include <math0x/SO3.h>

namespace math0x { 

	// some handy twist/wrench accessors
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
	

	// 3-dimensional euclidean group
	template<class U>
	struct SE<3, U> {
		typedef SE base;

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
		
		typedef translation_type domain;
		domain operator()(const domain& x) const {
			return rotation(x) + translation;
		}

		// TODO cleaner ? move to traits ?
		typedef vector<U, 6> algebra;
		static void ad_proj(algebra& u, algebra& vv, const algebra& x, const algebra& h) {
			auto omega = angular(x);
			auto v = linear(x);
						
			auto a = angular(h);
			auto b = linear(h);

			algebra res;
			
			typedef vector<U, 3> vec3;
			typedef SO<3, U> SO3;
			
			U theta2 = omega.squaredNorm();
			
			if( std::sqrt( theta2 ) < epsilon<U> ()) {
				// TODO
				throw error("not implemented");
				
			} else {
				vec3 v_x, v_y, tmp;
				SO3::ad_proj(tmp, v_x, omega, a);
				
				U alpha = omega.dot(a) / theta2;
				SO<3, U>::ad_proj(tmp, v_y, omega, b - v.cross(v_x) - alpha * v);
				
				angular(vv) = v_x;
				linear(vv) = v_y;
				
				// TODO
				u = h - typename lie::traits<SE>::ad(x)(vv);
			}
			
		}


		class push {
			func::push<rotation_type> rot;
		public:
			
			push(const SE& of, const domain& at) : rot(of.rotation, at) { }
			push(const SE& of) : rot(of.rotation) { }
			
			domain operator()(const domain& v) const { 
				return rot(v);
			}

		};

		class pull {
			func::pull<rotation_type> rotT;
		public:
			
			pull(const SE& of, const domain& at) : rotT(of.rotation, at) { }
			pull(const SE& of) : rotT(of.rotation) { }
			
			euclid::dual<domain> operator()(const euclid::dual<domain>& f) const { 
				return rotT(f);
			}

		};

	
		

	};


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
			G inv(const G& g) const { return g.inv(); }
			G prod( const G& a, const G& b) const { return a * b; }
			
			euclid::space< algebra > alg() const { return {}; }

			struct Ad { 
				typedef Ad base;
				G at;
				Ad(const G& at) : at(at) { }
				
				algebra operator()(const algebra& x) const {
					algebra res;
					
					angular(res) = at.rotation( angular(x));
					linear(res) = at.translation.cross( angular(res)) + at.rotation(linear(x));
					
					return res;
				}
				
			};

			class AdT {
				func::pull< SO3 > RT;
				vector<U, 3> t;
			public:
				typedef AdT base;
				
				AdT(const G& at)  : RT(at.rotation),
				                    t(at.translation) 
				{ }
				
				coalgebra operator()(const coalgebra& x) const {
					coalgebra res;
					
					linear(res) = RT(linear(x));
					angular(res) = RT(angular(x) - t.cross(linear(x).transpose()).transpose());

					return res;
				}
				
			};
			

			struct ad {
				algebra x;
				
				ad(const algebra& x) : x(x) { }
				
				algebra operator()(const algebra& y) const {
					algebra res;
					
					angular(res) = angular(x).cross( angular(y) );
					linear(res) = linear(x).cross(angular(y)) + angular(x).cross(linear(y));
					
					return res;
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
					algebra x;
					Ad Ad_exp_inv;
					
					push(const exp& of, const algebra& at) 
						: x(at),
						  Ad_exp_inv( of(-at) )
					{ }
					
					algebra operator()(const algebra& h) const {
						
						algebra u, v;
						G::ad_proj(u, v, x, h);
						return u + v - Ad_exp_inv(v);
					}

				};

				struct pull {

					pull(const exp& of, const algebra& at) { }
					
					coalgebra operator()(const coalgebra& h) const {
						throw error("not implemented");
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
						throw error("not implemented");						
					}
					
				};
				
				struct pull {

					pull(const log& of, const G& at) { }
					
					coalgebra operator()(const coalgebra& h) const {
						throw error("not implemented");						
					}
					
				};
				
			};


		};

	}

	namespace func {

		template<class U>
		struct apply< SE<3, U> > {
			
			typedef SE<3, U> G;
			
			typedef std::tuple<G, func::domain<G> > domain;
			typedef func::range<G> range;
		
			apply( const lie::group<G>& = {} ) { }

			range operator()(const domain& x) const {
				return std::get<0>(x)(std::get<1>(x));
			}
      

			struct push {

				domain at;
				
				push( const apply&, const domain& at) : at(at) { }
				
				lie::algebra<range> operator()(const lie::algebra<domain>& dx) const {
					return std::get<0>(at).rotation( angular( std::get<0>(dx)).cross( std::get<1>(at)) 
					                                 + linear( std::get<0>(dx) ) + std::get<1>(dx) );
				}
				
			};
			
			class pull {
				domain at;
				func::pull< SO<3, U> > RT;
			public:
				pull( const apply&, const domain& at) 
					: at(at), RT( std::get<0>(at).rotation ) { }
				
				lie::coalgebra<domain> operator()(const lie::coalgebra<range>& f) const {
					lie::coalgebra<range> RTf = RT(f);
					
					lie::coalgebra<domain> res;

					linear( std::get<0>(res) ) = RTf;
					std::get<1>(res) = RTf;
					
					angular( std::get<0>(res) ) = std::get<1>(at).cross(RTf.transpose()).transpose();
				
					return res;
				}
				
				
			};
      

		};


	}


}





#endif
