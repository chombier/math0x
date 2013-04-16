#ifndef MATH0X_EIGEN_H
#define MATH0X_EIGEN_H

#include <math0x/euclid.h>
#include <math0x/lie.h>
#include <math0x/error.h>

#include <math0x/func/array.h>
#include <math0x/each.h>
#include <Eigen/Core>

// adapts eigen vectors for use with this library.
namespace math0x {

	template<int K>
	NN ensure_static() {
		static_assert( K != -1, "should not be called for dynamic size" );
		return K;
	}
	
	

	// template<class Vector, class F>
	// Vector map(F&& f,
	//            int size = ensure_static< Vector::SizeAtComileTime >() ) {
	// 	Vector res;
	// 	res.resize( size );
		
	// 	each(res, [&](NN i ) {
	// 			res(i) = f(i);
	// 		});
		
	// 	return res;
	// }
	

	namespace euclid {
 
		template<class U, int M, int N>
		struct traits< Eigen::Matrix<U, M, N> > {

			static_assert( (M == 1) || (N == 1),
			               "sorry, unimplemented: at least one dimension must be 1");
    
			NN n;
			space<U> sub;

			typedef euclid::field<U> field;
			typedef Eigen::Matrix<U, M, N> E;

    
			traits(NN n = ensure_static< E::SizeAtCompileTime >(),
			       const space<U>& sub = space<U>())
				: n(n),
				  sub(sub) {
				assert( dim() );
			}
    
    
			traits( const E& x ) 
			: n( x.size() ),
			  sub( x(0)) {
				assert( dim() );
			}
    
			// dual space: swaps dimensions
			typedef Eigen::Matrix< euclid::dual<U>, N, M> dual;
			
			NN dim() const { return n * sub.dim(); }
    
			E zero() const {
				E res;
				res.resize( n );

				each(res, [&](NN i) {
						res(i) = sub.zero();
					});
					
				return res;
			}
    
			field& coord(NN i, E& x) const {
				return sub.coord(i % sub.dim(), x( i / sub.dim() ) );
			}
    
			const field& coord(NN i, const E& x) const {
				return sub.coord(i % sub.dim(), x( i / sub.dim() ) );
			}
    
			space< dual > operator*() const {
				return space<dual>(n, *sub);
			}
    

			static constexpr int static_dim = M == -1 ? -1 
				: N == -1 ? -1 
				: M * N;
    
			// // prevent nested dynamic vectors
			// static_assert( static_dim > 0 || euclid::space<U>::static_dim > 0,
			//                "can't have nested dynamic vectors (yet) for efficiency reasons" );


			class range {
				NN outer;
				const E* data;
				
				euclid::range<U> sub;
			public:
				
				range(const E& data) 
				: outer(0),
				  data(&data),
				  sub(data( outer )) { 
					
				}
				
				bool empty() const { return outer >= data->size(); }
				
				const field& front() const {
					return sub.front();
				}
				
				void pop() {
					sub.pop();
					if( sub.empty() ) {
						++outer;
						
						if( !empty() ) sub = euclid::range<U>((*data)(outer));
					}
				}
				
			};

		};

  
	}


	namespace lie {
  
		// TODO termwise operations ?
		template<class U, int M, int N>
		struct traits< Eigen::Matrix<U, M, N> > {

			static_assert( (M == 1) || (N == 1),
			               "sorry, unimplemented: at least one dimension must be 1");
    
			NN n;
			group<U> sub;

			typedef Eigen::Matrix< lie::algebra<U>, M, N > algebra;
			typedef euclid::dual< algebra > coalgebra;

			typedef Eigen::Matrix<U, M, N> G;
			
			traits(NN n = ensure_static< G::SizeAtCompileTime >(),
			       const group<U>& sub = group<U>())
				: n(n),
				  sub(sub) {
				assert( n );
			}

	  
			traits( const G& x ) 
			: n( x.size() ),
			  sub( x(0)) {
				assert( n );
			}
    
			
			G id() const {
				G res; res.resize(n);
				each(res, [&](NN i) { res(i) = sub.id(); } );
				return res;
			}
    
			G inv( const G& x ) const { 
				assert( x.size() == n );
				
				G res; res.resize( n );
				each(res, [&](NN i) { res(i) = sub.inv(x(i)); } );
				return res;
			}

			G inv( G&& x ) const { 
				assert( x.size() == n );
				
				each(x, [&](NN i) { x(i) = sub.inv( std::move(x(i))); } );
				return x;
			}
			
			G prod( const G& x, const G& y) const { 
				assert( x.size() == n && y.size() == n );
				G res; res.resize( n );
				each(res, [&](NN i) { res(i) = sub.prod(x(i), y(i)); } );
				return res;
			}

			G prod( G&& x, const G& y) const { 
				assert( x.size() == n && y.size() == n );

				each(x, [&](NN i) { x(i) = sub.prod( std::move(x(i)), y(i)); } );
				return x;
			}

			G prod( const G& x, G&& y) const { 
				assert( x.size() == n && y.size() == n );

				each(y, [&](NN i) { y(i) = sub.prod( x(i), std::move(y(i))); } );
				return y;
			}



			euclid::space< lie::algebra<G> > alg() const { 
				return euclid::space<lie::algebra<G> >(n, sub.alg() ); 
			}
			

			algebra bracket(const algebra& x, const algebra& y) const {
				algebra res; res.resize(n);
				each(res, [&](NN i) { res(i) = sub.bracket(x(i), y(i)); });
				return res;
			}

			coalgebra cobracket(const coalgebra& x, const coalgebra& y) const {
				coalgebra res; res.resize(n);
				each(res, [&](NN i) { res(i) = sub.cobracket(x(i), y(i)); });
				return res;
			}

			struct Ad : func::array< lie::Ad< U >, 
			                         algebra,
			                         algebra,
			                         G::SizeAtCompileTime> {
				Ad(const G& at) 
				: Ad::base(at.size(), [&](NN i) {
						return lie::Ad< U >(at(i));
					}) {
					
				}
				
			};

			struct AdT : func::array< lie::AdT<U> , 
			                          coalgebra,
			                          coalgebra,
			                          G::SizeAtCompileTime> {
				AdT(const G& at) 
				: AdT::base(at.size(), [&](NN i) {
						return lie::AdT< U >(at(i));
					}) {
					
				}
				
			};
    
		

			struct exp : func::array< lie::exp<U>,
			                          algebra,
			                          G,
			                          G::SizeAtCompileTime > {

				exp( const group<G>& g) 
					: exp::base(g.impl.n, [&](NN ) {
							return lie::exp<U>(g.impl.sub);
						}) { 
					
				}

			};

			struct log : func::array< lie::log<U>,
			                          G,
			                          algebra,
			                          G::SizeAtCompileTime > {

				log( const group<G>& g) 
					: log::base(g.impl.n, [&](NN ) {
							return lie::log<U>(g.impl.sub);
						}) { 
					
				}

			};

		
    
		};
		
	} 
}
#endif
