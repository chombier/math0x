#ifndef MATH0X_EUCLID_H
#define MATH0X_EUCLID_H

#include <math0x/types.h>
#include <math0x/meta.h>

#include <cmath>
#include <cassert>

namespace math0x {
	namespace euclid {
  
		template<class E>
		class space {
			typedef traits<E> impl_type;
			traits<E> impl;

			void constraints() {
				// TODO more
      
				E x = zero();
				space euclid(x);

				meta::noop(x, euclid);

				static_assert( (static_dim > 0) || (static_dim == -1),
				               "static_dim must be positive or -1" );

			}    
    
		public:

			~space() {
				meta::noop( &space::constraints );
			}
      
			template<class...Args>
			explicit space(Args&& ... args) 
				: impl(std::forward<Args>(args)...) { }

			space() { }
    
			space(const space& ) = default;
			space(space&& ) = default;
			space(const E& x) : impl(x) { }

			// dual geometry
			space< dual<E> > operator*() const { 
				return *impl;
			}

			// coordinate accessors 
			field<E>& coord(NN i, E& e) const { 
				return impl.coord(i, e); 
			}
    
			const field<E>& coord(NN i, const E& e) const { 
				return impl.coord(i, e); 
			}
    
			// space dimension
			NN dim() const { return impl.dim(); }

			// zero element
			E zero() const { return impl.zero(); }
    
			// dimension if known at compile-time
			static constexpr int static_dim = impl_type::static_dim;
    
			// coordinate iterators
			template<class F>
			void each(const E& x, const F& f) const {
				for(natural i = 0, n = dim(); i < n; ++i) {
					f( coord(i, x) );
				}
			};

			template<class F>
			void each(E& x, const F& f) const {
				for(natural i = 0, n = dim(); i < n; ++i) {
					f( coord(i, x) );
				}
			};


			// vector space operations : TODO this is useless with scal
			E minus(E&& x) const {
				each(x, [&](field<E>& xi) {
						xi = -xi;
					});
      
				return x;
			} 
    

			E minus(const E& res) const {
				return minus( E(res) );
			} 

			E scal(field<E> lambda, E&& x) const {
				each(x, [&](field<E>& xi) {
						xi = lambda * xi;
					});
				
				return x;
			}
    
			E scal(field<E> lambda, const E& res) const {
				return scal(lambda, E(res) );
			}
    

			E sum(E&& x, const E& y) const {
				NN i = 0;

				each(x, [&](field<E>& xi) {
						xi += this->coord(i, y);
						++i;
					});
      
				return x;
			}

			E sum(const E& x, E&& y) const {
				return sum(std::move(y), x);
			}

			E sum(E&& x, E&& y) const {
				return sum(std::move(x), y);
			}

			E sum(const E& x, const E& y) const {
				return sum(E(x), y);
			}

			// convenience
			E diff(const E& x, const E& y) const {
				return sum(x, minus(y) );
			}

			E diff(const E& x, E&& y) const {
				return sum(x, minus( std::move(y) ) );
			}
    
			// convenience
			template<class Vector>
			void get(Vector&& v, const E& x) const {
				assert(v.size() == dim());
      
				for(NN i = 0, n = dim(); i < n; ++i) {
					v(i) = coord(i, x);
				}
      
			}

			template<class Vector>
			void set(E& x, Vector&& v) const {
				assert(v.size() == dim());
      
				for(NN i = 0, n = dim(); i < n; ++i) {
					coord(i, x) = v(i);
				}
      
			}


			// canonical dot-product
			field<E> dot(const E& x, const E& y) const {
				field<E> res = 0;
				
				for(NN i = 0, n = dim(); i < n; ++i) {
					res += coord(i, x) * coord(i, y);
				}

				return res;
			}

			// canonical squared norm
			field<E> norm2(const E& x) const {
				field<E> res = 0;
				each(x, [&](const field<E>& xi) {
						res += xi * xi;
					});
				return res;
			}

			// canonical norm
			field<E> norm(const E& x) const {
				return std::sqrt( norm2(x) );
			}


			// canonical representation
			dual<E> transpose(const E& x) const {
				space< dual<E> > d = **this;
				
				dual<E> res = d.zero();
				for(NN i = 0, n = dim(); i < n; ++i) {
					d.coord(i, res) = coord(i, x);
				}

				return res;
			}
			
		};

		template<class E>
		space<E> space_of(const E& of) { return {of}; }
		
	}
}

#endif
