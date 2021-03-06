#ifndef MATH0X_LIE_H
#define MATH0X_LIE_H

#include <math0x/euclid.h>
#include <math0x/meta.h>

namespace math0x {
	namespace lie {
  
		template<class G>
		class group {
			// friend class traits<G>;
			typedef traits<G> impl_type;

		public:
			impl_type impl;
		private:
			void constraints() {
				// TODO more

				G x = id();
				group lie( x );
      
				// static_assert( std::is_same< G, func::range< lie::exp<G> > >::value, "exp range error" );
				// static_assert( std::is_same< lie::alg<G>, func::domain< lie::exp<G> > >::value, "exp domain error" );

				// static_assert( std::is_same< G, func::domain< lie::log<G> > >::value, "log range error" );
				// static_assert( std::is_same< lie::alg<G>, func::range< lie::log<G> > >::value, "exp domain error" );

				// static_assert( std::is_same< lie::alg<G>, func::range< lie::ad<G> > >::value, "ad range error" );
				// static_assert( std::is_same< lie::alg<G>, func::domain< lie::ad<G> > >::value, "ad domain error" );

				// static_assert( std::is_same< lie::coalg<G>, func::range< lie::adT<G> > >::value, "adT range error" );
				// static_assert( std::is_same< lie::coalg<G>, func::domain< lie::adT<G> > >::value, "adT domain error" );
      
				meta::noop(x, lie);
			}
    
		public:
    
			~group() {
				meta::noop( &group::constraints );
			}

			// identity
			G id() const { 
				return impl.id(); 
			}

			// inverse
			G inv(const G& x) const { 
				return impl.inv(x); 
			}

			// move inverse
			G inv(G&& x) const { 
				return impl.inv(std::move(x)); 
			}
    
			// product
			G prod(const G& x, const G& y) const { 
				return impl.prod(x, y); 
			}

			// move products
			G prod(G&& x, const G& y) const { 
				return impl.prod(std::move(x), y); 
			}

			G prod(const G& x, G&& y) const { 
				return impl.prod(x, std::move(y)); 
			}

			G prod(G&& x, G&& y) const { 
				return prod( std::move(x), y);
			}

			// algebra euclidean structure
			euclid::space< lie::algebra<G> > alg() const { return impl.alg(); }
			
			// coalgebra euclidean structure (convenience)
			euclid::space< euclid::dual<lie::algebra<G> > > coalg() const { return *impl.alg(); }
			
			// group adjoint
			lie::Ad<G> Ad(const G& g) const { return {g}; }

			// group adjoint transpose
			lie::AdT<G> AdT(const G& g) const { return {g}; }

			// // algebra adjoint
			// lie::ad<G> ad(const lie::algebra<G>& x) const { return {x, *this}; }
			
			// // algebra adjoint transpose
			// lie::adT<G> adT(const lie::algebra<G>& x) const { return {x, *this}; }
		
			lie::algebra<G> bracket(const lie::algebra<G>& x, 
			                        const lie::algebra<G>& y) const { 
				return impl.bracket(x, y); 
			}
			
	
			// exponential
			lie::exp<G> exp() const { return {*this}; }

			// logarithm
			lie::log<G> log() const { return {*this}; }
    
			// forward ctor
			template<class ... Args>
			explicit group(Args&& ... args) : impl(std::forward<Args>(args)...) {  }
    
			// default ctor in case it exists
			group() {} 

			// 
			group(const group& ) = default;
			group(group&& ) = default;
    
			// impl must provide this ctor
			group(const G& g) : impl(g) { }
		};


		template<class G>
		group<G> group_of(const G& of) { return {of}; }

	}
}

#endif
