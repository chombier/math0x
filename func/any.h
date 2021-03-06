#ifndef MATH0X_FUNC_ANY_H
#define MATH0X_FUNC_ANY_H

#include <math0x/lie.h>
#include <memory>

// #include <group/debug.h>


namespace math0x { 
	namespace func {

		// type erasure class
		template<class Domain, class Range>
		class any {

			typedef any< lie::algebra<Domain>, lie::algebra<Range> > push_type;
			typedef any< lie::coalgebra<Range>, lie::coalgebra<Domain> > pull_type;
    
    
			struct interface {
				virtual ~interface() { }

				virtual Range operator()(const Domain& x) const = 0;
				virtual Range operator()(Domain&& x) const = 0;
      
				virtual push_type push(const Domain& at) const = 0;
				virtual pull_type pull(const Domain& at) const = 0;
      
				virtual interface* copy() const = 0;
			};

			template<class F>
			struct impl : interface {
				F value;
      
				static_assert( std::is_same< Domain, func::domain<F> >::value, "bad domain" );
				static_assert( std::is_same< Range, func::range<F> >::value, "bad range" );
      
				impl(const F& value ) : value(value) { }
				impl(F&& value ) : value( std::move(value) ) { }
      
				Range operator()(const Domain& x) const { return value(x); }
				Range operator()(Domain&& x) const { return value( std::move(x) ); }
      
				push_type push(const Domain& at) const { return func::push<F>(value, at); }
				pull_type pull(const Domain& at) const { return func::pull<F>(value, at); }
      
				interface* copy() const { return new impl(value); }
			};

			typedef std::unique_ptr<interface> ptr_type;
			ptr_type ptr;
    
			const interface& get() const { 
				assert( ptr ); 
				return *ptr;
			}

			template<class F>
			static ptr_type make_ptr(F&& f) {
				return ptr_type(  new impl< meta::decay<F> >( std::forward<F>(f) ) );
			}

		
			
		public:
			typedef any base;

			any copy() const {
				// debug("func::any::copy");
				any res;
				res.ptr.reset( ptr ? ptr->copy() : 0 );
				return res;
			}
    
			Range operator()(const Domain& x) const { 
				return get()(x);
			}

			Range operator()(Domain&& x) const { 
				return get()( std::move(x) );
			}

    
			any( const any& other) :  any( other.copy() ) { }
			any& operator=(const any& other) {
				return (*this = other.copy());
			}

			any() { }

			template<class F>
			any(F&& f) : ptr( make_ptr( std::forward<F>(f) ) )  {}
    
			
			any( any&& ) = default;
    
			any& operator=(any&& ) = default;

			template<class F>
			any& operator=(F&& f) {
				ptr = make_ptr( std::forward<F>(f) );
				return *this;
			}
    
			explicit operator bool() const {
				return ptr.get();
			}
    
			void reset() { ptr.reset(); }

			struct push : push_type {
      
				push(const any& of, const Domain& at) 
					: push::base(of.get().push(at)) { }
      
			};

			struct pull : pull_type {
      
				pull(const any& of, const Domain& at) 
					: pull::base(of.get().pull(at)) { }
      
			};


		};

	}


}
#endif
