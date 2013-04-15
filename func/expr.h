#ifndef MATH0X_EXPR_H
#define MATH0X_EXPR_H

#include <math0x/macro.h>

#include <math0x/func/sum.h>
#include <math0x/func/minus.h>
#include <math0x/func/scal.h>
#include <math0x/func/tie.h>

#include <math0x/func/comp.h>

// an alternate mechanism for functional expressions using expression
// templates: this one is much more convenient when the range space
// needs contextual information (e.g. dimension)

// typical use:

// using namespace func::expr;
// euclid::space<vec> RRn(14);
// auto f = (a + b + c).in(RRn);

namespace math0x {

	namespace func {

		namespace expr {

			// crtp base for expressions
			template<class E>
			struct base {

				const E& derived() const { 
					return static_cast<const E&>(*this);
				}
				
			};
	
			template<class E>
			E check(const base<E>&& e) { return e.derived(); }

			// some handy shorthand
			template<class E>
			using range = typename E::range;

			// lvalues
			template<class F>
			struct lvalue : base < lvalue<F> > {
				const F& ref;
		
				lvalue(const F& ref) : ref(ref) { }
		
				typedef func::range<F> range;
		
				template<class ... Args>
				const F& in(Args&& ... ) const { 
					return ref;
				}
				
			};

			template<class F>
			lvalue< meta::decay<F> > check(const F& x,
			                               decltype( func::requires< meta::decay<F> >() )* = 0) { 
				return {x}; 
			}
 
			// rvalues
			template<class F>
			struct rvalue : base< rvalue<F> > {
				F value;
				
				rvalue(const F&& value) : value(std::move(value) ) {}
		
				typedef func::range<F> range;

				template<class ... Args>
				const F& in(Args&& ... ) const { 
					return value;
				}
				
			};

			template<class F>
			rvalue< meta::decay<F> > check(const F&& x,
			                               decltype( func::requires< meta::decay<F> >() )* = 0) { return { std::move(x)}; }
	


			// expressions 

			// sum
			template<class LHS, class RHS>
			struct sum : base< sum<LHS, RHS> > {
				
				LHS lhs;
				RHS rhs;

				sum(LHS&& lhs, RHS&& rhs) : lhs(std::move(lhs)),
				                            rhs(std::move(rhs)) { }

				typedef expr::range<RHS> range;
		
				static_assert( std::is_same< expr::range<LHS>, expr::range<RHS> >::value, 
				               "ranges must agree" );
		
				auto in(const euclid::space< expr::range<RHS> >& space = {} ) const 
					-> macro_returns( func::make_sum(space) << func::make_tie(lhs.in(space), 
					                                                          rhs.in(space)) );
		
			};
	
			template<class LHS, class RHS>
			sum< meta::decay<LHS>, meta::decay<RHS> > make_sum(LHS&& lhs, RHS&& rhs) {
				return { std::forward<LHS>(lhs), std::forward<RHS>(rhs) };
			}

		
	
			template<class LHS, class RHS>
			auto operator+(LHS&& lhs, RHS&& rhs) -> 
				macro_returns( make_sum( check( std::forward<LHS>(lhs)), 
				                         check( std::forward<RHS>(rhs)) ) );


			// difference
			template<class LHS, class RHS>
			struct diff : base< diff<LHS, RHS> > {
				
				LHS lhs;
				RHS rhs;

				diff(LHS&& lhs, RHS&& rhs) : lhs(std::move(lhs)),
				                             rhs(std::move(rhs)) { }
		
				typedef expr::range<RHS> range;
		
				static_assert( std::is_same< expr::range<LHS>, expr::range<RHS> >::value, 
				               "ranges must agree" );
		
				// TODO n-ary sums ?
				auto in(const euclid::space< expr::range<RHS> >& space = {} ) const 
					-> macro_returns( func::make_sum(space) << func::make_tie(lhs.in(space), 
					                                                          func::make_minus(space) << rhs.in(space)) );
			};
	
			template<class LHS, class RHS>
			diff< meta::decay<LHS>, meta::decay<RHS> > make_diff(LHS&& lhs, RHS&& rhs) {
				return { std::forward<LHS>(lhs), std::forward<RHS>(rhs) };
			}

			template<class LHS, class RHS>
			auto operator-(LHS&& lhs, RHS&& rhs) -> 
				macro_returns( make_diff( check( std::forward<LHS>(lhs)), 
				                          check( std::forward<RHS>(rhs)) ) );


			// unary minus
			template<class E>
			struct minus : base< minus<E> > {

				E e;

				minus(E&& e) : e(std::move(e)) { }
		
				typedef expr::range<E> range;
		
				auto in(const euclid::space< expr::range<E> >& space = {} ) const 
					-> macro_returns( func::make_minus(space) << e.in(space) );
			};
	
			template<class E>
			minus< meta::decay<E> >  make_minus(E&& e) {
				return { std::forward<E>(e) };
			}

			template<class E>
			auto operator-(E&& e) ->
				macro_returns( make_minus( check( std::forward<E>(e) ) ) );
	
	

			// scalar product
			template<class E>
			struct scal : base< scal<E> > {
		
				typedef expr::range<E> range;

				euclid::field<range> lambda;
				E e;
		
				scal(euclid::field<range> lambda,
				     E&& e) 
					: lambda(lambda),
					  e(std::move(e)) { }
		
				auto in(const euclid::space< expr::range<E> >& space = {} ) const 
					-> macro_returns( func::make_scal(lambda, space) << e.in(space) );
			};
	
			template<class E>
			scal< meta::decay<E> > make_scal(euclid::field< range< meta::decay<E> > > lambda,
			                                 E&& e) {
				return {lambda,  std::forward<E>(e) };
			}

			
				
			template<class Field, class E>
			auto operator*(Field lambda, E&& e) ->
				macro_returns( make_scal(lambda, check( std::forward<E>(e) ) ) );
			
			template<class E, class Field>
			auto operator*(E&& e,  Field lambda) ->
				macro_returns( make_scal(lambda, check( std::forward<E>(e) ) ) );

			template<class E, class Field>
			auto operator/(E&& e, Field lambda) ->
				macro_returns( make_scal(1.0 / lambda, check( std::forward<E>(e) ) ) );
	


			// TODO constants (needs domain ?)

			// TODO products, powers and friends

			
		}

	}

}


#endif
