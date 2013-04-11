#ifndef MATH0X_LEVMAR_H
#define MATH0X_LEVMAR_H

#include <math0x/func.h>
#include <math0x/iter.h>
#include <math0x/coords.h>

#include <math0x/func/coords.h>
#include <math0x/func/trans.h>
#include <math0x/func/jacobian.h>
#include <math0x/minres.h>

#include <math0x/func/push.h>
#include <math0x/func/pull.h>

#include <math0x/debug.h>

namespace math0x {

	// levenberg-marquardt
	struct levmar {
    
		iter outer;			// iterations for gauss-newton (outer loop)
		iter inner;			// iterations for minres (inner loop)
    
		RR lambda;			
		RR nu;
		
		// if lambda = 0 then gauss-newton
		// else if nu = 1 then fixed-damping levmar
		// else adaptive-damping levmar 
		levmar(RR lambda = 0, RR nu = 1) 
			: lambda(lambda), nu(nu) {
			
			assert( nu >= 1);
			
		}
    
		// return true when ok
		bool adaptive(real& llambda, real prev, real curr) const {
			if( curr < prev ) {
				llambda /= nu;
				return true;
			}
			else {
				llambda *= nu;
				return false;
			}
		}
		
		template<class F>
		struct data_type {

			typedef func::domain< F > domain;
			typedef func::range< F > range;
      
			typedef lie::algebra< domain > domain_algebra;
			typedef lie::algebra< range > range_algebra;
			
			typedef euclid::coords<domain_algebra> domain_coords;
			typedef euclid::coords<range_algebra> range_coords;
			
			// structures
			lie::group< domain > dmn;
			lie::group< range > rng;

			euclid::space< domain_algebra > dmn_alg;
			euclid::space< range_algebra > rng_alg;

			// functions
			typedef func::comp< func::comp< lie::log<range>, func::left<range> >, F > residual_type;
			residual_type residual;
			
			lie::exp<domain> exp;
			
			// vars
			data_type(const domain& x,
			          const F& f,
			          const range& y) 
				: dmn(x),
				  rng(y),
				  dmn_alg(dmn.alg()),
				  rng_alg(rng.alg()),
				  residual{ rng.log() << func::L(rng.inv(y)), f},
				  exp( dmn.exp() )
				  
			{

			}

			static constexpr int domain_dim = euclid::space<domain_algebra>::static_dim;
			
			// convenience
			void set(domain_algebra& delta, const domain_coords& dx) const {
				dmn_alg.set(delta, dx);
			}

			void get(range_coords& b, const range_algebra& r ) const {
				rng_alg.get(b, r);
			}

			void step(domain& x, const domain_algebra& delta) const {
				x = dmn.prod(x, exp(delta) );
			}
			
		};


		// assembling jacobian matrices, non-homogeneous damping
		// needs only dlog on range space
		template<class F>
		iter dense(func::domain< meta::decay<F> >& x,
		           const F& f,
		           const func::range< meta::decay<F> >& y) const {
			data_type<F> data(x, f, y);

			real llambda = lambda;
			
			// vars
			auto r = data.residual( x );

			auto Jf = func::J(data.residual, data.dmn, lie::group_of(r) );
			
			auto delta = data.dmn_alg.zero();
			
			real last = data.rng_alg.norm(r);
			
			// jacobian matrix and JTJ's diagonal
			typename func::jacobian<F>::range J;
			typename data_type<F>::domain_coords diag;
			
			// temporaries
			typedef typename data_type<F>::domain_coords domain_coords;
			domain_coords dmn_tmp, rhs, dx;
			
			typename data_type<F>::range_coords rng_tmp = vec::Zero( data.rng_alg.dim() );
			
			// tangent normal equations
			minres< data_type<F>::domain_dim > normal;
			normal.iter = inner;
			
			auto JTJ = [&](const domain_coords& u) -> const domain_coords& {
				rng_tmp.noalias() = J * u;
				dmn_tmp.noalias() = J.transpose() * rng_tmp;
				
				if( llambda ) dmn_tmp.array() += llambda * (diag.array() * u.array());
				
				return dmn_tmp;
			};

			return outer( [&] ()-> RR  {
	  
					J = Jf(x);
					if( llambda ) diag = J.colwise().squaredNorm().transpose();					
					
					data.get(rng_tmp, r);
					rhs.noalias() = -J.transpose() * rng_tmp;
	  
					// don't forget to clear dx !
					dx.setZero();
					normal.iter.epsilon = inner.epsilon * rhs.norm();
					normal.solve(dx, JTJ, rhs);
					
					data.set(delta, dx);
					
					// TODO wrap this in data as it's the same for sparse/dense
					auto old = x;
					data.step(x, delta );
					r = data.residual(x);
	  
					real norm = data.rng_alg.norm(r);
					
					if( !adaptive( llambda, last, norm ) ) {
						x = old;
					} else {
						last = norm;
					}

					if( rhs.norm() <= outer.epsilon ) return outer.epsilon;
					return norm;
				});
			
			
		}


		template<class F, class DF>
		static void update(typename data_type<F>::domain_coords& res,
		                   DF df,
		                   const data_type<F>& data) {
			typename data_type<F>::domain_algebra unit = data.dmn_alg.zero();

			// omp_each(res, [unit, df, &res, &data](NN i) mutable {
			// 		data.dmn_alg.coord(i, unit) = 1;
			// 		res(i) = data.rng_alg.norm2( df(unit) );
			// 		data.dmn_alg.coord(i, unit) = 0;
			// 	});

			auto task = [unit, df, &res, &data](NN i) mutable {
				data.dmn_alg.coord(i, unit) = 1;
				res(i) = data.rng_alg.norm2( df(unit) );
				data.dmn_alg.coord(i, unit) = 0;
			};
			
			unsigned n = res.size();

#pragma omp parallel for firstprivate(task)
			for(unsigned i = 0; i < n; ++i) task(i);
			
			
		};
		                 

		// without jacobian matrices assembly, but needs both dlog/dlogT on
		// range space.
		template<class F>
		iter sparse(func::domain<F>& x,
		            const F& f,
		            const func::range<F>& y) const {

			data_type<F> data(x, f, y);
				
			real llambda = lambda;
			
			auto r = data.residual(x);
			auto delta = data.dmn_alg.zero();
			
			real best = data.rng_alg.norm( r );
			real last = best;
			
			// temporaries
			typename data_type<F>::domain_coords dmn_tmp, rhs, dx, diag;
			typename data_type<F>::range_coords rng_tmp;
			
			dx.resize( data.dmn_alg.dim() );
			rhs.resize( data.dmn_alg.dim() );
			rng_tmp.resize( data.rng_alg.dim() );
			diag.setOnes( data.dmn_alg.dim() );
			
			// tangent normal equations
			minres< data_type<F>::domain_dim > normal;
			normal.iter = inner;
			
			typedef decltype( data.residual ) residual_type;
			typedef func::push< residual_type > residual_push;
			typedef func::pull< residual_type > residual_pull;

			// TODO lots of allocs !
			auto JTJ = [&](const func::domain<F>& at) {
				using namespace func;
				return  				
				func::make_coords( residual_pull(data.residual, at), *data.rng_alg, *data.dmn_alg) << 
				func::make_coords( residual_push(data.residual, at), data.dmn_alg, data.rng_alg);
			};
			
			typename data_type<F>::domain_algebra unit = data.dmn_alg.zero();
			
			return outer( [&] ()-> RR  {
					if(llambda) this->update(diag, residual_push(data.residual, x), data);
					
					auto A = JTJ(x);
					
					vec storage;
					
					// TODO alloc
					auto AA = [&] (const vec& x) -> const vec& {
						storage.noalias() = A(x) + llambda * diag.cwiseProduct(x);
						return storage;
					};
					
					data.get(rng_tmp, r);
					rhs = -func::make_coords( residual_pull(data.residual, x), *data.rng_alg, *data.dmn_alg)(rng_tmp);
					
					dx.setZero();
					normal.solve(dx, AA, rhs);
					
					data.set(delta, dx);

					auto old = x;
					data.step(x, delta );
					
					r = data.residual(x);
					
					real norm = data.rng_alg.norm(r);
					
					if( !adaptive(llambda, last, norm) ) {
						x = old;
					} else {
						last = norm;
					}
					
					if( rhs.norm() <= outer.epsilon ) return outer.epsilon;
					return norm;
				});
			
		}



		// TODO provide quick and dirty one (without dlog, fixed damping)

		



		
	};


}


#endif
