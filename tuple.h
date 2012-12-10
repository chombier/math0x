#ifndef MATH0X_TUPLE_H
#define MATH0X_TUPLE_H

#include <math0x/euclid.h>

#include <math0x/tuple/range.h>
#include <math0x/tuple/element.h>
#include <math0x/tuple/find.h>

#include <cassert>

namespace math0x { 
	namespace euclid {
  

		namespace impl {

			template<class ... Args> struct static_dim;
    
			template<class Head, class ... Tail> 
			struct static_dim<Head, Tail...> {
				static constexpr int head = euclid::space<Head>::static_dim;
				static constexpr int tail = static_dim<Tail...>::value;
      
				static constexpr int value = head == -1 ? -1
					: tail == -1 ? -1 
					: head + tail;
      
			}; 

			template<> struct static_dim<> { 
				static constexpr int value = 0;
			};
    

		}

	
		template<class ... Args>
		struct traits< std::tuple<Args...> > {
    
			// iterator type
			typedef tuple::range< sizeof...(Args) > each;
			
			static_assert( (each::N == sizeof...(Args)) , "derp" );
			
			typedef std::tuple< space<Args>... > impl_type;
			impl_type impl;
    
			static constexpr int size = sizeof...(Args);

			typedef tuple::repeat<NN, size > offset_type;
			offset_type offset;
    
			NN dimension;

			static constexpr int static_dim = impl::static_dim<Args...>::value;
    
			typedef std::tuple<Args...> E;

			template<int I>
			using type = math0x::tuple::element<I, E>;

			// choose first as field TODO check consistency
			typedef euclid::field< type<0> > field;
			typedef std::tuple< euclid::dual<Args>... > dual;
    
    
			// some helper functors for implementing operations
			struct helper {

				// term-wise euclidean structure
				struct impl {
	
					const E& value;

					template<int I>
					euclid::space< type<I> > operator()() const {
						return { std::get<I>(value) };
					}
	
				};

				// term-wise dual-space
				struct dual {
					const impl_type& impl;

					template<int I>
					euclid::space< euclid::dual< type<I> > > operator()() const {
						return *std::get<I>(impl);
					}
				};

				// computes successive offsets
				struct offset {
					NN& last;

					offset_type& offset;
					const impl_type& impl;
	
					template<int I>
					void operator()() const {
						std::get<I>(offset) = last;
						last += std::get<I>(impl).dim();
					}

				};

				// coordinate access
				struct coord {
					const E& value;
					const impl_type& impl;
					const offset_type& offset;
	
					const field*& res;
					NN i;
	
					template<int I>
					void operator()() const {
						const auto& impl_i = std::get<I>(impl);
						NN index = i - std::get<I>(offset);
						res = & impl_i.coord( index, std::get<I>(value) );
					}
	
				};

				// term-wise zero element
				struct zero {
	
					const impl_type& impl;
					
					template<int I>
					type<I> operator()() const {
						return std::get<I>(impl).zero();
					}
	
				};
      
			};
    
			// builds offset tuple and returns total dimension
			NN make_offset() {
				NN last = 0;
				each::apply( typename helper::offset{last, offset, impl} );
				return last;
			}
    
			// dual geometry
			space< dual > operator*() const {
				return { each::map( typename helper::dual{impl} ) };
			}
    
			// coordinate access
			const field& coord(NN i, const E& x) const {
				assert( i < dim() );

				const field* res = 0;
				tuple::find(offset, i, typename helper::coord{x, impl, offset, res, i} );
      
				return *res;
			}

			field& coord(NN i, E& x) const {
				assert( i < dim() );

				const field* res = 0;
				tuple::find(offset, i, typename helper::coord{x, impl, offset, res, i} );

				// uglyness i has it      
				return *const_cast<field*>(res);
			}
    
    
			E zero() const {
				return each::map( typename helper::zero{impl} );
			}

			NN dim() const { return dimension; }

			// ctors
			traits( const impl_type& impl = impl_type() ) 
				: impl(impl) {
				dimension = make_offset();
			}

			traits( const space<Args>&... args ) 
				: impl(args...) {
				dimension = make_offset();
			}

			traits( const E& value )
			: impl( each::map( typename helper::impl{value} ) ) {
				dimension = make_offset();
			}
    
		};

	}

	namespace lie {

		template<class ... Args>
		struct traits< std::tuple<Args...> > {
    
			typedef std::tuple<Args...> G;

			// iterator type
			typedef tuple::range< sizeof...(Args) > each;
    
			typedef std::tuple< group<Args>... > args_type;
			args_type args;
  
			typedef std::tuple< lie::algebra<Args>... > algebra;
			typedef euclid::dual<algebra> coalgebra;

			template<int I>
			using type = math0x::tuple::element<I, G>;


			struct get_alg {
				const args_type& args;

				template<int I>
				euclid::space< lie::algebra< type<I> > > operator()() const {
					return std::get<I>(args).alg();
				}

			};

			euclid::space< algebra > alg() const { 
				return { each::map( get_alg{args} ) };
			}

			struct get_id {
				const args_type& args;
      
				template<int I>
				type<I> operator()() const {
					return std::get<I>(args).id();
				}
      
			};

			struct get_inv {
				const args_type& args;
				const G& x;
      
				template<int I>
				type<I> operator()() const {
					return std::get<I>(args).inv( std::get<I>(x) );
				}
      
			};

			struct get_prod {
				const args_type& args;

				const G& lhs;
				const G& rhs;
      
				template<int I>
				type<I> operator()() const {
					return std::get<I>(args).prod( std::get<I>(lhs),
					                               std::get<I>(rhs));
				}
      
			};


			struct get_bracket {
				const args_type& args;
				
				const algebra& lhs;
				const algebra& rhs;
				
				template<int I>
				lie::algebra< type<I> > operator()() const {
					return std::get<I>(args).bracket( std::get<I>(lhs),
					                                  std::get<I>(rhs));
				}
				
			};

			struct get_cobracket {
				const args_type& args;
				
				const coalgebra& lhs;
				const coalgebra& rhs;
				
				template<int I>
				lie::coalgebra< type<I> > operator()() const {
					return std::get<I>(args).cobracket( std::get<I>(lhs),
					                                    std::get<I>(rhs));
				}
				
			};


	



			G id() const {
				return each::map( get_id{args} );
			}

			G inv(const G& x) const {
				return each::map( get_inv{args, x} );
			}

			G prod(const G& lhs, const G& rhs) const {
				return each::map( get_prod{args, lhs, rhs} );
			}


			algebra bracket(const algebra& x, const algebra& y) const {
				return each::map( get_bracket{args, x, y} );
			}

			coalgebra cobracket(const coalgebra& x, const coalgebra& y) const {
				return each::map( get_cobracket{args, x, y} );
			}


			struct get_group {
				const G& g;

				template<int I>
				lie::group< type<I> > operator()() const {
					return { std::get<I>(g) };
				}
      
			};

			traits( const G& g ) : args( each::map( get_group{g} ) ) { }
			traits( const args_type& args = args_type() ) : args(args) { }


			struct Ad : func::tuple< lie::Ad<Args>... > {
      
				struct get {
					const G& at;
	
					template<int I>
					lie::Ad< type<I> > operator()() const {
						return { std::get<I>(at) };
					}
				};
      
				Ad( const G& at) : Ad::base{ each::map( get{at} ) } { }
      
			};


			struct AdT : func::tuple< lie::AdT<Args>... > {
				
				struct get {
					const G& at;
	
					template<int I>
					lie::AdT< type<I> > operator()() const {
						return { std::get<I>(at) };
					}
				};
      
				AdT( const G& at) : AdT::base{ each::map( get{at} ) } { }
      
			};

    
			struct exp : func::tuple< lie::exp<Args>... >  {
				struct get {
					const args_type& args;

					template<int I>
					lie::exp< type<I> > operator()() const {
						return std::get<I>(args).exp();
					}
				};
      
				exp( const group<G>& g = group<G>() ) 
					: exp::base{ each::map( get{g.impl.args} ) }  { }
      
      
			};


			struct log : func::tuple< lie::log<Args>... >  {
				struct get {
					const args_type& args;

					template<int I>
					lie::log< type<I> > operator()() const {
						return std::get<I>(args).log();
					}
				};
      
				log( const group<G>& g = group<G>() ) 
					: log::base{ each::map( get{g.impl.args} ) }  { }
      
			};


		};
	}

}
#endif
