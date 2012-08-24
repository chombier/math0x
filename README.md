**Warning:** This is a complete work in progress, nowhere near a
releasable state.

**group** is a small c++11 template library for doing differential
geometry in Lie groups. It can be used to automatically compute
complicated derivatives/pullbacks, for instance when dealing with
robot kinematics.

*group* tries to provide a high degree of flexibility by providing a
collection of small, atomic functions that can be chained together to
form more complex functions. This composition happens at compile time,
thus the compiler will +beg for mercy+ do its best to produce
efficient code.

You can use your custom data types with *group* by implementing a few
/trait classes/ describing Euclidean and Lie group operations. Since
most of the functions in ~group~ only use these interfaces, you will
be able to use them right away in your code. For instance, you can use
*Eigen* vectors with *group*.

# TODO Examples

## Basic Types

```c++
// basic mathematical types
NN n = 1;   // natural numbers
ZZ p = -1;  // integers
RR x = 0.1; // real numbers
```
## Vector Types (from the Eigen library)

```c++
 // alias Eigen type
 typedef Eigen::Vector3d vec3;
  
 // some 3-vector, equals (1, 1, 1)
 vec3 u = vec3::Ones();
```

## Functions, Derivatives

```c++
// vector line with direction u: RR -> vec3
func::line<vec3> line(u);

// let v == (0.5, 0.5, 0.5)
vec3 v = line(0.5);

// construct a more complex function by composing line and
// squared norm functions
auto f = func::norm2<vec3>() << line;

// real variable, equals 3
RR alpha = f(1.0);

// derivative of f at 1.0
auto df = func::d(f)(1.0);

// derive f at 1.0, along tangent vector 2.0
RR beta = df(2.0);
```

## Euclidean Structure

Data types representing Euclidean spaces must provide a class
describing the Euclidean structure:

```c++
euclid::space<vec3> vec3_info;

u = vec3_info.zero();        // zero element
u = vec3_info.minus( u );    // additive inverse
u = vec3_info.scal(0.5, u);  // scalar multiplication
u = vec3_info.sum(u, v);     // sum

// coordinate access, field type
euclid::field<vec3>& u0 = vec3_info.coord(0, u);

u0 = 5;
// u is now (5, 1, 1)
```

For dynamic-sized types:

``` c++
// dynamic-sized vectors
typedef Eigen::VectorXd vec;
 
// we want 10-vectors
euclid::space<vec> vec_info(10);

// now we know that zero should be a 10-vector
vec w = vec_info.zero(); 
```

## Tuple Types

=std::tuple= can be used to construct mathematical direct
products. When possible, the Euclidean structure is automatically
infered:

```c++
typdedef std::tuple< vec3, vec3 > my_tuple_type;
euclid::space<my_tuple_type> my_tuple_info;
```

When using dynamic-sized types, we need to provide structure for each
of the types in the tuple:

```c++
// direct product type
typedef std::tuple<vec3, vec> my_type;
 
// Euclidean space structure for my_type: we indicate that dynamic-sized
// vectors should be 10-dimensional				
euclid::space<my_type> my_info( vec3_info, vec_info );

// squared norm: my_type -> RR
auto g = func::norm2<my_type>(info);
```

## Algorithms

The above Euclidean space abstraction can be used to build generic
algorithms that work for /any/ Euclidean space:

```c++
template<class E>
struct my_algo {
  euclid::space<E> info;
  my_algo( const euclid::space<E>& info = euclid::space<E>() ) : info(info) { }
 
  // ... fancy algorithm here, using only operations provided by info

};
```

This might seems a little overkill for now, but the same mechanism
applies for Lie groups as well, and enables us to write *one* single
spline interpolation algorithm that works for *any* classical Lie
group. Yes, even for that funny Lie group 

	\[ G = SE(3) \times SO(3) \times SL(2) \times \mathbb{R}^10 \]
	
should you want to do it.

# TODO Usage



# TODO Documentation

hahaha wat


