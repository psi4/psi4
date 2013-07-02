/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680

  $Id$
*/
#ifndef MADNESS_MRA_VMRA_H__INCLUDED
#define MADNESS_MRA_VMRA_H__INCLUDED

/*!
	\file vmra.h
	\brief Defines operations on vectors of Functions
	\ingroup mra

	This file defines a number of operations on vectors of functions.
	Assume v is a vector of NDIM-D functions of a certain type.


	Operations on array of functions

	*) copying: deep copying of vectors of functions to vector of functions
	\code
	vector2 = copy(world, vector1,fence);
	\endcode

	*) compress: convert multiwavelet representation to legendre representation
	\code
	compress(world, vector, fence);
	\endcode

	*) reconstruct: convert representation to multiwavelets
	\code
	reconstruct(world, vector, fence);
	\endcode

	*) nonstandard: convert to non-standard form
	\code
	nonstandard(world, v, fence);
	\endcode

	*) standard: convert to standard form
	\code
	standard(world, v, fence);
	\endcode

	*) truncate: truncating vectors of functions to desired precision
	\code
	truncate(world, v, tolerance, fence);
	\endcode


	*) zero function: create a vector of zero functions of length n
	\code
	v=zero(world, n);
	\endcode

	*) transform: transform a representation from one basis to another
	\code
	transform(world, vector, tensor, tolerance, fence )
	\endcode

	Setting thresh-hold for precision

	*) set_thresh: setting a finite thresh-hold for a vector of functions
	\code
	void set_thresh(World& world, std::vector< Function<T,NDIM> >& v, double thresh, bool fence=true);
	\endcode

	Arithmetic Operations on arrays of functions

	*) conjugation: conjugate a vector of complex functions

	*) add
	*) sub
	*) mul
	   - mul_sparse
	*) square
	*) gaxpy
	*) apply

	Norms, inner-products, blas-1 like operations on vectors of functions

	*) inner
	*) matrix_inner
	*) norm_tree
	*) normalize
	*) norm2
	    - norm2s
	*) scale(world, v, alpha);




*/

#include <mra/mra.h>
#include <mra/derivative.h>
#include <cstdio>

namespace madness {

    /// Compress a vector of functions
    template <typename T, std::size_t NDIM>
    void compress(World& world,
                  const std::vector< Function<T,NDIM> >& v,
                  bool fence=true) {

        PROFILE_BLOCK(Vcompress);
        bool must_fence = false;
        for (unsigned int i=0; i<v.size(); ++i) {
            if (!v[i].is_compressed()) {
                v[i].compress(false);
                must_fence = true;
            }
        }

        if (fence && must_fence) world.gop.fence();
    }


    /// Reconstruct a vector of functions
    template <typename T, std::size_t NDIM>
    void reconstruct(World& world,
                     const std::vector< Function<T,NDIM> >& v,
                     bool fence=true) {
        PROFILE_BLOCK(Vreconstruct);
        bool must_fence = false;
        for (unsigned int i=0; i<v.size(); ++i) {
            if (v[i].is_compressed()) {
                v[i].reconstruct(false);
                must_fence = true;
            }
        }

        if (fence && must_fence) world.gop.fence();
    }


    /// Generates non-standard form of a vector of functions
    template <typename T, std::size_t NDIM>
    void nonstandard(World& world,
                     std::vector< Function<T,NDIM> >& v,
                     bool fence=true) {
        PROFILE_BLOCK(Vnonstandard);
        reconstruct(world, v);
        for (unsigned int i=0; i<v.size(); ++i) {
            v[i].nonstandard(false,false);
        }
        if (fence) world.gop.fence();
    }


    /// Generates standard form of a vector of functions
    template <typename T, std::size_t NDIM>
    void standard(World& world,
                  std::vector< Function<T,NDIM> >& v,
                  bool fence=true) {
        PROFILE_BLOCK(Vstandard);
        for (unsigned int i=0; i<v.size(); ++i) {
            v[i].standard(false);
        }
        if (fence) world.gop.fence();
    }


    /// Truncates a vector of functions
    template <typename T, std::size_t NDIM>
    void truncate(World& world,
                  std::vector< Function<T,NDIM> >& v,
                  double tol=0.0,
                  bool fence=true) {
        PROFILE_BLOCK(Vtruncate);

        compress(world, v);

        for (unsigned int i=0; i<v.size(); ++i) {
            v[i].truncate(tol, false);
        }

        if (fence) world.gop.fence();
    }

    /// Applies a derivative operator to a vector of functions
    template <typename T, std::size_t NDIM>
    std::vector< Function<T,NDIM> >
    apply(World& world,
          const Derivative<T,NDIM>& D,
          const std::vector< Function<T,NDIM> >& v,
          bool fence=true)
    {
        reconstruct(world, v);
        std::vector< Function<T,NDIM> > df(v.size());
        for (unsigned int i=0; i<v.size(); ++i) {
            df[i] = D(v[i],false);
        }
        if (fence) world.gop.fence();
        return df;
    }

    /// Generates a vector of zero functions
    template <typename T, std::size_t NDIM>
    std::vector< Function<T,NDIM> >
    zero_functions(World& world, int n) {
        PROFILE_BLOCK(Vzero_functions);
        std::vector< Function<T,NDIM> > r(n);
        for (int i=0; i<n; ++i)
            r[i] = Function<T,NDIM>(FunctionFactory<T,NDIM>(world));

        return r;
    }


    /// Transforms a vector of functions according to new[i] = sum[j] old[j]*c[j,i]

    /// Uses sparsity in the transformation matrix --- set small elements to
    /// zero to take advantage of this.
    template <typename T, typename R, std::size_t NDIM>
    std::vector< Function<TENSOR_RESULT_TYPE(T,R),NDIM> >
    transform(World& world,
              const std::vector< Function<T,NDIM> >& v,
              const Tensor<R>& c,
              bool fence=true) {

        PROFILE_BLOCK(Vtransformsp);
        typedef TENSOR_RESULT_TYPE(T,R) resultT;
        int n = v.size();  // n is the old dimension
        int m = c.dim(1);  // m is the new dimension
        MADNESS_ASSERT(n==c.dim(0));

        std::vector< Function<resultT,NDIM> > vc = zero_functions<resultT,NDIM>(world, m);
        compress(world, v);
        compress(world, vc);

        for (int i=0; i<m; ++i) {
            for (int j=0; j<n; ++j) {
                if (c(j,i) != R(0.0)) vc[i].gaxpy(1.0,v[j],c(j,i),false);
            }
        }

        if (fence) world.gop.fence();
        return vc;
    }


    template <typename L, typename R, std::size_t NDIM>
    std::vector< Function<TENSOR_RESULT_TYPE(L,R),NDIM> >
    transform(World& world,  const std::vector< Function<L,NDIM> >& v, const Tensor<R>& c, double tol, bool fence) {
        PROFILE_BLOCK(Vtransform);
        MADNESS_ASSERT(v.size() == (unsigned int)(c.dim(0)));

        std::vector< Function<TENSOR_RESULT_TYPE(L,R),NDIM> > vresult(c.dim(1));
        for (int i=0; i<c.dim(1); ++i) {
            vresult[i] = Function<TENSOR_RESULT_TYPE(L,R),NDIM>(FunctionFactory<TENSOR_RESULT_TYPE(L,R),NDIM>(world));
        }
        compress(world, v, false);
        compress(world, vresult, false);
        world.gop.fence();
        vresult[0].vtransform(v, c, vresult, tol, fence);
        return vresult;
    }

    /// Scales inplace a vector of functions by distinct values
    template <typename T, typename Q, std::size_t NDIM>
    void scale(World& world,
               std::vector< Function<T,NDIM> >& v,
               const std::vector<Q>& factors,
               bool fence=true) {
        PROFILE_BLOCK(Vscale);
        for (unsigned int i=0; i<v.size(); ++i) v[i].scale(factors[i],false);
        if (fence) world.gop.fence();
    }

    /// Scales inplace a vector of functions by the same 
    template <typename T, typename Q, std::size_t NDIM>
    void scale(World& world,
               std::vector< Function<T,NDIM> >& v,
	       const Q factor,
               bool fence=true) {
        PROFILE_BLOCK(Vscale);
        for (unsigned int i=0; i<v.size(); ++i) v[i].scale(factor,false);
        if (fence) world.gop.fence();
    }

    /// Computes the 2-norms of a vector of functions
    template <typename T, std::size_t NDIM>
    std::vector<double> norm2s(World& world,
                              const std::vector< Function<T,NDIM> >& v) {
        PROFILE_BLOCK(Vnorm2);
        std::vector<double> norms(v.size());
        for (unsigned int i=0; i<v.size(); ++i) norms[i] = v[i].norm2sq_local();
        world.gop.sum(&norms[0], norms.size());
        for (unsigned int i=0; i<v.size(); ++i) norms[i] = sqrt(norms[i]);
        world.gop.fence();
        return norms;
    }

    /// Computes the 2-norm of a vector of functions
    template <typename T, std::size_t NDIM>
    double norm2(World& world,
                              const std::vector< Function<T,NDIM> >& v) {
        PROFILE_BLOCK(Vnorm2);
        std::vector<double> norms(v.size());
        for (unsigned int i=0; i<v.size(); ++i) norms[i] = v[i].norm2sq_local();
        world.gop.sum(&norms[0], norms.size());
        for (unsigned int i=1; i<v.size(); ++i) norms[0] += norms[i];
        world.gop.fence();
        return sqrt(norms[0]);
    }

    inline double conj(double x) {
        return x;
    }

    inline double conj(float x) {
        return x;
    }


    template <typename T, typename R, std::size_t NDIM>
    struct MatrixInnerTask : public TaskInterface {
        Tensor<TENSOR_RESULT_TYPE(T,R)> result; // Must be a copy
        const Function<T,NDIM>& f;
        const std::vector< Function<R,NDIM> >& g;
        long jtop;

        using PoolTaskInterface::run;

        MatrixInnerTask(const Tensor<TENSOR_RESULT_TYPE(T,R)>& result,
                        const Function<T,NDIM>& f,
                        const std::vector< Function<R,NDIM> >& g,
                        long jtop)
                : result(result), f(f), g(g), jtop(jtop) {}

        void run(World& world) {
            for (long j=0; j<jtop; ++j) {
                result(j) = f.inner_local(g[j]);
            }
        }
    };


    /// Computes the matrix inner product of two function vectors - q(i,j) = inner(f[i],g[j])

    /// For complex types symmetric is interpreted as Hermitian.
    ///
    /// The current parallel loop is non-optimal but functional.
    template <typename T, typename R, std::size_t NDIM>
    Tensor< TENSOR_RESULT_TYPE(T,R) > matrix_inner(World& world,
            const std::vector< Function<T,NDIM> >& f,
            const std::vector< Function<R,NDIM> >& g,
            bool sym=false) {
        PROFILE_BLOCK(Vmatrix_inner);
        long n=f.size(), m=g.size();
        Tensor< TENSOR_RESULT_TYPE(T,R) > r(n,m);
        if (sym) MADNESS_ASSERT(n==m);

        world.gop.fence();
        compress(world, f);
        if (&f != &g) compress(world, g);

//         for (long i=0; i<n; ++i) {
//             long jtop = m;
//             if (sym) jtop = i+1;
//             for (long j=0; j<jtop; ++j) {
//                 r(i,j) = f[i].inner_local(g[j]);
//                 if (sym) r(j,i) = conj(r(i,j));
//             }
//         }

        for (long i=n-1; i>=0; --i) {
            long jtop = m;
            if (sym) jtop = i+1;
            world.taskq.add(new MatrixInnerTask<T,R,NDIM>(r(i,_), f[i], g, jtop));
        }
        world.gop.fence();
        world.gop.sum(r.ptr(),n*m);

        if (sym) {
            for (int i=0; i<n; ++i) {
                for (int j=0; j<i; ++j) {
                    r(j,i) = conj(r(i,j));
                }
            }
        }
        return r;
    }

    /// Computes the element-wise inner product of two function vectors - q(i) = inner(f[i],g[i])
    template <typename T, typename R, std::size_t NDIM>
    Tensor< TENSOR_RESULT_TYPE(T,R) > inner(World& world,
                                            const std::vector< Function<T,NDIM> >& f,
                                            const std::vector< Function<R,NDIM> >& g) {
        PROFILE_BLOCK(Vinnervv);
        long n=f.size(), m=g.size();
        MADNESS_ASSERT(n==m);
        Tensor< TENSOR_RESULT_TYPE(T,R) > r(n);

        compress(world, f);
        compress(world, g);

        for (long i=0; i<n; ++i) {
            r(i) = f[i].inner_local(g[i]);
        }

        world.taskq.fence();
        world.gop.sum(r.ptr(),n);
        world.gop.fence();
        return r;
    }


    /// Computes the inner product of a function with a function vector - q(i) = inner(f,g[i])
    template <typename T, typename R, std::size_t NDIM>
    Tensor< TENSOR_RESULT_TYPE(T,R) > inner(World& world,
                                            const Function<T,NDIM>& f,
                                            const std::vector< Function<R,NDIM> >& g) {
        PROFILE_BLOCK(Vinner);
        long n=g.size();
        Tensor< TENSOR_RESULT_TYPE(T,R) > r(n);

        f.compress();
        compress(world, g);

        for (long i=0; i<n; ++i) {
            r(i) = f.inner_local(g[i]);
        }

        world.taskq.fence();
        world.gop.sum(r.ptr(),n);
        world.gop.fence();
        return r;
    }


    /// Multiplies a function against a vector of functions --- q[i] = a * v[i]
    template <typename T, typename R, std::size_t NDIM>
    std::vector< Function<TENSOR_RESULT_TYPE(T,R), NDIM> >
    mul(World& world,
        const Function<T,NDIM>& a,
        const std::vector< Function<R,NDIM> >& v,
        bool fence=true) {
        PROFILE_BLOCK(Vmul);
        a.reconstruct(false);
        reconstruct(world, v, false);
        world.gop.fence();
        return vmulXX(a, v, 0.0, fence);
    }

    /// Multiplies a function against a vector of functions using sparsity of a and v[i] --- q[i] = a * v[i]
    template <typename T, typename R, std::size_t NDIM>
    std::vector< Function<TENSOR_RESULT_TYPE(T,R), NDIM> >
    mul_sparse(World& world,
               const Function<T,NDIM>& a,
               const std::vector< Function<R,NDIM> >& v,
               double tol,
               bool fence=true) {
        PROFILE_BLOCK(Vmulsp);
        a.reconstruct(false);
        reconstruct(world, v, false);
        world.gop.fence();
        for (unsigned int i=0; i<v.size(); ++i) {
            v[i].norm_tree(false);
        }
        a.norm_tree();
        return vmulXX(a, v, tol, fence);
    }

    /// Makes the norm tree for all functions in a vector
    template <typename T, std::size_t NDIM>
    void norm_tree(World& world,
                   const std::vector< Function<T,NDIM> >& v,
                   bool fence=true)
    {
        PROFILE_BLOCK(Vnorm_tree);
        for (unsigned int i=0; i<v.size(); ++i) {
            v[i].norm_tree(false);
        }
        if (fence) world.gop.fence();
    }

    /// Multiplies two vectors of functions q[i] = a[i] * b[i]
    template <typename T, typename R, std::size_t NDIM>
    std::vector< Function<TENSOR_RESULT_TYPE(T,R), NDIM> >
    mul(World& world,
        const std::vector< Function<T,NDIM> >& a,
        const std::vector< Function<R,NDIM> >& b,
        bool fence=true) {
        PROFILE_BLOCK(Vmulvv);
        reconstruct(world, a, false);
        if (&a != &b) reconstruct(world, b, false);
        world.gop.fence();

        std::vector< Function<TENSOR_RESULT_TYPE(T,R),NDIM> > q(a.size());
        for (unsigned int i=0; i<a.size(); ++i) {
            q[i] = mul(a[i], b[i], false);
        }
        if (fence) world.gop.fence();
        return q;
    }


    /// Computes the square of a vector of functions --- q[i] = v[i]**2
    template <typename T, std::size_t NDIM>
    std::vector< Function<T,NDIM> >
    square(World& world,
           const std::vector< Function<T,NDIM> >& v,
           bool fence=true) {
        return mul<T,T,NDIM>(world, v, v, fence);
//         std::vector< Function<T,NDIM> > vsq(v.size());
//         for (unsigned int i=0; i<v.size(); ++i) {
//             vsq[i] = square(v[i], false);
//         }
//         if (fence) world.gop.fence();
//         return vsq;
    }


    /// Sets the threshold in a vector of functions
    template <typename T, std::size_t NDIM>
    void set_thresh(World& world, std::vector< Function<T,NDIM> >& v, double thresh, bool fence=true) {
        for (unsigned int j=0; j<v.size(); ++j) {
            v[j].set_thresh(thresh,false);
        }
        if (fence) world.gop.fence();
    }

    /// Returns the complex conjugate of the vector of functions
    template <typename T, std::size_t NDIM>
    std::vector< Function<T,NDIM> >
    conj(World& world,
         const std::vector< Function<T,NDIM> >& v,
         bool fence=true) {
        PROFILE_BLOCK(Vconj);
        std::vector< Function<T,NDIM> > r = copy(world, v); // Currently don't have oop conj
        for (unsigned int i=0; i<v.size(); ++i) {
            r[i].conj(false);
        }
        if (fence) world.gop.fence();
        return r;
    }


    /// Returns a deep copy of a vector of functions
    template <typename T, std::size_t NDIM>
    std::vector< Function<T,NDIM> >
    copy(World& world,
         const std::vector< Function<T,NDIM> >& v,
         bool fence=true) {
        PROFILE_BLOCK(Vcopy);
        std::vector< Function<T,NDIM> > r(v.size());
        for (unsigned int i=0; i<v.size(); ++i) {
            r[i] = copy(v[i], false);
        }
        if (fence) world.gop.fence();
        return r;
    }

   /// Returns a vector of deep copies of of a function
    template <typename T, std::size_t NDIM>
    std::vector< Function<T,NDIM> >
    copy(World& world,
         const Function<T,NDIM>& v,
         const unsigned int n,
         bool fence=true) {
        PROFILE_BLOCK(Vcopy1);
        std::vector< Function<T,NDIM> > r(n);
        for (unsigned int i=0; i<n; ++i) {
            r[i] = copy(v, false);
        }
        if (fence) world.gop.fence();
        return r;
    }

    /// Returns new vector of functions --- q[i] = a[i] + b[i]
    template <typename T, typename R, std::size_t NDIM>
    std::vector< Function<TENSOR_RESULT_TYPE(T,R), NDIM> >
    add(World& world,
        const std::vector< Function<T,NDIM> >& a,
        const std::vector< Function<R,NDIM> >& b,
        bool fence=true) {
        PROFILE_BLOCK(Vadd);
        MADNESS_ASSERT(a.size() == b.size());
        compress(world, a);
        compress(world, b);

        std::vector< Function<TENSOR_RESULT_TYPE(T,R),NDIM> > r(a.size());
        for (unsigned int i=0; i<a.size(); ++i) {
            r[i] = add(a[i], b[i], false);
        }
        if (fence) world.gop.fence();
        return r;
    }

     /// Returns new vector of functions --- q[i] = a + b[i]
    template <typename T, typename R, std::size_t NDIM>
    std::vector< Function<TENSOR_RESULT_TYPE(T,R), NDIM> >
    add(World& world,
        const Function<T,NDIM> & a,
        const std::vector< Function<R,NDIM> >& b,
        bool fence=true) {
        PROFILE_BLOCK(Vadd1);
        a.compress();
        compress(world, b);

        std::vector< Function<TENSOR_RESULT_TYPE(T,R),NDIM> > r(b.size());
        for (unsigned int i=0; i<b.size(); ++i) {
            r[i] = add(a, b[i], false);
        }
        if (fence) world.gop.fence();
        return r;
    }
    template <typename T, typename R, std::size_t NDIM>
    inline std::vector< Function<TENSOR_RESULT_TYPE(T,R), NDIM> >
    add(World& world,
        const std::vector< Function<R,NDIM> >& b,
        const Function<T,NDIM> & a,
        bool fence=true) {
        return add(world, a, b, fence);
    }

    /// Returns new vector of functions --- q[i] = a[i] - b[i]
    template <typename T, typename R, std::size_t NDIM>
    std::vector< Function<TENSOR_RESULT_TYPE(T,R), NDIM> >
    sub(World& world,
        const std::vector< Function<T,NDIM> >& a,
        const std::vector< Function<R,NDIM> >& b,
        bool fence=true) {
        PROFILE_BLOCK(Vsub);
        MADNESS_ASSERT(a.size() == b.size());
        compress(world, a);
        compress(world, b);

        std::vector< Function<TENSOR_RESULT_TYPE(T,R),NDIM> > r(a.size());
        for (unsigned int i=0; i<a.size(); ++i) {
            r[i] = sub(a[i], b[i], false);
        }
        if (fence) world.gop.fence();
        return r;
    }


    /// Generalized A*X+Y for vectors of functions ---- a[i] = alpha*a[i] + beta*b[i]
    template <typename T, typename Q, typename R, std::size_t NDIM>
    void gaxpy(World& world,
               Q alpha,
               std::vector< Function<T,NDIM> >& a,
               Q beta,
               const std::vector< Function<R,NDIM> >& b,
               bool fence=true) {
        PROFILE_BLOCK(Vgaxpy);
        MADNESS_ASSERT(a.size() == b.size());
        compress(world, a);
        compress(world, b);

        for (unsigned int i=0; i<a.size(); ++i) {
            a[i].gaxpy(alpha, b[i], beta, false);
        }
        if (fence) world.gop.fence();
    }


    /// Applies a vector of operators to a vector of functions --- q[i] = apply(op[i],f[i])
    template <typename opT, typename R, std::size_t NDIM>
    std::vector< Function<TENSOR_RESULT_TYPE(typename opT::opT,R), NDIM> >
    apply(World& world,
          const std::vector< std::shared_ptr<opT> >& op,
          const std::vector< Function<R,NDIM> > f) {

        PROFILE_BLOCK(Vapplyv);
        MADNESS_ASSERT(f.size()==op.size());

        std::vector< Function<R,NDIM> >& ncf = *const_cast< std::vector< Function<R,NDIM> >* >(&f);

        reconstruct(world, f);
        nonstandard(world, ncf);

        std::vector< Function<TENSOR_RESULT_TYPE(typename opT::opT,R), NDIM> > result(f.size());
        for (unsigned int i=0; i<f.size(); ++i) {
            result[i] = apply_only(*op[i], f[i], false);
        }

        world.gop.fence();

        standard(world, ncf, false);  // restores promise of logical constness
        world.gop.fence();
        reconstruct(world, result);

        return result;
    }


    /// Applies an operator to a vector of functions --- q[i] = apply(op,f[i])
    template <typename T, typename R, std::size_t NDIM>
    std::vector< Function<TENSOR_RESULT_TYPE(T,R), NDIM> >
    apply(World& world,
          const SeparatedConvolution<T,NDIM>& op,
          const std::vector< Function<R,NDIM> > f) {
        PROFILE_BLOCK(Vapply);

        std::vector< Function<R,NDIM> >& ncf = *const_cast< std::vector< Function<R,NDIM> >* >(&f);

        reconstruct(world, f);
        nonstandard(world, ncf);

        std::vector< Function<TENSOR_RESULT_TYPE(T,R), NDIM> > result(f.size());
        for (unsigned int i=0; i<f.size(); ++i) {
            result[i] = apply_only(op, f[i], false);
        }

        world.gop.fence();

        standard(world, ncf, false);  // restores promise of logical constness
        reconstruct(world, result);

        return result;
    }

    /// Normalizes a vector of functions --- v[i] = v[i].scale(1.0/v[i].norm2())
    template <typename T, std::size_t NDIM>
    void normalize(World& world, std::vector< Function<T,NDIM> >& v, bool fence=true) {
        PROFILE_BLOCK(Vnormalize);
        std::vector<double> nn = norm2s(world, v);
        for (unsigned int i=0; i<v.size(); ++i) v[i].scale(1.0/nn[i],false);
        if (fence) world.gop.fence();
    }

}
#endif // MADNESS_MRA_VMRA_H__INCLUDED
