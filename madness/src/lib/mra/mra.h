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


  $Id: mra.h 2371 2011-06-15 00:22:38Z rjharrison $

*/


#ifndef MADNESS_MRA_MRA_H__INCLUDED
#define MADNESS_MRA_MRA_H__INCLUDED

/*!
  \file mra/mra.h
  \brief Main include file for MADNESS and defines \c Function interface

  \addtogroup mra

*/


#include <world/world.h>
#include <misc/misc.h>
#include <tensor/tensor.h>

#define FUNCTION_INSTANTIATE_1
#define FUNCTION_INSTANTIATE_2
#define FUNCTION_INSTANTIATE_3
#ifndef HAVE_IBMBGP
#define FUNCTION_INSTANTIATE_4
#define FUNCTION_INSTANTIATE_5
#define FUNCTION_INSTANTIATE_6
#endif

static const bool VERIFY_TREE = false; //true;


namespace madness {
    void startup(World& world, int argc, char** argv);
}

#include <mra/key.h>
#include <mra/twoscale.h>
#include <mra/legendre.h>
#include <mra/indexit.h>
#include <mra/funcimpl.h>
#include <mra/loadbal.h>

namespace madness {
	/// \ingroup mra
    /// \addtogroup function

    /// A multiresolution adaptive numerical function
    template <typename T, std::size_t NDIM>
    class Function : public archive::ParallelSerializableObject {
        // We make all of the content of Function and FunctionImpl
        // public with the intent of avoiding the cumbersome forward
        // and friend declarations.  However, this open access should
        // not be abused.

    private:
        std::shared_ptr< FunctionImpl<T,NDIM> > impl;

    public:
        typedef FunctionImpl<T,NDIM> implT;
        typedef FunctionFactory<T,NDIM> factoryT;
        typedef typename implT::coordT coordT; ///< Type of vector holding coordinates

        /// Asserts that the function is initialized
        inline void verify() const {
            MADNESS_ASSERT(impl);
        }

        /// Returns true if the function is initialized
        bool is_initialized() const {
            return impl.get();
        }

        /// Default constructor makes uninitialized function.  No communication.

        /// An unitialized function can only be assigned to.  Any other operation will throw.
        Function() : impl() {}


        /// Constructor from FunctionFactory provides named parameter idiom.  Possible non-blocking communication.
        Function(const factoryT& factory)
                : impl(new FunctionImpl<T,NDIM>(factory)) {
            PROFILE_MEMBER_FUNC(Function);
        }


        /// Copy constructor is \em shallow.  No communication, works in either basis.
        Function(const Function<T,NDIM>& f)
                : impl(f.impl) {
        }


        /// Assignment is \em shallow.  No communication, works in either basis.
        Function<T,NDIM>& operator=(const Function<T,NDIM>& f) {
            PROFILE_MEMBER_FUNC(Function);
            if (this != &f) impl = f.impl;
            return *this;
        }

        /// Destruction of any underlying implementation is deferred to next global fence.
        ~Function() {}

        /// Evaluates the function at a point in user coordinates.  Possible non-blocking comm.

        /// Only the invoking process will receive the result via the future
        /// though other processes may be involved in the evaluation.
        ///
        /// Throws if function is not initialized.
        Future<T> eval(const coordT& xuser) const {
            PROFILE_MEMBER_FUNC(Function);
            const double eps=1e-15;
            verify();
            MADNESS_ASSERT(!is_compressed());
            coordT xsim;
            user_to_sim(xuser,xsim);
            // If on the boundary, move the point just inside the
            // volume so that the evaluation logic does not fail
            for (std::size_t d=0; d<NDIM; ++d) {
                if (xsim[d] < -eps) {
                    MADNESS_EXCEPTION("eval: coordinate lower-bound error in dimension", d);
                }
                else if (xsim[d] < eps) {
                    xsim[d] = eps;
                }

                if (xsim[d] > 1.0+eps) {
                    MADNESS_EXCEPTION("eval: coordinate upper-bound error in dimension", d);
                }
                else if (xsim[d] > 1.0-eps) {
                    xsim[d] = 1.0-eps;
                }
            }

            Future<T> result;
            impl->eval(xsim, impl->key0(), result.remote_ref(impl->world));
            return result;
        }

        /// Evaluate function only if point is local returning (true,value); otherwise return (false,0.0)

        /// maxlevel is the maximum depth to search down to --- the max local depth can be
        /// computed with max_local_depth();
        std::pair<bool,T> eval_local_only(const Vector<double,NDIM>& xuser, Level maxlevel) const {
            const double eps=1e-15;
            verify();
            MADNESS_ASSERT(!is_compressed());
            coordT xsim;
            user_to_sim(xuser,xsim);
            // If on the boundary, move the point just inside the
            // volume so that the evaluation logic does not fail
            for (std::size_t d=0; d<NDIM; ++d) {
                if (xsim[d] < -eps) {
                    MADNESS_EXCEPTION("eval: coordinate lower-bound error in dimension", d);
                }
                else if (xsim[d] < eps) {
                    xsim[d] = eps;
                }

                if (xsim[d] > 1.0+eps) {
                    MADNESS_EXCEPTION("eval: coordinate upper-bound error in dimension", d);
                }
                else if (xsim[d] > 1.0-eps) {
                    xsim[d] = 1.0-eps;
                }
            }
            return impl->eval_local_only(xsim,maxlevel);
        }

        /// Only the invoking process will receive the result via the future
        /// though other processes may be involved in the evaluation.
        ///
        /// Throws if function is not initialized.
        ///
        /// This function is a minimally-modified version of eval()
        Future<Level> evaldepthpt(const coordT& xuser) const {
            PROFILE_MEMBER_FUNC(Function);
            const double eps=1e-15;
            verify();
            MADNESS_ASSERT(!is_compressed());
            coordT xsim;
            user_to_sim(xuser,xsim);
            // If on the boundary, move the point just inside the
            // volume so that the evaluation logic does not fail
            for (std::size_t d=0; d<NDIM; ++d) {
                if (xsim[d] < -eps) {
                    MADNESS_EXCEPTION("eval: coordinate lower-bound error in dimension", d);
                }
                else if (xsim[d] < eps) {
                    xsim[d] = eps;
                }

                if (xsim[d] > 1.0+eps) {
                    MADNESS_EXCEPTION("eval: coordinate upper-bound error in dimension", d);
                }
                else if (xsim[d] > 1.0-eps) {
                    xsim[d] = 1.0-eps;
                }
            }

            Future<Level> result;
            impl->evaldepthpt(xsim, impl->key0(), result.remote_ref(impl->world));
            return result;
        }

        /// Evaluates a cube/slice of points (probably for plotting) ... collective but no fence necessary

        /// All processes recieve the entire result (which is a rather severe limit
        /// on the size of the cube that is possible).

        /// Set eval_refine=true to return the refinment levels of
        /// the given function.

        /// @param[in] cell A Tensor describe the cube where the function to be evaluated in
        /// @param[in] npt How many points to evaluate in each dimension
        /// @param[in] eval_refine Wether to return the refinment levels of the given function
        Tensor<T> eval_cube(const Tensor<double>& cell,
                            const std::vector<long>& npt,
                            bool eval_refine = false) const {
            MADNESS_ASSERT(static_cast<std::size_t>(cell.dim(0))>=NDIM && cell.dim(1)==2 && npt.size()>=NDIM);
            PROFILE_MEMBER_FUNC(Function);
            const double eps=1e-14;
            verify();
            reconstruct();
            coordT simlo, simhi;
            for (std::size_t d=0; d<NDIM; ++d) {
                simlo[d] = cell(d,0);
                simhi[d] = cell(d,1);
            }
            user_to_sim(simlo, simlo);
            user_to_sim(simhi, simhi);

            // Move the bounding box infintesimally inside dyadic
            // points so that the evaluation logic does not fail
            for (std::size_t d=0; d<NDIM; ++d) {
                MADNESS_ASSERT(simhi[d] >= simlo[d]);
                MADNESS_ASSERT(simlo[d] >= 0.0);
                MADNESS_ASSERT(simhi[d] <= 1.0);

                double delta = eps*(simhi[d]-simlo[d]);
                simlo[d] += delta;
                simhi[d] -= 2*delta;  // deliberate asymmetry
            }
            return impl->eval_plot_cube(simlo, simhi, npt, eval_refine);
        }


        /// Evaluates the function at a point in user coordinates.  Collective operation.

        /// Throws if function is not initialized.
        ///
        /// This function calls eval, blocks until the result is
        /// available and then broadcasts the result to everyone.
        /// Therefore, if you are evaluating many points in parallel
        /// it is \em vastly less efficient than calling eval
        /// directly, saving the futures, and then forcing all of the
        /// results.
        T operator()(const coordT& xuser) const {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            if (is_compressed()) reconstruct();
            T result;
            if (impl->world.rank() == 0) result = eval(xuser).get();
            impl->world.gop.broadcast(result);
            //impl->world.gop.fence();
            return result;
        }

        /// Evaluates the function at a point in user coordinates.  Collective operation.

        /// See "operator()(const coordT& xuser)" for more info
        T operator()(double x, double y=0, double z=0, double xx=0, double yy=0, double zz=0) const {
            coordT r;
            r[0] = x;
            if (NDIM>=2) r[1] = y;
            if (NDIM>=3) r[2] = z;
            if (NDIM>=4) r[3] = xx;
            if (NDIM>=5) r[4] = yy;
            if (NDIM>=6) r[5] = zz;
            return (*this)(r);
        }

        /// Throws if function is not initialized.
        ///
        /// This function mimics operator() by going through the
        /// tree looking for the depth of the tree at the point.
        /// It blocks until the result is
        /// available and then broadcasts the result to everyone.
        /// Therefore, if you are evaluating many points in parallel
        /// it is \em vastly less efficient than calling evaldepthpt
        /// directly, saving the futures, and then forcing all of the
        /// results.
        Level depthpt(const coordT& xuser) const {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            if (is_compressed()) reconstruct();
            Level result;
            if (impl->world.rank() == 0) result = evaldepthpt(xuser).get();
            impl->world.gop.broadcast(result);
            //impl->world.gop.fence();
            return result;
        }

        /// Returns an estimate of the difference ||this-func||^2 from local data

        /// No communication is performed.  If the function is not
        /// reconstructed, it throws an exception.  To get the global
        /// value either do a global sum of the local values or call
        /// errsq
        /// @param[in] func Templated interface to the a user specified function
        template <typename funcT>
        double errsq_local(const funcT& func) const {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            if (is_compressed()) MADNESS_EXCEPTION("Function:errsq_local:not reconstructed",0);
            return impl->errsq_local(func);
        }


        /// Returns an estimate of the difference ||this-func|| ... global sum performed

        /// If the function is compressed, it is reconstructed first.  For efficient use
        /// especially with many functions, reconstruct them all first, and use errsq_local
        /// instead so you can perform a global sum on all at the same time.
        /// @param[in] func Templated interface to the a user specified function
        template <typename funcT>
        double err(const funcT& func) const {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            if (VERIFY_TREE) verify_tree();
            if (is_compressed()) reconstruct();
            if (VERIFY_TREE) verify_tree();
            double local = impl->errsq_local(func);
            impl->world.gop.sum(local);
            impl->world.gop.fence();
            return sqrt(local);
        }

        /// Verifies the tree data structure ... global sync implied
        void verify_tree() const {
            PROFILE_MEMBER_FUNC(Function);
            if (impl) impl->verify_tree();
        }


        /// Returns true if compressed, false otherwise.  No communication.

        /// If the function is not initialized, returns false.
        bool is_compressed() const {
            PROFILE_MEMBER_FUNC(Function);
            if (impl)
                return impl->is_compressed();
            else
                return false;
        }


        /// Returns the number of nodes in the function tree ... collective global sum
        std::size_t tree_size() const {
            PROFILE_MEMBER_FUNC(Function);
            if (!impl) return 0;
            return impl->tree_size();
        }


        /// Returns the maximum depth of the function tree ... collective global sum
        std::size_t max_depth() const {
            PROFILE_MEMBER_FUNC(Function);
            if (!impl) return 0;
            return impl->max_depth();
        }


        /// Returns the maximum local depth of the function tree ... no communications
        std::size_t max_local_depth() const {
            PROFILE_MEMBER_FUNC(Function);
            if (!impl) return 0;
            return impl->max_local_depth();
        }


        /// Returns the max number of nodes on a processor
        std::size_t max_nodes() const {
            PROFILE_MEMBER_FUNC(Function);
            if (!impl) return 0;
            return impl->max_nodes();
        }

        /// Returns the min number of nodes on a processor
        std::size_t min_nodes() const {
            PROFILE_MEMBER_FUNC(Function);
            if (!impl) return 0;
            return impl->min_nodes();
        }


        /// Returns the number of coefficients in the function ... collective global sum
        std::size_t size() const {
            PROFILE_MEMBER_FUNC(Function);
            if (!impl) return 0;
            return impl->size();
        }


        /// Returns value of autorefine flag.  No communication.
        bool autorefine() const {
            PROFILE_MEMBER_FUNC(Function);
            if (!impl) return true;
            return impl->get_autorefine();
        }


        /// Sets the value of the autorefine flag.  Optional global fence.

        /// A fence is required to ensure consistent global state.
        void set_autorefine(bool value, bool fence = true) {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            impl->set_autorefine(value);
            if (fence) impl->world.gop.fence();
        }


        /// Returns value of truncation threshold.  No communication.
        double thresh() const {
            PROFILE_MEMBER_FUNC(Function);
            if (!impl) return 0.0;
            return impl->get_thresh();
        }


        /// Sets the vaule of the truncation threshold.  Optional global fence.

        /// A fence is required to ensure consistent global state.
        void set_thresh(double value, bool fence = true) {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            impl->set_thresh(value);
            if (fence) impl->world.gop.fence();
        }


        /// Returns the number of multiwavelets (k).  No communication.
        int k() const {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            return impl->get_k();
        }


        /// Truncate the function with optional fence.  Compresses with fence if not compressed.

        /// If the truncation threshold is less than or equal to zero the default value
        /// specified when the function was created is used.
        /// If the function is not initialized, it just returns.
        ///
        /// Returns this for chaining.
		/// @param[in] tol Tolerance for truncating the coefficients. Default 0.0 means use the implimentation's member value \c thresh instead.
		/// @param[in] fence Do fence
        Function<T,NDIM>& truncate(double tol = 0.0, bool fence = true) {
            PROFILE_MEMBER_FUNC(Function);
            if (!impl) return *this;
            verify();
            if (!is_compressed()) compress();
            impl->truncate(tol,fence);
            if (VERIFY_TREE) verify_tree();
            return *this;
        }


        /// Returns a shared-pointer to the implementation
        const std::shared_ptr< FunctionImpl<T,NDIM> >& get_impl() const {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            return impl;
        }

        /// Replace current FunctionImpl with provided new one
        void set_impl(const std::shared_ptr< FunctionImpl<T,NDIM> >& impl) {
            PROFILE_MEMBER_FUNC(Function);
            this->impl = impl;
        }


        /// Replace current FunctionImpl with a new one using the same parameters & map as f

        /// If zero is true the function is initialized to zero, otherwise it is empty
        template <typename R>
        void set_impl(const Function<R,NDIM>& f, bool zero = true) {
            impl = std::shared_ptr<implT>(new implT(*f.get_impl(), f.get_pmap(), zero));
        }

        /// Returns the world
        World& world() const {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            return  impl->world;
        }


        /// Returns a shared pointer to the process map
        const std::shared_ptr< WorldDCPmapInterface< Key<NDIM> > >& get_pmap() const {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            return impl->get_pmap();
        }


        /// Returns the square of the norm of the local function ... no communication

        /// Works in either basis
        double norm2sq_local() const {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            return impl->norm2sq_local();
        }


        /// Returns the 2-norm of the function ... global sum ... works in either basis

        /// See comments for err() w.r.t. applying to many functions.
        double norm2() const {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            if (VERIFY_TREE) verify_tree();
            double local = impl->norm2sq_local();

            impl->world.gop.sum(local);
            impl->world.gop.fence();
            return sqrt(local);
        }


        /// Initializes information about the function norm at all length scales
        void norm_tree(bool fence = true) const {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            if (VERIFY_TREE) verify_tree();
            if (is_compressed()) reconstruct();
            const_cast<Function<T,NDIM>*>(this)->impl->norm_tree(fence);
        }


        /// Compresses the function, transforming into wavelet basis.  Possible non-blocking comm.

        /// By default fence=true meaning that this operation completes before returning,
        /// otherwise if fence=false it returns without fencing and the user must invoke
        /// world.gop.fence() to assure global completion before using the function
        /// for other purposes.
        ///
        /// Noop if already compressed or if not initialized.
        ///
        /// Since reconstruction/compression do not discard information we define them
        /// as const ... "logical constness" not "bitwise constness".
        const Function<T,NDIM>& compress(bool fence = true) const {
            PROFILE_MEMBER_FUNC(Function);
            if (!impl || is_compressed()) return *this;
            if (VERIFY_TREE) verify_tree();
            const_cast<Function<T,NDIM>*>(this)->impl->compress(false, false, fence);
            return *this;
        }


        /// Compresses the function retaining scaling function coeffs.  Possible non-blocking comm.

        /// By default fence=true meaning that this operation completes before returning,
        /// otherwise if fence=false it returns without fencing and the user must invoke
        /// world.gop.fence() to assure global completion before using the function
        /// for other purposes.
        ///
        /// Noop if already compressed or if not initialized.
        void nonstandard(bool keepleaves, bool fence) {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            if (impl->is_nonstandard()) return;
            if (VERIFY_TREE) verify_tree();
            if (is_compressed()) reconstruct();
            impl->compress(true, keepleaves, fence);
        }

        /// Converts the function from nonstandard form to standard form.  Possible non-blocking comm.

        /// By default fence=true meaning that this operation completes before returning,
        /// otherwise if fence=false it returns without fencing and the user must invoke
        /// world.gop.fence() to assure global completion before using the function
        /// for other purposes.
        ///
        /// Must be already compressed.
        void standard(bool fence = true) {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            MADNESS_ASSERT(is_compressed());
            impl->standard(fence);
            if (fence && VERIFY_TREE) verify_tree();
        }

        /// Reconstructs the function, transforming into scaling function basis.  Possible non-blocking comm.

        /// By default fence=true meaning that this operation completes before returning,
        /// otherwise if fence=false it returns without fencing and the user must invoke
        /// world.gop.fence() to assure global completion before using the function
        /// for other purposes.
        ///
        /// Noop if already reconstructed or if not initialized.
        ///
        /// Since reconstruction/compression do not discard information we define them
        /// as const ... "logical constness" not "bitwise constness".
        void reconstruct(bool fence = true) const {
            PROFILE_MEMBER_FUNC(Function);
            if (!impl || !is_compressed()) return;
            const_cast<Function<T,NDIM>*>(this)->impl->reconstruct(fence);
            if (fence && VERIFY_TREE) verify_tree(); // Must be after in case nonstandard
        }


        /// Sums scaling coeffs down tree restoring state with coeffs only at leaves.  Optional fence.  Possible non-blocking comm.
        void sum_down(bool fence = true) const {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            MADNESS_ASSERT(!is_compressed());
            const_cast<Function<T,NDIM>*>(this)->impl->sum_down(fence);
            if (fence && VERIFY_TREE) verify_tree(); // Must be after in case nonstandard
        }


        /// Inplace autorefines the function.  Optional fence. Possible non-blocking comm.
        template <typename opT>
        void refine_general(const opT& op, bool fence = true) const {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            if (is_compressed()) reconstruct();
            impl->refine(op, fence);
        }


        struct autorefine_square_op {
            bool operator()(implT* impl, const Key<NDIM>& key, const Tensor<T>& t) const {
                return impl->autorefine_square_test(key, t);
            }

            template <typename Archive> void serialize (Archive& ar) {}
        };

        /// Inplace autorefines the function using same test as for squaring.
        void refine(bool fence = true) {
            refine_general(autorefine_square_op(), fence);
        }

        /// Inplace broadens support in scaling function basis
        void broaden(const BoundaryConditions<NDIM>& bc=FunctionDefaults<NDIM>::get_bc(),
                     bool fence = true) const {
            verify();
            reconstruct();
            impl->broaden(bc.is_periodic(), fence);
        }


        /// Get the scaling function coeffs at level n starting from NS form
        Tensor<T> coeffs_for_jun(Level n, long mode=0) {
            PROFILE_MEMBER_FUNC(Function);
            nonstandard(true,true);
            return impl->coeffs_for_jun(n,mode);
            //return impl->coeffs_for_jun(n);
        }


        /// Clears the function as if constructed uninitialized.  Optional fence.

        /// Any underlying data will not be freed until the next global fence.
        void clear(bool fence = true) {
            PROFILE_MEMBER_FUNC(Function);
            if (impl) {
                World& world = impl->world;
                impl.reset();
                if (fence) world.gop.fence();
            }
        }

        /// Process 0 prints a summary of all nodes in the tree (collective)
        void print_tree() const {
            PROFILE_MEMBER_FUNC(Function);
            if (impl) impl->print_tree();
        }

        /// Print a summary of the load balancing info

        /// This is serial and VERY expensive
        void print_info() const {
            PROFILE_MEMBER_FUNC(Function);
            if (impl) impl->print_info();
        }

        struct SimpleUnaryOpWrapper {
            T (*f)(T);
            SimpleUnaryOpWrapper(T (*f)(T)) : f(f) {}
            void operator()(const Key<NDIM>& key, Tensor<T>& t) const {
                UNARY_OPTIMIZED_ITERATOR(T, t, *_p0 = f(*_p0));
            }
            template <typename Archive> void serialize(Archive& ar) {}
        };

        /// Inplace unary operation on function values
        void unaryop(T (*f)(T)) {
            // Must fence here due to temporary object on stack
            // stopping us returning before complete
            this->unaryop(SimpleUnaryOpWrapper(f));
        }


        /// Inplace unary operation on function values
        template <typename opT>
        void unaryop(const opT& op, bool fence=true) {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            reconstruct();
            impl->unary_op_value_inplace(op, fence);
        }


        /// Unary operation applied inplace to the coefficients
        template <typename opT>
        void unaryop_coeff(const opT& op,
                           bool fence = true) {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            impl->unary_op_coeff_inplace(op, fence);
        }


        /// Unary operation applied inplace to the nodes
        template <typename opT>
        void unaryop_node(const opT& op,
                          bool fence = true) {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            impl->unary_op_node_inplace(op, fence);
        }


        static void doconj(const Key<NDIM>, Tensor<T>& t) {
            PROFILE_MEMBER_FUNC(Function);
            t.conj();
        }

        /// Inplace complex conjugate.  No communication except for optional fence.

        /// Returns this for chaining.  Works in either basis.
        Function<T,NDIM> conj(bool fence = true) {
            PROFILE_MEMBER_FUNC(Function);
            unaryop_coeff(&Function<T,NDIM>::doconj, fence);
            return *this;
        }


        /// Inplace, scale the function by a constant.  No communication except for optional fence.

        /// Works in either basis.  Returns reference to this for chaining.
        template <typename Q>
        Function<T,NDIM>& scale(const Q q, bool fence=true) {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            if (VERIFY_TREE) verify_tree();
            impl->scale_inplace(q,fence);
            return *this;
        }


        /// Inplace add scalar.  No communication except for optional fence.
        Function<T,NDIM>& add_scalar(T t, bool fence=true) {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            if (VERIFY_TREE) verify_tree();
            impl->add_scalar_inplace(t,fence);
            return *this;
        }


        /// Inplace, general bi-linear operation in wavelet basis.  No communication except for optional fence.

        /// If the functions are not in the wavelet basis an exception is thrown since this routine
        /// is intended to be fast and unexpected compression is assumed to be a performance bug.
        ///
        /// Returns this for chaining.
        ///
        /// this <-- this*alpha + other*beta
        template <typename Q, typename R>
        Function<T,NDIM>& gaxpy(const T& alpha,
                                const Function<Q,NDIM>& other, const R& beta, bool fence=true) {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            other.verify();
            MADNESS_ASSERT(is_compressed() && other.is_compressed());
            impl->gaxpy_inplace(alpha,*other.get_impl(),beta,fence);
            return *this;
        }


        /// Inplace addition of functions in the wavelet basis

        /// Using operator notation forces a global fence after every operation.
        /// Functions are compressed if not already so.
        template <typename Q>
        Function<T,NDIM>& operator+=(const Function<Q,NDIM>& other) {
            PROFILE_MEMBER_FUNC(Function);
            if (!is_compressed()) compress();
            if (!other.is_compressed()) other.compress();
            if (VERIFY_TREE) verify_tree();
            if (VERIFY_TREE) other.verify_tree();
            return gaxpy(T(1.0), other, Q(1.0), true);
        }


        /// Inplace subtraction of functions in the wavelet basis

        /// Using operator notation forces a global fence after every operation
        template <typename Q>
        Function<T,NDIM>& operator-=(const Function<Q,NDIM>& other) {
            PROFILE_MEMBER_FUNC(Function);
            if (!is_compressed()) compress();
            if (!other.is_compressed()) other.compress();
            if (VERIFY_TREE) verify_tree();
            if (VERIFY_TREE) other.verify_tree();
            return gaxpy(T(1.0), other, Q(-1.0), true);
        }


        /// Inplace scaling by a constant

        /// Using operator notation forces a global fence after every operation
        template <typename Q>
        typename IsSupported<TensorTypeData<Q>, Function<T,NDIM> >::type &
        operator*=(const Q q) {
            PROFILE_MEMBER_FUNC(Function);
            scale(q,true);
            return *this;
        }


        /// Inplace squaring of function ... global comm only if not reconstructed

        /// Returns *this for chaining.
        Function<T,NDIM>& square(bool fence = true) {
            PROFILE_MEMBER_FUNC(Function);
            if (is_compressed()) reconstruct();
            if (VERIFY_TREE) verify_tree();
            impl->square_inplace(fence);
            return *this;
        }

        /// Returns *this for chaining.
        Function<T,NDIM>& abs(bool fence = true) {
            PROFILE_MEMBER_FUNC(Function);
            if (is_compressed()) reconstruct();
            if (VERIFY_TREE) verify_tree();
            impl->abs_inplace(fence);
            return *this;
        }

        /// Returns *this for chaining.
        Function<T,NDIM>& abs_square(bool fence = true) {
            PROFILE_MEMBER_FUNC(Function);
            if (is_compressed()) reconstruct();
            if (VERIFY_TREE) verify_tree();
            impl->abs_square_inplace(fence);
            return *this;
        }

        /// Returns local contribution to \c int(f(x),x) ... no communication

        /// In the wavelet basis this is just the coefficient of the first scaling
        /// function which is a constant.  In the scaling function basis we
        /// must add up contributions from each box.
        T trace_local() const {
            PROFILE_MEMBER_FUNC(Function);
            if (!impl) return 0.0;
            if (VERIFY_TREE) verify_tree();
            return impl->trace_local();
        }


        /// Returns global value of \c int(f(x),x) ... global comm required
        T trace() const {
            PROFILE_MEMBER_FUNC(Function);
            if (!impl) return 0.0;
            T sum = impl->trace_local();
            impl->world.gop.sum(sum);
            impl->world.gop.fence();
            return sum;
        }


        /// Returns local part of inner product ... throws if both not compressed
        template <typename R>
        TENSOR_RESULT_TYPE(T,R) inner_local(const Function<R,NDIM>& g) const {
            PROFILE_MEMBER_FUNC(Function);
            MADNESS_ASSERT(is_compressed());
            MADNESS_ASSERT(g.is_compressed());
            if (VERIFY_TREE) verify_tree();
            if (VERIFY_TREE) g.verify_tree();
            return impl->inner_local(*(g.get_impl()));
        }


        /// Returns the inner product

        /// Not efficient for computing multiple inner products
        template <typename R>
        TENSOR_RESULT_TYPE(T,R) inner(const Function<R,NDIM>& g) const {
            PROFILE_MEMBER_FUNC(Function);
            if (!is_compressed()) compress();
            if (!g.is_compressed()) g.compress();
            if (VERIFY_TREE) verify_tree();
            if (VERIFY_TREE) g.verify_tree();

            TENSOR_RESULT_TYPE(T,R) local = impl->inner_local(*(g.get_impl()));
            impl->world.gop.sum(local);
            impl->world.gop.fence();
            return local;
        }


        /// Replaces this function with one loaded from an archive using the default processor map

        /// Archive can be sequential or parallel.
        ///
        /// The & operator for serializing will only work with parallel archives.
        template <typename Archive>
        void load(World& world, Archive& ar) {
            PROFILE_MEMBER_FUNC(Function);
            // Type checking since we are probably circumventing the archive's own type checking
            long magic, id, ndim, k;
            ar & magic & id & ndim & k;
            MADNESS_ASSERT(magic == 7776768); // Mellow Mushroom Pizza tel.# in Knoxville
            MADNESS_ASSERT(id == TensorTypeData<T>::id);
            MADNESS_ASSERT(ndim == NDIM);

            impl.reset(new implT(FunctionFactory<T,NDIM>(world).k(k).empty()));

            impl->load(ar);
        }


        /// Stores the function to an archive

        /// Archive can be sequential or parallel.
        ///
        /// The & operator for serializing will only work with parallel archives.
        template <typename Archive>
        void store(Archive& ar) const {
            PROFILE_MEMBER_FUNC(Function);
            verify();
            // For type checking, etc.
            ar & long(7776768) & long(TensorTypeData<T>::id) & long(NDIM) & long(k());

            impl->store(ar);
        }


        /// This is replaced with left*right ...  private
        template <typename Q, typename opT>
        Function<typename opT::resultT,NDIM>& unary_op_coeffs(const Function<Q,NDIM>& func,
                const opT& op, bool fence) {
            PROFILE_MEMBER_FUNC(Function);
            func.verify();
            MADNESS_ASSERT(!(func.is_compressed()));
            if (VERIFY_TREE) func.verify_tree();
            impl.reset(new implT(*func.get_impl(), func.get_pmap(), false));
            impl->unaryXX(func.get_impl().get(), op, fence);
            return *this;
        }

        /// Returns vector of FunctionImpl pointers corresponding to vector of functions
        template <typename Q, std::size_t D>
        static std::vector< std::shared_ptr< FunctionImpl<Q,D> > > vimpl(const std::vector< Function<Q,D> >& v) {
            PROFILE_MEMBER_FUNC(Function);
            std::vector< std::shared_ptr< FunctionImpl<Q,D> > > r(v.size());
            for (unsigned int i=0; i<v.size(); ++i) r[i] = v[i].get_impl();
            return r;
        }


        /// Multiplication of function * vector of functions using recursive algorithm of mulxx
        template <typename L, typename R>
        void vmulXX(const Function<L,NDIM>& left,
                    const std::vector< Function<R,NDIM> >& right,
                    std::vector< Function<T,NDIM> >& result,
                    double tol,
                    bool fence) {
            PROFILE_MEMBER_FUNC(Function);

            std::vector<FunctionImpl<T,NDIM>*> vresult(right.size());
            std::vector<const FunctionImpl<R,NDIM>*> vright(right.size());
            for (unsigned int i=0; i<right.size(); ++i) {
                result[i].set_impl(left,false);
                vresult[i] = result[i].impl.get();
                vright[i] = right[i].impl.get();
            }

            left.world().gop.fence(); // Is this still essential?  Yes.
            vresult[0]->mulXXvec(left.get_impl().get(), vright, vresult, tol, fence);
        }

        /// sparse transformation of a vector of functions ... private
        template <typename R, typename Q>
        void vtransform(const std::vector< Function<R,NDIM> >& v,
                        const Tensor<Q>& c,
                        std::vector< Function<T,NDIM> >& vresult,
                        double tol,
                        bool fence=true) {
            PROFILE_MEMBER_FUNC(Function);
            vresult[0].impl->vtransform(vimpl(v), c, vimpl(vresult), tol, fence);
        }

        /// This is replaced with alpha*left + beta*right ...  private
        template <typename L, typename R>
        Function<T,NDIM>& gaxpy_oop(T alpha, const Function<L,NDIM>& left,
                                    T beta,  const Function<R,NDIM>& right, bool fence) {
            PROFILE_MEMBER_FUNC(Function);
            left.verify();
            right.verify();
            MADNESS_ASSERT(left.is_compressed() && right.is_compressed());
            if (VERIFY_TREE) left.verify_tree();
            if (VERIFY_TREE) right.verify_tree();
            impl.reset(new implT(*left.get_impl(), left.get_pmap(), false));
            impl->gaxpy(alpha,*left.get_impl(),beta,*right.get_impl(),fence);
            return *this;
        }

        /// This is replaced with mapdim(f) ...  private
        Function<T,NDIM>& mapdim(const Function<T,NDIM>& f, const std::vector<long>& map, bool fence) {
            PROFILE_MEMBER_FUNC(Function);
            f.verify();
            if (VERIFY_TREE) f.verify_tree();
            for (std::size_t i=0; i<NDIM; ++i) MADNESS_ASSERT(map[i]>=0 && static_cast<std::size_t>(map[i])<NDIM);
            impl.reset(new implT(*f.impl, f.get_pmap(), false));
            impl->mapdim(*f.impl,map,fence);
            return *this;
        }
    };


    /// Returns new function equal to alpha*f(x) with optional fence
    template <typename Q, typename T, std::size_t NDIM>
    Function<TENSOR_RESULT_TYPE(Q,T),NDIM>
    mul(const Q alpha, const Function<T,NDIM>& f, bool fence=true) {
        PROFILE_FUNC;
        f.verify();
        if (VERIFY_TREE) f.verify_tree();
        Function<TENSOR_RESULT_TYPE(Q,T),NDIM> result;
        result.set_impl(f, false);
        result.get_impl()->scale_oop(alpha,*f.get_impl(),fence);
        return result;
    }


    /// Returns new function equal to f(x)*alpha with optional fence
    template <typename Q, typename T, std::size_t NDIM>
    Function<TENSOR_RESULT_TYPE(Q,T),NDIM>
    mul(const Function<T,NDIM>& f, const Q alpha, bool fence=true) {
        PROFILE_FUNC;
        return mul(alpha,f,fence);
    }


    /// Returns new function equal to f(x)*alpha

    /// Using operator notation forces a global fence after each operation
    template <typename Q, typename T, std::size_t NDIM>
    Function<TENSOR_RESULT_TYPE(Q,T),NDIM>
    operator*(const Function<T,NDIM>& f, const Q alpha) {
        return mul(alpha, f, true);
    }

    /// Returns new function equal to alpha*f(x)

    /// Using operator notation forces a global fence after each operation
    template <typename Q, typename T, std::size_t NDIM>
    Function<TENSOR_RESULT_TYPE(Q,T),NDIM>
    operator*(const Q alpha, const Function<T,NDIM>& f) {
        return mul(alpha, f, true);
    }

    /// Sparse multiplication --- left and right must be reconstructed and if tol!=0 have tree of norms already created
    template <typename L, typename R,std::size_t NDIM>
    Function<TENSOR_RESULT_TYPE(L,R),NDIM>
    mul_sparse(const Function<L,NDIM>& left, const Function<R,NDIM>& right, double tol, bool fence=true) {
        PROFILE_FUNC;
        left.verify();
        right.verify();
        MADNESS_ASSERT(!(left.is_compressed() || right.is_compressed()));
        if (VERIFY_TREE) left.verify_tree();
        if (VERIFY_TREE) right.verify_tree();

        Function<TENSOR_RESULT_TYPE(L,R),NDIM> result;
        result.set_impl(left, false);
        result.get_impl()->mulXX(left.get_impl().get(), right.get_impl().get(), tol, fence);
        return result;
    }

    /// Same as \c operator* but with optional fence and no automatic reconstruction
    template <typename L, typename R,std::size_t NDIM>
    Function<TENSOR_RESULT_TYPE(L,R),NDIM>
    mul(const Function<L,NDIM>& left, const Function<R,NDIM>& right, bool fence=true) {
        return mul_sparse(left,right,0.0,fence);
    }

    /// Generate new function = op(left,right) where op acts on the function values
    template <typename L, typename R, typename opT, std::size_t NDIM>
    Function<TENSOR_RESULT_TYPE(L,R),NDIM>
    binary_op(const Function<L,NDIM>& left, const Function<R,NDIM>& right, const opT& op, bool fence=true) {
        PROFILE_FUNC;
        if (left.is_compressed()) left.reconstruct();
        if (right.is_compressed()) right.reconstruct();

        Function<TENSOR_RESULT_TYPE(L,R),NDIM> result;
        result.set_impl(left, false);
        result.get_impl()->binaryXX(left.get_impl().get(), right.get_impl().get(), op, fence);
        return result;
    }

    /// Out of place application of unary operation to function values with optional fence
    template <typename Q, typename opT, std::size_t NDIM>
    Function<typename opT::resultT, NDIM>
    unary_op(const Function<Q,NDIM>& func, const opT& op, bool fence=true) {
        if (func.is_compressed()) func.reconstruct();
        Function<typename opT::resultT, NDIM> result;
        if (VERIFY_TREE) func.verify_tree();
        result.set_impl(func, false);
        result.get_impl()->unaryXXvalues(func.get_impl().get(), op, fence);
        return result;
    }


    /// Out of place application of unary operation to scaling function coefficients with optional fence
    template <typename Q, typename opT, std::size_t NDIM>
    Function<typename opT::resultT, NDIM>
    unary_op_coeffs(const Function<Q,NDIM>& func, const opT& op, bool fence=true) {
        if (func.is_compressed()) func.reconstruct();
        Function<typename opT::resultT, NDIM> result;
        return result.unary_op_coeffs(func,op,fence);
    }

    /// Use the vmra/mul(...) interface instead

    /// This so that we don't have to have friend functions in a different header.
    ///
    /// If using sparsity (tol != 0) you must have created the tree of norms
    /// already for both left and right.
    template <typename L, typename R, std::size_t D>
    std::vector< Function<TENSOR_RESULT_TYPE(L,R),D> >
    vmulXX(const Function<L,D>& left, const std::vector< Function<R,D> >& vright, double tol, bool fence=true) {
        if (vright.size() == 0) return std::vector< Function<TENSOR_RESULT_TYPE(L,R),D> >();
        std::vector< Function<TENSOR_RESULT_TYPE(L,R),D> > vresult(vright.size());
        vresult[0].vmulXX(left, vright, vresult, tol, fence);
        return vresult;
    }

    /// Multiplies two functions with the new result being of type TensorResultType<L,R>

    /// Using operator notation forces a global fence after each operation but also
    /// enables us to automatically reconstruct the input functions as required.
    template <typename L, typename R, std::size_t NDIM>
    Function<TENSOR_RESULT_TYPE(L,R), NDIM>
    operator*(const Function<L,NDIM>& left, const Function<R,NDIM>& right) {
        if (left.is_compressed())  left.reconstruct();
        if (right.is_compressed()) right.reconstruct();
        return mul(left,right,true);
    }


    /// Returns new function alpha*left + beta*right optional fence and no automatic compression
    template <typename L, typename R,std::size_t NDIM>
    Function<TENSOR_RESULT_TYPE(L,R),NDIM>
    gaxpy_oop(TENSOR_RESULT_TYPE(L,R) alpha, const Function<L,NDIM>& left,
              TENSOR_RESULT_TYPE(L,R) beta,  const Function<R,NDIM>& right, bool fence=true) {
        PROFILE_FUNC;
        Function<TENSOR_RESULT_TYPE(L,R),NDIM> result;
        return result.gaxpy_oop(alpha, left, beta, right, fence);
    }

    /// Same as \c operator+ but with optional fence and no automatic compression
    template <typename L, typename R,std::size_t NDIM>
    Function<TENSOR_RESULT_TYPE(L,R),NDIM>
    add(const Function<L,NDIM>& left, const Function<R,NDIM>& right, bool fence=true) {
        return gaxpy_oop(TENSOR_RESULT_TYPE(L,R)(1.0), left,
                         TENSOR_RESULT_TYPE(L,R)(1.0), right, fence);
    }


    /// Adds two functions with the new result being of type TensorResultType<L,R>

    /// Using operator notation forces a global fence after each operation
    template <typename L, typename R, std::size_t NDIM>
    Function<TENSOR_RESULT_TYPE(L,R), NDIM>
    operator+(const Function<L,NDIM>& left, const Function<R,NDIM>& right) {
        if (VERIFY_TREE) left.verify_tree();
        if (VERIFY_TREE) right.verify_tree();
        if (!left.is_compressed())  left.compress();
        if (!right.is_compressed()) right.compress();
        return add(left,right,true);
    }

    /// Same as \c operator- but with optional fence and no automatic compression
    template <typename L, typename R,std::size_t NDIM>
    Function<TENSOR_RESULT_TYPE(L,R),NDIM>
    sub(const Function<L,NDIM>& left, const Function<R,NDIM>& right, bool fence=true) {
        return gaxpy_oop(TENSOR_RESULT_TYPE(L,R)(1.0), left,
                         TENSOR_RESULT_TYPE(L,R)(-1.0), right, fence);
    }


    /// Subtracts two functions with the new result being of type TensorResultType<L,R>

    /// Using operator notation forces a global fence after each operation
    template <typename L, typename R, std::size_t NDIM>
    Function<TENSOR_RESULT_TYPE(L,R), NDIM>
    operator-(const Function<L,NDIM>& left, const Function<R,NDIM>& right) {
        PROFILE_FUNC;
        if (!left.is_compressed())  left.compress();
        if (!right.is_compressed()) right.compress();
        return sub(left,right,true);
    }


    /// Create a new function that is the square of f - global comm only if not reconstructed
    template <typename T, std::size_t NDIM>
    Function<T,NDIM> square(const Function<T,NDIM>& f, bool fence=true) {
        PROFILE_FUNC;
        Function<T,NDIM> result = copy(f,true);  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        return result.square(true); //fence);  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    }

    /// Create a new function that is the abs of f - global comm only if not reconstructed
    template <typename T, int NDIM>
    Function<T,NDIM> abs(const Function<T,NDIM>& f, bool fence=true) {
        PROFILE_FUNC;
        Function<T,NDIM> result = copy(f,true);  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        return result.abs(true); //fence);  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    }

    /// Create a new function that is the abs_square of f - global comm only if not reconstructed
    template <typename T, int NDIM>
    Function<T,NDIM> abs_square(const Function<T,NDIM>& f, bool fence=true) {
        PROFILE_FUNC;
        Function<T,NDIM> result = copy(f,true);  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        return result.abs_square(true); //fence);  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    }


    /// Create a new copy of the function with different distribution and optional fence

    /// Works in either basis.  Different distributions imply
    /// asynchronous communication and the optional fence is
    /// collective.
    template <typename T, std::size_t NDIM>
    Function<T,NDIM> copy(const Function<T,NDIM>& f,
                          const std::shared_ptr< WorldDCPmapInterface< Key<NDIM> > >& pmap,
                          bool fence = true) {
        PROFILE_FUNC;
        f.verify();
        Function<T,NDIM> result;
        typedef FunctionImpl<T,NDIM> implT;
        result.set_impl(std::shared_ptr<implT>(new implT(*f.get_impl(), pmap, false)));
        result.get_impl()->copy_coeffs(*f.get_impl(), fence);
        if (VERIFY_TREE) result.verify_tree();
        return result;
    }

    /// Create a new copy of the function with the same distribution and optional fence
    template <typename T, std::size_t NDIM>
    Function<T,NDIM> copy(const Function<T,NDIM>& f, bool fence = true) {
        PROFILE_FUNC;
        return copy(f, f.get_pmap(), fence);
    }

    /// Type conversion implies a deep copy.  No communication except for optional fence.

    /// Works in either basis but any loss of precision may result in different errors
    /// in applied in a different basis.
    ///
    /// The new function is formed with the options from the default constructor.
    ///
    /// There is no automatic type conversion since this is generally a rather dangerous
    /// thing and because there would be no way to make the fence optional.
    template <typename T, typename Q, std::size_t NDIM>
    Function<Q,NDIM> convert(const Function<T,NDIM>& f, bool fence = true) {
        PROFILE_FUNC;
        f.verify();
        Function<Q,NDIM> result;
        result.set_impl(f, false);
        result.get_impl()->copy_coeffs(*f.get_impl(), fence);
        return result;
    }


    /// Return the complex conjugate of the input function with the same distribution and optional fence

    /// !!! The fence is actually not optional in the current implementation !!!
    template <typename T, std::size_t NDIM>
    Function<T,NDIM> conj(const Function<T,NDIM>& f, bool fence = true) {
        PROFILE_FUNC;
        Function<T,NDIM> result = copy(f,true);
        return result.conj(fence);
    }

    /// Apply operator ONLY in non-standard form - required other steps missing !!
    template <typename opT, typename R, std::size_t NDIM>
    Function<TENSOR_RESULT_TYPE(typename opT::opT,R), NDIM>
    apply_only(const opT& op, const Function<R,NDIM>& f, bool fence=true) {
        PROFILE_FUNC;
        Function<TENSOR_RESULT_TYPE(typename opT::opT,R), NDIM> result;

        result.set_impl(f, true);
        result.get_impl()->apply(op, *f.get_impl(), op.get_bc().is_periodic(), fence);
        return result;
    }

    /// Apply operator in non-standard form

    /// Returns a new function with the same distribution
    ///
    /// !!! For the moment does NOT respect fence option ... always fences
    template <typename opT, typename R, std::size_t NDIM>
    Function<TENSOR_RESULT_TYPE(typename opT::opT,R), NDIM>
    apply(const opT& op, const Function<R,NDIM>& f, bool fence=true) {
        Function<R,NDIM>& ff = const_cast< Function<R,NDIM>& >(f);
        if (VERIFY_TREE) ff.verify_tree();
        ff.reconstruct();
        ff.nonstandard(op.doleaves, true);

        Function<TENSOR_RESULT_TYPE(typename opT::opT,R), NDIM> result = apply_only(op, ff, fence);

        ff.standard();
        result.reconstruct();
        return result;
    }


    template <typename opT, typename R, std::size_t NDIM>
    Function<TENSOR_RESULT_TYPE(typename opT::opT,R), NDIM>
    apply_1d_realspace_push(const opT& op, const Function<R,NDIM>& f, int axis, bool fence=true) {
        PROFILE_FUNC;
        Function<R,NDIM>& ff = const_cast< Function<R,NDIM>& >(f);
        if (VERIFY_TREE) ff.verify_tree();
        ff.reconstruct();

        Function<TENSOR_RESULT_TYPE(typename opT::opT,R), NDIM> result;

        result.set_impl(ff, false);
        result.get_impl()->apply_1d_realspace_push(op, ff.get_impl().get(), axis, fence);
        return result;
    }


    /// Generate a new function by reordering dimensions ... optional fence

    /// You provide an array of dimension NDIM that maps old to new dimensions
    /// according to
    /// \code
    ///    newdim = mapdim[olddim]
    /// \endcode
    /// Works in either scaling function or wavelet basis.
    ///
    /// Would be easy to modify this to also change the procmap here
    /// if desired but presently it uses the same procmap as f.
    template <typename T, std::size_t NDIM>
    Function<T,NDIM>
    mapdim(const Function<T,NDIM>& f, const std::vector<long>& map, bool fence=true) {
        PROFILE_FUNC;
        Function<T,NDIM> result;
        return result.mapdim(f,map,fence);
    }

    template <typename T, std::size_t NDIM>
    Function<T,NDIM>
    project(const Function<T,NDIM>& other,
            int k=FunctionDefaults<NDIM>::get_k(),
            double thresh=FunctionDefaults<NDIM>::get_thresh(),
            bool fence=true)
    {
        PROFILE_FUNC;
        Function<T,NDIM> result = FunctionFactory<T,NDIM>(other.world()).k(k).thresh(thresh).empty();
        other.reconstruct();
        result.get_impl()->project(*other.get_impl(),fence);
        return result;
    }


    /// Computes the scalar/inner product between two functions

    /// In Maple this would be \c int(conjugate(f(x))*g(x),x=-infinity..infinity)
    template <typename T, typename R, std::size_t NDIM>
    TENSOR_RESULT_TYPE(T,R) inner(const Function<T,NDIM>& f, const Function<R,NDIM>& g) {
        PROFILE_FUNC;
        return f.inner(g);
    }


    template <typename T, typename R, std::size_t NDIM>
    typename IsSupported<TensorTypeData<R>, Function<TENSOR_RESULT_TYPE(T,R),NDIM> >::type
    operator+(const Function<T,NDIM>& f, R r) {
        return (f*R(1.0)).add_scalar(r);
    }

    template <typename T, typename R, std::size_t NDIM>
    typename IsSupported<TensorTypeData<R>, Function<TENSOR_RESULT_TYPE(T,R),NDIM> >::type
    operator+(R r, const Function<T,NDIM>& f) {
        return (f*R(1.0)).add_scalar(r);
    }

    template <typename T, typename R, std::size_t NDIM>
    typename IsSupported<TensorTypeData<R>, Function<TENSOR_RESULT_TYPE(T,R),NDIM> >::type
    operator-(const Function<T,NDIM>& f, R r) {
        return (f*R(1.0)).add_scalar(-r);
    }

    template <typename T, typename R, std::size_t NDIM>
    typename IsSupported<TensorTypeData<R>, Function<TENSOR_RESULT_TYPE(T,R),NDIM> >::type
    operator-(R r, const Function<T,NDIM>& f) {
        return (f*R(-1.0)).add_scalar(r);
    }

    namespace detail {
        template <std::size_t NDIM>
        struct realop {
            typedef double resultT;
            Tensor<double> operator()(const Key<NDIM>& key, const Tensor<double_complex>& t) const {
                return real(t);
            }

            template <typename Archive> void serialize (Archive& ar) {}
        };

        template <std::size_t NDIM>
        struct imagop {
            typedef double resultT;
            Tensor<double> operator()(const Key<NDIM>& key, const Tensor<double_complex>& t) const {
                return imag(t);
            }

            template <typename Archive> void serialize (Archive& ar) {}
        };

        template <std::size_t NDIM>
        struct abssqop {
            typedef double resultT;
            Tensor<double> operator()(const Key<NDIM>& key, const Tensor<double_complex>& t) const {
                Tensor<double> r = abs(t);
                return r.emul(r);
            }

            template <typename Archive> void serialize (Archive& ar) {}
        };
    }

    /// Returns a new function that is the real part of the input
    template <std::size_t NDIM>
    Function<double,NDIM> real(const Function<double_complex,NDIM>& z, bool fence=true) {
        return unary_op_coeffs(z, detail::realop<NDIM>(), fence);
    }

    /// Returns a new function that is the imaginary part of the input
    template <std::size_t NDIM>
    Function<double,NDIM> imag(const Function<double_complex,NDIM>& z, bool fence=true) {
        return unary_op_coeffs(z, detail::imagop<NDIM>(), fence);
    }

    /// Returns a new function that is the square of the absolute value of the input
    template <std::size_t NDIM>
    Function<double,NDIM> abssq(const Function<double_complex,NDIM>& z, bool fence=true) {
        return unary_op(z, detail::abssqop<NDIM>(), fence);
    }


}

#include <mra/funcplot.h>

namespace madness {
    namespace archive {
        template <class T, std::size_t NDIM>
        struct ArchiveLoadImpl< ParallelInputArchive, Function<T,NDIM> > {
            static inline void load(const ParallelInputArchive& ar, Function<T,NDIM>& f) {
                f.load(*ar.get_world(), ar);
            }
        };

        template <class T, std::size_t NDIM>
        struct ArchiveStoreImpl< ParallelOutputArchive, Function<T,NDIM> > {
            static inline void store(const ParallelOutputArchive& ar, const Function<T,NDIM>& f) {
                f.store(ar);
            }
        };
    }


}
/* @} */

#include <mra/derivative.h>
#include <mra/operator.h>
#include <mra/functypedefs.h>
#include <mra/vmra.h>

#endif // MADNESS_MRA_MRA_H__INCLUDED
