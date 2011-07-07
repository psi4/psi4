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


  $Id: mraimpl.h 1602 2009-12-27 19:53:06Z rjharrison $
*/

#ifndef MADNESS_DERIVATIVE_H__INCLUDED
#define MADNESS_DERIVATIVE_H__INCLUDED

#include <world/world.h>
#include <iostream>
#include <world/print.h>
#include <misc/misc.h>
#include <tensor/tensor.h>
#include <mra/key.h>
#include <mra/funcdefaults.h>
#include <mra/funcimpl.h>
#include <mra/loadbal.h>

/// \file mra/derivative.h
/// \brief Declaration and initialization of tree traversal functions and generic derivative
/// \ingroup mra

namespace madness {

    /// Tri-diagonal operator traversing tree primarily for derivative operator

    /// \ingroup mra
    template <typename T, std::size_t NDIM>
    class DerivativeBase : public WorldObject< DerivativeBase<T, NDIM> > {
        typedef WorldObject< DerivativeBase<T, NDIM> > woT;
    protected:
        World& world;
        const std::size_t axis      ;  // Axis along which the operation is performed
        const int k         ;  // Number of wavelets of the function
        const BoundaryConditions<NDIM> bc;
        const std::vector<long> vk; ///< (k,...) used to initialize Tensors

    public:
        friend class FunctionImpl<T, NDIM>;

        typedef Tensor<T>               tensorT  ;
        typedef Key<NDIM>               keyT     ;
        typedef std::pair<keyT,tensorT> argT     ;
        typedef FunctionImpl<T,NDIM>    implT    ;
        typedef Function<T,NDIM>        functionT;
        typedef WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> > dcT;
        typedef FunctionNode<T,NDIM> nodeT;


        DerivativeBase(World& world, std::size_t axis, int k, BoundaryConditions<NDIM> bc)
            : WorldObject< DerivativeBase<T, NDIM> >(world)
            , world(world)
            , axis(axis)
            , k(k)
            , bc(bc)
            , vk(NDIM,k)
        {
            // No!  Cannot process incoming messages until the *derived* class is constructed.
            // this->process_pending();
        }

        virtual ~DerivativeBase() { }

        Void forward_do_diff1(const implT* f, implT* df, const keyT& key,
                              const std::pair<keyT,tensorT>& left,
                              const std::pair<keyT,tensorT>& center,
                              const std::pair<keyT,tensorT>& right)  const {

            const dcT& coeffs = f->get_coeffs();
            ProcessID owner = coeffs.owner(key);

            if (owner == world.rank()) {
                if (left.second.size() == 0) {
                    woT::task(owner, &madness::DerivativeBase<T,NDIM>::do_diff1,
                            f, df, key, find_neighbor(f, key,-1), center, right,
                            TaskAttributes::hipri());
                }
                else if (right.second.size() == 0) {
                    woT::task(owner, &madness::DerivativeBase<T,NDIM>::do_diff1,
                            f, df, key, left, center, find_neighbor(f, key,1),
                            TaskAttributes::hipri());
                }
                // Boundary node
                else if (left.first.is_invalid() || right.first.is_invalid()) {
                    woT::task(owner, &madness::DerivativeBase<T,NDIM>::do_diff2b,
                            f, df, key, left, center, right);
                }
                // Interior node
                else {
                    woT::task(owner, &madness::DerivativeBase<T,NDIM>::do_diff2i,
                            f, df, key, left, center, right);
                }
            }
            else {
                df->task(owner, &implT::do_diff1, this, f, key, left, center, right, TaskAttributes::hipri());
            }
            return None;
        }

        Void do_diff1(const implT* f, implT* df, const keyT& key,
                      const std::pair<keyT,tensorT>& left,
                      const std::pair<keyT,tensorT>& center,
                      const std::pair<keyT,tensorT>& right) const {
            MADNESS_ASSERT(axis<NDIM);

            if (left.second.size()==0 || right.second.size()==0) {
                // One of the neighbors is below us in the tree ... recur down
                df->get_coeffs().replace(key,nodeT(tensorT(),true));
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    const keyT& child = kit.key();
                    if ((child.translation()[axis]&1) == 0) {
                        // leftmost child automatically has right sibling
                        forward_do_diff1(f, df, child, left, center, center);
                    }
                    else {
                        // rightmost child automatically has left sibling
                        forward_do_diff1(f, df, child, center, center, right);
                    }
                }
            }
            else {
                forward_do_diff1(f, df, key, left, center, right);
            }
            return None;
        }

        virtual Void do_diff2b(const implT* f, implT* df, const keyT& key,
                               const std::pair<keyT,tensorT>& left,
                               const std::pair<keyT,tensorT>& center,
                               const std::pair<keyT,tensorT>& right) const = 0;

        virtual Void do_diff2i(const implT* f, implT* df, const keyT& key,
                               const std::pair<keyT,tensorT>& left,
                               const std::pair<keyT,tensorT>& center,
                               const std::pair<keyT,tensorT>& right) const = 0;


        /// Differentiate w.r.t. given coordinate (x=0, y=1, ...) with optional fence

        /// Returns a new function with the same distribution
        Function<T,NDIM>
        operator()(const functionT& f, bool fence=true) const {
            if (VERIFY_TREE) f.verify_tree();

            if (f.is_compressed()) {
                if (fence) {
                    f.reconstruct();
                }
                else {
                    MADNESS_EXCEPTION("diff: trying to diff a compressed function without fencing",0);
                }
            }

            functionT df;
            df.set_impl(f,false);

            df.get_impl()->diff(this, f.get_impl().get(), fence);
            return df;
        }


        static bool enforce_bc(int bc_left, int bc_right, Level n, Translation& l) {
            Translation two2n = 1ul << n;
            if (l < 0) {
                if (bc_left == BC_ZERO || bc_left == BC_FREE || bc_left == BC_DIRICHLET || bc_left == BC_ZERONEUMANN || bc_left == BC_NEUMANN) {
                    return false; // f=0 BC, or no BC, or nonzero f BC, or zero deriv BC, or nonzero deriv BC
                }
                else if (bc_left == BC_PERIODIC) {
                    l += two2n; // Periodic BC
                    MADNESS_ASSERT(bc_left == bc_right);   //check that both BCs are periodic
                }
                else {
                    MADNESS_EXCEPTION("enforce_bc: confused left BC?",bc_left);
                }
            }
            else if (l >= two2n) {
                if (bc_right == BC_ZERO || bc_right == BC_FREE || bc_right == BC_DIRICHLET || bc_right == BC_ZERONEUMANN || bc_right == BC_NEUMANN) {
                    return false; // f=0 BC, or no BC, or nonzero f BC, or zero deriv BC, or nonzero deriv BC
                }
                else if (bc_right == BC_PERIODIC) {
                    l -= two2n; // Periodic BC
                    MADNESS_ASSERT(bc_left == bc_right);   //check that both BCs are periodic
                }
                else {
                    MADNESS_EXCEPTION("enforce_bc: confused BC right?",bc_right);
                }
            }
            return true;
        }

        Key<NDIM> neighbor(const keyT& key, int step) const {
            Vector<Translation,NDIM> l = key.translation();
            l[axis] += step;
            if (!enforce_bc(bc(axis,0), bc(axis,1), key.level(), l[axis])) {
                return keyT::invalid();
            }
            else {
                return keyT(key.level(),l);
            }
        }

        Future< std::pair< Key<NDIM>,Tensor<T> > >
        find_neighbor(const implT* f, const Key<NDIM>& key, int step) const {
            keyT neigh = neighbor(key, step);
            if (neigh.is_invalid()) {
                return Future<argT>(argT(neigh,tensorT(vk))); // Zero bc
            }
            else {
                Future<argT> result;
                f->task(f->get_coeffs().owner(neigh), &implT::sock_it_to_me, neigh, result.remote_ref(world), TaskAttributes::hipri());
                return result;
            }
        }


        template <typename Archive> void serialize(const Archive& ar) const {
            throw "NOT IMPLEMENTED";
        }

    };  // End of the DerivativeBase class


    /// Implements derivatives operators with variety of boundary conditions on simulation domain
    template <typename T, std::size_t NDIM>
    class Derivative : public DerivativeBase<T, NDIM> {
    private:
        typedef DerivativeBase<T, NDIM> baseT;

    public:
        typedef Tensor<T>               tensorT  ;
        typedef Key<NDIM>               keyT     ;
        typedef std::pair<keyT,tensorT> argT     ;
        typedef FunctionImpl<T,NDIM>    implT    ;
        typedef Function<T,NDIM>        functionT;
        typedef WorldContainer< Key<NDIM> , FunctionNode<T, NDIM> > dcT;
        typedef FunctionNode<T,NDIM> nodeT;

    private:
        const functionT g1;  ///< Function describing the boundary condition on the right side
        const functionT g2;  ///< Function describing the boundary condition on the left side

        // Tensors for holding the modified coefficients
        Tensor<double> rm, r0, rp        ; ///< Blocks of the derivative operator
        Tensor<double> left_rm, left_r0  ; ///< Blocks of the derivative for the left boundary
        Tensor<double> right_r0, right_rp; ///< Blocks of the derivative for the right boundary
        Tensor<double> bv_left, bv_right ; ///< Blocks of the derivative operator for the boundary contribution

        Void do_diff2b(const implT* f, implT* df, const keyT& key,
                       const std::pair<keyT,tensorT>& left,
                       const std::pair<keyT,tensorT>& center,
                       const std::pair<keyT,tensorT>& right) const {
            Vector<Translation,NDIM> l = key.translation();
            double lev   = (double) key.level();

            tensorT d;

            //left boundary
            if (l[this->axis] == 0) {
                d = madness::inner(left_rm ,
                        df->parent_to_child(right.second, right.first,
                            baseT::neighbor(key,1)).swapdim(this->axis,0), 1, 0);
                inner_result(left_r0,
                             df->parent_to_child(center.second, center.first, key).swapdim(this->axis,0),
                             1, 0, d);
            }
            else {
                d = madness::inner(right_rp,
                                   df->parent_to_child(left.second, left.first, baseT::neighbor(key,-1)).swapdim(this->axis,0),
                                   1, 0);
                inner_result(right_r0,
                             df->parent_to_child(center.second, center.first, key).swapdim(this->axis,0),
                             1, 0, d);
            }
            if (this->axis) d = copy(d.swapdim(this->axis,0)); // make it contiguous
            d.scale(FunctionDefaults<NDIM>::get_rcell_width()[this->axis]*pow(2.0,lev));
            df->get_coeffs().replace(key,nodeT(d,false));


            // This is the boundary contribution (formally in BoundaryDerivative)
            int bc_left  = this->bc(this->axis,0);
            int bc_right = this->bc(this->axis,1);

            Future<argT> found_argT;
            tensorT bf, bdry_t;
            //left boundary
            if (l[this->axis] == 0) {
                if (bc_left != BC_PERIODIC && bc_left != BC_FREE && bc_left != BC_ZERO && bc_left != BC_ZERONEUMANN) {
                    bf = copy(bv_left);
                    found_argT = g1.get_impl()->find_me(key);
                }
                else {
                    return None;
                }
            }
            else { //right boundary
                if (bc_right != BC_PERIODIC && bc_right != BC_FREE && bc_right != BC_ZERO && bc_right != BC_ZERONEUMANN) {
                    bf = copy(bv_right);
                    found_argT = g2.get_impl()->find_me(key);
                }
                else {
                    return None;
                }
            }

            tensorT gcoeffs = df->parent_to_child(found_argT.get().second, found_argT.get().first,key);

            //if (this->bc.get_bc().dim(0) == 1) {
            if (NDIM == 1) {
                bdry_t = gcoeffs[0]*bf;
            }
            else {
                tensorT slice_aid(this->k);  //vector of zeros
                slice_aid[0] = 1;
                tensorT tmp = inner(slice_aid, gcoeffs, 0, this->axis);
                bdry_t = outer(bf,tmp);
                if (this->axis) bdry_t = copy(bdry_t.cycledim(this->axis,0,this->axis)); // make it contiguous
            }
            bdry_t.scale(FunctionDefaults<NDIM>::get_rcell_width()[this->axis]);

            if (l[this->axis]==0) {
                if (bc_left == BC_DIRICHLET)
                    bdry_t.scale( pow(2.0,lev));
				else if (bc_left ==BC_NEUMANN)
					bdry_t.scale(FunctionDefaults<NDIM>::get_cell_width()[this->axis]);
            }
            else {
                if (bc_right == BC_DIRICHLET)
                    bdry_t.scale( pow(2.0,lev));
				else if (bc_right ==BC_NEUMANN)
					bdry_t.scale(FunctionDefaults<NDIM>::get_cell_width()[this->axis]);
            }

            bdry_t = bdry_t + d;
            df->get_coeffs().replace(key,nodeT(bdry_t,false));

            return None;
        }

        Void do_diff2i(const implT* f, implT*df, const keyT& key,
                       const std::pair<keyT,tensorT>& left,
                       const std::pair<keyT,tensorT>& center,
                       const std::pair<keyT,tensorT>& right) const
        {
            tensorT d = madness::inner(rp,
                                       df->parent_to_child(left.second, left.first, baseT::neighbor(key,-1)).swapdim(this->axis,0),
                                       1, 0);
            inner_result(r0,
                         df->parent_to_child(center.second, center.first, key).swapdim(this->axis,0),
                         1, 0, d);
            inner_result(rm,
                         df->parent_to_child(right.second, right.first, baseT::neighbor(key,1)).swapdim(this->axis,0),
                         1, 0, d);
            if (this->axis) d = copy(d.swapdim(this->axis,0)); // make it contiguous
            d.scale(FunctionDefaults<NDIM>::get_rcell_width()[this->axis]*pow(2.0,(double) key.level()));
            df->get_coeffs().replace(key,nodeT(d,false));
            return None;
        }

        void initCoefficients()  {
            r0 = Tensor<double>(this->k,this->k);
            rp = Tensor<double>(this->k,this->k);
            rm = Tensor<double>(this->k,this->k);

            left_rm = Tensor<double>(this->k,this->k);
            left_r0 = Tensor<double>(this->k,this->k);

            right_r0 = Tensor<double>(this->k,this->k);
            right_rp = Tensor<double>(this->k,this->k);

            // These are the coefficients for the boundary contribution
            bv_left  = Tensor<double>(this->k);
            bv_right = Tensor<double>(this->k);


            int bc_left  = this->bc(this->axis,0);
            int bc_right = this->bc(this->axis,1);

            double kphase = -1.0;
            if (this->k%2 == 0) kphase = 1.0;
            double iphase = 1.0;
            for (int i=0; i<this->k; ++i) {
                double jphase = 1.0;
                for (int j=0; j<this->k; ++j) {
                    double gammaij = sqrt(double((2*i+1)*(2*j+1)));
                    double Kij;
                    if (((i-j)>0) && (((i-j)%2)==1))
                        Kij = 2.0;
                    else
                        Kij = 0.0;

                    r0(i,j) = 0.5*(1.0 - iphase*jphase - 2.0*Kij)*gammaij;
                    rm(i,j) = 0.5*jphase*gammaij;
                    rp(i,j) =-0.5*iphase*gammaij;

                    // Constraints on the derivative
                    if (bc_left == BC_ZERONEUMANN || bc_left == BC_NEUMANN) {
                        left_rm(i,j) = jphase*gammaij*0.5*(1.0 + iphase*kphase/this->k);

                        double phi_tmpj_left = 0;

                        for (int l=0; l<this->k; ++l) {
                            double gammalj = sqrt(double((2*l+1)*(2*j+1)));
                            double Klj;

                            if (((l-j)>0) && (((l-j)%2)==1))  Klj = 2.0;
                            else   Klj = 0.0;

                            phi_tmpj_left += sqrt(double(2*l+1))*Klj*gammalj;
                        }
                        phi_tmpj_left = -jphase*phi_tmpj_left;
                        left_r0(i,j) = (0.5*(1.0 + iphase*kphase/this->k) - Kij)*gammaij + iphase*sqrt(double(2*i+1))*phi_tmpj_left/pow(this->k,2.);
                    }
                    else if (bc_left == BC_ZERO || bc_left == BC_DIRICHLET || bc_left == BC_FREE) {
                        left_rm(i,j) = rm(i,j);

                        // B.C. with a function
                        if (bc_left == BC_ZERO || bc_left == BC_DIRICHLET)
                            left_r0(i,j) = (0.5 - Kij)*gammaij;

                        // No B.C.
                        else if (bc_left == BC_FREE)
                            left_r0(i,j) = (0.5 - iphase*jphase - Kij)*gammaij;
                    }

                    // Constraints on the derivative
                    if (bc_right == BC_ZERONEUMANN || bc_right == BC_NEUMANN) {
                        right_rp(i,j) = -0.5*(iphase + kphase / this->k)*gammaij;

                        double phi_tmpj_right = 0;
                        for (int l=0; l<this->k; ++l) {
                            double gammalj = sqrt(double((2*l+1)*(2*j+1)));
                            double Klj;
                            if (((l-j)>0) && (((l-j)%2)==1))  Klj = 2.0;
                            else   Klj = 0.0;
                            phi_tmpj_right += sqrt(double(2*l+1))*Klj*gammalj;
                        }
                        right_r0(i,j) = -(0.5*jphase*(iphase+ kphase/this->k) + Kij)*gammaij + sqrt(double(2*i+1))*phi_tmpj_right/pow(this->k,2.);
                    }
                    else if (bc_right == BC_ZERO || bc_right == BC_FREE || bc_right == BC_DIRICHLET) {
                        right_rp(i,j) = rp(i,j);

                        // Zero BC
                        if (bc_right == BC_ZERO || bc_right == BC_DIRICHLET)
                            right_r0(i,j) = -(0.5*iphase*jphase + Kij)*gammaij;

                        // No BC
                        else if (bc_right == BC_FREE)
                            right_r0(i,j) = (1.0 - 0.5*iphase*jphase - Kij)*gammaij;

                    }

                    jphase = -jphase;
                }
                iphase = -iphase;
            }

            // Coefficients for the boundary contributions
            iphase = 1.0;
            for (int i=0; i<this->k; ++i) {
                iphase = -iphase;

                if (bc_left == BC_DIRICHLET)
                    bv_left(i) = iphase*sqrt(double(2*i+1));            // vector for left dirichlet BC
                else if(bc_left == BC_NEUMANN)
                    bv_left(i) = -iphase*sqrt(double(2*i+1))/pow(this->k,2.);  // vector for left deriv BC
                else
                    bv_left(i) = 0.0;

                if (bc_right == BC_DIRICHLET)
                    bv_right(i) = sqrt(double(2*i+1));                  // vector for right dirichlet BC
                else if (bc_right == BC_NEUMANN)
                    bv_right(i) = sqrt(double(2*i+1))/pow(this->k,2.);         // vector for right deriv BC
                else
                    bv_right(i) = 0.0;
            }

            //print(rm.normf(),r0.normf(),rp.normf(),left_rm.normf(),left_r0.normf(),right_r0.normf(),right_rp.normf(),bv_left.normf(),bv_right.normf());
        }

    public:
        typedef T opT;

        /// Constructs a derivative operator

        /// @param world The world
        /// @param axis The direction to differentiate
        /// @param bc Boundary conditions (default from FunctionDefaults)
        /// @param g1 Function providing left boundary value (default empty)
        /// @param g2 Function providing right boundary value (default empty)
        /// @param k Wavelet order (default from FunctionDefaults)
        Derivative(World& world,
                   std::size_t axis,
                   const BoundaryConditions<NDIM>& bc=FunctionDefaults<NDIM>::get_bc(),
                   const functionT g1=functionT(),
                   const functionT g2=functionT(),
                   int k=FunctionDefaults<NDIM>::get_k())
            :  DerivativeBase<T, NDIM>(world, axis, k, bc)
            , g1(g1)
            , g2(g2)
        {
            MADNESS_ASSERT(axis<NDIM);
            initCoefficients();
            g1.reconstruct();
            g2.reconstruct();

            this->process_pending();
        }

        virtual ~Derivative() { }
    };


    /// Convenience function returning derivative operator with free-space boundary conditions
    template <typename T, std::size_t NDIM>
    Derivative<T,NDIM>
    free_space_derivative(World& world, int axis, int k=FunctionDefaults<NDIM>::get_k()) {
        return Derivative<T, NDIM>(world, axis, BoundaryConditions<NDIM>(BC_FREE), Function<T,NDIM>(), Function<T,NDIM>(), k);
    }


    /// Conveinence function returning derivative operator with periodic boundary conditions
    template <typename T, std::size_t NDIM>
    Derivative<T,NDIM>
    periodic_derivative(World& world, int axis, int k=FunctionDefaults<NDIM>::get_k()) {
        return Derivative<T, NDIM>(world, axis, BoundaryConditions<NDIM>(BC_PERIODIC), Function<T,NDIM>(), Function<T,NDIM>(), k);
    }

    /// Applies derivative operator to function (for syntactic equivalence to integral operator apply)
    template <typename T, std::size_t NDIM>
    Function<T,NDIM>
    apply(const Derivative<T,NDIM>& D, const Function<T,NDIM>& f, bool fence=true) {
        return D(f,fence);
    }

    /// Convenience function returning vector of derivative operators implementing grad (\f$ \nabla \f$)

    /// This will only work for BC_ZERO, BC_PERIODIC, BC_FREE and
    /// BC_ZERONEUMANN since we are not passing in any boundary
    /// functions.
    template <typename T, std::size_t NDIM>
    std::vector< std::shared_ptr< Derivative<T,NDIM> > >
    gradient_operator(World& world,
                      const BoundaryConditions<NDIM>& bc = FunctionDefaults<NDIM>::get_bc(),
                      int k = FunctionDefaults<NDIM>::get_k()) {
        std::vector< std::shared_ptr< Derivative<T,NDIM> > > r(NDIM);
        for (std::size_t d=0; d<NDIM; ++d) {
            MADNESS_ASSERT(bc(d,0)!=BC_DIRICHLET && bc(d,1)!=BC_DIRICHLET);
            MADNESS_ASSERT(bc(d,0)!=BC_NEUMANN   && bc(d,1)!=BC_NEUMANN);
            r[d].reset(new Derivative<T,NDIM>(world,d,bc,Function<T,NDIM>(),Function<T,NDIM>(),k));
        }
        return r;
    }


    namespace archive {
        template <class Archive, class T, std::size_t NDIM>
        struct ArchiveLoadImpl<Archive,const DerivativeBase<T,NDIM>*> {
            static void load(const Archive& ar, const DerivativeBase<T,NDIM>*& ptr) {
                WorldObject< DerivativeBase<T,NDIM> >* p;
                ar & p;
                ptr = static_cast< const DerivativeBase<T,NDIM>* >(p);
            }
        };

        template <class Archive, class T, std::size_t NDIM>
        struct ArchiveStoreImpl<Archive,const DerivativeBase<T,NDIM>*> {
            static void store(const Archive& ar, const DerivativeBase<T,NDIM>* const & ptr) {
                ar & ptr->id();
            }
        };
    }

}  // End of the madness namespace

#endif // MADNESS_MRA_DERIVATIVE_H_INCLUDED
