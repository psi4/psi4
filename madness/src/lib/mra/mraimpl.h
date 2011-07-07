
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


  $Id: mraimpl.h 2235 2011-03-22 18:25:06Z jakiej $
*/

#ifndef MADNESS_MRA_MRAIMPL_H__INCLUDED
#define MADNESS_MRA_MRAIMPL_H__INCLUDED

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <world/worldhashmap.h>
#include <math.h>

/// \file mra/mraimpl.h
/// \brief Declaration and initialization of static data, some implementation, some instantiation

namespace madness {

    // Definition and initialization of FunctionDefaults static members
    // It cannot be an instance of FunctionFactory since we want to
    // set the defaults independent of the data type.

    template <typename T, std::size_t NDIM>
    void FunctionCommonData<T,NDIM>::_init_twoscale() {
        if (! two_scale_hg(k, &hg)) throw "failed to get twoscale coefficients";
        hgT = transpose(hg);
        hgsonly = copy(hg(Slice(0,k-1),_));
    }

    template <typename T, std::size_t NDIM>
    void FunctionCommonData<T,NDIM>::_init_quadrature
    (int k, int npt, Tensor<double>& quad_x, Tensor<double>& quad_w,
     Tensor<double>& quad_phi, Tensor<double>& quad_phiw, Tensor<double>& quad_phit) {
        quad_x = Tensor<double>(npt);
        quad_w = Tensor<double>(npt);
        quad_phi = Tensor<double>(npt,k);
        quad_phiw = Tensor<double>(npt,k);

        gauss_legendre(npt,0.0,1.0,quad_x.ptr(),quad_w.ptr());
        for (int mu=0; mu<npt; ++mu) {
            double phi[200];
            legendre_scaling_functions(quad_x(mu),k,phi);
            for (int j=0; j<k; ++j) {
                quad_phi(mu,j) = phi[j];
                quad_phiw(mu,j) = quad_w(mu)*phi[j];
            }
        }
        quad_phit = transpose(quad_phi);
    }


    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::verify_tree() const {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        world.gop.fence();  // Make sure nothing is going on

        // Verify consistency of compression status, existence and size of coefficients,
        // and has_children() flag.
        for (typename dcT::const_iterator it=coeffs.begin(); it!=coeffs.end(); ++it) {
            const keyT& key = it->first;
            const nodeT& node = it->second;
            bool bad;

            if (is_compressed()) {
                if (node.has_children()) {
                    bad = node.coeff().dim(0) != 2*cdata.k;
                }
                else {
                    bad = node.coeff().size() != 0;
                }
            }
            else {
                if (node.has_children()) {
                    bad = node.coeff().size() != 0;
                }
                else {
                    bad = node.coeff().dim(0) != cdata.k;
                }
            }

            if (bad) {
                print(world.rank(), "FunctionImpl: verify: INCONSISTENT TREE NODE, key =", key, ", node =", node,
                      ", dim[0] =",node.coeff().dim(0),", compressed =",is_compressed());
                std::cout.flush();
                MADNESS_EXCEPTION("FunctionImpl: verify: INCONSISTENT TREE NODE", 0);
            }
        }

        // Ensure that parents and children exist appropriately
        for (typename dcT::const_iterator it=coeffs.begin(); it!=coeffs.end(); ++it) {
            const keyT& key = it->first;
            const nodeT& node = it->second;

            if (key.level() > 0) {
                const keyT parent = key.parent();
                typename dcT::const_iterator pit = coeffs.find(parent).get();
                if (pit == coeffs.end()) {
                    print(world.rank(), "FunctionImpl: verify: MISSING PARENT",key,parent);
                    std::cout.flush();
                    MADNESS_EXCEPTION("FunctionImpl: verify: MISSING PARENT", 0);
                }
                const nodeT& pnode = pit->second;
                if (!pnode.has_children()) {
                    print(world.rank(), "FunctionImpl: verify: PARENT THINKS IT HAS NO CHILDREN",key,parent);
                    std::cout.flush();
                    MADNESS_EXCEPTION("FunctionImpl: verify: PARENT THINKS IT HAS NO CHILDREN", 0);
                }
            }

            for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                typename dcT::const_iterator cit = coeffs.find(kit.key()).get();
                if (cit == coeffs.end()) {
                    if (node.has_children()) {
                        print(world.rank(), "FunctionImpl: verify: MISSING CHILD",key,kit.key());
                        std::cout.flush();
                        MADNESS_EXCEPTION("FunctionImpl: verify: MISSING CHILD", 0);
                    }
                }
                else {
                    if (! node.has_children()) {
                        print(world.rank(), "FunctionImpl: verify: UNEXPECTED CHILD",key,kit.key());
                        std::cout.flush();
                        MADNESS_EXCEPTION("FunctionImpl: verify: UNEXPECTED CHILD", 0);
                    }
                }
            }
        }

        world.gop.fence();
    }

    template <typename T, std::size_t NDIM>
    T FunctionImpl<T,NDIM>::eval_cube(Level n, coordT& x, const tensorT& c) const {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        const int k = cdata.k;
        double px[NDIM][k];
        T sum = T(0.0);

        for (std::size_t i=0; i<NDIM; ++i) legendre_scaling_functions(x[i],k,px[i]);

        if (NDIM == 1) {
            for (int p=0; p<k; ++p)
                sum += c(p)*px[0][p];
        }
        else if (NDIM == 2) {
            for (int p=0; p<k; ++p)
                for (int q=0; q<k; ++q)
                    sum += c(p,q)*px[0][p]*px[1][q];
        }
        else if (NDIM == 3) {
            for (int p=0; p<k; ++p)
                for (int q=0; q<k; ++q)
                    for (int r=0; r<k; ++r)
                        sum += c(p,q,r)*px[0][p]*px[1][q]*px[2][r];
        }
        else if (NDIM == 4) {
            for (int p=0; p<k; ++p)
                for (int q=0; q<k; ++q)
                    for (int r=0; r<k; ++r)
                        for (int s=0; s<k; ++s)
                            sum += c(p,q,r,s)*px[0][p]*px[1][q]*px[2][r]*px[3][s];
        }
        else if (NDIM == 5) {
            for (int p=0; p<k; ++p)
                for (int q=0; q<k; ++q)
                    for (int r=0; r<k; ++r)
                        for (int s=0; s<k; ++s)
                            for (int t=0; t<k; ++t)
                                sum += c(p,q,r,s,t)*px[0][p]*px[1][q]*px[2][r]*px[3][s]*px[4][t];
        }
        else if (NDIM == 6) {
            for (int p=0; p<k; ++p)
                for (int q=0; q<k; ++q)
                    for (int r=0; r<k; ++r)
                        for (int s=0; s<k; ++s)
                            for (int t=0; t<k; ++t)
                                for (int u=0; u<k; ++u)
                                    sum += c(p,q,r,s,t,u)*px[0][p]*px[1][q]*px[2][r]*px[3][s]*px[4][t]*px[5][u];
        }
        else {
            MADNESS_EXCEPTION("FunctionImpl:eval_cube:NDIM?",NDIM);
        }
        return sum*pow(2.0,0.5*NDIM*n)/sqrt(FunctionDefaults<NDIM>::get_cell_volume());
    }

    template <typename T, std::size_t NDIM>
    Void FunctionImpl<T,NDIM>::reconstruct_op(const keyT& key, const tensorT& s) {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        // Note that after application of an integral operator not all
        // siblings may be present so it is necessary to check existence
        // and if absent insert an empty leaf node.
        //
        // If summing the result of an integral operator (i.e., from
        // non-standard form) there will be significant scaling function
        // coefficients at all levels and possibly difference coefficients
        // in leaves, hence the tree may refine as a result.
        typename dcT::iterator it = coeffs.find(key).get();
        if (it == coeffs.end()) {
            coeffs.replace(key,nodeT(tensorT(),false));
            it = coeffs.find(key).get();
        }
        nodeT& node = it->second;

        // The integral operator will correctly connect interior nodes
        // to children but may leave interior nodes without coefficients
        // ... but they still need to sum down so just give them zeros
        if (node.has_children() && !node.has_coeff()) {
            node.set_coeff(tensorT(cdata.v2k));
        }

        if (node.has_children() || node.has_coeff()) { // Must allow for inconsistent state from transform, etc.
            tensorT d = node.coeff();
            if (d.size() == 0) d = tensorT(cdata.v2k);
            if (key.level() > 0) d(cdata.s0) += s; // -- note accumulate for NS summation
            d = unfilter(d);
            node.clear_coeff();
            node.set_has_children(true);
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                const keyT& child = kit.key();
                tensorT ss = copy(d(child_patch(child)));
                PROFILE_BLOCK(recon_send);
                woT::task(coeffs.owner(child), &implT::reconstruct_op, child, ss);
            }
        }
        else {
            if (key.level()) node.set_coeff(copy(s));
            else node.set_coeff(s);
        }
        return None;
    }

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::fcube(const keyT& key, T (*f)(const coordT&), const Tensor<double>& qx, tensorT& fval) const {
		fcube(key,typename FunctionFactory<T,NDIM>::FunctorInterfaceWrapper(f) , qx, fval);
	}

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::fcube(const keyT& key, const FunctionFunctorInterface<T,NDIM>& f, const Tensor<double>& qx, tensorT& fval) const {
    //~ template <typename T, std::size_t NDIM> template< typename FF>
    //~ void FunctionImpl<T,NDIM>::fcube(const keyT& key, const FF& f, const Tensor<double>& qx, tensorT& fval) const {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        const Vector<Translation,NDIM>& l = key.translation();
        const Level n = key.level();
        const double h = std::pow(0.5,double(n));
        coordT c; // will hold the point in user coordinates
        const int npt = qx.dim(0);

        const Tensor<double>& cell_width = FunctionDefaults<NDIM>::get_cell_width();
        const Tensor<double>& cell = FunctionDefaults<NDIM>::get_cell();

        if (NDIM == 1) {
            for (int i=0; i<npt; ++i) {
                c[0] = cell(0,0) + h*cell_width[0]*(l[0] + qx(i)); // x
                fval(i) = f(c);
            }
        }
        else if (NDIM == 2) {
            for (int i=0; i<npt; ++i) {
                c[0] = cell(0,0) + h*cell_width[0]*(l[0] + qx(i)); // x
                for (int j=0; j<npt; ++j) {
                    c[1] = cell(1,0) + h*cell_width[1]*(l[1] + qx(j)); // y
                    fval(i,j) = f(c);
                }
            }
        }
        else if (NDIM == 3) {
            for (int i=0; i<npt; ++i) {
                c[0] = cell(0,0) + h*cell_width[0]*(l[0] + qx(i)); // x
                for (int j=0; j<npt; ++j) {
                    c[1] = cell(1,0) + h*cell_width[1]*(l[1] + qx(j)); // y
                    for (int k=0; k<npt; ++k) {
                        c[2] = cell(2,0) + h*cell_width[2]*(l[2] + qx(k)); // z
                        fval(i,j,k) = f(c);
                    }
                }
            }
        }
        else if (NDIM == 4) {
            for (int i=0; i<npt; ++i) {
                c[0] = cell(0,0) + h*cell_width[0]*(l[0] + qx(i)); // x
                for (int j=0; j<npt; ++j) {
                    c[1] = cell(1,0) + h*cell_width[1]*(l[1] + qx(j)); // y
                    for (int k=0; k<npt; ++k) {
                        c[2] = cell(2,0) + h*cell_width[2]*(l[2] + qx(k)); // z
                        for (int m=0; m<npt; ++m) {
                            c[3] = cell(3,0) + h*cell_width[3]*(l[3] + qx(m)); // xx
                            fval(i,j,k,m) = f(c);
                        }
                    }
                }
            }
        }
        else if (NDIM == 5) {
            for (int i=0; i<npt; ++i) {
                c[0] = cell(0,0) + h*cell_width[0]*(l[0] + qx(i)); // x
                for (int j=0; j<npt; ++j) {
                    c[1] = cell(1,0) + h*cell_width[1]*(l[1] + qx(j)); // y
                    for (int k=0; k<npt; ++k) {
                        c[2] = cell(2,0) + h*cell_width[2]*(l[2] + qx(k)); // z
                        for (int m=0; m<npt; ++m) {
                            c[3] = cell(3,0) + h*cell_width[3]*(l[3] + qx(m)); // xx
                            for (int n=0; n<npt; ++n) {
                                c[4] = cell(4,0) + h*cell_width[4]*(l[4] + qx(n)); // yy
                                fval(i,j,k,m,n) = f(c);
                            }
                        }
                    }
                }
            }
        }
        else if (NDIM == 6) {
            for (int i=0; i<npt; ++i) {
                c[0] = cell(0,0) + h*cell_width[0]*(l[0] + qx(i)); // x
                for (int j=0; j<npt; ++j) {
                    c[1] = cell(1,0) + h*cell_width[1]*(l[1] + qx(j)); // y
                    for (int k=0; k<npt; ++k) {
                        c[2] = cell(2,0) + h*cell_width[2]*(l[2] + qx(k)); // z
                        for (int m=0; m<npt; ++m) {
                            c[3] = cell(3,0) + h*cell_width[3]*(l[3] + qx(m)); // xx
                            for (int n=0; n<npt; ++n) {
                                c[4] = cell(4,0) + h*cell_width[4]*(l[4] + qx(n)); // yy
                                for (int p=0; p<npt; ++p) {
                                    c[5] = cell(5,0) + h*cell_width[5]*(l[5] + qx(p)); // zz
                                    fval(i,j,k,m,n,p) = f(c);
                                }
                            }
                        }
                    }
                }
            }
        }
        else {
            MADNESS_EXCEPTION("FunctionImpl: fcube: confused about NDIM?",NDIM);
        }
    }

    template <typename T, std::size_t NDIM>
    Void FunctionImpl<T,NDIM>::project_refine_op(const keyT& key,
                                                 bool do_refine,
                                                 const std::vector<Vector<double,NDIM> >& specialpts) {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        if (do_refine && key.level() < max_refine_level) {

            // Restrict special points to this box
            std::vector<Vector<double,NDIM> > newspecialpts;
            if (key.level() < functor->special_level()) {
                for (unsigned int i = 0; i < specialpts.size(); ++i) {
                    coordT simpt;
                    user_to_sim(specialpts[i], simpt);
                    Key<NDIM> specialkey = simpt2key(simpt, key.level());
                    if (specialkey.is_neighbor_of(key)) {
                        newspecialpts.push_back(specialpts[i]);
                    }
                }
            }

            // If refining compute scaling function coefficients and
            // norm of difference coefficients
            tensorT r, s0;
            double dnorm = 0.0;
            if (newspecialpts.size() == 0) {
                // Make in r child scaling function coeffs at level n+1
                r = tensorT(cdata.v2k);
                for (KeyChildIterator<NDIM> it(key); it; ++it) {
                    const keyT& child = it.key();
                    r(child_patch(child)) = project(child);
                }
                // Filter then test difference coeffs at level n
                tensorT d = filter(r);
                if (truncate_on_project) s0 = copy(d(cdata.s0));
                d(cdata.s0) = T(0);
                dnorm = d.normf();
            }

            // If have special points always refine.  If don't have special points
            // refine if difference norm is big
            if (newspecialpts.size() > 0 || dnorm >=truncate_tol(thresh,key.level())) {
                coeffs.replace(key,nodeT(tensorT(),true)); // Insert empty node for parent
                for (KeyChildIterator<NDIM> it(key); it; ++it) {
                    const keyT& child = it.key();
                    ProcessID p;
                    if (FunctionDefaults<NDIM>::get_project_randomize()) {
                        p = world.random_proc();
                    }
                    else {
                        p = coeffs.owner(child);
                    }
                    PROFILE_BLOCK(proj_refine_send);
                    woT::task(p, &implT::project_refine_op, child, do_refine, newspecialpts);
                }
            }
            else {
                if (truncate_on_project) {
                    coeffs.replace(key,nodeT(s0,false));
                }
                else {
                    coeffs.replace(key,nodeT(tensorT(),true)); // Insert empty node for parent
                    for (KeyChildIterator<NDIM> it(key); it; ++it) {
                        const keyT& child = it.key();
                        coeffs.replace(child,nodeT(copy(r(child_patch(child))),false));
                    }
                }
            }
        }
        else {
            coeffs.replace(key,nodeT(project(key),false));
        }
        return None;
    }

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::add_scalar_inplace(T t, bool fence) {
        std::vector<long> v0(NDIM,0L);
        if (is_compressed()) {
            if (world.rank() == coeffs.owner(cdata.key0)) {
                typename dcT::iterator it = coeffs.find(cdata.key0).get();
                MADNESS_ASSERT(it != coeffs.end());
                nodeT& node = it->second;
                MADNESS_ASSERT(node.has_coeff());
                node.coeff()(v0) += t*sqrt(FunctionDefaults<NDIM>::get_cell_volume());
            }
        }
        else {
            for (typename dcT::iterator it=coeffs.begin(); it!=coeffs.end(); ++it) {
                Level n = it->first.level();
                nodeT& node = it->second;
                if (node.has_coeff()) {
                    node.coeff()(v0) += t*sqrt(FunctionDefaults<NDIM>::get_cell_volume()*pow(0.5,double(NDIM*n)));
                }
            }
        }
        if (fence) world.gop.fence();
    }

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::insert_zero_down_to_initial_level(const keyT& key) {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        if (compressed) initial_level = std::max(initial_level,1); // Otherwise zero function is confused
        if (coeffs.is_local(key)) {
            if (compressed) {
                if (key.level() == initial_level) {
                    coeffs.replace(key, nodeT(tensorT(), false));
                }
                else {
                    coeffs.replace(key, nodeT(tensorT(cdata.v2k), true));
                }
            }
            else {
                if (key.level()<initial_level) {
                    coeffs.replace(key, nodeT(tensorT(), true));
                }
                else {
                    coeffs.replace(key, nodeT(tensorT(cdata.vk), false));
                }
            }
        }
        if (key.level() < initial_level) {
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                insert_zero_down_to_initial_level(kit.key());
            }
        }

    }


    template <typename T, std::size_t NDIM>
    Future<bool> FunctionImpl<T,NDIM>::truncate_spawn(const keyT& key, double tol) {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        typename dcT::iterator it = coeffs.find(key).get();
        if (it == coeffs.end()) {
            // In a standard tree all children would exist but some ops (transform)
            // can leave the tree in a messy state.  Just make the missing node as an
            // empty leaf.
            coeffs.replace(key,nodeT());
            it = coeffs.find(key).get();
        }
        nodeT& node = it->second;
        if (node.has_children()) {
            std::vector< Future<bool> > v = future_vector_factory<bool>(1<<NDIM);
            int i=0;
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit,++i) {
                v[i] = woT::task(coeffs.owner(kit.key()), &implT::truncate_spawn, kit.key(), tol, TaskAttributes::generator());
            }
            return woT::task(world.rank(),&implT::truncate_op, key, tol, v);
        }
        else {
            // In compressed form leaves should not have coeffs ... however the
            // transform op could leave the tree with leaves that do have coeffs
            // in which case we want something sensible to happen
            //MADNESS_ASSERT(!node.has_coeff());
            if (node.has_coeff() && key.level()>1) {
                double dnorm = node.coeff().normf();
                if (dnorm < truncate_tol(tol,key)) {
                    node.clear_coeff();
                }
            }
            return Future<bool>(node.has_coeff());
        }
    }


    template <typename T, std::size_t NDIM>
    bool FunctionImpl<T,NDIM>::truncate_op(const keyT& key, double tol, const std::vector< Future<bool> >& v) {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        // If any child has coefficients, a parent cannot truncate
        for (int i=0; i<(1<<NDIM); ++i) if (v[i].get()) return true;
        nodeT& node = coeffs.find(key).get()->second;

        // Interior nodes should always have coeffs but transform might
        // leave empty interior nodes ... hence just force no coeffs to
        // be zero coeff unless it is a leaf.
        if (node.has_children() && !node.has_coeff()) node.set_coeff(tensorT(cdata.v2k));

        if (key.level() > 1) { // >1 rather >0 otherwise reconstruct might get confused
            double dnorm = node.coeff().normf();
            if (dnorm < truncate_tol(tol,key)) {
                node.clear_coeff();
                if (node.has_children()) {
                    node.set_has_children(false);
                    for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                        coeffs.erase(kit.key());
                    }
                }
            }
        }
        return node.has_coeff();
    }


    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::print_tree(Level maxlevel) const {
        if (world.rank() == 0) do_print_tree(cdata.key0, maxlevel);
        world.gop.fence();
        if (world.rank() == 0) std::cout.flush();
        world.gop.fence();
    }


    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::do_print_tree(const keyT& key, Level maxlevel) const {
        typename dcT::const_iterator it = coeffs.find(key).get();
        if (it == coeffs.end()) {
            //MADNESS_EXCEPTION("FunctionImpl: do_print_tree: null node pointer",0);
            for (int i=0; i<key.level(); ++i) std::cout << "  ";
            std::cout << key << "  missing --> " << coeffs.owner(key) << "\n";
        }
        else {
            const nodeT& node = it->second;
            for (int i=0; i<key.level(); ++i) std::cout << "  ";
            std::cout << key << "  " << node << " --> " << coeffs.owner(key) << "\n";
            if (key.level() < maxlevel  &&  node.has_children()) {
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    do_print_tree(kit.key(),maxlevel);
                }
            }
        }
    }

    template <typename T, std::size_t NDIM>
    Tensor<T> FunctionImpl<T,NDIM>::project(const keyT& key) const {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        MADNESS_ASSERT(cdata.npt == cdata.k); // only necessary due to use of fast transform
        tensorT fval(cdata.vq,false); // this will be the returned result
        tensorT work(cdata.vk,false); // initially evaluate the function in here
        tensorT workq(cdata.vq,false); // initially evaluate the function in here

        if (functor) {
            fcube(key,*functor,cdata.quad_x,work);
        }
        else {
            MADNESS_EXCEPTION("FunctionImpl: project: confusion about function?",0);
        }

        work.scale(sqrt(FunctionDefaults<NDIM>::get_cell_volume()*pow(0.5,double(NDIM*key.level()))));
        //return transform(work,cdata.quad_phiw);
        return fast_transform(work,cdata.quad_phiw,fval,workq);
    }

    template <typename T, std::size_t NDIM>
    Future<double> FunctionImpl<T,NDIM>::get_norm_tree_recursive(const keyT& key) const {
        if (coeffs.probe(key)) {
            return Future<double>(coeffs.find(key).get()->second.get_norm_tree());
        }
        MADNESS_ASSERT(key.level());
        keyT parent = key.parent();
        return woT::task(coeffs.owner(parent), &implT::get_norm_tree_recursive, parent, TaskAttributes::hipri());
    }


    template <typename T, std::size_t NDIM>
    Void FunctionImpl<T,NDIM>::sock_it_to_me(const keyT& key,
            const RemoteReference< FutureImpl< std::pair<keyT,tensorT> > >& ref) const {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        if (coeffs.probe(key)) {
            const nodeT& node = coeffs.find(key).get()->second;
            Future< std::pair<keyT,tensorT> > result(ref);
            if (node.has_coeff()) {
                //madness::print("sock found it with coeff",key);
                result.set(std::pair<keyT,tensorT>(key,node.coeff()));
            }
            else {
                //madness::print("sock found it without coeff",key);
                result.set(std::pair<keyT,tensorT>(key,tensorT()));
            }
        }
        else {
            keyT parent = key.parent();
            //madness::print("sock forwarding to parent",key,parent);
            PROFILE_BLOCK(sitome_send);
            woT::task(coeffs.owner(parent), &FunctionImpl<T,NDIM>::sock_it_to_me, parent, ref, TaskAttributes::hipri());
        }
        return None;
    }

// like sock_it_to_me, but it replaces empty node with averaged coeffs from further down the tree
    template <typename T, std::size_t NDIM>
    Void FunctionImpl<T,NDIM>::sock_it_to_me_too(const keyT& key,
            const RemoteReference< FutureImpl< std::pair<keyT,tensorT> > >& ref) const {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        if (coeffs.probe(key)) {
            const nodeT& node = coeffs.find(key).get()->second;
            Future< std::pair<keyT,tensorT> > result(ref);
            if (node.has_coeff()) {
               result.set(std::pair<keyT,tensorT>(key,node.coeff()));
            }
            else {
                result.set(std::pair<keyT,tensorT>(key,nodeT(project(key),false).coeff()));
            }
        }
        else {
            keyT parent = key.parent();
            PROFILE_BLOCK(sitome2_send);
            woT::task(coeffs.owner(parent), &FunctionImpl<T,NDIM>::sock_it_to_me_too, parent, ref, TaskAttributes::hipri());
        }
        return None;
    }


    template <typename T, std::size_t NDIM>
    Void FunctionImpl<T,NDIM>::eval(const Vector<double,NDIM>& xin,
                                    const keyT& keyin,
                                    const typename Future<T>::remote_refT& ref) {

        PROFILE_MEMBER_FUNC(FunctionImpl);
        // This is ugly.  We must figure out a clean way to use
        // owner computes rule from the container.
        Vector<double,NDIM> x = xin;
        keyT key = keyin;
        Vector<Translation,NDIM> l = key.translation();
        ProcessID me = world.rank();
        while (1) {
            ProcessID owner = coeffs.owner(key);
            if (owner != me) {
                PROFILE_BLOCK(eval_send);
                woT::task(owner, &implT::eval, x, key, ref, TaskAttributes::hipri());
                return None;
            }
            else {
                typename dcT::futureT fut = coeffs.find(key);
                typename dcT::iterator it = fut.get();
                nodeT& node = it->second;
                if (node.has_coeff()) {
                    Future<T>(ref).set(eval_cube(key.level(), x, node.coeff()));
                    return None;
                }
                else {
                    for (std::size_t i=0; i<NDIM; ++i) {
                        double xi = x[i]*2.0;
                        int li = int(xi);
                        if (li == 2) li = 1;
                        x[i] = xi - li;
                        l[i] = 2*l[i] + li;
                    }
                    key = keyT(key.level()+1,l);
                }
            }
        }
        //MADNESS_EXCEPTION("should not be here",0);
    }


    template <typename T, std::size_t NDIM>
    std::pair<bool,T>
    FunctionImpl<T,NDIM>::eval_local_only(const Vector<double,NDIM>& xin, Level maxlevel) {
        Vector<double,NDIM> x = xin;
        keyT key(0);
        Vector<Translation,NDIM> l = key.translation();
        const ProcessID me = world.rank();
        while (key.level() <= maxlevel) {
            if (coeffs.owner(key) == me) {
                typename dcT::futureT fut = coeffs.find(key);
                typename dcT::iterator it = fut.get();
                if (it != coeffs.end()) {
                    nodeT& node = it->second;
                    if (node.has_coeff()) {
                        return std::pair<bool,T>(true,eval_cube(key.level(), x, node.coeff()));
                    }
                }
            }
            for (std::size_t i=0; i<NDIM; ++i) {
                double xi = x[i]*2.0;
                int li = int(xi);
                if (li == 2) li = 1;
                x[i] = xi - li;
                l[i] = 2*l[i] + li;
            }
            key = keyT(key.level()+1,l);
        }
        return std::pair<bool,T>(false,0.0);
    }

    template <typename T, std::size_t NDIM>
    Void FunctionImpl<T,NDIM>::evaldepthpt(const Vector<double,NDIM>& xin,
                                const keyT& keyin,
                                const typename Future<Level>::remote_refT& ref) {

        PROFILE_MEMBER_FUNC(FunctionImpl);
        // This is ugly.  We must figure out a clean way to use
        // owner computes rule from the container.
        Vector<double,NDIM> x = xin;
        keyT key = keyin;
        Vector<Translation,NDIM> l = key.translation();
        ProcessID me = world.rank();
        while (1) {
            ProcessID owner = coeffs.owner(key);
            if (owner != me) {
                PROFILE_BLOCK(eval_send);
                woT::task(owner, &implT::evaldepthpt, x, key, ref, TaskAttributes::hipri());
                return None;
            }
            else {
                typename dcT::futureT fut = coeffs.find(key);
                typename dcT::iterator it = fut.get();
                nodeT& node = it->second;
                if (node.has_coeff()) {
                    Future<Level>(ref).set(key.level());
                    return None;
                }
                else {
                    for (std::size_t i=0; i<NDIM; ++i) {
                        double xi = x[i]*2.0;
                        int li = int(xi);
                        if (li == 2) li = 1;
                        x[i] = xi - li;
                        l[i] = 2*l[i] + li;
                    }
                    key = keyT(key.level()+1,l);
                }
            }
        }
        //MADNESS_EXCEPTION("should not be here",0);
    }

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::tnorm(const tensorT& t, double* lo, double* hi) const {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        // Chosen approach looks stupid but it is more accurate
        // than the simple approach of summing everything and
        // subtracting off the low-order stuff to get the high
        // order (assuming the high-order stuff is small relative
        // to the low-order)
        tensorT work = copy(t);
        tensorT tlo = work(cdata.sh);
        *lo = tlo.normf();
        tlo.fill(0.0);
        *hi = work.normf();
    }

    namespace detail {
        template <typename A, typename B>
        struct noop {
            void operator()(const A& a, const B& b) const {};

            template <typename Archive> void serialize(Archive& ar) {}
        };

        template <typename T, std::size_t NDIM>
        struct scaleinplace {
            T q;
            scaleinplace() {}
	    // G++ 4.1.2 ICEs on BGP ... scaleinplace(T q) : q(q) {}
            scaleinplace(T q) {this->q = q;}
            void operator()(const Key<NDIM>& key, Tensor<T>& t) const {
                t.scale(q);
            }
            template <typename Archive> void serialize(Archive& ar) {
                ar & q;
            }
        };

        template <typename T, std::size_t NDIM>
        struct squareinplace {
            void operator()(const Key<NDIM>& key, Tensor<T>& t) const {
                t.emul(t);
            }
            template <typename Archive> void serialize(Archive& ar) {}
        };
          
        template <typename T, std::size_t NDIM>
        struct absinplace {
            void operator()(const Key<NDIM>& key, Tensor<T>& t) const {abs(t);}
            template <typename Archive> void serialize(Archive& ar) {}
        };      

        template <typename T, std::size_t NDIM>
        struct abssquareinplace {
            void operator()(const Key<NDIM>& key, Tensor<T>& t) const {abs(t.emul(t));}
            template <typename Archive> void serialize(Archive& ar) {}
        };      
    }

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::scale_inplace(const T q, bool fence) {
        unary_op_coeff_inplace(detail::scaleinplace<T,NDIM>(q), fence);
    }

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::square_inplace(bool fence) {
        //unary_op_value_inplace(&implT::autorefine_square_test, detail::squareinplace<T,NDIM>(), fence);
        unary_op_value_inplace(detail::squareinplace<T,NDIM>(), fence);
    }

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::abs_inplace(bool fence) {
        unary_op_value_inplace(detail::absinplace<T,NDIM>(), fence);
    }

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::abs_square_inplace(bool fence) {
        unary_op_value_inplace(detail::abssquareinplace<T,NDIM>(), fence);
    }

    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::phi_for_mul(Level np, Translation lp, Level nc, Translation lc, Tensor<double>& phi) const {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        double p[200];
        double scale = pow(2.0,double(np-nc));
        for (int mu=0; mu<cdata.npt; ++mu) {
            double xmu = scale*(cdata.quad_x(mu)+lc) - lp;
            MADNESS_ASSERT(xmu>-1e-15 && xmu<(1+1e-15));
            legendre_scaling_functions(xmu,cdata.k,p);
            for (int i=0; i<k; ++i) phi(i,mu) = p[i];
        }
        phi.scale(pow(2.0,0.5*np));
    }

    template <typename T, std::size_t NDIM>
    const Tensor<T> FunctionImpl<T,NDIM>::parent_to_child(const tensorT& s, const keyT& parent, const keyT& child) const {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        // An invalid parent/child means that they are out of the box
        // and it is the responsibility of the caller to worry about that
        // ... most likely the coefficients (s) are zero to reflect
        // zero B.C. so returning s makes handling this easy.
        if (parent == child || parent.is_invalid() || child.is_invalid()) return s;

        tensorT result = fcube_for_mul<T>(child, parent, s);
        result.scale(sqrt(FunctionDefaults<NDIM>::get_cell_volume()*pow(0.5,double(NDIM*child.level()))));
        result = transform(result,cdata.quad_phiw);

        return result;
    }


    template <typename T, std::size_t NDIM>
    T FunctionImpl<T,NDIM>::trace_local() const {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        std::vector<long> v0(NDIM,0);
        T sum = 0.0;
        if (compressed) {
            if (world.rank() == coeffs.owner(cdata.key0)) {
                typename dcT::const_iterator it = coeffs.find(cdata.key0).get();
                if (it != coeffs.end()) {
                    const nodeT& node = it->second;
                    if (node.has_coeff()) sum = node.coeff()(v0);
                }
            }
        }
        else {
            for (typename dcT::const_iterator it=coeffs.begin(); it!=coeffs.end(); ++it) {
                const keyT& key = it->first;
                const nodeT& node = it->second;
                if (node.has_coeff()) sum += node.coeff()(v0)*pow(0.5,NDIM*key.level()*0.5);
            }
        }
        return sum*sqrt(FunctionDefaults<NDIM>::get_cell_volume());
    }


    static inline bool enforce_bc(bool is_periodic, Level n, Translation& l) {
        Translation two2n = 1ul << n;
        if (l < 0) {
            if (is_periodic)
                l += two2n; // Periodic BC
            else
                return false; // Zero BC
        }
        else if (l >= two2n) {
            if (is_periodic)
                l -= two2n; // Periodic BC
            else
                return false; // Zero BC
        }
        return true;
    }

    template <typename T, std::size_t NDIM>
    Key<NDIM> FunctionImpl<T,NDIM>::neighbor(const keyT& key, const Key<NDIM>& disp, const std::vector<bool>& is_periodic) const {
        Vector<Translation,NDIM> l = key.translation();

        for (std::size_t axis=0; axis<NDIM; ++axis) {
            l[axis] += disp.translation()[axis];

            //if (!enforce_bc(bc(axis,0), bc(axis,1), key.level(), l[axis])) {
            if (!enforce_bc(is_periodic[axis], key.level(), l[axis])) {
                return keyT::invalid();
            }
        }
        return keyT(key.level(),l);
    }


   template <typename T, std::size_t NDIM>
    Future< std::pair< Key<NDIM>,Tensor<T> > >
    FunctionImpl<T,NDIM>::find_me(const Key<NDIM>& key) const {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        typedef std::pair< Key<NDIM>,Tensor<T> > argT;
        Future<argT> result;
        PROFILE_BLOCK(find_me_send);
        woT::task(coeffs.owner(key), &implT::sock_it_to_me_too, key, result.remote_ref(world), TaskAttributes::hipri());
        return result;
    }


    template <typename T, std::size_t NDIM>
    void FunctionImpl<T,NDIM>::mapdim(const implT& f, const std::vector<long>& map, bool fence) {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        for (typename dcT::const_iterator it=f.coeffs.begin(); it!=f.coeffs.end(); ++it) {
            const keyT& key = it->first;
            const nodeT& node = it->second;

            Vector<Translation,NDIM> l;
            for (std::size_t i=0; i<NDIM; ++i) l[map[i]] = key.translation()[i];
            tensorT c = node.coeff();
            if (c.size()) c = copy(c.mapdim(map));

            coeffs.replace(keyT(key.level(),l), nodeT(c,node.has_children()));
        }
        if (fence) world.gop.fence();
    }

    template <typename T, std::size_t NDIM>
    Future< Tensor<T> > FunctionImpl<T,NDIM>::compress_spawn(const Key<NDIM>& key, bool nonstandard, bool keepleaves) {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        MADNESS_ASSERT(coeffs.probe(key));
        nodeT& node = coeffs.find(key).get()->second;
        if (node.has_children()) {
            std::vector< Future<tensorT> > v = future_vector_factory<tensorT>(1<<NDIM);
            int i=0;
            for (KeyChildIterator<NDIM> kit(key); kit; ++kit,++i) {
                PROFILE_BLOCK(compress_send);
                v[i] = woT::task(coeffs.owner(kit.key()), &implT::compress_spawn, kit.key(), nonstandard, keepleaves, TaskAttributes::hipri());
            }
            return woT::task(world.rank(),&implT::compress_op, key, v, nonstandard);
        }
        else {
            Future<tensorT> result(node.coeff());
            if (!keepleaves) node.clear_coeff();
            return result;
        }
    }

    template <typename T, std::size_t NDIM>
    Void FunctionImpl<T,NDIM>::plot_cube_kernel(archive::archive_ptr< Tensor<T> > ptr,
                                                const keyT& key,
                                                const coordT& plotlo, const coordT& plothi, const std::vector<long>& npt,
                                                bool eval_refine) const {

        Tensor<T>& r = *ptr;

        coordT h; // Increment between points in each dimension
        for (std::size_t i=0; i<NDIM; ++i) {
            if (npt[i] > 1) {
                h[i] = (plothi[i]-plotlo[i])/(npt[i]-1);
            }
            else {
                MADNESS_ASSERT(plotlo[i] == plothi[i]);
                h[i] = 0.0;
            }
        }

        const Level n = key.level();
        const Vector<Translation,NDIM>& l = key.translation();
        const double twon = pow(2.0,double(n));
        const tensorT& coeff = coeffs.find(key).get()->second.coeff(); // Ugh!
        long ind[NDIM];
        coordT x;

        coordT boxlo, boxhi;
        Vector<int,NDIM> boxnpt;
        double fac = pow(0.5,double(key.level()));
        int npttotal = 1;
        for (std::size_t d=0; d<NDIM; ++d) {
            // Coords of box
            boxlo[d] = fac*key.translation()[d];
            boxhi[d] = boxlo[d]+fac;

            if (boxlo[d] > plothi[d] || boxhi[d] < plotlo[d]) {
                // Discard boxes out of the plot range
                npttotal = boxnpt[d] = 0;
                //print("OO range?");
                break;
            }
            else if (npt[d] == 1) {
                // This dimension is only a single point
                boxlo[d] = boxhi[d] = plotlo[d];
                boxnpt[d] = 1;
            }
            else {
                // Restrict to plot range
                boxlo[d] = std::max(boxlo[d],plotlo[d]);
                boxhi[d] = std::min(boxhi[d],plothi[d]);

                // Round lo up to next plot point; round hi down
                double xlo = long((boxlo[d]-plotlo[d])/h[d])*h[d] + plotlo[d];
                if (xlo < boxlo[d]) xlo += h[d];
                boxlo[d] =  xlo;
                double xhi = long((boxhi[d]-plotlo[d])/h[d])*h[d] + plotlo[d];
                if (xhi > boxhi[d]) xhi -= h[d];
                // MADNESS_ASSERT(xhi >= xlo);  // nope
                boxhi[d] = xhi;
                boxnpt[d] = long(round((boxhi[d] - boxlo[d])/h[d])) + 1;
            }
            npttotal *= boxnpt[d];
        }
        //print("    box", boxlo, boxhi, boxnpt, npttotal);
        if (npttotal > 0) {
            for (IndexIterator it(boxnpt); it; ++it) {
                for (std::size_t d=0; d<NDIM; ++d) {
                    double xd = boxlo[d] + it[d]*h[d]; // Sim. coords of point
                    x[d] = twon*xd - l[d]; // Offset within box
                    MADNESS_ASSERT(x[d]>=0.0 && x[d] <=1.0);  // sanity
                    if (npt[d] > 1) {
                        ind[d] = long(round((xd-plotlo[d])/h[d])); // Index of plot point
                    }
                    else {
                        ind[d] = 0;
                    }
                    MADNESS_ASSERT(ind[d]>=0 && ind[d]<npt[d]); // sanity
                }
                if (eval_refine) {
                    r(ind) = n;
                }
                else {
                    T tmp = eval_cube(n, x, coeff);
                    r(ind) = tmp;
                    //print("    eval", ind, tmp, r(ind));
                }
            }
        }

        return None;
    }

    /// Set plot_refine=true to get a plot of the refinment levels of
    /// the given function (defaulted to false in prototype).
    template <typename T, std::size_t NDIM>
    Tensor<T> FunctionImpl<T,NDIM>::eval_plot_cube(const coordT& plotlo,
                                                   const coordT& plothi,
                                                   const std::vector<long>& npt,
                                                   const bool eval_refine) const {
        PROFILE_MEMBER_FUNC(FunctionImpl);
        Tensor<T> r(NDIM, &npt[0]);
        //r(___) = 99.0;
        MADNESS_ASSERT(!compressed);

        for (typename dcT::const_iterator it=coeffs.begin(); it!=coeffs.end(); ++it) {
            const keyT& key = it->first;
            const nodeT& node = it->second;
            if (node.has_coeff()) {
                woT::task(world.rank(), &implT::plot_cube_kernel,
                    archive::archive_ptr< Tensor<T> >(&r), key, plotlo, plothi, npt, eval_refine);
            }
        }

        //        ITERATOR(r, if (r(IND) == 99.0) {print("BAD", IND); error("bad",0);});

        world.taskq.fence();
        world.gop.sum(r.ptr(), r.size());
        world.gop.fence();

        return r;
    }

    static inline void dxprintvalue(FILE* f, const double t) {
        fprintf(f,"%.6e\n",t);
    }

    static inline void dxprintvalue(FILE* f, const double_complex& t) {
        fprintf(f,"%.6e %.6e\n", t.real(), t.imag());
    }

    template <typename T, std::size_t NDIM>
    void plotdx(const Function<T,NDIM>& function,
                const char* filename,
                const Tensor<double>& cell,
                const std::vector<long>& npt,
                bool binary) {
        PROFILE_FUNC;
        MADNESS_ASSERT(NDIM<=6);
        const char* element[6] = {"lines","quads","cubes","cubes4D","cubes5D","cubes6D"};

        function.verify();
        World& world = const_cast< Function<T,NDIM>& >(function).world();
        FILE *f=0;
        if (world.rank() == 0) {
            f = fopen(filename, "w");
            if (!f) MADNESS_EXCEPTION("plotdx: failed to open the plot file", 0);

            fprintf(f,"object 1 class gridpositions counts ");
            for (std::size_t d=0; d<NDIM; ++d) fprintf(f," %ld",npt[d]);
            fprintf(f,"\n");

            fprintf(f,"origin ");
            for (std::size_t d=0; d<NDIM; ++d) fprintf(f, " %.6e", cell(d,0));
            fprintf(f,"\n");

            for (std::size_t d=0; d<NDIM; ++d) {
                fprintf(f,"delta ");
                for (std::size_t c=0; c<d; ++c) fprintf(f, " 0");
                double h = 0.0;
                if (npt[d]>1) h = (cell(d,1)-cell(d,0))/(npt[d]-1);
                fprintf(f," %.6e", h);
                for (std::size_t c=d+1; c<NDIM; ++c) fprintf(f, " 0");
                fprintf(f,"\n");
            }
            fprintf(f,"\n");

            fprintf(f,"object 2 class gridconnections counts ");
            for (std::size_t d=0; d<NDIM; ++d) fprintf(f," %ld",npt[d]);
            fprintf(f,"\n");
            fprintf(f, "attribute \"element type\" string \"%s\"\n", element[NDIM-1]);
            fprintf(f, "attribute \"ref\" string \"positions\"\n");
            fprintf(f,"\n");

            int npoint = 1;
            for (std::size_t d=0; d<NDIM; ++d) npoint *= npt[d];
            const char* iscomplex = "";
            if (TensorTypeData<T>::iscomplex) iscomplex = "category complex";
            const char* isbinary = "";
            if (binary) isbinary = "binary";
            fprintf(f,"object 3 class array type double %s rank 0 items %d %s data follows\n",
                    iscomplex, npoint, isbinary);
        }

        world.gop.fence();
        Tensor<T> r = function.eval_cube(cell, npt);

        if (world.rank() == 0) {
            if (binary) {
                // This assumes that the values are double precision
                fflush(f);
                fwrite((void *) r.ptr(), sizeof(T), r.size(), f);
                fflush(f);
            }
            else {
                for (IndexIterator it(npt); it; ++it) {
                    //fprintf(f,"%.6e\n",r(*it));
                    dxprintvalue(f,r(*it));
                }
            }
            fprintf(f,"\n");

            fprintf(f,"object \"%s\" class field\n",filename);
            fprintf(f,"component \"positions\" value 1\n");
            fprintf(f,"component \"connections\" value 2\n");
            fprintf(f,"component \"data\" value 3\n");
            fprintf(f,"\nend\n");
            fclose(f);
        }
        world.gop.fence();
    }

    template <std::size_t NDIM>
    void FunctionDefaults<NDIM>::set_defaults(World& world) {
        k = 6;
        thresh = 1e-4;
        initial_level = 2;
        max_refine_level = 30;
        truncate_mode = 0;
        refine = true;
        autorefine = true;
        debug = false;
        truncate_on_project = false;
        apply_randomize = false;
        project_randomize = false;
        bc = BoundaryConditions<NDIM>(BC_FREE);
        cell = Tensor<double>(NDIM,2);
        cell(_,1) = 1.0;
        recompute_cell_info();

        //pmap = std::shared_ptr< WorldDCPmapInterface< Key<NDIM> > >(new WorldDCDefaultPmap< Key<NDIM> >(world));
        pmap = std::shared_ptr< WorldDCPmapInterface< Key<NDIM> > >(new MyPmap<NDIM>(world));
        //pmap = std::shared_ptr< WorldDCPmapInterface< Key<NDIM> > >(new SimpleMap< Key<NDIM> >(world));
    }

    template <typename T, std::size_t NDIM>
    const FunctionCommonData<T,NDIM>* FunctionCommonData<T,NDIM>::data[MAXK] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    template <std::size_t NDIM> int FunctionDefaults<NDIM>::k;
    template <std::size_t NDIM> double FunctionDefaults<NDIM>::thresh;
    template <std::size_t NDIM> int FunctionDefaults<NDIM>::initial_level;
    template <std::size_t NDIM> int FunctionDefaults<NDIM>::max_refine_level;
    template <std::size_t NDIM> int FunctionDefaults<NDIM>::truncate_mode;
    template <std::size_t NDIM> bool FunctionDefaults<NDIM>::refine;
    template <std::size_t NDIM> bool FunctionDefaults<NDIM>::autorefine;
    template <std::size_t NDIM> bool FunctionDefaults<NDIM>::debug;
    template <std::size_t NDIM> bool FunctionDefaults<NDIM>::truncate_on_project;
    template <std::size_t NDIM> bool FunctionDefaults<NDIM>::apply_randomize;
    template <std::size_t NDIM> bool FunctionDefaults<NDIM>::project_randomize;
    template <std::size_t NDIM> BoundaryConditions<NDIM> FunctionDefaults<NDIM>::bc;
    template <std::size_t NDIM> Tensor<double> FunctionDefaults<NDIM>::cell;
    template <std::size_t NDIM> Tensor<double> FunctionDefaults<NDIM>::cell_width;
    template <std::size_t NDIM> Tensor<double> FunctionDefaults<NDIM>::rcell_width;
    template <std::size_t NDIM> double FunctionDefaults<NDIM>::cell_volume;
    template <std::size_t NDIM> double FunctionDefaults<NDIM>::cell_min_width;
    template <std::size_t NDIM> std::shared_ptr< WorldDCPmapInterface< Key<NDIM> > > FunctionDefaults<NDIM>::pmap;

    template <std::size_t NDIM> std::vector< Key<NDIM> > Displacements<NDIM>::disp;
    template <std::size_t NDIM> std::vector< Key<NDIM> > Displacements<NDIM>::disp_periodicsum[64];

}

#endif // MADNESS_MRA_MRAIMPL_H__INCLUDED
