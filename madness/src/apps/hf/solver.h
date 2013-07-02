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

  $Id: solver.h 2350 2011-06-06 17:31:35Z wsttiger@gmail.com $
*/
#include <mra/mra.h>
//#include <mra/lbdeux.h>
#include <world/world.h>
#include <linalg/solvers.h>
#include <vector>
#include <fortran_ctypes.h>
#include <cmath>

//#include "poperator.h"
#include "libxc.h"
#include "electronicstructureparams.h"
#include "complexfun.h"
#include "esolver.h"

#ifndef SOLVER_H_

using std::endl;
using std::setfill;
using std::setw;
using std::ostream;

/*!
  \ingroup applications
  \defgroup periodic_solver Periodic Solver
  \brief The Periodic Solver group is a group that contains the software
  objects that are needed to solve a periodic Kohn-Sham hamiltonian.
*/

//*****************************************************************************
/* static double onesfunc(const coordT& x) */
/* { */
/*   return 1.0; */
/* } */
//*****************************************************************************

namespace madness
{
  //***************************************************************************
  /*!
   \ingroup periodic_solver
   \Fermi-Dirac distribution function for fixing the occupation numbers
   */
  template <typename T>
  T stheta_fd(const T& x)
  {
    if (x > 50.0)
    {
      return 1.0;
    }
    else if (x < -50.0)
    {
      return 0.0;
    }
    else
    {
      return 1.0/(1.0 + exp(-x));
    }
  }
  //***************************************************************************

  // NOTE: so this is totally hacked because of the requirement of using a box
  //       of L*L*L rather than the convention of using a combination of
  //       lattice vectors and reciprocal lattice vectors (FOR NOW)
  // On second thought, it might be ok for general use if used "right".
  //***************************************************************************
  template <std::size_t NDIM>
  class ComplexExp : public FunctionFunctorInterface<double_complex,NDIM> {
  public:
      typedef Vector<double,NDIM> coordT;
      typedef Vector<double,NDIM> vec3dT;
      const double_complex coeff;
      const vec3dT exponent;

      ComplexExp(vec3dT exponent, double_complex coeff)
              : exponent(exponent), coeff(coeff) {}

      double_complex operator()(const coordT& x) const {
          double sum = 0.0;
          for (std::size_t i=0; i<NDIM; ++i) {
              sum += x[i]*exponent[i];
          };
          return coeff*exp(double_complex(0.0,sum));
      }
  };
  //***************************************************************************


  /*!
   \ingroup periodic_solver

   \brief The SubspaceK class is a container class holding previous orbitals
   and residuals.
   \par
   The Solver class uses the Krylov Accelerated Inexact Newton Solver (KAIN)
   accelerate the convergence a given calculation. The KAIN solver needs to
   store a subspace of previous orbitals and residuals. In the case is this
   implementation, the orbitals are store according to which k-point to which
   they belong.
  */

  //***************************************************************************
  template <typename T, int NDIM>
  class SubspaceK
  {
    // Typedef's
    typedef std::complex<T> valueT;
    typedef Function<valueT,NDIM> functionT;
    typedef std::vector<functionT> vecfuncT;
    typedef std::pair<vecfuncT,vecfuncT> pairvecfuncT;
    typedef std::vector<pairvecfuncT> subspaceT;
    typedef Tensor<valueT> tensorT;
    typedef std::vector<tensorT> vectensorT;
    typedef std::vector<subspaceT> vecsubspaceT;

    //*************************************************************************
    World& _world;
    //*************************************************************************

    //*************************************************************************
    vectensorT _Q;
    //*************************************************************************

    //*************************************************************************
    vecsubspaceT _subspace;
    //*************************************************************************

    //*************************************************************************
    std::vector<KPoint> _kpoints;
    //*************************************************************************

    //*************************************************************************
    ElectronicStructureParams _params;
    //*************************************************************************

    //*************************************************************************
    double _residual;
    //*************************************************************************

  public:

    //*************************************************************************
    SubspaceK(World& world, const ElectronicStructureParams& params,
            const std::vector<KPoint>& kpoints) : _world(world), _kpoints(kpoints),
            _params(params)
    {
      _residual = 1e6;
      for (unsigned int kp = 0; kp < _kpoints.size(); kp++)
      {
        _Q.push_back(tensorT(1,1));
        _subspace.push_back(subspaceT());
      }
    }
    //*************************************************************************

    //*************************************************************************
    void update_subspace(vecfuncT& awfs_new,
                         vecfuncT& bwfs_new,
                         const vecfuncT& awfs_old,
                         const vecfuncT& bwfs_old)
    {
      // compute residuals (total)
      vecfuncT t_rm = sub(_world, awfs_old, awfs_new);
      if (_params.spinpol)
      {
        vecfuncT br = sub(_world, bwfs_old, bwfs_new);
        t_rm.insert(t_rm.end(), br.begin(), br.end());
      }
      std::vector<double> rnvec = norm2<valueT,NDIM>(_world, t_rm);
      if (_world.rank() == 0)
      {
        double rnorm = 0.0;
        for (unsigned int i = 0; i < rnvec.size(); i++) rnorm += rnvec[i];
        if (_world.rank() == 0) print("residual = ", rnorm);
        _residual = rnorm;
      }
      _world.gop.broadcast(_residual, 0);

      for (unsigned int kp = 0; kp < _kpoints.size(); kp++)
      {
        KPoint kpoint = _kpoints[kp];

        vecfuncT k_phisa(awfs_old.begin() + kpoint.begin, awfs_old.begin() + kpoint.end);
        vecfuncT k_phisb(bwfs_old.begin() + kpoint.begin, bwfs_old.begin() + kpoint.end);
        vecfuncT k_awfs(awfs_new.begin() + kpoint.begin, awfs_new.begin() + kpoint.end);
        vecfuncT k_bwfs(bwfs_new.begin() + kpoint.begin, bwfs_new.begin() + kpoint.end);

        // compute residuals for k-point
        // concatentate up and down spins
        vecfuncT k_rm = sub(_world, k_phisa, k_awfs);
        vecfuncT k_vm(k_phisa);
        if (_params.spinpol)
        {
          vecfuncT k_br = sub(_world, k_phisb, k_bwfs);
          k_rm.insert(k_rm.end(), k_br.begin(), k_br.end());
          k_vm.insert(k_vm.end(), k_phisb.begin(), k_phisb.end());
        }

        // Update subspace and matrix Q
        compress(_world, k_vm, false);
        compress(_world, k_rm, false);
        _world.gop.fence();
        subspaceT k_subspace = _subspace[kp];
        k_subspace.push_back(pairvecfuncT(k_vm,k_rm));

        int m = k_subspace.size();
        tensorT ms(m);
        tensorT sm(m);
        for (int s = 0; s < m; s++)
        {
            const vecfuncT& k_vs = k_subspace[s].first;
            const vecfuncT& k_rs = k_subspace[s].second;
            for (unsigned int i = 0; i < k_vm.size(); i++)
            {
                ms[s] += k_vm[i].inner_local(k_rs[i]);
                sm[s] += k_vs[i].inner_local(k_rm[i]);
            }
        }
        _world.gop.sum(ms.ptr(),m);
        _world.gop.sum(sm.ptr(),m);

        tensorT newQ(m,m);
        if (m > 1) newQ(Slice(0,-2),Slice(0,-2)) = _Q[kp];
        newQ(m-1,_) = ms;
        newQ(_,m-1) = sm;

        _Q[kp] = newQ;
        if (_world.rank() == 0) print(_Q[kp]);

        // Solve the subspace equations
        tensorT c;
        if (_world.rank() == 0) {
            double rcond = 1e-12;
            while (1) {
                c = KAIN(_Q[kp],rcond);
                if (abs(c[m-1]) < 3.0) {
                    break;
                }
                else if (rcond < 0.01) {
                    if (_world.rank() == 0)
                      print("Increasing subspace singular value threshold ", c[m-1], rcond);
                    rcond *= 100;
                }
                else {
                    if (_world.rank() == 0)
                      print("Forcing full step due to subspace malfunction");
                    c = 0.0;
                    c[m-1] = 1.0;
                    break;
                }
            }
        }

        _world.gop.broadcast_serializable(c, 0);
        if (_world.rank() == 0) {
            //print("Subspace matrix");
            //print(Q);
            print("Subspace solution", c);
        }

        // Form linear combination for new solution
        vecfuncT k_phisa_new = zero_functions<valueT,NDIM>(_world, k_phisa.size());
        vecfuncT k_phisb_new = zero_functions<valueT,NDIM>(_world, k_phisb.size());
        compress(_world, k_phisa_new, false);
        compress(_world, k_phisb_new, false);
        _world.gop.fence();
        std::complex<double> one = std::complex<double>(1.0,0.0);
        unsigned int norbitals = awfs_old.size() / _kpoints.size();
        for (unsigned int m = 0; m < k_subspace.size(); m++)
        {
            const vecfuncT& k_vm = k_subspace[m].first;
            const vecfuncT& k_rm = k_subspace[m].second;
            // WSTHORNTON Stopped here!
            const vecfuncT  vma(k_vm.begin(),k_vm.begin() + norbitals);
            const vecfuncT  rma(k_rm.begin(),k_rm.begin() + norbitals);
            const vecfuncT  vmb(k_vm.end() - norbitals, k_vm.end());
            const vecfuncT  rmb(k_rm.end() - norbitals, k_rm.end());

            gaxpy(_world, one, k_phisa_new, c(m), vma, false);
            gaxpy(_world, one, k_phisa_new,-c(m), rma, false);
            gaxpy(_world, one, k_phisb_new, c(m), vmb, false);
            gaxpy(_world, one, k_phisb_new,-c(m), rmb, false);
        }
        _world.gop.fence();

        if (_params.maxsub <= 1) {
            // Clear subspace if it is not being used
            k_subspace.clear();
        }
        else if (k_subspace.size() == _params.maxsub) {
            // Truncate subspace in preparation for next iteration
            k_subspace.erase(k_subspace.begin());
            _Q[kp] = _Q[kp](Slice(1,-1),Slice(1,-1));
        }
        // Save subspace
        _subspace[kp] = k_subspace;

        for (unsigned int wi = kpoint.begin, fi = 0; wi < kpoint.end;
          wi++, fi++)
        {
          awfs_new[wi] = k_phisa_new[fi];
          bwfs_new[wi] = k_phisb_new[fi];
        }
      }
    }
    //*************************************************************************

  };

  /*!
   \ingroup periodic_solver

   \brief The SubspaceK class is a container class holding previous orbitals
   and residuals.
   \par
   The Solver class uses the Krylov Accelerated Inexact Newton Solver (KAIN)
   accelerate the convergence a given calculation. The KAIN solver needs to
   store a subspace of previous orbitals and residuals.
  */

  //***************************************************************************
  template <typename T, int NDIM>
  class Subspace
  {
    // Typedef's
    typedef std::complex<T> valueT;
    typedef Function<valueT,NDIM> functionT;
    typedef std::vector<functionT> vecfuncT;
    typedef std::pair<vecfuncT,vecfuncT> pairvecfuncT;
    typedef std::vector<pairvecfuncT> subspaceT;
    typedef Tensor<valueT> tensorT;

    //*************************************************************************
    World& _world;
    //*************************************************************************

    //*************************************************************************
    tensorT _Q;
    //*************************************************************************

    //*************************************************************************
    subspaceT _subspace;
    //*************************************************************************

    //*************************************************************************
    std::vector<KPoint> _kpoints;
    //*************************************************************************

    //*************************************************************************
    ElectronicStructureParams _params;
    //*************************************************************************

  public:

    //*************************************************************************
    Subspace(World& world, const ElectronicStructureParams& params)
      : _world(world), _params(params)
    {
    }
    //*************************************************************************

    //*************************************************************************
    void update_subspace(vecfuncT& awfs_new,
                         vecfuncT& bwfs_new,
                         const vecfuncT& awfs_old,
                         const vecfuncT& bwfs_old,
                         const vecfuncT& rm)
    {
      // concatentate up and down spins
      vecfuncT vm = awfs_old;
      if (_params.spinpol)
      {
        vm.insert(vm.end(), bwfs_old.begin(), bwfs_old.end());
      }

      // Update subspace and matrix Q
      compress(_world, vm, false);
      compress(_world, rm, false);
      _world.gop.fence();
      _subspace.push_back(pairvecfuncT(vm,rm));

      int m = _subspace.size();
      tensorT ms(m);
      tensorT sm(m);
      for (int s=0; s<m; s++)
      {
          const vecfuncT& vs = _subspace[s].first;
          const vecfuncT& rs = _subspace[s].second;
          for (unsigned int i=0; i<vm.size(); i++)
          {
              ms[s] += vm[i].inner_local(rs[i]);
              sm[s] += vs[i].inner_local(rm[i]);
          }
      }
      _world.gop.sum(ms.ptr(),m);
      _world.gop.sum(sm.ptr(),m);

      tensorT newQ(m,m);
      if (m > 1) newQ(Slice(0,-2),Slice(0,-2)) = _Q;
      newQ(m-1,_) = ms;
      newQ(_,m-1) = sm;

      _Q = newQ;
      if (_world.rank() == 0) print(_Q);

      // Solve the subspace equations
      tensorT c;
      if (_world.rank() == 0) {
          double rcond = 1e-12;
          while (1) {
              c = KAIN(_Q,rcond);
              if (abs(c[m-1]) < 3.0) {
                  break;
              }
              else if (rcond < 0.01) {
                  if (_world.rank() == 0)
                    print("Increasing subspace singular value threshold ", c[m-1], rcond);
                  rcond *= 100;
              }
              else {
                  if (_world.rank() == 0)
                    print("Forcing full step due to subspace malfunction");
                  c = 0.0;
                  c[m-1] = 1.0;
                  break;
              }
          }
      }

      _world.gop.broadcast_serializable(c, 0);
      if (_world.rank() == 0) {
          //print("Subspace matrix");
          //print(Q);
          print("Subspace solution", c);
      }

      // Form linear combination for new solution
      vecfuncT phisa_new = zero_functions<valueT,NDIM>(_world, awfs_old.size());
      vecfuncT phisb_new = zero_functions<valueT,NDIM>(_world, bwfs_old.size());
      compress(_world, phisa_new, false);
      compress(_world, phisb_new, false);
      _world.gop.fence();
      std::complex<double> one = std::complex<double>(1.0,0.0);
      for (unsigned int m=0; m<_subspace.size(); m++) {
          const vecfuncT& vm = _subspace[m].first;
          const vecfuncT& rm = _subspace[m].second;
          const vecfuncT  vma(vm.begin(),vm.begin()+awfs_old.size());
          const vecfuncT  rma(rm.begin(),rm.begin()+awfs_old.size());
          const vecfuncT  vmb(vm.end()-bwfs_old.size(), vm.end());
          const vecfuncT  rmb(rm.end()-bwfs_old.size(), rm.end());

          gaxpy(_world, one, phisa_new, c(m), vma, false);
          gaxpy(_world, one, phisa_new,-c(m), rma, false);
          gaxpy(_world, one, phisb_new, c(m), vmb, false);
          gaxpy(_world, one, phisb_new,-c(m), rmb, false);
      }
      _world.gop.fence();

      if (_params.maxsub <= 1) {
          // Clear subspace if it is not being used
          _subspace.clear();
      }
      else if (_subspace.size() == _params.maxsub) {
          // Truncate subspace in preparation for next iteration
          _subspace.erase(_subspace.begin());
          _Q = _Q(Slice(1,-1),Slice(1,-1));
      }
      awfs_new = phisa_new;
      bwfs_new = phisb_new;
    }
    //*************************************************************************

    //*************************************************************************
    void reproject()
    {
      //  //if (_world.rank() == 0)
      //    //printf("\n\nreprojecting subspace to wavelet order: %d and thresh: %.5e\n\n",
      //    //FunctionDefaults<3>::get_k(), FunctionDefaults<3>::get_thresh());
      //
      //  unsigned int m = _subspace.size();
      //  for (unsigned int s = 0; s < m; s++)
      //  {
      //      vecfuncT& vs = _subspace[s].first;
      //      vecfuncT& rs = _subspace[s].second;
      //      reconstruct(_world, vs);
      //      reconstruct(_world, rs);
      //      unsigned int vm = vs.size();
      //      for (unsigned int i = 0; i < vm; i++)
      //      {
      //        vs[i] = madness::project(vs[i], FunctionDefaults<3>::get_k(),
      //          FunctionDefaults<3>::get_thresh(), false);
      //        rs[i] = madness::project(rs[i], FunctionDefaults<3>::get_k(),
      //          FunctionDefaults<3>::get_thresh(), false);
      //      }
      //      _world.gop.fence();
      //      truncate(_world, vs);
      //      truncate(_world, rs);
      //      normalize(_world, vs);
      //  }
      //  _world.gop.fence();

    }
    //*************************************************************************

  };

//  template <typename T, int NDIM>
//  struct lbcost {
//      double leaf_value;
//      double parent_value;
//      lbcost(double leaf_value=1.0, double parent_value=0.0) : leaf_value(leaf_value), parent_value(parent_value) {}
//      double operator()(const Key<NDIM>& key, const FunctionNode<T,NDIM>& node) const {
//          if (key.level() <= 1) {
//              return 100.0*(leaf_value+parent_value);
//          }
//          else if (node.is_leaf()) {
//              return leaf_value;
//          }
//          else {
//              return parent_value;
//          }
//      }
//  };

  //***************************************************************************

  /*! \ingroup periodic_solver
      \brief The main class of the periodic DFT solver
      \f[
      z = frac{x}{1 - y^2}
      \f]
  */

  //***************************************************************************
  template <typename T, int NDIM>
  class Solver
  {
    // Typedef's
    typedef std::complex<T> valueT;
    typedef Function<T,NDIM> rfuntionT;
    typedef FunctionFactory<T,NDIM> rfactoryT;
    typedef Function<valueT,NDIM> functionT;
    typedef std::vector<functionT> vecfuncT;
    typedef FunctionFactory<valueT,NDIM> factoryT;
    typedef Vector<double,NDIM> kvecT;
    typedef SeparatedConvolution<T,3> operatorT;
    typedef std::shared_ptr<operatorT> poperatorT;
    typedef Tensor<double> rtensorT;
    typedef Tensor<std::complex<double> > ctensorT;
    typedef Tensor<valueT> tensorT;
    typedef std::pair<vecfuncT,vecfuncT> pairvecfuncT;
    typedef std::vector<pairvecfuncT> subspaceT;
    typedef std::vector<tensorT> vectensorT;
    typedef std::vector<subspaceT> vecsubspaceT;

    //*************************************************************************
    World& _world;
    //*************************************************************************

    //*************************************************************************
    // This variable could either be a nuclear potiential or a nuclear charge
    // density depending on the "ispotential" boolean variable in the
    // ElectronicStructureParams class.
    rfuntionT _vnucrhon;
    //*************************************************************************

    //*************************************************************************
    vecfuncT _phisa;
    //*************************************************************************

    //*************************************************************************
    vecfuncT _phisb;
    //*************************************************************************

    //*************************************************************************
    std::vector<T> _eigsa;
    //*************************************************************************

    //*************************************************************************
    std::vector<T> _eigsb;
    //*************************************************************************

    //*************************************************************************
    ElectronicStructureParams _params;
    //*************************************************************************

    //*************************************************************************
    std::vector<KPoint> _kpoints;
    //*************************************************************************

    //*************************************************************************
    std::vector<double> _occsa;
    //*************************************************************************

    //*************************************************************************
    std::vector<double> _occsb;
    //*************************************************************************

    //*************************************************************************
    rfuntionT _rhoa;
    //*************************************************************************

    //*************************************************************************
    rfuntionT _rhob;
    //*************************************************************************

    //*************************************************************************
    rfuntionT _rho;
    //*************************************************************************

    //*************************************************************************
    rfuntionT _vnuc;
    //*************************************************************************

    //*************************************************************************
    SeparatedConvolution<T,NDIM>* _cop;
    //*************************************************************************

//    //*************************************************************************
//    vecsubspaceT _subspace;
//    //*************************************************************************
//
//    //*************************************************************************
//    vectensorT _Q;
//    //*************************************************************************

    //*************************************************************************
    Subspace<T,NDIM>* _subspace;
    //*************************************************************************

    //*************************************************************************
    MolecularEntity _mentity;
    //*************************************************************************

    //*************************************************************************
    double _residual;
    //*************************************************************************

    //*************************************************************************
    AtomicBasisSet _aobasis;
    //*************************************************************************

    //*************************************************************************
    double _maxthresh;
    //*************************************************************************

    //*************************************************************************
    std::ofstream _outputF;
    //*************************************************************************

    //*************************************************************************
    std::ofstream _matF;
    //*************************************************************************

    //*************************************************************************
    std::ofstream _eigF;
    //*************************************************************************

    //*************************************************************************
    std::ofstream _kmeshF;
    //*************************************************************************

    //*************************************************************************
    int _it;
    //*************************************************************************

  public:

    //*************************************************************************
    double ttt, sss;
    void START_TIMER(World& world) {
        world.gop.fence(); ttt=wall_time(); sss=cpu_time();
    }

    void END_TIMER(World& world, const char* msg) {
    ttt=wall_time()-ttt; sss=cpu_time()-sss; if (world.rank()==0) printf("timer: %20.20s %8.2fs %8.2fs\n", msg, sss, ttt);
    }
    //*************************************************************************

    //*************************************************************************
    Solver(World& world, const std::string& filename) : _world(world)
    {
      init(filename);
      if (_params.periodic) FunctionDefaults<3>::set_bc(BC_PERIODIC);
      _residual = 1e5;
      make_nuclear_potential();

      if (_params.restart==0) initial_guess();
//      for (unsigned int kp = 0; kp < _kpoints.size(); kp++)
//      {
//        _Q.push_back(tensorT(1,1));
//        _subspace.push_back(subspaceT());
//      }
      _subspace = new Subspace<T,NDIM>(world, _params);
    }
    //*************************************************************************

    //*************************************************************************
    void init(const std::string& filename)
    {
      // params
      if (_world.rank() == 0)
      {
        _params.read_file(filename);
        //_params.fractional = false;
      }
      // Send params
      _world.gop.broadcast_serializable(_params, 0);
      if (_params.fractional)
        //FunctionDefaults<3>::set_cubic_cell(0,_params.L);
        FunctionDefaults<3>::set_cubic_cell(-_params.L/2,_params.L/2);
      else
        FunctionDefaults<3>::set_cubic_cell(-_params.L/2,_params.L/2);
      FunctionDefaults<3>::set_thresh(_params.thresh);
      FunctionDefaults<3>::set_k(_params.waveorder);

      // mentity and aobasis
      if (_world.rank() == 0)
      {
        _aobasis.read_file(_params.basis);
        _mentity.read_file(filename, _params.fractional);
        _mentity.center();
      }
      // Send mentity and aobasis
      _world.gop.broadcast_serializable(_mentity, 0);
      _world.gop.broadcast_serializable(_aobasis, 0);
      // set number of electrons to the total nuclear charge of the mentity
      _params.nelec = _mentity.total_nuclear_charge();
      // total number of bands include empty
      _params.nbands = (_params.nelec/2) + _params.nempty;
      _params.ncharge = _mentity.total_nuclear_charge();
      if ((_params.nelec % 2) == 1) _params.nelec++;

      // Open output files
      _outputF.open("output.txt");
      _matF.open("matrix.txt");
      _eigF.open("eigs.txt");

      // kmesh
      if (_params.restart == 0)
      {
        if (_params.periodic) // PERIODIC
        {
          // GAMMA POINT
          if ((_params.ngridk0 == 1) && (_params.ngridk1 == 1) && (_params.ngridk2 == 1))
          {
            _kpoints.push_back(KPoint(coordT(0.0), 1.0));
          }
//          if ((_params.ngridk0 == 0) && (_params.ngridk1 == 0) && (_params.ngridk2 == 0))
//          {
//            _kpoints.push_back(KPoint(coordT(0.0), 0.5));
//            _kpoints.push_back(KPoint(coordT(0.0), 0.5));
//          }
          else // NORMAL BANDSTRUCTURE
          {
            _kpoints = genkmesh(_params.ngridk0, _params.ngridk1,
                                _params.ngridk2, _params.koffset0,
                                _params.koffset1, _params.koffset2,
                                _params.L);
          }
        }
        else // NOT-PERIODIC
        {
          _kpoints.push_back(KPoint(coordT(0.0), 1.0));
        }
        if (_world.rank() == 0)
        {
          _kmeshF.open("kpoints.txt");
          _kmeshF << "kpts: " << _kpoints.size() << endl;
          _kmeshF << "ik" << setw(10) << "kpt" << setw(30) << "weight" << endl;
          _kmeshF << "--" << setw(10) << "---" << setw(30) << "------" << endl;

          //_kmeshF << setfill('-') << setw(55) << "-" << endl;
          //_kmeshF << setfill(' ') << endl;
          _kmeshF << endl;
          for (unsigned int i = 0; i < _kpoints.size(); i++)
          {
            KPoint kpoint = _kpoints[i];
            _kmeshF << i << setw(10) << kpoint.k[0];
            _kmeshF << setw(10) << kpoint.k[1];
            _kmeshF << setw(10) << kpoint.k[2];
            _kmeshF << setw(10) << kpoint.weight << endl;
          }
          _kmeshF.close();
        }
      }
      else
      {
        load_orbitals();
      }
    }
    //*************************************************************************

    //*************************************************************************
    std::vector<KPoint> genkmesh(unsigned int ngridk0, unsigned ngridk1, unsigned int ngridk2,
                                 double koffset0, double koffset1, double koffset2, double R)
    {
      std::vector<KPoint> kmesh;
      double step0 = 1.0/ngridk0;
      double step1 = 1.0/ngridk1;
      double step2 = 1.0/ngridk2;
      double weight = 1.0/(ngridk0*ngridk1*ngridk2);
      double TWO_PI = 2.0 * madness::constants::pi;
      for (unsigned int i = 0; i < ngridk0; i++)
      {
        for (unsigned int j = 0; j < ngridk1; j++)
        {
          for (unsigned int k = 0; k < ngridk2; k++)
          {
            //double k0 = (i*step0 - step0/2) * TWO_PI/R;
            //double k1 = (j*step1 - step1/2) * TWO_PI/R;
            //double k2 = (k*step2 - step2/2) * TWO_PI/R;
            double k0 = ((i*step0)+koffset0) * TWO_PI/R;
            double k1 = ((j*step1)+koffset1) * TWO_PI/R;
            double k2 = ((k*step2)+koffset2) * TWO_PI/R;
            KPoint kpoint(k0, k1, k2, weight);
            kmesh.push_back(kpoint);
          }
        }
      }
      return kmesh;
    }
    //*************************************************************************

    //*************************************************************************
    void save_orbitals()
    {
      archive::ParallelOutputArchive ar(_world, "orbitals", _params.nio);
      ar & _params.spinpol;
      ar & (unsigned int)(_kpoints.size());
      for (unsigned int i = 0; i < _kpoints.size(); i++) ar & _kpoints[i];
//      ar & (unsigned int)(_occsa.size());
//      for (unsigned int i = 0; i < _occsa.size(); i++) ar & _occsa[i];
      ar & (unsigned int)(_phisa.size());
      for (unsigned int i = 0; i < _phisa.size(); i++) ar & _phisa[i];
      if (_params.spinpol)
      {
//        ar & (unsigned int)(_occsb.size());
//        for (unsigned int i = 0; i < _occsb.size(); i++) ar & _occsb[i];
        ar & (unsigned int)(_phisb.size());
        for (unsigned int i = 0; i < _phisb.size(); i++) ar & _phisb[i];
      }
    }
    //*************************************************************************

    //*************************************************************************
    void load_orbitals()
    {
      const double thresh = FunctionDefaults<3>::get_thresh();
      const int k = FunctionDefaults<3>::get_k();

      archive::ParallelInputArchive ar(_world, "orbitals");

      // spin-polarized
      bool spinrest;
      ar & spinrest;
      // kpoints
      unsigned int nkpts;
      ar & nkpts;
      _kpoints.clear();
      for (unsigned int i = 0; i < nkpts; i++)
      {
        KPoint tkpt;
        ar & tkpt;
        _kpoints.push_back(tkpt);
      }
      // occs
//      unsigned int noccs;
//      ar & noccs;
//      _occs.clear();
//      for (unsigned int i = 0; i < noccs; i++)
//      {
//        double tocc;
//        ar & tocc;
//        _occs.push_back(tocc);
//      }
      // orbitals
      unsigned int norbs;
      ar & norbs;
      _phisa.clear();
      _eigsa.clear();
      for (unsigned int i = 0; i < norbs; i++)
      {
        functionT tfunc;
        ar & tfunc;
        _phisa.push_back(tfunc);
        _eigsa.push_back(-0.1);
      }
      // check for k mismatch
      if (_phisa[0].k() != k)
      {
        reconstruct(_world,_phisa);
        for (unsigned int i = 0; i < _phisa.size(); i++)
          _phisa[i] = madness::project(_phisa[i], k, thresh, false);
        _world.gop.fence();
      }
      // orthonormalize
      for (unsigned int i = 0; i < _kpoints.size(); i++)
        gram_schmidt(_phisa, _kpoints[i]);

      if (_params.spinpol)
      {
        _phisb.clear();
        _eigsb.clear();
        for (unsigned int i = 0; i < norbs; i++)
        {
          functionT tfunc;
          ar & tfunc;
          _phisb.push_back(tfunc);
          _eigsb.push_back(-0.1);
        }
        // check for k mismatch
        if (_phisb[0].k() != k)
        {
          reconstruct(_world,_phisb);
          for (unsigned int i = 0; i < _phisb.size(); i++)
            _phisb[i] = madness::project(_phisb[i], k, thresh, false);
          _world.gop.fence();
        }
        // orthonormalize
        for (unsigned int i = 0; i < _kpoints.size(); i++)
          gram_schmidt(_phisb, _kpoints[i]);
      }
      else
      {
        for (unsigned int i = 0; i < norbs; i++)
        {
          _phisb.push_back(_phisa[i]);
          _eigsb.push_back(_eigsa[i]);
        }
      }
      // create vector for occupation numbers
      _occsa = std::vector<double>(norbs, 0.0);
      _occsb = std::vector<double>(norbs, 0.0);
    }
    //*************************************************************************

    //*************************************************************************
    void make_nuclear_potential()
    {

      if (_params.periodic) // periodic
      {
        Tensor<double> cellsize = FunctionDefaults<3>::get_cell_width();
        //_cop = PeriodicCoulombOpPtr<double,3>(_world, _params.waveorder,_params.lo, _params.thresh * 0.1, cellsize);
        _cop = CoulombOperatorPtr(_world, _params.lo, FunctionDefaults<3>::get_thresh() * 0.1);
      }
      else // not periodic
      {
        _cop = CoulombOperatorPtr(_world, _params.lo, FunctionDefaults<3>::get_thresh() * 0.1);
        //_cop = CoulombOperatorPtr(_world, _params.lo, _params.thresh * 0.1);
      }

      if (_world.rank() == 0)
        _outputF << "Making nuclear potential .." << endl << endl;
      Tensor<double> csize = FunctionDefaults<3>::get_cell();
      if (_world.rank() == 0)
      {
        print("cell(x) is ",csize(0,0), csize(0,1));
        print("cell(y) is ",csize(1,0), csize(1,1));
        print("cell(z) is ",csize(2,0), csize(2,1));
      }
      if (_params.ispotential) // potential
      {
        _vnucrhon = rfactoryT(_world).functor(rfunctorT(new MolecularPotentialFunctor(_mentity))).thresh(_params.thresh * 0.1).truncate_on_project();
        _vnuc = copy(_vnucrhon);
        double vnuctrace = _vnuc.trace();
        //if (_world.rank == 0) print("vnuc trace:     ", vnuctrace);
        _vnuc.reconstruct();
        vnuctrace = _vnuc.trace();
        //if (_world.rank == 0) print("vnuc trace:     ", vnuctrace);
      }
      else // charge density
      {
        std::vector<coordT> specialpts;
        for (int i = 0; i < _mentity.natom(); i++)
        {
          coordT pt(0.0);
          Atom atom = _mentity.get_atom(i);
          pt[0] = atom.x; pt[1] = atom.y; pt[2] = atom.z;
          specialpts.push_back(pt);
          if (_world.rank() == 0) print("Special point: ", pt);
        }
        double now = wall_time();
        _vnucrhon = rfactoryT(_world).functor(
            rfunctorT(new MolecularNuclearChargeDensityFunctor(_mentity, _params.L, _params.periodic, specialpts))).
            thresh(_params.thresh).initial_level(6).truncate_on_project();

        if (_world.rank() == 0) printf("%f\n", wall_time() - now);
        if (_world.rank() == 0) print("calculating trace of rhon ..\n\n");
        double rtrace = _vnucrhon.trace();
        if (_world.rank() == 0) print("rhon trace = ", rtrace);
        now = wall_time();
        _vnucrhon.truncate();
        _vnuc = apply(*_cop, _vnucrhon);
        if (_world.rank() == 0) printf("%f\n", wall_time() - now);
        if (_world.rank() == 0) print("Done creating nuclear potential ..\n");
      }

      std::vector<long> npt(3,101);
      plotdx(_vnuc, "vnuc.dx", FunctionDefaults<3>::get_cell(), npt);
    }
    //*************************************************************************

    //*************************************************************************
    struct GuessDensity : public FunctionFunctorInterface<double,3> {
        const MolecularEntity& mentity;
        const AtomicBasisSet& aobasis;
        double R;
        const bool periodic;

        double operator()(const coordT& x) const
        {
          double value = 0.0;
          if (periodic)
          {
            for (int xr = -1; xr <= 1; xr += 1)
            {
              for (int yr = -1; yr <= 1; yr += 1)
              {
                for (int zr = -1; zr <= 1; zr += 1)
                {
                  value += aobasis.eval_guess_density(mentity,
                      x[0]+xr*R, x[1]+yr*R, x[2]+zr*R);
                }
              }
            }
          }
          else
          {
            value = aobasis.eval_guess_density(mentity, x[0], x[1], x[2]);
          }
          return value;
        }

        GuessDensity(const MolecularEntity& mentity, const AtomicBasisSet& aobasis,
            const double& R, const bool& periodic)
        : mentity(mentity), aobasis(aobasis), R(R), periodic(periodic) {}
    };
    //*************************************************************************

    //*************************************************************************
    rfunctionT
    make_lda_potential(World& world,
                       const rfunctionT& arho,
                       const rfunctionT& brho,
                       const rfunctionT& adelrhosq,
                       const rfunctionT& bdelrhosq)
    {
  //      MADNESS_ASSERT(!_params.spinpol);
        rfunctionT vlda = copy(arho);
//        vlda.unaryop(&::libxc_ldaop<double>);
        vlda.unaryop(&::ldaop<double>);
        return vlda;
    }

    vecfuncT project_ao_basis(World& world, KPoint kpt) {
        vecfuncT ao(_aobasis.nbf(_mentity));

        Level initial_level = 3;
        for (int i=0; i < _aobasis.nbf(_mentity); i++) {
            functorT aofunc(new AtomicBasisFunctor(_aobasis.get_atomic_basis_function(_mentity,i),
                _params.L, _params.periodic, kpt));
            ao[i] = factoryT(world).functor(aofunc).initial_level(initial_level).truncate_on_project().nofence();
            // WSTHORNTON
            // if at k-point == (0.0,0.0,0.0) basis function should be real
//            if ((abs(kpt.k[0]) < 1e-12) && (abs(kpt.k[1]) < 1e-12) &&
//                (abs(kpt.k[2]) < 1e-12))
//            {
//              if (!::madness::is_real<T,NDIM>(ao[i]))
//              {
//                double t1 = (imag(ao[i])).norm2();
//                if (_world.rank() == 0)
//                  printf("ERROR: ao basis function imaginary: %5.3d     %15.8f\n\n", i, t1);
//
//              }
//            }
        }
        world.gop.fence();

        std::vector<double> norms;

        norms = norm2s(world, ao);

        for (int i=0; i<_aobasis.nbf(_mentity); i++) {
            if (world.rank() == 0 && fabs(norms[i]-1.0)>1e-3) print(i," bad ao norm?", norms[i]);
            norms[i] = 1.0/norms[i];
        }

        scale(world, ao, norms);

        norms = norm2s(world, ao);

        for (int i=0; i<_aobasis.nbf(_mentity); i++) {
            if (world.rank() == 0 && fabs(norms[i]-1.0)>1e-3) print(i," bad ao norm?", norms[i]);
            norms[i] = 1.0/norms[i];
        }

        scale(world, ao, norms);

        return ao;
    }
    //*************************************************************************

    //*************************************************************************
    // Constructor
    Solver(World& world,
           rfuntionT vnucrhon,
           vecfuncT phisa,
           vecfuncT phisb,
           std::vector<T> eigsa,
           std::vector<T> eigsb,
           ElectronicStructureParams params,
           MolecularEntity mentity)
       : _world(world), _vnucrhon(vnucrhon), _phisa(phisa), _phisb(phisb),
       _eigsa(eigsa), _eigsb(eigsb), _params(params), _mentity(mentity)
    {
      _residual = 1e5;
      if (params.periodic)
      {
        Tensor<double> box = FunctionDefaults<NDIM>::get_cell_width();
        //_cop = PeriodicCoulombOpPtr<T,NDIM>(const_cast<World&>(world),
        //    FunctionDefaults<NDIM>::get_k(), params.lo, params.thresh * 0.1, box);
        _cop = CoulombOperatorPtr(const_cast<World&>(world), params.lo, params.thresh * 0.1);
      }
      else
      {
        _cop = CoulombOperatorPtr(const_cast<World&>(world), params.lo, params.thresh * 0.1);
      }

      if (params.ispotential)
      {
        _vnuc = copy(_vnucrhon);
      }
      else
      {
        _vnuc = apply(*_cop, _vnucrhon);
      }
    }
    //*************************************************************************

    /// Initializes alpha and beta mos, occupation numbers, eigenvalues
    //*************************************************************************
    void initial_guess()
    {
//      // Get initial guess for the electronic density
//      if (_world.rank() == 0) print("Guessing rho ...\n\n");
//      rfunctionT rho = (_params.restart == 0) ? rfactoryT(_world).functor(rfunctorT(
//          new GuessDensity(_mentity, _aobasis, _params.L, _params.periodic))).initial_level(3) :
//            compute_rho(_phisa, _kpoints);
//      if (_params.restart == 1)
//      {
//        if (_params.spinpol)
//          rho += compute_rho(_phisb, _kpoints);
//        else
//          rho.scale(2.0);
//      }

      // Get initial guess for the electronic density
      if (_world.rank() == 0) print("Guessing rho ...\n\n");

      rfunctionT rho = rfactoryT(_world);
      if (_params.restart == 0)
      {
        rho = rfactoryT(_world).functor(rfunctorT(new GuessDensity(_mentity,
              _aobasis, _params.L, _params.periodic))).initial_level(3);
      }
      else
      {
        MADNESS_EXCEPTION("restart not working right now!",0);
      }

      // This is a cheat
      double rtrace = rho.trace();
      if (_world.rank() == 0) print("trace of rho = ", rtrace);
      rho.scale(_params.ncharge/rho.trace());

      // load balance
      //if(_world.size() > 1) {
      //    START_TIMER(_world);
      //    LoadBalanceDeux<3> lb(_world);
      //    lb.add_tree(_vnuc, lbcost<double,3>(1.0, 0.0), false);
      //    lb.add_tree(rho, lbcost<double,3>(1.0, 1.0), true);

      //    FunctionDefaults<3>::redistribute(_world, lb.load_balance(6.0));
      //}

      if (_params.restart != 1)
      {
        // build effective potential
        rfunctionT vlocal;
        // Is this a many-body system?
        //int rank = _world.rank();
        if (_params.ncharge > 1.0)
        {
          if (_world.rank() == 0) print("Creating Coulomb op ...\n\n");
          SeparatedConvolution<double, 3>* op = 0;
          // Is this system periodic?
          if (_params.periodic)
          {
            Tensor<double> cellsize = FunctionDefaults<3>::get_cell_width();
            //op = PeriodicCoulombOpPtr<double, 3> (_world, _params.waveorder,
            //    _params.lo, _params.thresh * 0.1, cellsize);
            op = CoulombOperatorPtr(_world, _params.lo, _params.thresh * 0.1);
          }
          else
          {
            op = CoulombOperatorPtr(_world, _params.lo, _params.thresh * 0.1);
          }
          if (_world.rank() == 0) print("Building effective potential ...\n\n");
          rfunctionT vc = apply(*op, rho);
          vlocal = _vnuc + vc; //.scale(1.0-1.0/nel); // Reduce coulomb to increase binding
          rho.scale(0.5);
          // Do the LDA
          rfunctionT vlda = make_lda_potential(_world, rho, rho, rfunctionT(), rfunctionT());
          vlocal = vlocal + vlda;
          double vctrace = vc.trace();
          double vldatrace = vlda.trace();
          double vlocaltrace = vlocal.trace();
          if (_world.rank() == 0) print("vctrace:     ", vctrace);
          if (_world.rank() == 0) print("vldatrace:     ", vldatrace);
          if (_world.rank() == 0) print("vlocaltrace:     ", vlocaltrace);
          delete op;
    //      vector<long> npt(3,101);
    //      plotdx(vc, "vc.dx", FunctionDefaults<3>::get_cell(), npt);
        }
        else
        {
          vlocal = _vnuc;
        }

        // Clear these functions
        rho.clear();
        vlocal.reconstruct();

        // Get size information from k-points and ao_basis so that we can correctly size
        // the _orbitals data structure and the eigs tensor
        // number of orbitals in the basis set
        int nao = _aobasis.nbf(_mentity);
        // number of kpoints
        int nkpts = _kpoints.size();
        // total number of orbitals to be processed (no symmetry)
        int norbs = _params.nbands * nkpts;
        // Check to see if the basis set can accomodate the number of bands
        if (_params.nbands > nao)
          MADNESS_EXCEPTION("Error: basis not large enough to accomodate number of bands", nao);
        // set the number of orbitals
        _eigsa = std::vector<double>(norbs, 0.0);
        _eigsb = std::vector<double>(norbs, 0.0);
        _occsa = std::vector<double>(norbs, 0.0);
        _occsb = std::vector<double>(norbs, 0.0);
        int kp = 0;
        if (_world.rank() == 0) print("Building kinetic energy matrix ...\n\n");
        // Need to do kinetic piece for every k-point
        for (int ki = 0; ki < nkpts; ki++, kp += _params.nbands)
        {
          // These are our initial basis functions
          if (_world.rank() == 0) print("Projecting atomic orbitals ...\n\n");
          START_TIMER(_world);
          cvecfuncT ao = project_ao_basis(_world, _kpoints[ki]);
          END_TIMER(_world, "projecting atomic orbital basis");

      //    for (unsigned int ai = 0; ai < ao.size(); ai++)
      //    {
      //      std::ostringstream strm;
      //      strm << "aod" << ai << ".dx";
      //      std::string fname = strm.str();
      //      vector<long> npt(3,101);
      //      plotdx(ao[ai], fname.c_str(), FunctionDefaults<3>::get_cell(), npt);
      //    }

          // load balancing
          //if(_world.size() > 1)
          //{
          //  LoadBalanceDeux<3> lb(_world);
          //  lb.add_tree(_vnuc, lbcost<double,3>(1.0, 1.0), false);
          //  for(unsigned int i = 0; i < ao.size(); i++)
          //  {
          //      lb.add_tree(ao[i], lbcost<valueT,3>(1.0, 1.0), false);
          //  }

          //  FunctionDefaults<3>::redistribute(_world, lb.load_balance(6.0));
          //}

          // Get k-point from list
          KPoint& kpt = _kpoints[ki];
          // Build kinetic matrx
          ctensorT kinetic = ::kinetic_energy_matrix<T,NDIM>(_world, ao, _params.periodic, kpt);
          // Build the overlap matrix
          if (_world.rank() == 0) print("Building overlap matrix ...\n\n");
          ctensorT overlap = matrix_inner(_world, ao, ao, true);
          // Build the potential matrix
          reconstruct(_world, ao);
          if (_world.rank() == 0) print("Building potential energy matrix ...\n\n");
          cvecfuncT vpsi = mul_sparse(_world, vlocal, ao, _params.thresh);
          // I don't know why fence is called twice here
          _world.gop.fence();
          _world.gop.fence();
          compress(_world, vpsi);
          truncate(_world, vpsi);
          compress(_world, ao);
          // Build the potential matrix
          ctensorT potential = matrix_inner(_world, vpsi, ao, true);
          _world.gop.fence();
          // free memory
          vpsi.clear();
          _world.gop.fence();

          // Set occupation numbers
//          if (_params.spinpol)
//          {
//            MADNESS_EXCEPTION("spin polarized not implemented", 0);
//          }
//          else
//          {
//            int filledbands = _params.nelec / 2;
//            int occstart = kp;
//            int occend = kp + filledbands;
////            for (int i = occstart; i < occend; i++) _occs[i] = _params.maxocc;
////            if (_params.nelec > _params.ncharge)
////              _occs[occend] = 1.0;
//            double availcharge = _params.ncharge;
//            //double tol = 1e-6;
//            for (int i = occstart; i < occend; i++)
//            {
//              if (_world.rank() == 0) printf("availcharge = %.8f\n", availcharge);
//              // do we have charge to give? can we give a full dose?
//              if ((int) round(availcharge - _params.maxocc) >= 0)
//              {
//                if (_world.rank() == 0) printf("%d\tcase #1: availcharge - maxocc = %.8f\n", i, availcharge - _params.maxocc);
//                _occs[i] = _params.maxocc;
//                availcharge -= _params.maxocc;
//              }
//              // can we give 1 electron?
//              else if ((int) round(availcharge - 1.0) > 0)
//              {
//                if (_world.rank() == 0) printf("%d\tcase #2: availcharge - 1.0 = %.8f\n", i, availcharge - 1.0);
//                _occs[i] = 1.0;
//                availcharge -= 1.0;
//              }
//              else
//              {
//                if (_world.rank() == 0) printf("case #2: availcharge = %.8f\n", availcharge);
//                _occs[i] = 0.0;
//              }
//            }
//          }
          // Construct and diagonlize Fock matrix
          ctensorT fock = potential + kinetic;
          fock = 0.5 * (fock + transpose(fock));
          ctensorT c; rtensorT e;
          sygv(fock, overlap, 1, c, e);

          // diagonlize kinetic
          ctensorT ck; rtensorT ek;
          sygv(kinetic, overlap, 1, ck, ek);
          // diagonalize potential
          ctensorT cp; rtensorT ep;
          sygv(potential, overlap, 1, cp, ep);
          // diagonalize overlap
          ctensorT co; rtensorT eo;
          syev(overlap, co, eo);

          if (_world.rank() == 0)
          {
            print("fock eigenvectors dims:",c.dim(0),c.dim(1));
            print("fock eigenvectors:");
            print(c);
          }

          if (_world.rank() == 0)
          {
            print("kinetic eigenvalues");
            print(ek);
          }

          if (_world.rank() == 0)
          {
            print("potential eigenvalues");
            print(ep);
          }

          if (_world.rank() == 0)
          {
            print("overlap eigenvalues");
            print(eo);
          }
          if (_world.rank() == 0)
          {
            print("fock eigenvalues");
            print(e);
          }

          compress(_world, ao);
          _world.gop.fence();
          // Take linear combinations of the gaussian basis orbitals as the starting
          // orbitals for solver
          vecfuncT tmp_orbitals = transform(_world, ao, c(_, Slice(0, nao - 1)));
          _world.gop.fence();
          truncate(_world, tmp_orbitals);
          normalize(_world, tmp_orbitals);

          // WSTHORNTON
          // if at k-point == (0.0,0.0,0.0) basis function should be real
//          if ((abs(kpt.k[0]) < 1e-12) && (abs(kpt.k[1]) < 1e-12) &&
//              (abs(kpt.k[2]) < 1e-12))
//          {
//            if (_world.rank() == 0) printf("checking gamma point orbitals \n\n");
//            for (unsigned int io = 0; io < tmp_orbitals.size(); io++)
//            {
//              if (!is_real<T,NDIM>(tmp_orbitals[io]))
//              {
//                double t1 = (imag(tmp_orbitals[io])).norm2();
//                if (_world.rank() == 0)
//                  printf("ERROR: mo function imaginary: %5.3d     %15.8f\n\n", io, t1);
//              }
//            }
//          }

          // Build the overlap matrix
          if (_world.rank() == 0) print("Building overlap matrix ...\n\n");
          ctensorT overlap2 = matrix_inner(_world, tmp_orbitals, tmp_orbitals, true);

          rtensorT tmp_eigs = e(Slice(0, nao - 1));

          if (_world.rank() == 0) printf("(%8.4f,%8.4f,%8.4f)\n",kpt.k[0], kpt.k[1], kpt.k[2]);
          if (_world.rank() == 0) print(tmp_eigs);
          if (_world.rank() == 0) print("\n");

          //if (_world.rank() == 0) print("kinetic energy for kp = ", kp);
          //if (_world.rank() == 0) print(kinetic);
          //if (_world.rank() == 0) print("\n");

          // DEBUG
          if (_params.print_matrices && _world.rank() == 0) {
              printf("Overlap: \n");
              for (int i = 0; i < kinetic.dim(0); i++)
                  {
                      for (int j = 0; j < kinetic.dim(1); j++)
                          {
                              printf("%10.5f", real(overlap(i,j)));
                          }
                      printf("\n");
                  }
              printf("\n");
              printf("\n");

              printf("Kinetic: \n");
              for (int i = 0; i < kinetic.dim(0); i++)
                  {
                      for (int j = 0; j < kinetic.dim(1); j++)
                          {
                              printf("%10.5f", real(kinetic(i,j)));
                          }
                      printf("\n");
                  }
              printf("\n");
              printf("\n");

              printf("V: \n");
              for (int i = 0; i < potential.dim(0); i++)
                  {
                      for (int j = 0; j < potential.dim(1); j++)
                          {
                              printf("%10.5f", real(potential(i,j)));
                          }
                      printf("\n");
                  }
              printf("\n");
              printf("\n");

              printf("Fock: \n");
              for (int i = 0; i < fock.dim(0); i++)
                  {
                      for (int j = 0; j < fock.dim(1); j++)
                          {
                              printf("%10.5f", real(fock(i,j)));
                          }
                      printf("\n");
                  }
              printf("\n");
              printf("\n");

              printf("New overlap: \n");
              for (int i = 0; i < overlap2.dim(0); i++)
                  {
                      for (int j = 0; j < overlap2.dim(1); j++)
                          {
                              printf("%10.5f", real(overlap2(i,j)));
                          }
                      printf("\n");
                  }
              printf("\n");
              printf("\n");
          }

          // Fill in orbitals and eigenvalues
          int kend = kp + _params.nbands;
          kpt.begin = kp;
          kpt.end = kend;
          for (int oi = kp, ti = 0; oi < kend; oi++, ti++)
          {
            //if (_world.rank() == 0) print(oi, ti, kpt.begin, kpt.end);
            // normalize the orbitals
            //tmp_orbitals[ti].scale(1.0/tmp_orbitals[ti].norm2());
            _phisa.push_back(tmp_orbitals[ti]);
            _phisb.push_back(tmp_orbitals[ti]);
            _eigsa[oi] = tmp_eigs[ti];
            _eigsb[oi] = tmp_eigs[ti];
          }
        }
      }
    }
    //*************************************************************************

    //*************************************************************************
    // Constructor
    Solver(World& world,
           const rfuntionT& vnucrhon,
           const vecfuncT& phis,
           const std::vector<T>& eigs,
           const ElectronicStructureParams& params,
           MolecularEntity mentity)
       : _world(world), _vnucrhon(vnucrhon), _phisa(phis), _phisb(phis),
       _eigsa(eigs), _eigsb(eigs), _params(params), _mentity(mentity)
    {
      _residual = 1e5;
      if (params.periodic)
      {
        Tensor<double> box = FunctionDefaults<NDIM>::get_cell_width();
        //_cop = PeriodicCoulombOpPtr<T,NDIM>(const_cast<World&>(world),
        //    FunctionDefaults<NDIM>::get_k(), params.lo, params.thresh * 0.1, box);
        _cop = CoulombOperatorPtr(const_cast<World&>(world), params.lo, params.thresh * 0.1);
      }
      else
      {
        _cop = CoulombOperatorPtr(const_cast<World&>(world), params.lo, params.thresh * 0.1);
      }

      if (params.ispotential)
      {
        _vnuc = copy(_vnucrhon);
      }
      else
      {
        _vnuc = apply(*_cop, _vnucrhon);
      }
    }
    //*************************************************************************

    //*************************************************************************
    // Constructor
    Solver(World& world,
           rfuntionT vnucrhon,
           vecfuncT phis,
           std::vector<T> eigs,
           std::vector<KPoint> kpoints,
           std::vector<double> occs,
           ElectronicStructureParams params,
           MolecularEntity mentity)
       : _world(world), _vnucrhon(vnucrhon), _phisa(phis), _phisb(phis),
         _eigsa(eigs), _eigsb(eigs), _params(params),
         _kpoints(kpoints), _occsa(occs), _occsb(occs), _mentity(mentity)
    {
      _residual = 1e5;
      if (params.periodic)
      {
        Tensor<double> box = FunctionDefaults<NDIM>::get_cell_width();
        _cop = CoulombOperatorPtr(const_cast<World&>(world), params.lo, params.thresh * 0.1);
      }
      else
      {
        _cop = CoulombOperatorPtr(const_cast<World&>(world), params.lo, params.thresh * 0.1);
      }

      if (params.ispotential)
      {
        _vnuc = copy(_vnucrhon);
      }
      else
      {
        _vnuc = apply(*_cop, _vnucrhon);
      }
    }
    //*************************************************************************

    //*************************************************************************
    virtual ~Solver()
    {
      _outputF.close();
      _matF.close();
      _eigF.close();
      delete _subspace;
    }
    //*************************************************************************

    //***************************************************************************
    // set occupation numbers (only for insulators ... no smearing)
    void set_occs2(const std::vector<KPoint>& kpoints,
                   const std::vector<double>& eigsa,
                   const std::vector<double>& eigsb,
                   std::vector<double>& occsa,
                   std::vector<double>& occsb)
    {
      for (unsigned int ik = 0; ik < kpoints.size(); ik++)
      {
        // get k-point
        KPoint k = kpoints[ik];

        // pull subset of data that corresponds to k
        const std::vector<double> k_eigsa(eigsa.begin() + k.begin, eigsa.begin() + k.end);
        const std::vector<double> k_eigsb(eigsb.begin() + k.begin, eigsb.begin() + k.end);
        std::vector<double> k_occsa(occsa.begin() + k.begin, occsa.begin() + k.end);
        std::vector<double> k_occsb(occsb.begin() + k.begin, occsb.begin() + k.end);

        // demand all vectors have the same size
        unsigned int sz = k_eigsa.size();
        MADNESS_ASSERT(k_eigsb.size() == sz);
        MADNESS_ASSERT(k_occsa.size() == sz);
        MADNESS_ASSERT(k_occsb.size() == sz);
        // concatenate eigenvalues
        std::vector<double> teigs;
        //std::copy(k_eigsa.begin(), k_eigsa.end(), teigs.begin());
        //teigs.insert(teigs.end(), k_eigsb.begin(), k_eigsb.end());
        //std::copy(k_eigsb.begin(), k_eigsb.end(), back_inserter(teigs));
        for (unsigned int ist = 0; ist < k_eigsa.size(); ist++) teigs.push_back(k_eigsa[ist]);
        for (unsigned int ist = 0; ist < k_eigsb.size(); ist++) teigs.push_back(k_eigsb[ist]);
        // indicies
        std::vector<unsigned int> inds(2*sz);
        for (unsigned int i = 0; i < 2*sz; i++) inds[i] = i;
        // sort by eigenvalue
        for (unsigned int i = 0; i < 2*sz; i++)
        {
          for (unsigned int j = i+1; j < 2*sz; j++)
          {
            if (teigs[j] < teigs[i])
            {
              double t1 = teigs[i];
              teigs[i] = teigs[j];
              teigs[j] = t1;
              int it1 = inds[i];
              inds[i] = inds[j];
              inds[j] = it1;
            }
          }
        }
        // assign occupation numbers
        double availablecharge = _params.ncharge;
        for (unsigned int i = 0; (i < 2*sz) && (availablecharge > 0.0) ; i++)
        {
          unsigned int current = inds[i];
          if (current >= sz)
          {
            current -= sz;
            k_occsb[current] = 1.0;
            availablecharge -= 1.0;
          }
          else
          {
            k_occsa[current] = 1.0;
            availablecharge -= 1.0;
          }
        }

        for (unsigned int ik1 = k.begin, ik2 = 0; ik1 < k.end; ik1++,ik2++)
        {
          occsa[ik1] = k_occsa[ik2];
          occsb[ik1] = k_occsb[ik2];
        }
      }

      for (unsigned int ik = 0; ik < kpoints.size(); ik++)
      {
        KPoint k = kpoints[ik];
        if (_world.rank() == 0)
        {
          printf("k-point is: %d: \n",ik);
        }
        for (unsigned int ist = k.begin; ist < k.end; ist++)
        {
          if (_world.rank() == 0)
          {
            printf("occa:    %12.5f          occb:    %12.5f \n",occsa[ist],occsb[ist]);
          }
        }
      }

    }
    //***************************************************************************

    //***************************************************************************
    /*!
     \ingroup periodic_solver
     \brief Compute the electronic density for either a molecular or periodic
            system.
     */
    rfuntionT compute_rho_slow(const vecfuncT& phis, std::vector<KPoint> kpoints,
                               std::vector<double> occs)
    {
      // Electron density
      rfuntionT rho = rfactoryT(_world);
      _world.gop.fence();
      if (_world.rank() == 0) _outputF << "computing rho ..." << endl;
      // Loop over k-points
      for (unsigned int kp = 0; kp < kpoints.size(); kp++)
      {
        // get k-point
        KPoint kpoint = kpoints[kp];
        // loop through bands
        for (unsigned int j = kpoint.begin; j < kpoint.end; j++)
        {
          // Get phi(j) from iterator
          const functionT& phij = phis[j];
          // Compute the j-th density
          rfuntionT prod = abs_square(phij);
          double rnrm = prod.trace();
          prod.scale(1/rnrm);
//          rho += 0.5 * _occs[j] * kpoint.weight * prod;
          rho += occs[j] * kpoint.weight * prod;
        }
      }
      rho.truncate();

      // WSTHORNTON
      // Test for George
      cfunctionT rho_complex = function_real2complex(rho);

      return rho;
    }


    rfuntionT compute_rho(const vecfuncT& phis, std::vector<KPoint> kpoints,
                          std::vector<double> occs)
    {
      if (_world.rank() == 0) _outputF << "computing rho ..." << endl;
      rfunctionT rho = rfactoryT(_world);       // Electron density

      reconstruct(_world, phis); // For max parallelism
      std::vector<rfunctionT> phisq(phis.size());
      for (unsigned int i=0; i<phis.size(); i++) {
          phisq[i] = abssq(phis[i],false);
      }
      _world.gop.fence();
      std::vector<double> phinorm = norm2s(_world, phis);

      compress(_world,phisq); // since will be using gaxpy for accumulation
      rho.compress();

      // Loop over k-points
      for (unsigned int kp = 0; kp < kpoints.size(); kp++)  {
          // get k-point
          KPoint kpoint = kpoints[kp];
          // loop through bands
          for (unsigned int j = kpoint.begin; j < kpoint.end; j++)  {
//            rho.gaxpy(1.0, phisq[j], 0.5 * _occs[j] * kpoint.weight / (phinorm[j]*phinorm[j]), false);
            rho.gaxpy(1.0, phisq[j], occs[j] * kpoint.weight / (phinorm[j]*phinorm[j]), false);
        }
      }
      _world.gop.fence();
      phisq.clear();
      rho.truncate();

      return rho;
    }


    //***************************************************************************

    //***************************************************************************
    std::vector<poperatorT> make_bsh_operators(const std::vector<T>& eigs)
    {
      // Make BSH vector
      std::vector<poperatorT> bops;
      // Get defaults
      double tol = FunctionDefaults<NDIM>::get_thresh();
      // Loop through eigenvalues, adding a BSH operator to bops
      // for each eigenvalue
      int sz = eigs.size();
      for (int i = 0; i < sz; i++)
      {
          T eps = eigs[i];
          if (eps > 0)
          {
              if (_world.rank() == 0)
              {
                  std::cout << "bsh: warning: positive eigenvalue" << i << eps << endl;
              }
              eps = -0.1;
          }
          if (_params.periodic)
          {
            Tensor<double> cellsize = FunctionDefaults<3>::get_cell_width();
            //bops.push_back(poperatorT(PeriodicBSHOpPtr<T,NDIM>(_world, sqrt(-2.0*eps), k, _params.lo, tol * 0.1,
            //    cellsize)));
            bops.push_back(poperatorT(BSHOperatorPtr3D(_world, sqrt(-2.0*eps), _params.lo, tol * 0.1)));
          }
          else
          {
            bops.push_back(poperatorT(BSHOperatorPtr3D(_world, sqrt(-2.0*eps), _params.lo, tol * 0.1)));
          }
      }
      return bops;
    }
    //*************************************************************************

    //*************************************************************************
//    void loadbal()
//    {
//        if(_world.size() == 1)
//            return;
//
//        LoadBalanceDeux<3> lb(_world);
//        lb.add_tree(_vnuc, lbcost<double,3>(1.0, 0.0), false);
//        lb.add_tree(_rhoa, lbcost<double,3>(1.0, 1.0), false);
//        for(unsigned int i = 0;i < _phisa.size();i++)
//        {
//            lb.add_tree(_phisa[i], lbcost<valueT,3>(1.0, 1.0), false);
//        }
//        if(_params.spinpol)
//        {
//            lb.add_tree(_rhob, lbcost<double,3>(1.0, 1.0), false);
//            for(unsigned int i = 0;i < _phisb.size();i++)
//            {
//                lb.add_tree(_phisa[i], lbcost<valueT,3>(1.0, 1.0), false);
//            }
//        }
//
//        FunctionDefaults<3>::redistribute(_world, lb.load_balance(6.0));
//   }
   //*************************************************************************

    //*************************************************************************
    double calculate_kinetic_energy()
    {
      double ke = 0.0;
      if (!_params.periodic)
      {
        complex_derivative_3d Dx(_world,0);
        complex_derivative_3d Dy(_world,1);
        complex_derivative_3d Dz(_world,2);
        for (unsigned int i = 0; i < _phisa.size(); i++)
        {
          functionT dpdx = Dx(_phisa[i]);
          functionT dpdy = Dy(_phisa[i]);
          functionT dpdz = Dz(_phisa[i]);
          ke += 0.5 * (real(inner(dpdx,dpdx)) + real(inner(dpdy,dpdy))
              + real(inner(dpdz,dpdz)));
        }
        if (_params.spinpol)
        {
          for (unsigned int i = 0; i < _phisb.size(); i++)
          {
            functionT dpdx = Dx(_phisb[i]);
            functionT dpdy = Dy(_phisb[i]);
            functionT dpdz = Dz(_phisb[i]);
            ke += 0.5 * (real(inner(dpdx,dpdx)) + real(inner(dpdy,dpdy))
                + real(inner(dpdz,dpdz)));
          }
        }
        else
        {
          ke *= 2.0;
        }
      }
      return ke;
    }
    //*************************************************************************

    //*************************************************************************
    /*!
     \ingroup periodic_solver
     \brief Applies the LDA effective potential to each orbital. Currently only
            lda and spin-polarized is not implemented.
     */
    void apply_potential(vecfuncT& pfuncsa,
        vecfuncT& pfuncsb, const vecfuncT& phisa,
        const vecfuncT& phisb, const rfuntionT& rhoa, const rfuntionT& rhob,
        const rfuntionT& rho)
    {
      // Nuclear and coulomb potentials
      rfuntionT vc = apply(*_cop, rho);
      rfuntionT vlocal = _vnuc + vc;
      // Calculate energies for Coulomb and nuclear
      double ce = 0.5*inner(vc,rho);
      double pe = inner(_vnuc,rho);
      double xc = 0.0;
      double ke = calculate_kinetic_energy();
      // Exchange
      if (_params.functional == 1)
      {
        // LDA, is calculation spin-polarized?
        if (_params.spinpol)
        {
        	MADNESS_EXCEPTION("Spin polarized not implemented!",0);
//          // potential
//          rfuntionT vxca = binary_op(rhoa, rhob, &::libxc_ldaop_sp<double>);
//          rfuntionT vxcb = binary_op(rhob, rhoa, &::libxc_ldaop_sp<double>);
//          pfuncsa = mul_sparse(_world, vlocal + vxca, phisa, _params.thresh * 0.1);
//          pfuncsb = mul_sparse(_world, vlocal + vxcb, phisb, _params.thresh * 0.1);
//          // energy
//          rfuntionT fca = binary_op(rhoa, rhob, &::libxc_ldaeop_sp<double>);
//          rfuntionT fcb = binary_op(rhob, rhoa, &::libxc_ldaeop_sp<double>);
//          xc = fca.trace() + fcb.trace();
        }
        else
        {
          // potential
          rfuntionT vxc = copy(rhoa);
//          vxc.unaryop(&::libxc_ldaop<double>);
          START_TIMER(_world);
          vxc.unaryop(&::ldaop<double>);
          pfuncsa = mul_sparse(_world, vlocal + vxc, phisa, _params.thresh * 0.1);
//          rfuntionT vxc2 = binary_op(rhoa, rhoa, &::libxc_ldaop_sp<double>);
//          pfuncsa = mul_sparse(_world, vlocal + vxc2, phisa, _params.thresh * 0.1);
          END_TIMER(_world,"Applying LDA potential");
          // energy
          rfuntionT fc = copy(rhoa);
          fc.unaryop(&::ldaeop<double>);
          xc = fc.trace();
        }
      }
      else if (_params.functional == 2)
      {
        START_TIMER(_world);
        pfuncsa = mul_sparse(_world, vlocal, phisa, _params.thresh * 0.1);
        END_TIMER(_world,"Applying local potential");
        START_TIMER(_world);
        // gamma-point?
        if ((_params.ngridk0 == 1) && (_params.ngridk1 == 1) && (_params.ngridk2 == 1))
        {
          apply_hf_exchange3(_phisa, _phisb, pfuncsa, pfuncsb, xc);
        }
        else
        {
          apply_hf_exchange4(_phisa, _phisb, pfuncsa, pfuncsb, xc);
        }
        END_TIMER(_world, "Applying HF exchange");
      }
      std::cout.precision(8);
      if (_world.rank() == 0)
      {
        printf("Energies:\n");
        printf("Kinetic energy:     %20.10f\n", ke);
        printf("Potential energy:   %20.10f\n", pe);
        printf("Coulomb energy:     %20.10f\n", ce);
        printf("Exchange energy:    %20.10f\n", xc);
        printf("Total energy:       %20.10f\n\n", ke + pe + ce + xc);
      }
    }
    //*************************************************************************

    // working for a molecule / atom ... need to test with periodic boundary
    // conditions as a gamma point only calculation
    //*************************************************************************
//    void apply_hf_exchange2(vecfuncT& phisa, vecfuncT& phisb,
//                           vecfuncT& funcsa, vecfuncT& funcsb,
//                           double& xc)
//    {
//      Vector<double,3> q = vec(0.0,0.0,0.0);
//      SeparatedConvolution<double_complex,3> hfexop =
//          PeriodicHFExchangeOperator(_world, q, _params.lo,
//              FunctionDefaults<3>::get_thresh() * 0.1);
//      for (unsigned int i = 0; i < phisa.size(); i++)
//      {
//        bool isreal1 = is_real<T,NDIM>(phisa[i]);
//        MADNESS_ASSERT(isreal1);
////        rfunctionT phi_i = real(phisa[i]);
//        for (unsigned int j = 0; j < phisa.size(); j++)
//        {
//          bool isreal2 = is_real<T,NDIM>(phisa[j]);
//          MADNESS_ASSERT(isreal2);
////          rfunctionT phi_j = real(phisa[j]);
////          rfunctionT prod = phi_i*phi_j;
//          functionT prod = phisa[i]*phisa[j];
//          prod.truncate();
////          rfunctionT Vex = apply(*_cop,prod);
//          rfunctionT Vex = real(apply(hfexop,prod));
//          functionT tf1 = Vex*phisa[j];
//          funcsa[i] -= tf1;
//          xc -= real(inner(phisa[i],tf1));
//        }
//      }
//    }
    //*************************************************************************

    // working for a molecule / atom ... works for a gamma point calculation
    //*************************************************************************
    void apply_hf_exchange3(vecfuncT& phisa, vecfuncT& phisb,
                            vecfuncT& funcsa, vecfuncT& funcsb,
                            double& xc)
    {
      for (unsigned int j = 0; j < phisa.size(); j++)
      {
        rfunctionT phi_j = real(phisa[j]);
        // do diagonal piece first
        rfunctionT dprod = phi_j*phi_j;
        dprod.truncate();
        rfunctionT dVex = apply(*_cop,dprod);
        functionT tf_jjj = dVex*phisa[j];
        funcsa[j] -= tf_jjj;
        xc -= real(inner(phisa[j],tf_jjj));
        for (unsigned int i = j+1; i < phisa.size(); i++)
        {
          rfunctionT phi_i = real(phisa[i]);
          rfunctionT prod = phi_i*phi_j;
          prod.truncate();
          rfunctionT Vex = apply(*_cop,prod);
          // do the jij-th term
          functionT tf_jij = Vex*phisa[j];
          funcsa[i] -= tf_jij;
          xc -= real(inner(phisa[i],tf_jij));
          // do the iji-th term (in the complex case, use the complex
          // conjugate)
          functionT tf_iji = Vex*phisa[i];
          funcsa[j] -= tf_iji;
          xc -= real(inner(phisa[j],tf_iji));
        }
      }
    }
    //*************************************************************************

    //*************************************************************************
    KPoint find_kpt_from_orb(unsigned int idx)
    {
      for (unsigned int i = 0; i < _kpoints.size(); i++)
      {
        KPoint k1 = _kpoints[i];
        if (k1.is_orb_in_kpoint(idx)) return k1;
      }
      MADNESS_EXCEPTION("Error: find_kpt_from_orb: didn't find kpoint\n", 0);
    }
    //************************************************************************

    // hf exchange with k-points for periodic solid
    //*************************************************************************
    void apply_hf_exchange4(vecfuncT& phisa, vecfuncT& phisb,
                           vecfuncT& funcsa, vecfuncT& funcsb,
                           double& xc)
    {
      for (unsigned int i = 0; i < phisa.size(); i++)
      {
        functionT phi_i = phisa[i];
        KPoint k_i = find_kpt_from_orb(i);
        for (unsigned int j = 0; j < phisa.size(); j++)
        {
          functionT phi_j = phisa[j];
          KPoint k_j = find_kpt_from_orb(j);
          Vector<double,3> q = vec((k_i.k[0]-k_j.k[0])*_params.L,
                                   (k_i.k[1]-k_j.k[1])*_params.L,
                                   (k_i.k[2]-k_j.k[2])*_params.L);
          Vector<double,3> q2 = vec(k_i.k[0]-k_j.k[0],
                                    k_i.k[1]-k_j.k[1],
                                    k_i.k[2]-k_j.k[2]);
          functionT cexp =
              factoryT(_world).functor(functorT(new
                  ComplexExp<3>(q2, double_complex(1.0,0.0))));
          functionT cexp2 = conj(cexp);
          cexp.truncate();
          cexp2.truncate();
          functionT prod = phi_i*phi_j*cexp;
          SeparatedConvolution<double_complex,3> hfexop =
              PeriodicHFExchangeOperator(_world, q, _params.lo,
                  FunctionDefaults<3>::get_thresh() * 0.1);
          functionT Vex = apply(hfexop,prod);
          functionT tf1 = Vex*phisa[j]*cexp2;
          funcsa[i] -= tf1;
          xc -= real(inner(phisa[i],tf1));
        }
      }
    }
    //*************************************************************************

    // not working and not tested .... supposed to be for the periodic boundary
    // conditions with k-points
    //*************************************************************************
    void apply_hf_exchange(vecfuncT& phisa, vecfuncT& phisb,
                           vecfuncT& funcsa, vecfuncT& funcsb)
    {
      for (unsigned int ink1 = 0, ik1 = 0; ink1 < _phisa.size(); ++ink1)
      {
        for (unsigned int ink2 = 0, ik2 = 0; ink2 <= ink1; ink2++)
        {
          KPoint k1 = _kpoints[ik1];
          KPoint k2 = _kpoints[ik2];

          if (ink1 == k1.end) ik1++;
          if (ink2 == k2.end) ik2++;

          MADNESS_ASSERT(ik1 == 0);
          MADNESS_ASSERT(ik2 == 0);

          // no phase factor
          if (ik1 == ik2)
          {
//            // same state (diagonal piece)
//            if (ink1 == ink2)
//            {
//              rfunctionT prod = abs_square(phisa[ink1]);
//              rfunctionT fr = apply(*_cop,prod);
//              funcsa[ink1] -= funcsa[ink1]*fr;
//            }
//            else
            {
              Vector<double,3> q = 0.0;
              functionT f = phisa[ink1]*conj(phisa[ink2]);
              SeparatedConvolution<double_complex,3> hfexop =
                  PeriodicHFExchangeOperator(_world, q, _params.lo, FunctionDefaults<3>::get_thresh() * 0.1);
              functionT fr = apply(hfexop,f);
              funcsa[ink1] -= funcsa[ink2]*fr;
              funcsa[ink2] -= funcsa[ink1]*conj(fr);

              // test code
              rfunctionT g1 = abs(phisa[ink1]);
              rfunctionT g2 = abs(phisa[ink2]);
              rfunctionT ff = g1*g2;

              printf("norm diff: %20.8f\n", abs(ff.norm2()-f.norm2()));
              MADNESS_ASSERT(abs(ff.norm2()-f.norm2()) <= 1e-5);
              rfunctionT fr2 = apply(*_cop,ff);
              MADNESS_ASSERT(abs(fr.norm2()-fr2.norm2()) <= 1e-5);
            }
          }
//          else
//          {
//            Vector<double,3> q = VectorFactory(k1.k[0]-k2.k[0],
//                                               k1.k[1]-k2.k[1],
//                                               k1.k[2]-k2.k[2]);
//            functionT cexp = factoryT(_world).functor(functorT(new ComplexExp<3>(q, double_complex(1.0,0.0))));
//            cexp.truncate();
//
//            functionT f = phisa[ink1]*conj(phisa[ink2])*cexp;
//            SeparatedConvolution<double_complex,3> hfexop =
//                PeriodicHFExchangeOperator(_world, q, _params.lo, FunctionDefaults<3>::get_thresh() * 0.1);
//            functionT fr = apply(hfexop,f);
//            funcsa[ink1] -= phisa[ink1]*fr*conj(cexp);
//            funcsa[ink2] -= phisa[ink2]*conj(fr)*cexp;
//          }
        }
      }
    }
    //*************************************************************************

    //*************************************************************************
    void reproject()
    {
      _params.waveorder += 2;
      _params.thresh /= 100;
      FunctionDefaults<3>::set_thresh(_params.thresh);
      FunctionDefaults<3>::set_k(_params.waveorder);
      if (_world.rank() == 0) _outputF << "reprojecting to wavelet order "
          << _params.waveorder << endl;
      reconstruct(_world, _phisa);
      for(unsigned int i = 0; i < _phisa.size(); i++)
      {
        _phisa[i] = madness::project(_phisa[i], FunctionDefaults<3>::get_k(),
          FunctionDefaults<3>::get_thresh(), false);
      }
      _world.gop.fence();
      truncate(_world, _phisa);
      normalize(_world, _phisa);
      if(_params.spinpol)
      {
        reconstruct(_world, _phisb);
        for(unsigned int i = 0; i < _phisb.size(); i++)
        {
            _phisb[i] = madness::project(_phisb[i], FunctionDefaults<3>::get_k(),
              FunctionDefaults<3>::get_thresh(), false);
        }
        _world.gop.fence();
        truncate(_world, _phisb);
        normalize(_world, _phisb);
      }

      delete _cop;
      make_nuclear_potential();
      //_subspace->reproject();
      delete _subspace;
      _subspace = new Subspace<T,NDIM>(_world, _params);
    }
    //*************************************************************************

    //*************************************************************************
    void solve()
    {

      // WSTHORNTON
      // multiply gamma-point WF by random phase
//      double_complex t1 = RandomValue<double_complex>();
//      double_complex t2 = exp(t1);
//      std::vector<double_complex> phase(_phisa.size(), t2);
//      scale(_world, _phisa, phase);

      // keep track of how many iterations have gone by without reprojecting
      int rit = 0;
      int rpthresh = 20;
      for (_it = 0; _it < _params.maxits && _residual > _params.rcriterion; _it++, rit++)
      {
        // should we reproject?
        if ((_it > 0) && ((_residual < 20*_params.thresh) || (rit == rpthresh)))
        {
          // reproject orbitals and reset threshold
          reproject();
          rit = 0;
        }

        if (_world.rank() == 0) _outputF << "_it = " << _it << endl;

        // Set occupation numbers
        set_occs2(_kpoints,_eigsa,_eigsb,_occsa,_occsb);

        // Compute density
        _rhoa = compute_rho(_phisa, _kpoints, _occsa);
        _rhob = (_params.spinpol) ? compute_rho(_phisb, _kpoints, _occsb) : _rhoa;
        _rho = _rhoa + _rhob;
        double rtrace = _rho.trace();
        if (_world.rank() == 0) _outputF << "trace of rho" << rtrace << endl;

//        if(_it < 2 || (_it % 10) == 0)
//        {
//          START_TIMER(_world);
//          //loadbal();
//          END_TIMER(_world, "Load balancing");
//        }

        std::vector<functionT> pfuncsa =
                zero_functions<valueT,NDIM>(_world, _phisa.size());
        std::vector<functionT> pfuncsb =
                zero_functions<valueT,NDIM>(_world, _phisb.size());

        // Apply the potentials to the orbitals
        if (_world.rank() == 0) _outputF << "applying potential ...\n" << endl;
        //START_TIMER(_world);
        apply_potential(pfuncsa, pfuncsb, _phisa, _phisb, _rhoa, _rhob, _rho);
        //END_TIMER(_world,"apply potential");

        // Do right hand side for all k-points
        std::vector<double> alpha(_phisa.size(), 0.0);
        do_rhs(_phisa, pfuncsa, _kpoints, alpha, _eigsa);

        // Make BSH Green's function
        std::vector<poperatorT> bopsa = make_bsh_operators(alpha);
        std::vector<T> sfactor(pfuncsa.size(), -2.0);
        scale(_world, pfuncsa, sfactor);

        // Apply Green's function to orbitals
        if (_world.rank() == 0) _outputF << "applying BSH operator ...\n" << endl;
        truncate<valueT,NDIM>(_world, pfuncsa);
        START_TIMER(_world);
        std::vector<functionT> tmpa = apply(_world, bopsa, pfuncsa);
        END_TIMER(_world,"apply BSH");
        bopsa.clear();

        // Do other spin
        vecfuncT tmpb = zero_functions<valueT,NDIM>(_world, _phisb.size());
        if (_params.spinpol)
        {
          alpha = std::vector<double>(_phisb.size(), 0.0);
          do_rhs(_phisb, pfuncsb,  _kpoints, alpha, _eigsb);
          std::vector<poperatorT> bopsb = make_bsh_operators(alpha);
          scale(_world, pfuncsb, sfactor);
          truncate<valueT,NDIM>(_world, pfuncsb);
          tmpb = apply(_world, bopsb, pfuncsb);
          bopsb.clear();
        }

        // Update orbitals
        update_orbitals(tmpa, tmpb, _kpoints);
        save_orbitals();
      }

      if (_params.plotorbs)
      {
        std::vector<long> npt(3,101);
        for (unsigned int ik = 0; ik < _kpoints.size(); ik++)
        {
          KPoint kpoint = _kpoints[ik];
          int ist = 0;
          for (unsigned int kst = kpoint.begin; kst < kpoint.end; kst++, ist++)
          {
            std::ostringstream strm;
            strm << "unk_" << ik << "_" << ist << ".dx";
            std::string fname = strm.str();
            plotdx(_phisa[kst], fname.c_str(), FunctionDefaults<3>::get_cell(), npt);
          }
        }
      }
      save_orbitals();
    }
    //*************************************************************************

    //*************************************************************************
    ctensorT matrix_exponential(const ctensorT& A) {
        const double tol = 1e-13;
        MADNESS_ASSERT(A.dim(0) == A.dim(1));

        // Scale A by a power of 2 until it is "small"
        double anorm = A.normf();
        int n = 0;
        double scale = 1.0;
        while (anorm*scale > 0.1)
        {
            n++;
            scale *= 0.5;
        }
        tensorT B = scale*A;    // B = A*2^-n

        // Compute exp(B) using Taylor series
        ctensorT expB = ctensorT(2, B.dims());
        for (int i = 0; i < expB.dim(0); i++) expB(i,i) = std::complex<T>(1.0,0.0);

        int k = 1;
        ctensorT term = B;
        while (term.normf() > tol)
        {
            expB += term;
            term = inner(term,B);
            k++;
            term.scale(1.0/k);
        }

        // Repeatedly square to recover exp(A)
        while (n--)
        {
            expB = inner(expB,expB);
        }

        return expB;
    }
    //*************************************************************************

    //*************************************************************************
    template<typename Q>
    void print_tensor2d(ostream& os, Tensor<Q> t)
    {
      os.precision(5);
      for (int i = 0; i < t.dim(0); i++)
      {
        for (int j = 0; j < t.dim(1); j++)
        {
          os << t(i,j) << setw(12);
        }
        os << endl;
      }
      os << endl;
    }
    //*************************************************************************

    //*************************************************************************
    void do_rhs(vecfuncT& wf,
                vecfuncT& vwf,
                std::vector<KPoint> kpoints,
                std::vector<T>& alpha,
                std::vector<double>& eigs)
    {
      // tolerance
      double trantol = 0.1*_params.thresh/std::min(30.0,double(wf.size()));
      double thresh = 1e-4;

      if (_world.rank() == 0) _eigF << "Iteration: " << _it << endl;
      for (unsigned int kp = 0; kp < kpoints.size(); kp++)
      {
        // Get k-point and orbitals for this k-point
        KPoint kpoint = kpoints[kp];
        double k0 = kpoint.k[0];
        double k1 = kpoint.k[1];
        double k2 = kpoint.k[2];
        // Extract the relevant portion of the list of orbitals and the list of the
        // V times the orbitals
        vecfuncT k_wf(wf.begin() + kpoint.begin, wf.begin() + kpoint.end);
        vecfuncT k_vwf(vwf.begin() + kpoint.begin, vwf.begin() + kpoint.end);
        // Build fock matrix
        tensorT fock = build_fock_matrix(k_wf, k_vwf, kpoint);
        tensorT overlap = matrix_inner(_world, k_wf, k_wf, true);
        // Do right hand side stuff for kpoint
        bool isgamma = (is_equal(k0,0.0,1e-5) &&
                        is_equal(k1,0.0,1e-5) &&
                        is_equal(k2,0.0,1e-5));
        if (_params.periodic && !isgamma) // Non-zero k-point
        {
          // Do the gradient term and k^2/2
          vecfuncT d_wf = zero_functions<valueT,NDIM>(_world, k_wf.size());
          complex_derivative_3d Dx(_world,0);
          complex_derivative_3d Dy(_world,1);
          complex_derivative_3d Dz(_world,2);
          for (unsigned int i = 0; i < k_wf.size(); i++)
          {
            // gradient
            functionT dx_wf = Dx(k_wf[i]);
            functionT dy_wf = Dy(k_wf[i]);
            functionT dz_wf = Dz(k_wf[i]);
            d_wf[i] = std::complex<T>(0.0,k0)*dx_wf +
                      std::complex<T>(0.0,k1)*dy_wf +
                      std::complex<T>(0.0,k2)*dz_wf;
            // k^/2
            double ksq = k0*k0 + k1*k1 + k2*k2;
            k_vwf[i] += 0.5 * ksq * k_wf[i];
          }
          gaxpy(_world, 1.0, k_vwf, -1.0, d_wf);
        }

        if (_params.canon) // canonical orbitals
        {
          ctensorT U; rtensorT e;
          sygv(fock, overlap, 1, U, e);

          unsigned int nmo = k_wf.size();
          // Fix phases.
          long imax;
          for (long j = 0; j < nmo; j++)
          {
              // Get index of largest value in column
              U(_,j).absmax(&imax);
              T ang = arg(U(imax,j));
              std::complex<T> phase = exp(std::complex<T>(0.0,-ang));
              // Loop through the rest of the column and divide by the phase
              for (long i = 0; i < nmo; i++)
              {
                U(i,j) *= phase;
              }
          }

          // Within blocks with the same occupation number attempt to
          // keep orbitals in the same order (to avoid confusing the
          // non-linear solver).  Have to run the reordering multiple
          // times to handle multiple degeneracies.
          for (int pass = 0; pass < 5; pass++)
          {
              long j;
              for (long i = 0; i < nmo; i++)
              {
                U(_, i).absmax(&j);
                if (i != j)
                {
                  tensorT tmp = copy(U(_, i));
                  U(_, i) = U(_, j);
                  U(_, j) = tmp;
                  //swap(e[i], e[j]);
                  T ti = e[i];
                  T tj = e[j];
                  e[i] = tj; e[j] = ti;
                }
              }
          }

          // Rotations between effectively degenerate states confound
          // the non-linear equation solver ... undo these rotations
          long ilo = 0; // first element of cluster
          while (ilo < nmo-1) {
              long ihi = ilo;
              while (fabs(real(e[ilo]-e[ihi+1])) < thresh*10.0*max(fabs(real(e[ilo])),1.0)) {
                  ihi++;
                  if (ihi == nmo-1) break;
              }
              long nclus = ihi - ilo + 1;
              if (nclus > 1) {
                  if (_world.rank() == 0) print("   found cluster", ilo, ihi);
                  tensorT q = copy(U(Slice(ilo,ihi),Slice(ilo,ihi)));
                  //print(q);
                  // Special code just for nclus=2
                  // double c = 0.5*(q(0,0) + q(1,1));
                  // double s = 0.5*(q(0,1) - q(1,0));
                  // double r = sqrt(c*c + s*s);
                  // c /= r;
                  // s /= r;
                  // q(0,0) = q(1,1) = c;
                  // q(0,1) = -s;
                  // q(1,0) = s;

                  // Iteratively construct unitary rotation by
                  // exponentiating the antisymmetric part of the matrix
                  // ... is quadratically convergent so just do 3
                  // iterations
                  ctensorT rot = matrix_exponential(-0.5*(q - conj_transpose(q)));
                  q = inner(q,rot);
                  ctensorT rot2 = matrix_exponential(-0.5*(q - conj_transpose(q)));
                  q = inner(q,rot2);
                  ctensorT rot3 = matrix_exponential(-0.5*(q - conj_transpose(q)));
                  q = inner(rot,inner(rot2,rot3));
                  U(_,Slice(ilo,ihi)) = inner(U(_,Slice(ilo,ihi)),q);
              }
              ilo = ihi+1;
          }

          // Debug output
          if (_params.print_matrices && _world.rank() == 0)
          {
              print("Overlap matrix:");
              print(overlap);

              print("Fock matrix:");
              print(fock);

              print("U matrix: (eigenvectors)");
              print(U);
          }

          if (_params.solver == 0)
          {
            // transform orbitals and V * (orbitals)
            k_vwf = transform(_world, k_vwf, U, 1e-5 / std::min(30.0, double(k_wf.size())), false);
            k_wf = transform(_world, k_wf, U, FunctionDefaults<3>::get_thresh() / std::min(30.0, double(k_wf.size())), true);
          }

          for (unsigned int ei = kpoint.begin, fi = 0; ei < kpoint.end;
            ei++, fi++)
          {
            // Save the latest eigenvalues
            eigs[ei] = abs(e(fi,fi));
            valueT t1 = (_params.solver == 0) ? e(fi,fi) : fock(fi,fi);
            if (real(t1) > -0.1)
            {
              alpha[ei] = -0.5;
              k_vwf[fi] += (alpha[ei]-real(t1))*k_wf[fi];
            }
            else
            {
              alpha[ei] = real(t1);
            }
          }
          if (_world.rank() == 0)
          {
            _eigF << "kpt: " << kp << endl;
            _eigF << setfill('-') << setw(20) << " " << endl;
            for (unsigned int ei = 0; ei < e.dim(0); ei++)
            {
              _eigF << ei << setfill(' ') << setw(12) << real(e(ei,ei)) << endl;
            }
            _eigF << "\n\n" << endl;
          }
        }
        else // non-canonical orbitals
        {
          // diagonlize just to print eigenvalues
          tensorT overlap = matrix_inner(_world, k_wf, k_wf, true);
          ctensorT c; rtensorT e;
          sygv(fock, overlap, 1, c, e);
          for (unsigned int ei = 0; ei < e.dim(0); ei++)
          {
            double diffe = (ei == 0) ? 0.0 : real(e(ei,ei))-real(e(ei-1,ei-1));
            if (_world.rank() == 0)
              print("kpoint ", kp, "ei ", ei, "eps ", real(e(ei,ei)), "\tdiff\t", diffe);
          }

          for (unsigned int ei = kpoint.begin, fi = 0;
            ei < kpoint.end; ei++, fi++)
          {
            alpha[ei] = std::min(-0.1, real(fock(fi,fi)));
            fock(fi,fi) -= std::complex<T>(alpha[ei], 0.0);
          }

          std::vector<functionT> fwf = transform(_world, k_wf, fock, trantol);
          gaxpy(_world, 1.0, k_vwf, -1.0, fwf);
          fwf.clear();
        }
        for (unsigned int wi = kpoint.begin, fi = 0; wi < kpoint.end;
          wi++, fi++)
        {
          wf[wi] = k_wf[fi];
          vwf[wi] = k_vwf[fi];
        }
      }
    }
    //*************************************************************************

    //*************************************************************************
    tensorT build_fock_matrix(vecfuncT& psi,
                              vecfuncT& vpsi,
                              KPoint kpoint)
    {
      // Build the potential matrix
      START_TIMER(_world);
      tensorT potential = matrix_inner(_world, psi, vpsi, true);
      _world.gop.fence();
      END_TIMER(_world,"kinetic energy matrix");

      START_TIMER(_world);
      if (_world.rank() == 0) _outputF << "Building kinetic energy matrix ...\n\n" << endl;
        tensorT kinetic = ::kinetic_energy_matrix<T,NDIM>(_world, psi,
                                                  _params.periodic,
                                                  kpoint);
      END_TIMER(_world,"potential energy matrix");
      if (_world.rank() == 0) printf("\n");

      if (_world.rank() == 0) _outputF << "Constructing Fock matrix ...\n\n" << endl;
      tensorT fock = potential + kinetic;
      fock = 0.5 * (fock + transpose(fock));
      _world.gop.fence();

      return fock;
    }
    //*************************************************************************

    //*************************************************************************
    void gram_schmidt(vecfuncT& f, KPoint kpoint)
    {
      for (unsigned int fi = kpoint.begin; fi < kpoint.end; ++fi)
      {
        // Project out the lower states
        for (unsigned int fj = kpoint.begin; fj < fi; ++fj)
        {
          valueT overlap = inner(f[fj], f[fi]);
          f[fi] -= overlap*f[fj];
        }
        f[fi].scale(1.0/f[fi].norm2());
      }
    }
    //*************************************************************************

    vecfuncT compute_residual(const vecfuncT& awfs,
                              const vecfuncT& bwfs)
    {
      // vector of residual functions
      vecfuncT rm = sub(_world, _phisa, awfs);
      // if spin-polarized insert beta spin-orbital residual functions
      if (_params.spinpol)
      {
        vecfuncT br = sub(_world, _phisb, bwfs);
        rm.insert(rm.end(), br.begin(), br.end());
      }
      // scalar residual
      std::vector<double> rnvec = norm2s<valueT,NDIM>(_world, rm);
      double rnorm = 0.0;
      for (unsigned int i = 0; i < rnvec.size(); i++) rnorm += rnvec[i];
      // renormalize and print
      _residual = rnorm / rnvec.size();
      if (_world.rank() == 0) _outputF << "\nResiduals\n---------" << endl;
      if (_world.rank() == 0) _outputF << std::setiosflags(std::ios::scientific) << "residual = " << _residual << endl;
      if (_world.rank() == 0)
      {
        _outputF << endl;
        for (unsigned int i = 0; i < rnvec.size(); i++)
        {
          _outputF << "residual" << i << "\t" << rnvec[i] << endl;
        }
        _outputF << endl;
      }

//      if (_world.rank() == 0)
//      {
//        _outputF << endl;
//        printf("Occupations\n-----------\n");
//        for (unsigned int i = 0; i < rnvec.size(); i++)
//        {
//          _outputF << "occ " << i << _occs[i] << endl;
//        }
//        _outputF << endl;
//      }

      return rm;
    }

    //*************************************************************************
    void update_orbitals(vecfuncT& awfs,
                         vecfuncT& bwfs,
                         std::vector<KPoint> kpoints)
    {
      // truncate before we do anyting
      truncate<valueT,NDIM> (_world, awfs);
      truncate<valueT,NDIM> (_world, _phisa);
      if (_params.spinpol)
      {
        truncate<valueT,NDIM> (_world, bwfs);
        truncate<valueT,NDIM> (_world, _phisb);
      }
      // compute residual
      vecfuncT rm = compute_residual(awfs, bwfs);
      if (_params.solver > 0 && _params.maxsub > 1)
      {
        // nonlinear solver
        _subspace->update_subspace(awfs, bwfs, _phisa, _phisb, rm);
      }

      // do step restriction
      step_restriction(_phisa, awfs, 0);
      if (_params.spinpol)
      {
        step_restriction(_phisb, bwfs, 1);
      }
      // do gram-schmidt
      for (unsigned int kp = 0; kp < kpoints.size(); kp++)
      {
        gram_schmidt(awfs, kpoints[kp]);
        if (_params.spinpol)
        {
          gram_schmidt(bwfs, kpoints[kp]);
        }
      }

      // fix occupation numbers


      // update alpha and beta orbitals
      truncate<valueT,NDIM>(_world, awfs);
      for (unsigned int ai = 0; ai < awfs.size(); ai++) {
          _phisa[ai] = awfs[ai].scale(1.0 / awfs[ai].norm2());
      }
      if (_params.spinpol)
      {
        truncate<valueT,NDIM>(_world, bwfs);
        for (unsigned int bi = 0; bi < bwfs.size(); bi++) {
            _phisb[bi] = bwfs[bi].scale(1.0 / bwfs[bi].norm2());
        }
      }
    }
    //*************************************************************************

    //*************************************************************************
    void step_restriction(vecfuncT& owfs,
                          vecfuncT& nwfs,
                          int aorb)
    {
      std::vector<double> rnorm = norm2s(_world, sub(_world, owfs, nwfs));
      // Step restriction
      int nres = 0;
      for (unsigned int i = 0; i < owfs.size(); i++)
      {
        if (rnorm[i] > _params.maxrotn)
        {
          double s = _params.maxrotn / rnorm[i];
          nres++;
          if (_world.rank() == 0)
          {
            if (!aorb && nres == 1) _outputF << "  restricting step for alpha orbitals:" << endl;
            if (aorb && nres == 1) _outputF << "  restricting step for beta orbitals:" << endl;
            _outputF << i;
          }
          nwfs[i].gaxpy(s, owfs[i], 1.0 - s, false);
        }
      }
      if (nres > 0 && _world.rank() == 0) printf("\n");
      _world.gop.fence();
    }
    //*************************************************************************

    //*************************************************************************
    void fix_occupations(const std::vector<T>& eps,
                         std::vector<double>& occs)
    {
      // Find max/min eigenvalues
      double emax = eps[0];
      double emin = emax;
      for (int i = 0; i < eps.size(); i++)
      {
        emax = (eps[i] > emax) ? eps[i] : emax;
        emin = (eps[i] < emin) ? eps[i] : emin;
      }

      int maxits = 1000;
      // This is hardcoded to 2.0 (non-spinpolarized case) for now.
      double occmax = 2.0;
      // Fermi energy
      double efermi = 0.0;
      // Use bisection method to find the fermi energy and update occupation numbers
      bool bstop = false;
      // Some smoothing parameter
      double t1 = 1.0/_params.swidth;
      for (int it = 0; (it < maxits)&&(!bstop); it++)
      {
        // Proposed fermi energy
        efermi = 0.5 * (emax + emin);
        // Accumulated charge
        double charge = 0.0;
        // Loop over all eigenvalues and count the charge
        for (int i = 0; i < eps.size(); i++)
        {
          double x = (efermi-eps[i]) * t1;
          // need to add some smearing function here
          occs[i] = occmax*stheta_fd(x);
          charge += _kpoints[i].weight() * occs[i];
        }
        if (fabs(emax-emin) < 1e-5)
          bstop = true;
        else if (charge < _params.ncharge)
          emin = efermi;
        else
          emax = efermi;
      }
    }
    //*************************************************************************

//    //*************************************************************************
//    void update_eigenvalues(const vecfuncT& wavefs,
//        const vecfuncT& pfuncs, const vecfuncT& phis,
//        std::vector<T>& eigs)
//    {
//      // Update e
//      if (_world.rank() == 0) printf("Updating e ...\n\n");
//      for (unsigned int ei = 0; ei < eigs.size(); ei++)
//      {
//        functionT r = wavefs[ei] - phis[ei];
//        double tnorm = wavefs[ei].norm2();
//        // Compute correction to the eigenvalues
//        T ecorrection = -0.5*real(inner(pfuncs[ei], r)) / (tnorm*tnorm);
//        T eps_old = eigs[ei];
//        T eps_new = eps_old + ecorrection;
////        if (_world.rank() == 0) printf("ecorrection = %.8f\n\n", ecorrection);
////        if (_world.rank() == 0) printf("eps_old = %.8f eps_new = %.8f\n\n", eps_old, eps_new);
//        // Sometimes eps_new can go positive, THIS WILL CAUSE THE ALGORITHM TO CRASH. So,
//        // I bounce the new eigenvalue back into the negative side of the real axis. I
//        // keep doing this until it's good or I've already done it 10 times.
//        int counter = 50;
//        while (eps_new >= 0.0 && counter < 20)
//        {
//          // Split the difference between the new and old estimates of the
//          // pi-th eigenvalue.
//          eps_new = eps_old + 0.5 * (eps_new - eps_old);
//          counter++;
//        }
//        // Still no go, forget about it. (1$ to Donnie Brasco)
//        if (eps_new >= 0.0)
//        {
//          if (_world.rank() == 0) printf("FAILURE OF WST: exiting!!\n\n");
//          _exit(0);
//        }
//        // Set new eigenvalue
//        eigs[ei] = eps_new;
//      }
//    }
//    //*************************************************************************

//    //*************************************************************************
//    double get_eig(int indx)
//    {
//      return _solver->get_eig(indx);
//    }
//    //*************************************************************************
//
//    //*************************************************************************
//    functionT get_phi(int indx)
//    {
//      return _solver->get_phi(indx);
//    }
//    //*************************************************************************
//
//    //*************************************************************************
//    const std::vector<double>& eigs()
//    {
//      return _solver->eigs();
//    }
//    //*************************************************************************
//
//    //*************************************************************************
//    const vecfuncT& phis()
//    {
//      return _solver->phis();
//    }
//    //*************************************************************************

  };
  //***************************************************************************

}
#define SOLVER_H_

#endif /* SOLVER_H_ */
