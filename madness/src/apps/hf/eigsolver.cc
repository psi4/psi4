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
#define WORLD_INSTANTIATE_STATIC_TEMPLATES

#include "eigsolver.h"
#include "util.h"
#include "poperator.h"
#include "outputwriter.h"
#include <mra/operator.h>

//#define DEBUG_STREAM *(OutputWriter::instance()->debug_stream())
//#define LOG_STREAM *(OutputWriter::instance()->log_stream())
//#define EIGV_STREAM *(OutputWriter::instance()->eigv_stream())

#define DEBUG_STREAM cout
#define LOG_STREAM cout
#define EIGV_STREAM cout

using std::cout;
using std::endl;

namespace madness
{
  //***************************************************************************
  template <typename T, int NDIM>
  EigSolver<T,NDIM>::EigSolver(World& world, std::vector<funcT> phis,
      std::vector<double> eigs, std::vector< EigSolverOp<T,NDIM>* > ops,
      std::vector<kvecT> kpoints, ElectronicStructureParams params)
  : _phis(phis), _eigs(eigs), _ops(ops), _kpoints(kpoints),
    _world(world), _params(params)
  {
    // fill the occupation numbers
    int size = eigs.size();
    for (int i = 0; i < size; i++) _occs.push_back(2.0);
    _rho = EigSolver::compute_rho(phis, _occs, world);
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  EigSolver<T,NDIM>::EigSolver(World& world, std::vector<funcT> phis,
      std::vector<double> eigs, std::vector< EigSolverOp<T,NDIM>* > ops,
      ElectronicStructureParams params)
  : _phis(phis), _eigs(eigs), _ops(ops), _world(world), _params(params)
  {
    if (params.periodic)
    {
      kvecT gammap(0.0); _kpoints.push_back(gammap);
    }
    // fill the occupation numbers
    int size = eigs.size();
    for (int i = 0; i < size; i++)  _occs.push_back(2.0);
    _rho = EigSolver::compute_rho(phis, _occs, world);
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  EigSolver<T,NDIM>::~EigSolver()
  {
    // Eigsolver is responsible for deleting the ops
    for (typename std::vector< EigSolverOp<T,NDIM>* >::iterator it = _ops.begin(); it != _ops.end();
      it++) delete (*it);
    _ops.clear();
    // Clear eigenvectors
    _phis.clear();
    // Clear eigenvalues
    _eigs.clear();
    // Clear observers
//    for (typename std::vector< IEigSolverObserver<T,NDIM>* >::iterator it = _obs.begin();
//      it != _obs.end(); it++) delete (*it);
    _obs.clear();
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  Function<T, NDIM> EigSolver<T,NDIM>::compute_rho(typename std::vector<funcT> phis,
      std::vector<double> occs, const World& world)
  {
    // Electron density
    funcT rho = FunctionFactory<double,NDIM>(const_cast<World&>(world));
    // Loop over all wavefunctions to compute density
    for (unsigned int j = 0; j < phis.size(); j++)
    {
      // Get phi(j) from iterator
      const funcT& phij = phis[j];
      // Compute the j-th density
      funcT prod = square(phij);
      rho += occs[j]*prod;
    }
    rho.truncate();
    return rho;
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  void EigSolver<T,NDIM>::prepare_ops()
  {
    // Loop through all of the density-dependent ops and prepare them, i.e.
    // build the rho-dependent potentials.
    for (unsigned int oi = 0; oi < _ops.size(); oi++)
    {
      EigSolverOp<T,NDIM>* op = _ops[oi];
      // Prepare density-dependent operator
      if (op->is_rd()) op->prepare_op(_rho);
    }
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  T EigSolver<T,NDIM>::matrix_element(const funcT& phii, const funcT& phij)
  {
    double value = 0.0;
    // Kinetic energy operator
    for (int axis = 0; axis < NDIM; axis++)
    {
      funcT dpsi_j = diff(phij, axis);
      funcT dpsi_i = diff(phii, axis);
      value += 0.5 * inner(dpsi_i, dpsi_j);
    }
    // Loop through all ops
    for (unsigned int oi = 0; oi < _ops.size(); oi++)
    {
      EigSolverOp<T,NDIM>* op = _ops[oi];
      // Operate with density-dependent operator
      if (op->is_rd()) value += op->coeff() * phii.inner(op->op_r(_rho, phij));
      // Operate with orbital-dependent operator
      if (op->is_od()) value += op->coeff() * phii.inner(op->op_o(_phis, phij));
    }
    return value;
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  void EigSolver<T,NDIM>::update_occupation()
  {
//    // Find max/min eigenvalues
//    double emax = -1.0e12;
//    double emin = 1.0e12;
//    for (int i = 0; i < _eigs.size(); i++)
//    {
//      emax = (_eigs[i] > emax) ? _eigs[i] : emax;
//      emin = (_eigs[i] < emin) ? _eigs[i] : emin;
//    }
//
//    int maxits = 100;
//    // This is hardcoded to 2.0 (non-spinpolarized case) for now.
//    double occmax = 2.0;
//    // Fermi energy
//    double efermi = 0.0;
//    // Use bisection method to find the fermi energy and update occupation numbers
//    bool bstop = false;
//    for (int it = 0; (it < maxits)&&(!bstop); it++)
//    {
//      // Proposed fermi energy
//      double efermi = 0.5 * (emax + emin);
//      // Accumulated charge
//      double charge = 0.0;
//      // Some smoothing parameter
//      double t1 = 0.1;
//      // Loop over all orbitals and count the charge
//      for (int i = 0; i < _phis.size(); i++)
//      {
//        double x = (efermi-_eigs[i]) * t1;
//        // need to add some smearing function here
//        _occs[i] = occmax;
//        //charge += _kpoints[i].weight() * _occs[i];
//      }
//      if (fabs(emax-emin) < 1e-5)
//        bstop = true;
//      else if (charge < _ncharge)
//        emin = efermi;
//      else
//        emax = efermi;
//    }
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  void EigSolver<T,NDIM>::make_bsh_operators()
  {
    // Clear BSH vector
    _bops.clear();
    // Get defaults
    int k = FunctionDefaults<NDIM>::get_k();
    double tol = FunctionDefaults<NDIM>::get_thresh();
    // Loop through eigenvalues, adding a BSH operator to _bops
    // for each eigenvalue
    int sz = _phis.size();
    for (int i = 0; i < sz; i++)
    {
        double eps = _eigs[i];
        if (eps > 0)
        {
            if (_world.rank() == 0)
            {
                DEBUG_STREAM << "bsh: warning: positive eigenvalue" << i << eps << endl;
            }
            eps = -0.1;
        }
        _bops.push_back(poperatorT(BSHOperatorPtr3D<double,NDIM>(_world, sqrt(-2.0*eps), k, 1e-4, tol)));
    }
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  void EigSolver<T,NDIM>::print_matrix_elements(const funcT& phii, const funcT& phij)
  {
    T value = 0.0;
    prepare_ops();
    // Kinetic energy operator
    for (int axis = 0; axis < 3; axis++)
    {
      funcT dpsi_j = diff(phij, axis);
      funcT dpsi_i = diff(phii, axis);
      value += 0.5 * inner(dpsi_i, dpsi_j);
    }
    if (_world.rank() == 0)
    {
      DEBUG_STREAM << "***** Evaluation of matrix elements *****" << endl;
      DEBUG_STREAM << "KineticEnergyOp:\t\t\t" << value << endl;
    }

    // Loop through all ops
    for (unsigned int oi = 0; oi < _ops.size(); oi++)
    {
      value = 0.0;
      EigSolverOp<T,NDIM>* op = _ops[oi];
      // Operate with density-dependent operator
      if (op->is_rd()) value += op->coeff() * phii.inner(op->op_r(_rho, phij));
      // Operate with orbital-dependent operator
      if (op->is_od()) value += op->coeff() * phii.inner(op->op_o(_phis, phij));
      if (_world.rank() == 0)
      {
        DEBUG_STREAM << op->messsageME() << ":\t\t\t" << value << endl;
      }
    }
    if (_world.rank() == 0) printf("\n\n");
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  void EigSolver<T,NDIM>::solve(int maxits)
  {
	std::cout.setf(std::ios::fixed,std::ios::floatfield);
	std::cout.precision(8);
  double thresh = FunctionDefaults<NDIM>::get_thresh();
	for (int it = 0; it < maxits; it++)
    {
      // Since, the density has already been computed (it's fresh outta the
      // oven), go ahead and build all of the density-dependent potentials that
      // we can.
      prepare_ops();
      if (_world.rank() == 0) DEBUG_STREAM << "Iteration #" << it
        << endl << endl;
      for (unsigned int pi = 0; pi < _phis.size(); pi++)
      {
        // Get psi from collection
        funcT psi = _phis[pi];
        funcT pfunc = FunctionFactory<T,NDIM>(_world);
        // Loop through all ops
        if (_world.rank() == 0) DEBUG_STREAM << "Looping through the ops ..."
          << endl << endl;
        for (unsigned int oi = 0; oi < _ops.size(); oi++)
        {
          EigSolverOp<T,NDIM>* op = _ops[oi];
          // Operate with density-dependent operator
          if (op->is_rd()) pfunc += op->coeff() * op->op_r(_rho, psi);
          // Operate with orbital-dependent operator
          if (op->is_od()) pfunc += op->coeff() * op->op_o(_phis, psi);
        }
        if (_world.rank() == 0) DEBUG_STREAM << "Creating BSH operator ..."
          << endl << endl;
        SeparatedConvolution<T,NDIM>* op = 0;
        if (_params.periodic)
        {
          // Subtract the k dot nabla part
          kvecT k = _kpoints[pi];
          pfunc -= k[0] * diff(psi, 0);
          pfunc -= k[1] * diff(psi, 1);
          pfunc -= k[2] * diff(psi, 2);
          pfunc.scale(-2.0).truncate();
          Tensor<double> L = FunctionDefaults<NDIM>::get_cell_width();
          op = PeriodicBSHOpPtr<T,NDIM>(_world, sqrt(-2.0*_eigs[pi]),
              FunctionDefaults<NDIM>::get_k(), 1e-4, thresh, L);
        }
        else
        {
          pfunc.scale(-2.0).truncate();
          op = BSHOperatorPtr3D<T>(_world, sqrt(-2.0*_eigs[pi]),
                                   FunctionDefaults<NDIM>::get_k(), 1e-4, thresh);
        }
        // Apply the Green's function operator (stubbed)
        if (_world.rank() == 0) DEBUG_STREAM << "Applying BSH operator ..."
          << endl << endl;
        pfunc.truncate();
        funcT tmp = apply(*op, pfunc);
        tmp.truncate();
        // delete op
        delete op;
        // (Not sure whether we have to do this mask thing or not!)
        // WSTHORNTON DEBUG
        double ttnorm = tmp.norm2();
        if (_world.rank() == 0) DEBUG_STREAM << "pi = " << pi
          << "\tttnorm = " << ttnorm << endl << endl;
        if (_world.rank() == 0) DEBUG_STREAM << "Gram-Schmidt ..."
          << endl << endl;
        for (unsigned int pj = 0; pj < pi; ++pj)
        {
//          // Make sure that pi != pj
//          MADNESS_ASSERT(pi != pj);
          // Project out the lower states
          // Get other wavefunction
          funcT psij = _phis[pj];
          double overlap = inner(tmp, psij);
          tmp -= overlap*psij;
        }
        // WSTHORNTON DEBUG
        double tttnorm = tmp.norm2();
        if (_world.rank() == 0) DEBUG_STREAM << "pi = " << pi << "ttnorm = "
          << tttnorm << endl;
        // Update e
        if (_world.rank() == 0) DEBUG_STREAM << "Updating e ..."
          << endl << endl;
        funcT r = tmp - psi;
        double tnorm = tmp.norm2();
        double eps_old = _eigs[pi];
        double ecorrection = -0.5*inner(pfunc, r) / (tnorm*tnorm);
        double eps_new = eps_old + ecorrection;
        double rnorm = r.norm2();
        if (_world.rank() == 0) DEBUG_STREAM << "pi = " << pi
          << "rnorm = " << rnorm << endl << endl;
        if (_world.rank() == 0) EIGV_STREAM <<  "pi = " << pi
          << " enew = " << eps_new << " eps_old = " << eps_old << endl;
        // Sometimes eps_new can go positive, THIS WILL CAUSE THE ALGORITHM TO CRASH. So,
        // I bounce the new eigenvalue back into the negative side of the real axis. I
        // keep doing this until it's good or I've already done it 10 times.
        // WSTHORNTON DEBUG
//        double rnorm = r.norm2();
//        if (_world.rank() == 0) printf("pi = %d\trnorm = %.5f\ttnorm = %.5f\n\n", pi, rnorm, tnorm);
        int counter = 0;
        while (eps_new >= 0.0 && counter < 10)
        {
          // Split the difference between the new and old estimates of the
          // pi-th eigenvalue.
          eps_new = eps_old + 0.5*(eps_new - eps_old);
          counter++;
        }
        // Still no go, forget about it. (1$ to Donnie Brasco)
        if (eps_new >= 0.0)
        {
          LOG_STREAM << "FAILURE OF WST: exiting!!\n" << endl;
          _exit(0);
        }
        // Update the eigenvalue estimates and wavefunctions.
        tmp.truncate();
        _eigs[pi] = eps_new;
        _phis[pi] = tmp.scale(1.0/tmp.norm2());
      }
      // Update rho
//      if (_world.rank() == 0) printf("Computing new density for it == #%d\n\n", it);
      _rho = EigSolver::compute_rho(_phis, _occs, _world);
      // Trace of rho
      double rhotrace = _rho.trace();
      // Output to observables
      for (typename std::vector<IEigSolverObserver<T,NDIM>*>::iterator itr = _obs.begin(); itr
        != _obs.end(); ++itr)
      {
        (*itr)->iterateOutput(_phis, _eigs, _rho, it, _params.periodic);
      }
    }
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  void EigSolver<T,NDIM>::multi_solve(int maxits)
  {
    for (int it = 0; it < maxits; it++)
    {
      // Since, the density has already been computed (it's fresh outta the
      // oven), go ahead and build all of the density-dependent potentials that
      // we can.
      prepare_ops();
      if (_world.rank() == 0) printf("Iteration #%d\n\n", it);
      // Create empty functions for calculations
      std::vector<funcT> pfuncs(_phis.size());
      for (unsigned int pi = 0; pi < _phis.size(); pi++)
        pfuncs[pi] = FunctionFactory<T, NDIM>(_world);
      // Loop through all ops to work on a vector of functions
      if (_world.rank() == 0) printf("Looping through the ops ...\n\n");
      for (unsigned int oi = 0; oi < _ops.size(); oi++)
      {
        EigSolverOp<T,NDIM>* op = _ops[oi];
        // Operate with density-dependent operator
        if (op->is_rd()) gaxpy(_world, 1.0, pfuncs, op->coeff(), op->multi_op_r(_rho, _phis));
        // Operate with orbital-dependent operator
        if (op->is_od()) gaxpy(_world, 1.0, pfuncs, op->coeff(), op->multi_op_o(_phis));
      }
      // Make BSH operators
      if (_world.rank() == 0) printf("Creating BSH operator ...\n\n");
      make_bsh_operators();
      // Apply the Green's function operator (stubbed)
      if (_world.rank() == 0) printf("Applying BSH operator ...\n\n");
      std::vector<double> sfactor(pfuncs.size());
      for (unsigned int si = 0; si < sfactor.size(); si++) sfactor[si] = -2.0;
      scale(_world, pfuncs, sfactor);
      std::vector<funcT> tmp = apply(_world, _bops, pfuncs);
      // WSTHORNTON DEBUG
      for (unsigned int ti = 0; ti < tmp.size(); ti++)
      {
        double ttnorm = tmp[ti].norm2();
        if (_world.rank() == 0) printf("ti = %d\tttnorm = %.5f\n\n", ti, ttnorm);
      }
      // Do Gram-Schmidt
      if (_world.rank() == 0) printf("Gram-Schmidt ...\n\n");
      for (unsigned int ti = 0; ti < tmp.size(); ++ti)
      {
        // Project out the lower states
        for (unsigned int pj = 0; pj < ti; ++pj)
        {
          double overlap = inner(tmp[ti], _phis[pj]);
          tmp[ti] -= overlap*_phis[pj];
        }
      }
      _world.gop.fence();
      // WSTHORNTON DEBUG
      for (unsigned int ti = 0; ti < tmp.size(); ti++)
      {
        double ttnorm = tmp[ti].norm2();
        if (_world.rank() == 0) printf("ti = %d\tttnorm = %.5f\n\n", ti, ttnorm);
      }
      // Update e
      if (_world.rank() == 0) printf("Updating e ...\n\n");
      for (unsigned int ei = 0; ei < _eigs.size(); ei++)
      {
        funcT r = tmp[ei] - _phis[ei];
        double tnorm = tmp[ei].norm2();
        double rnorm = r.norm2();
        if (_world.rank() == 0) printf("ei = %d\trnorm = %.5f\ttnorm = %.5f\n\n", ei, rnorm, tnorm);
        // Compute correction to the eigenvalues
        double ecorrection = -0.5*inner(pfuncs[ei], r) / (tnorm*tnorm);
        double eps_old = _eigs[ei];
        double eps_new = eps_old + ecorrection;
        // Sometimes eps_new can go positive, THIS WILL CAUSE THE ALGORITHM TO CRASH. So,
        // I bounce the new eigenvalue back into the negative side of the real axis. I
        // keep doing this until it's good or I've already done it 10 times.
        int counter = 10;
        while (eps_new >= 0.0 && counter < 20)
        {
          // Split the difference between the new and old estimates of the
          // pi-th eigenvalue.
//          if (_world.rank() == 0) EIGV_STREAM  << "ei = %d\teps_new = %.5f\teps_old = %.5f\n\n"
//            << ei << eps_new, eps_old);
          eps_new = eps_old + 0.5*(eps_new - eps_old);
          counter++;
        }
        // Still no go, forget about it. (1$ to Donnie Brasco)
        if (eps_new >= 0.0)
        {
          if (_world.rank() == 0) printf("FAILURE OF WST: exiting!!\n\n");
          _exit(0);
        }
        // Set new eigenvalue
        _eigs[ei] = eps_new;
        if (_world.rank() == 0) printf("ei = %d\teps = %.5f\n\n", ei, eps_new);
      }
      // Update the eigenvalue estimates and wavefunctions.
      truncate(_world, tmp);
      for (unsigned int ti = 0; ti < tmp.size(); ti++)
      {
        _phis[ti] = tmp[ti].scale(1.0/tmp[ti].norm2());
      }
      // Update rho
      if (_world.rank() == 0) printf("Computing new density for it == #%d\n\n", it);
      _rho = EigSolver::compute_rho(_phis, _occs, _world);
      // Output to observables
      for (typename std::vector<IEigSolverObserver<T,NDIM>*>::iterator itr = _obs.begin(); itr
        != _obs.end(); ++itr)
      {
        (*itr)->iterateOutput(_phis, _eigs, _rho, it, _params.periodic);
      }
    }
  }
  //***************************************************************************

  //***************************************************************************
  template class EigSolver<double,1>;
  template class EigSolver<double,2>;
  template class EigSolver<double,3>;
  //***************************************************************************
}


