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
#include "dft.h"
#include "util.h"
//#include <moldft/xc/f2c.h>
#include <vector>
#include "poperator.h"
#include "libxc.h"

typedef madness::Vector<double,3> coordT;

namespace madness
{
  //***************************************************************************
  template <typename T, int NDIM>
  DFTNuclearChargeDensityOp<T,NDIM>::DFTNuclearChargeDensityOp(World& world, funcT rhon,
      double coeff, double thresh, bool periodic) : EigSolverOp<T,NDIM>(world, coeff, thresh)
  {
    // Message for the matrix element output
    this->messageME("NuclearChargeDensityOp");
    _rhon = rhon;
    SeparatedConvolution<T,NDIM>* cop;
    if (periodic)
    {
      Tensor<double> L = FunctionDefaults<NDIM>::get_cell_width();
      cop = PeriodicCoulombOpPtr<T,NDIM>(world, FunctionDefaults<NDIM>::get_k(),
          1e-8, thresh * 0.1, L);
    }
    else
    {
      cop =
        CoulombOperatorPtr<T,NDIM>(world,
            FunctionDefaults<NDIM>::get_k(), 1e-8, thresh * 0.1);
    }
    // Apply operator to get potential
    _Vnuc = apply(*cop, rhon);
    delete cop;
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  DFTNuclearPotentialOp<T,NDIM>::DFTNuclearPotentialOp(World& world, funcT V,
      double coeff, double thresh) : EigSolverOp<T,NDIM>(world, coeff, thresh)
  {
    // Message for the matrix element output
    this->messageME("NuclearPotentialOp");
    _V = copy(V);
  }
  //***************************************************************************

  //*************************************************************************
  template <typename T, int NDIM>
  DFTCoulombOp<T,NDIM>::DFTCoulombOp(World& world, double coeff,
      double thresh) : EigSolverOp<T,NDIM>(world, coeff, thresh)
  {
    // Message for the matrix element output
    this->messageME("CoulombOp");
    // For now, no spin polarized
    _spinpol = false;
    // Create Coulomb operator
    _cop = CoulombOperatorPtr<T,NDIM>(world, FunctionDefaults<NDIM>::get_k(),
        1e-4, thresh);
    // Initialize potential
    _Vc = FunctionFactory<T,NDIM>(world);
  }
  //*************************************************************************

  //*************************************************************************
  template <typename T, int NDIM>
  DFTCoulombPeriodicOp<T,NDIM>::DFTCoulombPeriodicOp(World& world, funcT rhon, double coeff,
      double thresh) : EigSolverOp<T,NDIM>(world, coeff, thresh)
  {
    // Nuclear charge density
    _rhon = rhon;
    // Message for the matrix element output
    this->messageME("PeriodicCoulombOp");
    // For now, no spin polarized
    _spinpol = false;
    // Create Coulomb operator
    Tensor<double> L = FunctionDefaults<NDIM>::get_cell_width();
    _cop = PeriodicCoulombOpPtr<T,NDIM>(world, FunctionDefaults<NDIM>::get_k(),
        1e-8, thresh * 0.1, L);
  }
  //*************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  void DFTCoulombOp<T,NDIM>::prepare_op(Function<double,NDIM> rho)
  {
    _Vc = apply(*_cop, rho);
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  void DFTCoulombPeriodicOp<T,NDIM>::prepare_op(Function<double,NDIM> rho)
  {
    _Vc = apply(*_cop, rho + _rhon);
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  Function<T,NDIM> DFTNuclearChargeDensityOp<T,NDIM>::op_r(const funcT& rho, const funcT& psi)
  {
    funcT rfunc = _Vnuc * psi;
    return rfunc;
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  Function<T,NDIM> DFTNuclearPotentialOp<T,NDIM>::op_r(const funcT& rho, const funcT& psi)
  {
    funcT rfunc = _V * psi;
    return rfunc;
  }
  //***************************************************************************

  //*************************************************************************
  template <typename T, int NDIM>
  Function<T,NDIM> DFTCoulombOp<T,NDIM>::op_r(const funcT& rho, const funcT& psi)
  {
    funcT rfunc = _Vc * psi;
    return  rfunc;
  }
  //*************************************************************************

  //*************************************************************************
  template <typename T, int NDIM>
  Function<T,NDIM> DFTCoulombPeriodicOp<T,NDIM>::op_r(const funcT& rho, const funcT& psi)
  {
    funcT rfunc = _Vc * psi;
    return  rfunc;
  }
  //*************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  XCFunctionalLDA<T,NDIM>::XCFunctionalLDA(World& world, double coeff, double thresh)
    : EigSolverOp<T,NDIM>(world, coeff, thresh)
  {
    // Message for the matrix element output
    this->messageME("XCFunctionalLDA");
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  Function<T,NDIM> XCFunctionalLDA<T,NDIM>::op_r(const funcT& rho, const funcT& psi)
  {
    Function<T,NDIM> V_rho = copy(rho);
    V_rho.scale(0.5);
    V_rho.unaryop(&::libxc_ldaop);
    funcT rfunc = V_rho * psi;
    return rfunc;
  }
  //***************************************************************************
}

namespace madness
{
  //***************************************************************************
  template <typename T, int NDIM>
  DFT<T,NDIM>::DFT(World& world, funcT vnucrhon, std::vector<funcT> phis,
      std::vector<double> eigs, ElectronicStructureParams params)
  : _world(world), _vnucrhon(vnucrhon), _params(params)
  {

    if (world.rank() == 0 && !params.periodic) printf("DFT constructor (non-peridic) ...\n\n");
    if (world.rank() == 0 && params.periodic) printf("DFT constructor (periodic) ...\n\n");

    // Create ops list
    std::vector<EigSolverOp<T,NDIM>*> ops;
    // Add nuclear potential to ops list
//    ops.push_back(new DFTNuclearPotentialOp<T,NDIM>(world, V, 1.0, thresh));
    if (params.periodic)
    {
      ops.push_back(new DFTCoulombPeriodicOp<T,NDIM>(world, vnucrhon, 1.0, params.thresh));
    }
    else
    {
      if (params.ispotential)
      {
        ops.push_back(new DFTNuclearPotentialOp<T,NDIM>(world, vnucrhon, 1.0, params.thresh));
      }
      else
      {
        ops.push_back(new DFTNuclearChargeDensityOp<T,NDIM>(world, vnucrhon, 1.0, params.thresh, false));
      }
      ops.push_back(new DFTCoulombOp<T,NDIM>(world, 1.0, params.thresh));
    }
    _xcfunc = new XCFunctionalLDA<T,NDIM>(world, 1.0, params.thresh);
    ops.push_back(_xcfunc);

    // Create solver
    if (params.periodic)
    {
      std::vector<kvecT> kpoints;
      kvecT gammap(0.0);
      kpoints.push_back(gammap);
      _solver = new EigSolver<T,NDIM>(world, phis, eigs, ops, kpoints, params);
    }
    else
    {
      _solver = new EigSolver<T,NDIM>(world, phis, eigs, ops, params);
    }
    _solver->addObserver(this);
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  DFT<T,NDIM>::~DFT()
  {
    delete _solver;
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  void DFT<T,NDIM>::solve(int maxits)
  {
    _solver->multi_solve(maxits);
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  double DFT<T,NDIM>::calculate_ke_sp(funcT psi, bool periodic)
  {
    // Do calculation
    double kenergy = 0.0;
    for (int axis = 0; axis < 3; axis++)
    {
      funcT dpsi = diff(psi, axis);
      kenergy += 0.5 * inner(dpsi, dpsi);
    }
    return kenergy;
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  double DFT<T,NDIM>::calculate_tot_ke_sp(const std::vector<funcT>& phis, bool spinpol,
      bool periodic)
  {
    double tot_ke = 0.0;
    for (unsigned int pi = 0; pi < phis.size(); pi++)
    {
      // Get psi from collection
      const funcT psi = phis[pi];
      // Calculate kinetic energy contribution from psi
      tot_ke += calculate_ke_sp(psi);
    }
    if (!spinpol) tot_ke *= 2.0;
    return tot_ke;
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  double DFT<T,NDIM>::calculate_tot_pe_sp(const World& world, const Function<double,NDIM>& rho,
      const Function<double,NDIM>& vnucrhon, bool spinpol, const double thresh, bool periodic,
      bool ispotential)
  {
    funcT vnuc = copy(vnucrhon);
    if (!ispotential)
    {
      // Create Coulomb operator
      SeparatedConvolution<T,NDIM>* op = 0;
      if (periodic)
      {
        Tensor<double> L = FunctionDefaults<NDIM>::get_cell_width();
        op = CoulombOperatorPtr<T,NDIM>(const_cast<World&>(world),
            FunctionDefaults<NDIM>::get_k(), 1e-8, thresh * 0.1);
      }
      else
      {
        op = CoulombOperatorPtr<T,NDIM>(const_cast<World&>(world),
            FunctionDefaults<NDIM>::get_k(), 1e-8, thresh * 0.1);
      }
      // Apply Coulomb operator and trace with the density
      vnuc = apply(*op, vnucrhon);
      delete op;
    }
    double tot_pe = inner(vnuc, rho);
    return tot_pe;
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  double DFT<T,NDIM>::calculate_tot_coulomb_energy(const World& world, const Function<double, NDIM>& rho,
      bool spinpol, const double thresh, bool periodic)
  {
    // Create Coulomb operator
    SeparatedConvolution<T,NDIM>* op;
    if (periodic)
    {
      Tensor<double> L = FunctionDefaults<NDIM>::get_cell_width();
      op = PeriodicCoulombOpPtr<T,NDIM>(const_cast<World&>(world),
        FunctionDefaults<NDIM>::get_k(), 1e-4, thresh, L);
    }
    else
    {
      op = CoulombOperatorPtr<T,NDIM>(const_cast<World&>(world),
        FunctionDefaults<NDIM>::get_k(), 1e-4, thresh);
    }
    // Apply Coulomb operator and trace with the density
    funcT Vc = apply(*op, rho);

    double tot_ce = 0.5 * Vc.inner(rho);
    delete op;
    return tot_ce;
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  double DFT<T,NDIM>::calculate_tot_xc_energy(const Function<double, NDIM>& rho)
  {
    funcT enefunc = copy(rho);
    enefunc.scale(0.5);
    enefunc.unaryop(&::ldaeop);
    return enefunc.trace();
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  void DFT<T,NDIM>::iterateOutput(const std::vector<funcT>& phis,
      const std::vector<double>& eigs, const Function<double, NDIM>& rho,
      const int& iter, bool periodic)
  {
    if (iter%3 == 0)
    {
      if (world().rank() == 0) printf("Calculating energies ...\n");
      if (world().rank() == 0) printf("Calculating KE ...\n");
      double ke = DFT::calculate_tot_ke_sp(phis, false, periodic);
//      if (world().rank() == 0) printf("Calculating PE and CE...\n");
      double pe = DFT::calculate_tot_pe_sp(_world, rho, _vnucrhon, _params.spinpol, _params.thresh, _params.periodic, _params.ispotential);
      if (world().rank() == 0) printf("Calculating CE ...\n");
      double ce = DFT::calculate_tot_coulomb_energy(_world, rho, _params.spinpol, _params.thresh, _params.periodic);
      if (world().rank() == 0) printf("Calculating EE ...\n");
      double xce = DFT::calculate_tot_xc_energy(rho);
      if (world().rank() == 0) printf("Calculating NE ...\n");
      double ne = 0.0;
      if (world().rank() == 0) printf("Kinetic energy:\t\t\t\t %.8f\n", ke);
      if (world().rank() == 0) printf("Potential energy:\t\t\t %.8f\n", pe);
      if (world().rank() == 0) printf("Coulomb energy:\t\t\t\t %.8f\n", ce);
      if (world().rank() == 0) printf("XC energy:\t\t\t\t %.8f\n", xce);
      if (world().rank() == 0) printf("Total energy:\t\t\t\t %.8f\n", ke + pe + ce + xce + ne);
      if (world().rank() == 0) printf("gs ene\t\t\t\t\t%.4f\n", eigs[0]);
      if (world().rank() == 0) printf("1st es ene\t\t\t\t%.4f\n", eigs[1]);
      T mtxe = matrix_element(phis[0], phis[0]);
      if (world().rank() == 0) printf("\nKS matrix element:\t\t\t%.8f\n\n", mtxe);
      print_matrix_elements(phis[0], phis[0]);
    }
  }
  //***************************************************************************

  //***************************************************************************
//  template class DFT<double, 1>;
//  template class DFT<double, 2>;
  template class DFT<double, 3>;


//  template class DFT< std::complex<double>, 3>;
  //***************************************************************************
}
