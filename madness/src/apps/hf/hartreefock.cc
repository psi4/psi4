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
#include "hartreefock.h"

namespace madness
{

  //*************************************************************************
  template <typename T, int NDIM>
  HartreeFockNuclearPotentialOp<T,NDIM>::HartreeFockNuclearPotentialOp(World& world,
    funcT V, double coeff, double thresh) :
  EigSolverOp<T,NDIM>(world, coeff, thresh)
  {
    // Message for the matrix element output
    this->messageME("HartreeFockNuclearPotentialOp");
    _V = V;
  }
  //*************************************************************************

  //*************************************************************************
  template <typename T, int NDIM>
  HartreeFockCoulombOp<T,NDIM>::HartreeFockCoulombOp(World& world, double coeff,
      double thresh) : EigSolverOp<T,NDIM>(world, coeff, thresh)
  {
    // Message for the matrix element output
    this->messageME("HartreeFockCoulombOp");
  }
  //*************************************************************************

  //*************************************************************************
  template <typename T, int NDIM>
  HartreeFockExchangeOp<T,NDIM>::HartreeFockExchangeOp(World& world, double coeff,
      double thresh) : EigSolverOp<T,NDIM>(world, coeff, thresh)
  {
    // Message for the matrix element output
    this->messageME("HartreeFockExchangeOp");
  }
  //*************************************************************************

  //*************************************************************************
  template <typename T, int NDIM>
  Function<T,NDIM> HartreeFockNuclearPotentialOp<T,NDIM>::op_r(const funcT& rho, const funcT& rhon, const funcT& psi)
  {
    funcT rfunc = _V*psi;
    return rfunc;
  }
  //*************************************************************************

  //*************************************************************************
  template <typename T, int NDIM>
  Function<T,NDIM> HartreeFockCoulombOp<T,NDIM>::op_r(const funcT& rho, const funcT& rhon, const funcT& psi)
  {
    // Create Coulomb operator
    SeparatedConvolution<T,3> cop =
      CoulombOperator<T,3>(this->world(), FunctionDefaults<3>::get_k(), 1e-4, this->thresh());
    // Apply the Coulomb operator
    funcT Vc = apply(cop, rho);
    funcT rfunc = Vc*psi;
    return  rfunc;
  }
  //*************************************************************************

  //*************************************************************************
  template <typename T, int NDIM>
  Function<T,NDIM> HartreeFockExchangeOp<T,NDIM>::op_o(const std::vector<funcT>& phis, const funcT& rhon, const funcT& psi)
  {
    // Return value
    funcT rfunc = FunctionFactory<T,3>(this->world());
    // Create Coulomb operator
    SeparatedConvolution<T,NDIM> cop = CoulombOperator<T,NDIM>(this->world(),
        FunctionDefaults<3>::get_k(), 1e-4, this->thresh());
    // Use the psi and pj wavefunctions to build a product so that the K
    // operator can be applied to the wavefunction indexed by pj, NOT PSI.
    for (typename std::vector<funcT>::const_iterator pj = phis.begin(); pj != phis.end(); ++pj)
    {
      // Get phi(j) from iterator
      const funcT& phij = (*pj);
      // NOTE that psi is involved in this calculation
      funcT prod = phij*psi;
      // Transform Coulomb operator into a function (stubbed)
      prod.truncate(this->thresh());
      funcT Vex = apply(cop, prod);
      // NOTE that the index is j.
      rfunc += Vex*phij;
    }
    return rfunc;

  }
  //*************************************************************************

  //*************************************************************************
  template <typename T, int NDIM>
  HartreeFock<T,NDIM>::HartreeFock(World& world, funcT V, std::vector<funcT> phis,
    std::vector<double> eigs, bool bCoulomb, bool bExchange, double thresh)
   : _world(world), _V(V), _thresh(thresh)
  {
    _bCoulomb = bCoulomb;
    _bExchange = bExchange;
    // Create ops list
    std::vector<EigSolverOp<T,NDIM>*> ops;
    // Add nuclear potential to ops list
    ops.push_back(new HartreeFockNuclearPotentialOp<T,NDIM>(world, V, 1.0, thresh));
    // Check for coulomb and exchange, and add as appropriate
    if (bCoulomb)
    {
      ops.push_back(new HartreeFockCoulombOp<T,NDIM>(world, 2.0, thresh));
    }
    if (bExchange)
    {
      ops.push_back(new HartreeFockExchangeOp<T,NDIM>(world, -1.0, thresh));
    }
    // Create solver
    _solver = new EigSolver<T,NDIM>(world, phis, eigs, ops, thresh, false);
    _solver->addObserver(this);
  }
  //***************************************************************************

  //***************************************************************************
    template <typename T, int NDIM>
  HartreeFock<T,NDIM>::~HartreeFock()
  {
    delete _solver;
  }
  //***************************************************************************

  //***************************************************************************
    template <typename T, int NDIM>
  void HartreeFock<T,NDIM>::hartree_fock(int maxits)
  {
    _solver->solve(maxits);
  }
  //***************************************************************************

  //***************************************************************************
    template <typename T, int NDIM>
  double HartreeFock<T,NDIM>::calculate_ke_sp(funcT psi)
  {
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
  double HartreeFock<T,NDIM>::calculate_pe_sp(funcT psi)
  {
    funcT vpsi = _V*psi;
    return vpsi.inner(psi);
  }
  //***************************************************************************

  //***************************************************************************
    template <typename T, int NDIM>
  double HartreeFock<T,NDIM>::calculate_coulomb_energy(const std::vector<funcT>& phis, const funcT& psi)
  {
    if (include_coulomb())
    {
      // Electron density
      funcT density = FunctionFactory<T,3>(_world);
      // Create Coulomb operator
      SeparatedConvolution<T,NDIM> op =
        CoulombOperator<T,3>(_world, FunctionDefaults<3>::get_k(), 1e-4, _thresh);
      for (typename std::vector<funcT>::const_iterator pj = phis.begin(); pj != phis.end(); ++pj)
      {
        // Get phi(j) from iterator
        const funcT& phij = (*pj);
        // Compute the j-th density
        funcT prod = phij*phij;
        density += prod;
      }
      // Transform Coulomb operator into a function (stubbed)
      density.truncate(_thresh);
      funcT Vc = apply(op, density);
      // Note that we are not using psi
      // The density is built from all of the wavefunctions. The contribution
      // psi will be subtracted out later during the exchange.
      funcT vpsi = Vc*psi;
      return inner(vpsi, psi);
    }
    return 0.0;
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  double HartreeFock<T,NDIM>::calculate_exchange_energy(const std::vector<funcT>& phis,
      const funcT& psi)
  {
    // Return value
    funcT rfunc = FunctionFactory<double,NDIM>(world());
    if (include_exchange())
    {
      // Create Coulomb operator
      SeparatedConvolution<T,NDIM> op = CoulombOperator<T,NDIM>(world(),
          FunctionDefaults<NDIM>::get_k(), 1e-4, thresh());
      // Use the psi and pj wavefunctions to build a product so that the K
      // operator can be applied to the wavefunction indexed by pj, NOT PSI.
      for (typename std::vector<funcT>::const_iterator pj = phis.begin(); pj != phis.end(); ++pj)
      {
        // Get phi(j) from iterator
        const funcT& phij = (*pj);
        // NOTE that psi is involved in this calculation
        funcT prod = phij*psi;
        // Transform Coulomb operator into a function (stubbed)
        funcT Vex = apply(op, prod);
        // NOTE that the index is j.
        rfunc += Vex*phij;
      }
    }
    return inner(rfunc, psi);
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  double HartreeFock<T,NDIM>::calculate_tot_ke_sp(const std::vector<funcT>& phis)
  {
    double tot_ke = 0.0;
    for (unsigned int pi = 0; pi < phis.size(); pi++)
    {
      // Get psi from collection
      const funcT psi = phis[pi];
      // Calculate kinetic energy contribution from psi
      tot_ke += calculate_ke_sp(psi);
    }
    return tot_ke;
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  double HartreeFock<T,NDIM>::calculate_tot_pe_sp(const std::vector<funcT>& phis)
  {
    double tot_pe = 0.0;
    for (unsigned int pi = 0; pi < phis.size(); pi++)
    {
      // Get psi from collection
      const funcT psi = phis[pi];
      // Calculate potential energy contribution from psi
      tot_pe += calculate_pe_sp(psi);
    }
    return tot_pe;
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  double HartreeFock<T,NDIM>::calculate_tot_coulomb_energy(const std::vector<funcT>& phis)
  {
    double tot_ce = 0.0;
    for (unsigned int pi = 0; pi < phis.size(); pi++)
    {
      // Get psi from collection
      const funcT psi = phis[pi];
      // Calculate coulomb energy contribution from psi
      tot_ce += calculate_coulomb_energy(phis, psi);
    }
    return tot_ce;
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  double HartreeFock<T,NDIM>::calculate_tot_exchange_energy(const std::vector<funcT>& phis)
  {
    double tot_ee = 0.0;
    for (unsigned int pi = 0; pi < phis.size(); pi++)
    {
      // Get psi from collection
      const funcT psi = phis[pi];
      // Calculate exchange energy contribution from psi
      tot_ee += calculate_exchange_energy(phis, psi);
    }
    return tot_ee;
  }
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  void HartreeFock<T,NDIM>::iterateOutput(const std::vector<funcT>& phis,
      const std::vector<double>& eigs,const funcT& rho, const int& iter)
  {
    if (iter%3 == 0)
    {
      if (world().rank() == 0) printf("Calculating energies ...\n");
      if (world().rank() == 0) printf("Calculating KE ...\n");
      double ke = 2.0 * calculate_tot_ke_sp(phis);
      if (world().rank() == 0) printf("Calculating PE ...\n");
      double pe = 2.0 * calculate_tot_pe_sp(phis);
      if (world().rank() == 0) printf("Calculating CE ...\n");
      double ce = calculate_tot_coulomb_energy(phis);
      if (world().rank() == 0) printf("Calculating EE ...\n");
      double ee = calculate_tot_exchange_energy(phis);
      if (world().rank() == 0) printf("Calculating NE ...\n");
      double ne = 0.0;
      if (world().rank() == 0) printf("Kinetic energy:\t\t\t %.8f\n", ke);
      if (world().rank() == 0) printf("Potential energy:\t\t %.8f\n", pe);
      if (world().rank() == 0) printf("Two-electron energy:\t\t %.8f\n", 2.0*ce - ee);
      if (world().rank() == 0) printf("Total energy:\t\t\t %.8f\n", ke + pe + 2.0*ce - ee + ne);
      if (world().rank() == 0) printf("gs eigv: = \t\t\t %.4f\n", eigs[0]);
    }
  }
  //***************************************************************************
}
//*****************************************************************************

