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
#ifndef HARTREEFOCK_H_
#define HARTREEFOCK_H_

#include <mra/mra.h>
#include <world/world.h>
#include <vector>

#include "eigsolver.h"

namespace madness
{

  //***************************************************************************
  template <typename T, int NDIM>
  // Typedef's
  class HartreeFockNuclearPotentialOp : public EigSolverOp<T,NDIM>
  {
  public:
    typedef Function<T,NDIM> funcT;
    //*************************************************************************
    // Constructor
    HartreeFockNuclearPotentialOp(World& world, funcT V, double coeff, double thresh);
    //*************************************************************************

    //*************************************************************************
    // Is there an orbitally-dependent term?
    virtual bool is_od() {return false;}
    //*************************************************************************

    //*************************************************************************
    // Is there a density-dependent term?
    virtual bool is_rd() {return true;}
    //*************************************************************************

    //*************************************************************************
    virtual funcT op_r(const funcT& rho, const funcT& rhon, const funcT& psi);
    //*************************************************************************

  private:
    //*************************************************************************
    funcT _V;
    //*************************************************************************
  };
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  // Typedef's
  class HartreeFockCoulombOp : public EigSolverOp<T,NDIM>
  {
  public:
    typedef Function<T,NDIM> funcT;
    //*************************************************************************
    // Constructor
    HartreeFockCoulombOp(World& world, double coeff, double thresh);
    //*************************************************************************

    //*************************************************************************
    // Is there an orbitally-dependent term?
    virtual bool is_od() {return false;}
    //*************************************************************************

    //*************************************************************************
    // Is there a density-dependent term?
    virtual bool is_rd() {return true;}
    //*************************************************************************

    //*************************************************************************
    virtual funcT op_r(const funcT& rho, const funcT& rhon, const funcT& psi);
    //*************************************************************************

  };
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  // Typedef's
  class HartreeFockExchangeOp : public EigSolverOp<T,NDIM>
  {
  public:
    typedef Function<T,NDIM> funcT;
    //*************************************************************************
    // Constructor
    HartreeFockExchangeOp(World& world, double coeff, double thresh);
    //*************************************************************************

    //*************************************************************************
    // Is there an orbitally-dependent term?
    virtual bool is_od() {return true;}
    //*************************************************************************

    //*************************************************************************
    // Is there a density-dependent term?
    virtual bool is_rd() {return false;}
    //*************************************************************************

    //*************************************************************************
    virtual funcT op_o(const std::vector<funcT>& phis, const funcT& psi);
    //*************************************************************************

  };
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  class HartreeFock : public IEigSolverObserver<T,NDIM>
  {
  public:
    typedef Function<T,NDIM> funcT;
    //*************************************************************************
    // Constructor
    HartreeFock(World& world, funcT V, std::vector<funcT> phis,
      std::vector<double> eigs, bool bCoulomb, bool bExchange, double thresh);
    //*************************************************************************

    //*************************************************************************
    // Destructor
  	virtual ~HartreeFock();
    //*************************************************************************

    //*************************************************************************
  	void hartree_fock(int maxits);
    //*************************************************************************

    //*************************************************************************
    funcT calculate_coulomb(funcT psi);
    //*************************************************************************

    //*************************************************************************
    funcT calculate_exchange(funcT psi);
    //*************************************************************************

    //*************************************************************************
    double calculate_ke_sp(funcT psi);
    //*************************************************************************

    //*************************************************************************
    double calculate_pe_sp(funcT psi);
    //*************************************************************************

    //*************************************************************************
    double calculate_coulomb_energy(const std::vector<funcT>& phis,
        const funcT& psi);
    //*************************************************************************

    //*************************************************************************
    double calculate_exchange_energy(const std::vector<funcT>& phis,
        const funcT& psi);
    //*************************************************************************

    //*************************************************************************
    T matrix_element(const funcT& phii, const funcT& phij)
    {
      return _solver->matrix_element(phii, phij);
    }
    //*************************************************************************

    //*************************************************************************
    virtual void iterateOutput(const std::vector<funcT>& phis,
        const std::vector<double>& eigs, const funcT& rho, const int& iter);
    //*************************************************************************

    //*************************************************************************
    bool include_coulomb()
    {
      return _bCoulomb;
    }
    //*************************************************************************

    //*************************************************************************
    bool include_exchange()
    {
      return _bExchange;
    }
    //*************************************************************************

    //*************************************************************************
    double get_eig(int indx)
    {
      return _solver->get_eig(indx);
    }
    //*************************************************************************

    //*************************************************************************
    funcT get_phi(int indx)
    {
      return _solver->get_phi(indx);
    }
    //*************************************************************************

    //*************************************************************************
    const std::vector<double>& eigs()
    {
      return _solver->eigs();
    }
    //*************************************************************************

    //*************************************************************************
    const std::vector<funcT>& phis()
    {
      return _solver->phis();
    }
    //*************************************************************************

    //*************************************************************************
    double calculate_tot_ke_sp(const std::vector<funcT>& phis);
    //*************************************************************************

    //*************************************************************************
    double calculate_tot_pe_sp(const std::vector<funcT>& phis);
    //*************************************************************************

    //*************************************************************************
    double calculate_tot_coulomb_energy(const std::vector<funcT>& phis);
    //*************************************************************************

    //*************************************************************************
    double calculate_tot_exchange_energy(const std::vector<funcT>& phis);
    //*************************************************************************

private:

    //*************************************************************************
    // Eigenvalue solver
    EigSolver<T,NDIM>* _solver;
    //*************************************************************************

    //*************************************************************************
    // Flags for whether to include the coulomb and exchange
    bool _bCoulomb;
    bool _bExchange;
    //*************************************************************************

    //*************************************************************************
    World& _world;
    //*************************************************************************

    //*************************************************************************
    funcT _V;
    //*************************************************************************

    //*************************************************************************
    double _thresh;
    //*************************************************************************

    //*************************************************************************
    World& world() {return _world;}
    //*************************************************************************

    //*************************************************************************
    double thresh() {return _thresh;}
    //*************************************************************************

  };
  //***************************************************************************

}

#endif /*HARTREEFOCK_H_*/
