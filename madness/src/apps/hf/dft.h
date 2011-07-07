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
#ifndef DFT_H_
#define DFT_H_

#include <mra/mra.h>
#include <world/world.h>
#include <vector>

#include "eigsolver.h"
#include "electronicstructureparams.h"

class ElectronicStructureAppParams;

namespace madness
{

////*****************************************************************************
//void xc_rks_generic_lda(Tensor<double> rho_alpha,           ///< Alpha-spin density at each grid point
//                        Tensor<double> f,                   ///< Value of functional at each grid point
//                        Tensor<double> df_drho);            ///< Derivative of functional w.r.t. rho_alpha
////*****************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  class DFTNuclearPotentialOp : public EigSolverOp<T,NDIM>
  {
  public:
    typedef Function<T,NDIM> funcT;
    //*************************************************************************
    // Constructor
    DFTNuclearPotentialOp(World& world, funcT V, double coeff, double thresh);
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
    virtual funcT op_r(const funcT& rho, const funcT& psi);
    //*************************************************************************

  private:
    //*************************************************************************
    funcT _V;
    //*************************************************************************
  };
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  class DFTCoulombOp : public EigSolverOp<T,NDIM>
  {
    // Typedef's
    typedef Function<T,NDIM> funcT;
  public:
    //*************************************************************************
    // Constructor
    DFTCoulombOp(World& world, double coeff, double thresh);
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
    virtual funcT op_r(const funcT& rho, const funcT& psi);
    //*************************************************************************

    //*************************************************************************
    // Build the potential from a density if a density-dependent operator.
    virtual void prepare_op(Function<double,NDIM> rho);
    //*************************************************************************

    //*************************************************************************
    SeparatedConvolution<T,NDIM>* _cop;
    //*************************************************************************

  private:
    //*************************************************************************
    funcT _Vc;
    //*************************************************************************

    //*************************************************************************
    bool _spinpol;
    //*************************************************************************
  };
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  class DFTCoulombPeriodicOp : public EigSolverOp<T,NDIM>
  {
    // Typedef's
    typedef Function<T,NDIM> funcT;
  public:
    //*************************************************************************
    // Constructor
    DFTCoulombPeriodicOp(World& world, funcT rhon, double coeff, double thresh);
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
    virtual funcT op_r(const funcT& rho, const funcT& psi);
    //*************************************************************************

    //*************************************************************************
    // Build the potential from a density if a density-dependent operator.
    virtual void prepare_op(Function<double,NDIM> rho);
    //*************************************************************************

    //*************************************************************************
    SeparatedConvolution<T,NDIM>* _cop;
    //*************************************************************************

  private:
    //*************************************************************************
    funcT _Vc;
    //*************************************************************************

    //*************************************************************************
    bool _spinpol;
    //*************************************************************************

    //*************************************************************************
    funcT _rhon;
    //*************************************************************************

  };
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  class XCFunctionalLDA : public EigSolverOp<T,NDIM>
  {
    // Typedef's
    typedef Function<T,NDIM> funcT;
  public:
    //*************************************************************************
    // Constructor
    XCFunctionalLDA(World& world, double coeff, double thresh);
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
    virtual funcT op_r(const funcT& rho, const funcT& psi);
    //*************************************************************************
  };
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  class DFTNuclearChargeDensityOp : public EigSolverOp<T,NDIM>
  {
  public:
    typedef Function<T,NDIM> funcT;
    //*************************************************************************
    // Constructor
    DFTNuclearChargeDensityOp(World& world, funcT rhon, double coeff,
        double thresh, bool periodic);
    //*************************************************************************

    //*************************************************************************
    ~DFTNuclearChargeDensityOp()
    {
    }
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
    void prepare_op(Function<double,NDIM> rho) {}
    //*************************************************************************

    //*************************************************************************
    virtual funcT op_r(const funcT& rho, const funcT& psi);
    //*************************************************************************

  private:
    //*************************************************************************
    funcT _rhon;
    //*************************************************************************

    //*************************************************************************
    funcT _Vnuc;
    //*************************************************************************
  };
  //***************************************************************************

  //***************************************************************************
  template <typename T, int NDIM>
  class DFT : public IEigSolverObserver<T,NDIM>
  {
    // Typedef's
    typedef Function<T,NDIM> funcT;
    typedef Vector<double,NDIM> kvecT;
  public:
    //*************************************************************************
    // Constructor
    DFT(World& world, funcT vnucrhon, std::vector<funcT> phis,
      std::vector<double> eigs, ElectronicStructureParams params);
    //*************************************************************************

  	//*************************************************************************
    DFT();
    //*************************************************************************

    //*************************************************************************
    virtual ~DFT();
    //*************************************************************************

    //*************************************************************************
     void solve(int maxits);
     //*************************************************************************

     //***************************************************************************
     static double calculate_ke_sp(funcT psi, bool periodic = false);
     //***************************************************************************

     //***************************************************************************
     static double calculate_tot_ke_sp(const std::vector<funcT>& phis,
         bool spinpol, bool periodic = false);
     //***************************************************************************

     //***************************************************************************
     static double calculate_tot_pe_sp(const World& world,
         const Function<double,NDIM>& rho, const Function<double,NDIM>& vnucrhon,
         bool spinpol, const double thresh, bool periodic, bool ispotential);
     //***************************************************************************

     //***************************************************************************
     static double calculate_tot_coulomb_energy(const World& world,
         const Function<double,NDIM>& rho, bool spinpol, const double thresh,
         bool periodic = false);
     //***************************************************************************

     //***************************************************************************
     static double calculate_tot_xc_energy(const Function<double,NDIM>& rho);
     //***************************************************************************

     //*************************************************************************
     T matrix_element(const funcT& phii, const funcT& phij)
     {
       return _solver->matrix_element(phii, phij);
     }
     //*************************************************************************

     //*************************************************************************
     void print_matrix_elements(const funcT& phii, const funcT& phij)
     {
       _solver->print_matrix_elements(phii, phij);
     }
     //*************************************************************************

     //*************************************************************************
     virtual void iterateOutput(const std::vector<funcT>& phis,
         const std::vector<double>& eigs, const Function<double,NDIM>& rho,
         const int& iter, bool periodic = false);
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

  private:

      //*************************************************************************
      // Eigenvalue solver
      EigSolver<T,NDIM>* _solver;
      //*************************************************************************

      //*************************************************************************
      World& _world;
      //*************************************************************************

      //*************************************************************************
      // This variable could either be a nuclear potiential or a nuclear charge
      // density depending on the "ispotential" variable in the
      // ElectronicStructureParams class.
      Function<double,NDIM> _vnucrhon;
      //*************************************************************************

      //*************************************************************************
      // Exchange-correlation functional. Needed to compute the energy Exc[rho]
      // Gets deleted my the EigSolver class during the EigSolver destructor
      EigSolverOp<T,NDIM>* _xcfunc;
      //*************************************************************************

      //*************************************************************************
      ElectronicStructureParams _params;
      //*************************************************************************

      //*************************************************************************
      World& world() {return _world;}
      //*************************************************************************

  };
  //***************************************************************************

}
#endif /*DFT_H_*/
