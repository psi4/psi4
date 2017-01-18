/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

// usapt.h, to be included in dftsapt.h if needed.

#ifndef USAPT_H
#define USAPT_H

#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/typedefs.h"
#include "psi4/libqt/qt.h"
#include <map>

#ifdef _OPENMP
  #include <omp.h>
#endif

namespace psi {

class JK;
class Options;
class Tensor;
class AtomicDensity;
class PSIO;

namespace sapt {

/*- Open-shell generalization of SAPT0
    Will be the basis for DFTSAPT
    Will need to create better class hierarchy -*/

class USAPT0 {

private:
    Options & options_;

protected:

    // SAPT type (until properly subclassed)
    std::string type_;

    // Print flag
    int print_;
    // Debug flag
    int debug_;
    // Bench flag
    int bench_;
    
    // CPKS maximum iterations
    int cpks_maxiter_;
    // CPKS convergence threshold
    double cpks_delta_;

    // Scaling factor, only used to set the corresponding PsiVar for now
    double exch_scale_alpha_;

    // Do coupled induction ?
    bool coupled_ind_;
    
    // Memory in doubles
    unsigned long int memory_;

    // Energies table
    std::map<std::string, double> energies_;

    // Dimer primary basis set
    std::shared_ptr<BasisSet> primary_;
    // Monomer A primary basis set 
    std::shared_ptr<BasisSet> primary_A_;
    // Monomer B primary basis set
    std::shared_ptr<BasisSet> primary_B_;
    // Dimer -RI or -MP2FIT auxiliary basis set 
    std::shared_ptr<BasisSet> mp2fit_;
    // Dimer -JKFIT auxiliary basis set 
    std::shared_ptr<BasisSet> jkfit_;

    // Dimer SCF energy
    double E_dimer_;
    // Monomer A SCF energy
    double E_monomer_A_;
    // Monomer B SCF energy
    double E_monomer_B_;
    
    // Dimer geometry
    std::shared_ptr<Molecule> dimer_;
    // Monomer A geometry
    std::shared_ptr<Molecule> monomer_A_;
    // Monomer B geometry
    std::shared_ptr<Molecule> monomer_B_;

    // Monomer A C matrix (full occ), alpha spin
    std::shared_ptr<Matrix> Cocca_A_;
    // Monomer B C matrix (full occ), alpha spin
    std::shared_ptr<Matrix> Cocca_B_;
    // Monomer A C matrix (full vir), alpha spin
    std::shared_ptr<Matrix> Cvira_A_;
    // Monomer B C matrix (full vir), alpha spin
    std::shared_ptr<Matrix> Cvira_B_;
    
    // Monomer A C matrix (full occ), beta spin
    std::shared_ptr<Matrix> Coccb_A_;
    // Monomer B C matrix (full occ), beta spin
    std::shared_ptr<Matrix> Coccb_B_;
    // Monomer A C matrix (full vir), beta spin
    std::shared_ptr<Matrix> Cvirb_A_;
    // Monomer B C matrix (full vir), beta spin
    std::shared_ptr<Matrix> Cvirb_B_;

    // Monomer A eps vector (full occ), alpha spin
    std::shared_ptr<Vector> eps_occa_A_;
    // Monomer B eps vector (full occ), alpha spin
    std::shared_ptr<Vector> eps_occa_B_;
    // Monomer A eps vector (full vir), alpha spin
    std::shared_ptr<Vector> eps_vira_A_;
    // Monomer B eps vector (full vir), alpha spin
    std::shared_ptr<Vector> eps_vira_B_;

    // Monomer A eps vector (full occ), beta spin
    std::shared_ptr<Vector> eps_occb_A_;
    // Monomer B eps vector (full occ), beta spin
    std::shared_ptr<Vector> eps_occb_B_;
    // Monomer A eps vector (full vir), beta spin
    std::shared_ptr<Vector> eps_virb_A_;
    // Monomer B eps vector (full vir), beta spin
    std::shared_ptr<Vector> eps_virb_B_;

    // Monomer A C matrix (active occ), alpha spin
    std::shared_ptr<Matrix> Caocca_A_;
    // Monomer B C matrix (active occ), alpha spin
    std::shared_ptr<Matrix> Caocca_B_;
    // Monomer A C matrix (active vir), alpha spin
    std::shared_ptr<Matrix> Cavira_A_;
    // Monomer B C matrix (active vir), alpha spin
    std::shared_ptr<Matrix> Cavira_B_;

    // Monomer A C matrix (active occ), beta spin
    std::shared_ptr<Matrix> Caoccb_A_;
    // Monomer B C matrix (active occ), beta spin
    std::shared_ptr<Matrix> Caoccb_B_;
    // Monomer A C matrix (active vir), beta spin
    std::shared_ptr<Matrix> Cavirb_A_;
    // Monomer B C matrix (active vir), beta spin
    std::shared_ptr<Matrix> Cavirb_B_;

    // Monomer A C matrix (frozen occ), alpha spin
    std::shared_ptr<Matrix> Cfocca_A_;
    // Monomer B C matrix (frozen occ), alpha spin
    std::shared_ptr<Matrix> Cfocca_B_;
    // Monomer A C matrix (frozen vir), alpha spin
    std::shared_ptr<Matrix> Cfvira_A_;
    // Monomer B C matrix (frozen vir), alpha spin
    std::shared_ptr<Matrix> Cfvira_B_;

    // Monomer A C matrix (frozen occ), beta spin
    std::shared_ptr<Matrix> Cfoccb_A_;
    // Monomer B C matrix (frozen occ), beta spin
    std::shared_ptr<Matrix> Cfoccb_B_;
    // Monomer A C matrix (frozen vir), beta spin
    std::shared_ptr<Matrix> Cfvirb_A_;
    // Monomer B C matrix (frozen vir), beta spin
    std::shared_ptr<Matrix> Cfvirb_B_;

    // Monomer A eps vector (active occ), alpha spin
    std::shared_ptr<Vector> eps_aocca_A_;
    // Monomer B eps vector (active occ), alpha spin
    std::shared_ptr<Vector> eps_aocca_B_;
    // Monomer A eps vector (active vir), alpha spin
    std::shared_ptr<Vector> eps_avira_A_;
    // Monomer B eps vector (active vir), alpha spin
    std::shared_ptr<Vector> eps_avira_B_;

    // Monomer A eps vector (active occ), beta spin
    std::shared_ptr<Vector> eps_aoccb_A_;
    // Monomer B eps vector (active occ), beta spin
    std::shared_ptr<Vector> eps_aoccb_B_;
    // Monomer A eps vector (active vir), beta spin
    std::shared_ptr<Vector> eps_avirb_A_;
    // Monomer B eps vector (active vir), beta spin
    std::shared_ptr<Vector> eps_avirb_B_;

    // Monomer A eps vector (frozen occ), alpha spin
    std::shared_ptr<Vector> eps_focca_A_;
    // Monomer B eps vector (frozen occ), alpha spin
    std::shared_ptr<Vector> eps_focca_B_;
    // Monomer A eps vector (frozen vir), alpha spin
    std::shared_ptr<Vector> eps_fvira_A_;
    // Monomer B eps vector (frozen vir), alpha spin
    std::shared_ptr<Vector> eps_fvira_B_;

    // Monomer A eps vector (frozen occ), beta spin
    std::shared_ptr<Vector> eps_foccb_A_;
    // Monomer B eps vector (frozen occ), beta spin
    std::shared_ptr<Vector> eps_foccb_B_;
    // Monomer A eps vector (frozen vir), beta spin
    std::shared_ptr<Vector> eps_fvirb_A_;
    // Monomer B eps vector (frozen vir), beta spin
    std::shared_ptr<Vector> eps_fvirb_B_;

    // Shared matrices (Fock-like)
    std::map<std::string, std::shared_ptr<Matrix> > vars_;


    // Print author/sizing/spec info
    virtual void print_header() const;
    // Obligatory
    virtual void print_trailer();

    // Hartree-Fock-like terms (Elst, Exch, Ind)
    virtual void fock_terms();
    // MP2-like terms (Disp)
    virtual void mp2_terms();

    // => Helper Methods <= //

    // Build the AO-basis dimer overlap matrix
    std::shared_ptr<Matrix> build_S(std::shared_ptr<BasisSet> basis);
    // Build the potential integral matrix
    std::shared_ptr<Matrix> build_V(std::shared_ptr<BasisSet> basis);

    // Build the alpha and beta S_ij matrices in the dimer occupied space
    std::shared_ptr<Matrix> build_Sija(std::shared_ptr<Matrix> S);
    std::shared_ptr<Matrix> build_Sijb(std::shared_ptr<Matrix> S);
    // Build the S^\infty expansion in the dimer occupied space
    std::shared_ptr<Matrix> build_Sij_n(std::shared_ptr<Matrix> Sij);
    // Build the Cbar matrices from S^\infty
    std::map<std::string, std::shared_ptr<Matrix> > build_Cbar(std::shared_ptr<Matrix> Sa, std::shared_ptr<Matrix> Sb);

    // Compute the CPKS solution
    std::map< std::string, std::shared_ptr<Matrix> > compute_x(std::shared_ptr<JK> jk,
                                                                 std::shared_ptr<Matrix> wa_B,
                                                                 std::shared_ptr<Matrix> wb_B,
                                                                 std::shared_ptr<Matrix> wa_A,
                                                                 std::shared_ptr<Matrix> wb_A);

    // Build the ExchInd20 potential in the monomer A ov space
    std::shared_ptr<Matrix> build_exch_ind_pot(std::map<std::string, std::shared_ptr<Matrix> >& vars);
    // Build the Ind20 potential in the monomer A ov space
    std::shared_ptr<Matrix> build_ind_pot(std::map<std::string, std::shared_ptr<Matrix> >& vars);

    void initialize(SharedWavefunction mA, SharedWavefunction mB);

//    // ==> DFTSAPT <==
//
//    // Number of frequency points in Casimir-Poldar
//    int freq_points_;
//    // Frequency scale in Casimir-Poldar
//    double freq_scale_;
//    // Maximum number of terms in Casimir-Poldar susceptibility coupling
//    int freq_max_k_;
//
//    // Grab an uncoupled susceptibility in the RI basis
//    std::shared_ptr<Matrix> uncoupled_susceptibility(
//        double omega,
//        std::shared_ptr<Vector> ea,
//        std::shared_ptr<Vector> er,
//        std::shared_ptr<Tensor> Bar);
//    // Grab a coupled susceptibility in the RI basis
//    std::shared_ptr<Matrix> coupled_susceptibility(
//        double omega,
//        std::shared_ptr<Vector> ea,
//        std::shared_ptr<Vector> er,
//        std::map<std::string, std::shared_ptr<Tensor> >& vars,
//        int nmax);
//    // Grab a coupled susceptibility in the RI basis (N^6)
//    std::shared_ptr<Matrix> coupled_susceptibility_debug(
//        double omega,
//        std::shared_ptr<Vector> ea,
//        std::shared_ptr<Vector> er,
//        std::shared_ptr<Tensor> AaaT,
//        std::shared_ptr<Tensor> AarT,
//        std::shared_ptr<Tensor> ArrT,
//        std::shared_ptr<Tensor> DarT);
//
//    // => Utility Routines <= //
//
//    // Inner product LT' * lambda * RT => Result
//    std::shared_ptr<Matrix> inner(
//        std::shared_ptr<Tensor> LT,
//        std::shared_ptr<Tensor> RT,
//        std::shared_ptr<Matrix> lambda = std::shared_ptr<Matrix>());
//    // Fitting product RT * metric => Result
//    std::shared_ptr<Tensor> fitting(
//        const std::string& name,
//        std::shared_ptr<Tensor> RT,
//        std::shared_ptr<Matrix> metric);
//    // DAXPY, alpha L + beta R => R
//    void axpy(
//        std::shared_ptr<Tensor> LT,
//        std::shared_ptr<Tensor> RT,
//        double alpha = 1.0,
//        double beta = 1.0);

public:
    // Constructor, call this with 3 converged SCF jobs (dimer, monomer A, monomer B)
    USAPT0(SharedWavefunction d, SharedWavefunction mA, SharedWavefunction mB,
           Options& options, std::shared_ptr<PSIO> psio);
    virtual ~USAPT0();

    // Compute the USAPT0 analysis
    virtual double compute_energy();

//    void fd(int nA, int nB, double** CAp, double** Sp, int nso, int no, double** CBp, double** Cp, std::shared_ptr<Matrix> Sa);
};

class CPKS_USAPT0 {

friend class USAPT0;

protected:

    // => Global Data <= //

    // Convergence tolerance
    double delta_;
    // Maximum allowed iterations
    int maxiter_;
    // JK Object
    std::shared_ptr<JK> jk_;

    // => Monomer A Problem <= //

    // Perturbations applied to A
    std::shared_ptr<Matrix> wa_B_;
    std::shared_ptr<Matrix> wb_B_;
    // Response of A
    std::shared_ptr<Matrix> xa_A_;
    std::shared_ptr<Matrix> xb_A_;
    // Active occ orbital coefficients of A
    std::shared_ptr<Matrix> Cocca_A_;
    std::shared_ptr<Matrix> Coccb_A_;
    // Active vir orbital coefficients of A
    std::shared_ptr<Matrix> Cvira_A_;
    std::shared_ptr<Matrix> Cvirb_A_;
    // Active occ orbital eigenvalues of A
    std::shared_ptr<Vector> eps_occa_A_;
    std::shared_ptr<Vector> eps_occb_A_;
    // Active vir orbital eigenvalues of A
    std::shared_ptr<Vector> eps_vira_A_;
    std::shared_ptr<Vector> eps_virb_A_;

    // => Monomer B Problem <= //

    // Perturbations applied to B
    std::shared_ptr<Matrix> wa_A_;
    std::shared_ptr<Matrix> wb_A_;
    // Response of B
    std::shared_ptr<Matrix> xa_B_;
    std::shared_ptr<Matrix> xb_B_;
    // Active occ orbital coefficients of B
    std::shared_ptr<Matrix> Cocca_B_;
    std::shared_ptr<Matrix> Coccb_B_;
    // Active vir orbital coefficients of B
    std::shared_ptr<Matrix> Cvira_B_;
    std::shared_ptr<Matrix> Cvirb_B_;
    // Active occ orbital eigenvalues of B
    std::shared_ptr<Vector> eps_occa_B_;
    std::shared_ptr<Vector> eps_occb_B_;
    // Active vir orbital eigenvalues of B
    std::shared_ptr<Vector> eps_vira_B_;
    std::shared_ptr<Vector> eps_virb_B_;

    // Form the s = Ab product for the provided vectors b (may or may not need more iterations)
    std::map<std::string, std::shared_ptr<Matrix> > product(std::map<std::string, std::shared_ptr<Matrix> >& b);
    // Apply the denominator from r into zs
    void preconditioner(std::shared_ptr<Matrix> r,
                        std::shared_ptr<Matrix> z,
                        std::shared_ptr<Vector> o,
                        std::shared_ptr<Vector> v);

public:
    CPKS_USAPT0();
    virtual ~CPKS_USAPT0();

    void compute_cpks();
};

}} // end namespaces

#endif
