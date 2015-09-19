/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

// usapt.h, to be included in dftsapt.h if needed.

#ifndef USAPT_H
#define USAPT_H

#include <libmints/typedefs.h>
#include <libmints/wavefunction.h>
#include <map>

namespace psi {

class JK;
class Options;
class Tensor;
class AtomicDensity;

namespace sapt {

/*- Open-shell generalization of SAPT0
    Will be the basis for DFTSAPT
    Will need to create better class hierarchy -*/

class USAPT0 {

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

    // Do coupled induction ?
    bool coupled_ind_;
    
    // Memory in doubles
    unsigned long int memory_;

    // Energies table
    std::map<std::string, double> energies_;

    // Dimer primary basis set
    boost::shared_ptr<BasisSet> primary_;
    // Monomer A primary basis set 
    boost::shared_ptr<BasisSet> primary_A_;
    // Monomer B primary basis set
    boost::shared_ptr<BasisSet> primary_B_;
    // Dimer -RI or -MP2FIT auxiliary basis set 
    boost::shared_ptr<BasisSet> mp2fit_;

    // Dimer SCF energy
    double E_dimer_;
    // Monomer A SCF energy
    double E_monomer_A_;
    // Monomer B SCF energy
    double E_monomer_B_;
    
    // Dimer geometry
    boost::shared_ptr<Molecule> dimer_;
    // Monomer A geometry
    boost::shared_ptr<Molecule> monomer_A_;
    // Monomer B geometry
    boost::shared_ptr<Molecule> monomer_B_;

    // Monomer A C matrix (full occ), alpha spin
    boost::shared_ptr<Matrix> Cocca_A_;
    // Monomer B C matrix (full occ), alpha spin
    boost::shared_ptr<Matrix> Cocca_B_;
    // Monomer A C matrix (full vir), alpha spin
    boost::shared_ptr<Matrix> Cvira_A_;
    // Monomer B C matrix (full vir), alpha spin
    boost::shared_ptr<Matrix> Cvira_B_;
    
    // Monomer A C matrix (full occ), beta spin
    boost::shared_ptr<Matrix> Coccb_A_;
    // Monomer B C matrix (full occ), beta spin
    boost::shared_ptr<Matrix> Coccb_B_;
    // Monomer A C matrix (full vir), beta spin
    boost::shared_ptr<Matrix> Cvirb_A_;
    // Monomer B C matrix (full vir), beta spin
    boost::shared_ptr<Matrix> Cvirb_B_;

    // Monomer A eps vector (full occ), alpha spin
    boost::shared_ptr<Vector> eps_occa_A_;
    // Monomer B eps vector (full occ), alpha spin
    boost::shared_ptr<Vector> eps_occa_B_;
    // Monomer A eps vector (full vir), alpha spin
    boost::shared_ptr<Vector> eps_vira_A_;
    // Monomer B eps vector (full vir), alpha spin
    boost::shared_ptr<Vector> eps_vira_B_;

    // Monomer A eps vector (full occ), beta spin
    boost::shared_ptr<Vector> eps_occb_A_;
    // Monomer B eps vector (full occ), beta spin
    boost::shared_ptr<Vector> eps_occb_B_;
    // Monomer A eps vector (full vir), beta spin
    boost::shared_ptr<Vector> eps_virb_A_;
    // Monomer B eps vector (full vir), beta spin
    boost::shared_ptr<Vector> eps_virb_B_;

    // Monomer A C matrix (active occ), alpha spin
    boost::shared_ptr<Matrix> Caocca_A_;
    // Monomer B C matrix (active occ), alpha spin
    boost::shared_ptr<Matrix> Caocca_B_;
    // Monomer A C matrix (active vir), alpha spin
    boost::shared_ptr<Matrix> Cavira_A_;
    // Monomer B C matrix (active vir), alpha spin
    boost::shared_ptr<Matrix> Cavira_B_;

    // Monomer A C matrix (active occ), beta spin
    boost::shared_ptr<Matrix> Caoccb_A_;
    // Monomer B C matrix (active occ), beta spin
    boost::shared_ptr<Matrix> Caoccb_B_;
    // Monomer A C matrix (active vir), beta spin
    boost::shared_ptr<Matrix> Cavirb_A_;
    // Monomer B C matrix (active vir), beta spin
    boost::shared_ptr<Matrix> Cavirb_B_;

    // Monomer A C matrix (frozen occ), alpha spin
    boost::shared_ptr<Matrix> Cfocca_A_;
    // Monomer B C matrix (frozen occ), alpha spin
    boost::shared_ptr<Matrix> Cfocca_B_;
    // Monomer A C matrix (frozen vir), alpha spin
    boost::shared_ptr<Matrix> Cfvira_A_;
    // Monomer B C matrix (frozen vir), alpha spin
    boost::shared_ptr<Matrix> Cfvira_B_;

    // Monomer A C matrix (frozen occ), beta spin
    boost::shared_ptr<Matrix> Cfoccb_A_;
    // Monomer B C matrix (frozen occ), beta spin
    boost::shared_ptr<Matrix> Cfoccb_B_;
    // Monomer A C matrix (frozen vir), beta spin
    boost::shared_ptr<Matrix> Cfvirb_A_;
    // Monomer B C matrix (frozen vir), beta spin
    boost::shared_ptr<Matrix> Cfvirb_B_;

    // Monomer A eps vector (active occ), alpha spin
    boost::shared_ptr<Vector> eps_aocca_A_;
    // Monomer B eps vector (active occ), alpha spin
    boost::shared_ptr<Vector> eps_aocca_B_;
    // Monomer A eps vector (active vir), alpha spin
    boost::shared_ptr<Vector> eps_avira_A_;
    // Monomer B eps vector (active vir), alpha spin
    boost::shared_ptr<Vector> eps_avira_B_;

    // Monomer A eps vector (active occ), beta spin
    boost::shared_ptr<Vector> eps_aoccb_A_;
    // Monomer B eps vector (active occ), beta spin
    boost::shared_ptr<Vector> eps_aoccb_B_;
    // Monomer A eps vector (active vir), beta spin
    boost::shared_ptr<Vector> eps_avirb_A_;
    // Monomer B eps vector (active vir), beta spin
    boost::shared_ptr<Vector> eps_avirb_B_;

    // Monomer A eps vector (frozen occ), alpha spin
    boost::shared_ptr<Vector> eps_focca_A_;
    // Monomer B eps vector (frozen occ), alpha spin
    boost::shared_ptr<Vector> eps_focca_B_;
    // Monomer A eps vector (frozen vir), alpha spin
    boost::shared_ptr<Vector> eps_fvira_A_;
    // Monomer B eps vector (frozen vir), alpha spin
    boost::shared_ptr<Vector> eps_fvira_B_;

    // Monomer A eps vector (frozen occ), beta spin
    boost::shared_ptr<Vector> eps_foccb_A_;
    // Monomer B eps vector (frozen occ), beta spin
    boost::shared_ptr<Vector> eps_foccb_B_;
    // Monomer A eps vector (frozen vir), beta spin
    boost::shared_ptr<Vector> eps_fvirb_A_;
    // Monomer B eps vector (frozen vir), beta spin
    boost::shared_ptr<Vector> eps_fvirb_B_;

    // Shared matrices (Fock-like)
    std::map<std::string, boost::shared_ptr<Matrix> > vars_;


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
    boost::shared_ptr<Matrix> build_S(boost::shared_ptr<BasisSet> basis);
    // Build the potential integral matrix
    boost::shared_ptr<Matrix> build_V(boost::shared_ptr<BasisSet> basis);

    // Build the alpha and beta S_ij matrices in the dimer occupied space
    boost::shared_ptr<Matrix> build_Sija(boost::shared_ptr<Matrix> S);
    boost::shared_ptr<Matrix> build_Sijb(boost::shared_ptr<Matrix> S);
    // Build the S^\infty expansion in the dimer occupied space
    boost::shared_ptr<Matrix> build_Sij_n(boost::shared_ptr<Matrix> Sij);
    // Build the Cbar matrices from S^\infty
    std::map<std::string, boost::shared_ptr<Matrix> > build_Cbar(boost::shared_ptr<Matrix> Sa, boost::shared_ptr<Matrix> Sb);

    // Compute the CPKS solution
    std::map< std::string, boost::shared_ptr<Matrix> > compute_x(boost::shared_ptr<JK> jk,
                                                                 boost::shared_ptr<Matrix> wa_B,
                                                                 boost::shared_ptr<Matrix> wb_B,
                                                                 boost::shared_ptr<Matrix> wa_A,
                                                                 boost::shared_ptr<Matrix> wb_A);

    // Build the ExchInd20 potential in the monomer A ov space
    boost::shared_ptr<Matrix> build_exch_ind_pot(std::map<std::string, boost::shared_ptr<Matrix> >& vars);
    // Build the Ind20 potential in the monomer A ov space
    boost::shared_ptr<Matrix> build_ind_pot(std::map<std::string, boost::shared_ptr<Matrix> >& vars);

    void common_init();

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
//    boost::shared_ptr<Matrix> uncoupled_susceptibility(
//        double omega,
//        boost::shared_ptr<Vector> ea,
//        boost::shared_ptr<Vector> er,
//        boost::shared_ptr<Tensor> Bar);
//    // Grab a coupled susceptibility in the RI basis
//    boost::shared_ptr<Matrix> coupled_susceptibility(
//        double omega,
//        boost::shared_ptr<Vector> ea,
//        boost::shared_ptr<Vector> er,
//        std::map<std::string, boost::shared_ptr<Tensor> >& vars,
//        int nmax);
//    // Grab a coupled susceptibility in the RI basis (N^6)
//    boost::shared_ptr<Matrix> coupled_susceptibility_debug(
//        double omega,
//        boost::shared_ptr<Vector> ea,
//        boost::shared_ptr<Vector> er,
//        boost::shared_ptr<Tensor> AaaT,
//        boost::shared_ptr<Tensor> AarT,
//        boost::shared_ptr<Tensor> ArrT,
//        boost::shared_ptr<Tensor> DarT);
//
//    // => Utility Routines <= //
//
//    // Inner product LT' * lambda * RT => Result
//    boost::shared_ptr<Matrix> inner(
//        boost::shared_ptr<Tensor> LT,
//        boost::shared_ptr<Tensor> RT,
//        boost::shared_ptr<Matrix> lambda = boost::shared_ptr<Matrix>());
//    // Fitting product RT * metric => Result
//    boost::shared_ptr<Tensor> fitting(
//        const std::string& name,
//        boost::shared_ptr<Tensor> RT,
//        boost::shared_ptr<Matrix> metric);
//    // DAXPY, alpha L + beta R => R
//    void axpy(
//        boost::shared_ptr<Tensor> LT,
//        boost::shared_ptr<Tensor> RT,
//        double alpha = 1.0,
//        double beta = 1.0);

public:
    USAPT0();
    virtual ~USAPT0();

    // Factory constructor, call this with 3 converged SCF jobs (dimer, monomer A, monomer B)
    static boost::shared_ptr<USAPT0> build(boost::shared_ptr<Wavefunction> d,
                                            boost::shared_ptr<Wavefunction> mA,
                                            boost::shared_ptr<Wavefunction> mB);

    // Compute the USAPT0 analysis
    virtual double compute_energy();

//    void fd(int nA, int nB, double** CAp, double** Sp, int nso, int no, double** CBp, double** Cp, boost::shared_ptr<Matrix> Sa);
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
    boost::shared_ptr<JK> jk_;

    // => Monomer A Problem <= //

    // Perturbations applied to A
    boost::shared_ptr<Matrix> wa_B_;
    boost::shared_ptr<Matrix> wb_B_;
    // Response of A
    boost::shared_ptr<Matrix> xa_A_;
    boost::shared_ptr<Matrix> xb_A_;
    // Active occ orbital coefficients of A
    boost::shared_ptr<Matrix> Cocca_A_;
    boost::shared_ptr<Matrix> Coccb_A_;
    // Active vir orbital coefficients of A
    boost::shared_ptr<Matrix> Cvira_A_;
    boost::shared_ptr<Matrix> Cvirb_A_;
    // Active occ orbital eigenvalues of A
    boost::shared_ptr<Vector> eps_occa_A_;
    boost::shared_ptr<Vector> eps_occb_A_;
    // Active vir orbital eigenvalues of A
    boost::shared_ptr<Vector> eps_vira_A_;
    boost::shared_ptr<Vector> eps_virb_A_;

    // => Monomer B Problem <= //

    // Perturbations applied to B
    boost::shared_ptr<Matrix> wa_A_;
    boost::shared_ptr<Matrix> wb_A_;
    // Response of B
    boost::shared_ptr<Matrix> xa_B_;
    boost::shared_ptr<Matrix> xb_B_;
    // Active occ orbital coefficients of B
    boost::shared_ptr<Matrix> Cocca_B_;
    boost::shared_ptr<Matrix> Coccb_B_;
    // Active vir orbital coefficients of B
    boost::shared_ptr<Matrix> Cvira_B_;
    boost::shared_ptr<Matrix> Cvirb_B_;
    // Active occ orbital eigenvalues of B
    boost::shared_ptr<Vector> eps_occa_B_;
    boost::shared_ptr<Vector> eps_occb_B_;
    // Active vir orbital eigenvalues of B
    boost::shared_ptr<Vector> eps_vira_B_;
    boost::shared_ptr<Vector> eps_virb_B_;

    // Form the s = Ab product for the provided vectors b (may or may not need more iterations)
    std::map<std::string, boost::shared_ptr<Matrix> > product(std::map<std::string, boost::shared_ptr<Matrix> >& b);
    // Apply the denominator from r into zs
    void preconditioner(boost::shared_ptr<Matrix> r,
                        boost::shared_ptr<Matrix> z,
                        boost::shared_ptr<Vector> o,
                        boost::shared_ptr<Vector> v);

public:
    CPKS_USAPT0();
    virtual ~CPKS_USAPT0();

    void compute_cpks();
};

}} // end namespaces

#endif
