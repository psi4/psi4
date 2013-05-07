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

#ifndef DFTSAPT_H
#define DFTSAPT_H

#include <libmints/typedefs.h>
#include <libmints/wavefunction.h>
#include <map>

namespace psi {

class JK;
class Options;
class Tensor;

namespace dftsapt {

class DFTSAPT {

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

    // Monomer A C matrix (full occ)
    boost::shared_ptr<Matrix> Cocc_A_;
    // Monomer B C matrix (full occ)
    boost::shared_ptr<Matrix> Cocc_B_;
    // Monomer A C matrix (full vir)
    boost::shared_ptr<Matrix> Cvir_A_;
    // Monomer B C matrix (full vir)
    boost::shared_ptr<Matrix> Cvir_B_;
    
    // Monomer A eps vector (full occ)
    boost::shared_ptr<Vector> eps_occ_A_;
    // Monomer B eps vector (full occ)
    boost::shared_ptr<Vector> eps_occ_B_;
    // Monomer A eps vector (full vir)
    boost::shared_ptr<Vector> eps_vir_A_;
    // Monomer B eps vector (full vir)
    boost::shared_ptr<Vector> eps_vir_B_;


    // Monomer A C matrix (active occ)
    boost::shared_ptr<Matrix> Caocc_A_;
    // Monomer B C matrix (active occ)
    boost::shared_ptr<Matrix> Caocc_B_;
    // Monomer A C matrix (active vir)
    boost::shared_ptr<Matrix> Cavir_A_;
    // Monomer B C matrix (active vir)
    boost::shared_ptr<Matrix> Cavir_B_;

    // Monomer A eps vector (active occ)
    boost::shared_ptr<Vector> eps_aocc_A_;
    // Monomer B eps vector (active occ)
    boost::shared_ptr<Vector> eps_aocc_B_;
    // Monomer A eps vector (active vir)
    boost::shared_ptr<Vector> eps_avir_A_;
    // Monomer B eps vector (active vir)
    boost::shared_ptr<Vector> eps_avir_B_;

    // Monomer A eps vector (frozen occ)
    boost::shared_ptr<Vector> eps_focc_A_;
    // Monomer B eps vector (frozen occ)
    boost::shared_ptr<Vector> eps_focc_B_;
    // Monomer A eps vector (frozen vir)
    boost::shared_ptr<Vector> eps_fvir_A_;
    // Monomer B eps vector (frozen vir)
    boost::shared_ptr<Vector> eps_fvir_B_;

    // Shared matrices (Fock-like)
    std::map<std::string, boost::shared_ptr<Matrix> > vars_;

    // Number of frequency points in Casimir-Poldar
    int freq_points_;   
    // Frequency scale in Casimir-Poldar
    double freq_scale_;
    // Maximum number of terms in Casimir-Poldar susceptibility coupling
    int freq_max_k_;
    
    // Print author/sizing/spec info
    virtual void print_header() const;
    // Obligatory
    virtual void print_trailer();

    // Hartree-Fock-like terms (Elst, Exch, Ind)
    virtual void fock_terms();
    // MP2-like terms (Disp)
    virtual void mp2_terms();

    // Hartree-Fock-like terms (Elst, Exch, Ind)
    virtual void local_fock_terms();
    // MP2-like terms (Disp)
    virtual void local_mp2_terms();

    // Build the AO-basis dimer overlap matrix
    boost::shared_ptr<Matrix> build_S(boost::shared_ptr<BasisSet> basis);
    // Build the potential integral matrix
    boost::shared_ptr<Matrix> build_V(boost::shared_ptr<BasisSet> basis);

    // Build the S_ij matrix in the dimer occupied space
    boost::shared_ptr<Matrix> build_Sij(boost::shared_ptr<Matrix> S);
    // Build the S^\infty expansion in the dimer occupied space
    boost::shared_ptr<Matrix> build_Sij_n(boost::shared_ptr<Matrix> Sij);
    // Build the Cbar matrices from S^\infty
    std::map<std::string, boost::shared_ptr<Matrix> > build_Cbar(boost::shared_ptr<Matrix> S);

    // Compute the CPKS solution
    std::pair<boost::shared_ptr<Matrix>, boost::shared_ptr<Matrix> > compute_x(boost::shared_ptr<JK> jk, boost::shared_ptr<Matrix> w_B, boost::shared_ptr<Matrix> w_A);

    // Try out some TDHF Disp2
    void tdhf_demo();
    // Grab an uncoupled susceptibility in the RI basis
    boost::shared_ptr<Matrix> uncoupled_susceptibility(
        double omega, 
        boost::shared_ptr<Vector> ea, 
        boost::shared_ptr<Vector> er, 
        boost::shared_ptr<Tensor> Bar);
    // Grab a coupled susceptibility in the RI basis
    boost::shared_ptr<Matrix> coupled_susceptibility(
        double omega, 
        boost::shared_ptr<Vector> ea, 
        boost::shared_ptr<Vector> er, 
        std::map<std::string, boost::shared_ptr<Tensor> >& vars,
        int nmax);
    // Grab a coupled susceptibility in the RI basis (N^6)
    boost::shared_ptr<Matrix> coupled_susceptibility_debug(
        double omega, 
        boost::shared_ptr<Vector> ea, 
        boost::shared_ptr<Vector> er, 
        boost::shared_ptr<Tensor> AaaT,
        boost::shared_ptr<Tensor> AarT,
        boost::shared_ptr<Tensor> ArrT,
        boost::shared_ptr<Tensor> DarT);

    // => Utility Routines <= //

    // Inner product LT' * lambda * RT => Result
    boost::shared_ptr<Matrix> inner(
        boost::shared_ptr<Tensor> LT, 
        boost::shared_ptr<Tensor> RT, 
        boost::shared_ptr<Matrix> lambda = boost::shared_ptr<Matrix>());
    // Fitting product RT * metric => Result
    boost::shared_ptr<Tensor> fitting(
        const std::string& name, 
        boost::shared_ptr<Tensor> RT, 
        boost::shared_ptr<Matrix> metric);
    // DAXPY, alpha L + beta R => R
    void axpy(
        boost::shared_ptr<Tensor> LT, 
        boost::shared_ptr<Tensor> RT, 
        double alpha = 1.0,
        double beta = 1.0);

    // Double GEMM
    boost::shared_ptr<Matrix> doublet(boost::shared_ptr<Matrix> A, boost::shared_ptr<Matrix> B, bool tA = false, bool tB = false);
    // Triple GEMM
    boost::shared_ptr<Matrix> triplet(boost::shared_ptr<Matrix> A, boost::shared_ptr<Matrix> B, boost::shared_ptr<Matrix> C, bool tA = false, bool tB = false, bool tC = false);

    // Protected constructor (use factory below)
    DFTSAPT();
    // Common initialization
    void common_init();
public:
    // Destructor, frees memory
    virtual ~DFTSAPT();

    // Factory constructor, call this with 3 converged SCF jobs (dimer, monomer A, monomer B)
    static boost::shared_ptr<DFTSAPT> build(boost::shared_ptr<Wavefunction> d,
                                            boost::shared_ptr<Wavefunction> mB,
                                            boost::shared_ptr<Wavefunction> mA);

    // Compute the DFT-SAPT analysis
    virtual double compute_energy();

};

class CPKS_SAPT {

friend class DFTSAPT;

protected:

    // => Global Data <= //

    // Convergence tolerance
    double delta_;
    // Maximum allowed iterations
    int maxiter_;
    // JK Object 
    boost::shared_ptr<JK> jk_;
    
    // => Monomer A Problem <= //

    // Perturbation applied to A
    boost::shared_ptr<Matrix> w_A_;
    // Response of A
    boost::shared_ptr<Matrix> x_A_;
    // Active occ orbital coefficients of A
    boost::shared_ptr<Matrix> Cocc_A_;
    // Active vir orbital coefficients of A
    boost::shared_ptr<Matrix> Cvir_A_;
    // Active occ orbital eigenvalues of A
    boost::shared_ptr<Vector> eps_occ_A_;
    // Active vir orbital eigenvalues of A
    boost::shared_ptr<Vector> eps_vir_A_;

    // => Monomer B Problem <= //

    // Perturbation applied to B
    boost::shared_ptr<Matrix> w_B_;
    // Response of B
    boost::shared_ptr<Matrix> x_B_;
    // Active occ orbital coefficients of B
    boost::shared_ptr<Matrix> Cocc_B_;
    // Active vir orbital coefficients of B
    boost::shared_ptr<Matrix> Cvir_B_;
    // Active occ orbital eigenvalues of B
    boost::shared_ptr<Vector> eps_occ_B_;
    // Active vir orbital eigenvalues of B
    boost::shared_ptr<Vector> eps_vir_B_;

    // Form the s = Ab product for the provided vectors b (may or may not need more iterations)
    std::map<std::string, boost::shared_ptr<Matrix> > product(std::map<std::string, boost::shared_ptr<Matrix> > b);
    // Apply the denominator from r into z
    void preconditioner(boost::shared_ptr<Matrix> r,
                        boost::shared_ptr<Matrix> z,
                        boost::shared_ptr<Vector> o,
                        boost::shared_ptr<Vector> v);

public:
    CPKS_SAPT();
    virtual ~CPKS_SAPT();

    void compute_cpks();
};

class GaussChebyshev {

protected:
    int npoint_;
    double scale_;

    std::vector<double> nodes_;
    std::vector<double> weights_;

public:
    GaussChebyshev(int npoint, double scale);
    ~GaussChebyshev();
    
    void print_header();
    void compute();

    std::vector<double>& nodes() { return nodes_; }
    std::vector<double>& weights() { return weights_; }

};

}} // End namespace

#endif

