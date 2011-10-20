#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <utility>
#include <string>
#include <cstring>

#include <psifiles.h>
#include <physconst.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>

#include <libmints/mints.h>
#include <libfunctional/superfunctional.h>
#include <libscf_solver/ks.h>
#include <libscf_solver/dft.h>
#include <libscf_solver/integralfunctors.h>
#include <libscf_solver/omegafunctors.h>

#include "omega.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace psi;
using namespace psi::functional;
using namespace boost;

namespace psi{ namespace scf {

OmegaKS::OmegaKS(Options & options, boost::shared_ptr<PSIO> psio) :
    options_(options), psio_(psio)
{
    common_init();
}
OmegaKS::~OmegaKS()
{
}
void OmegaKS::common_init()
{
    memory_ = Process::environment.get_memory();
    nthread_ = 1;
    #ifdef _OPENMP
        nthread_ = omp_get_max_threads();
    #endif

    print_ = options_.get_int("PRINT");
    print_ = (print_ ? print_ : 1);

    debug_ = options_.get_int("DEBUG");

    // Read in energy convergence threshold
    int thresh = options_.get_int("E_CONVERGE");
    energy_threshold_ = pow(10.0, (double)-thresh);

    // Read in density convergence threshold
    thresh = options_.get_int("D_CONVERGE");
    density_threshold_ = pow(10.0, (double)-thresh);

    block_size_ = options_.get_int("DFT_MAX_POINTS");
    
    //Build the superfunctional
    functional_ = SuperFunctional::createSuperFunctional(options_.get_str("DFT_FUNCTIONAL"),block_size_,1);

    if (!functional_->isRangeCorrected())
        throw PSIEXCEPTION("No omega in this functional");
}
void OmegaKS::run_procedure()
{
    // What are we doing?
    print_header();

    // Pick a starting omega
    guess_omega();

    // Form H
    form_H();

    // Form X according to the usual rules (always canonical herein)
    form_X();

    // Build the DF Integrals objects
    form_DF();

    // Build the KS Integrators
    form_KS();

    // Populate the necessary N, N+1, N-1, etc wavefunctions 
    // Provide hooks to DF and V objects, as well as X, H, C, D, etc
    populate();

    // Burn-in the wavefunctions
    fprintf(outfile, "  ==> Burn-In Procedure <==\n\n");
    fprintf(outfile, "  *** Step %4d: Burn-In Step ***\n\n", 1);

    for (std::map<std::string, boost::shared_ptr<OmegaWavefunction> >::iterator it = wfns_.begin(); it != wfns_.end(); it++) {

        std::string name = (*it).first;
        boost::shared_ptr<OmegaWavefunction> wfn = (*it).second;
      
        // Reset DIIS and whatever 
        wfn->reset();

        fprintf(outfile, "   => %s SCF <=\n\n", name.c_str());
        fprintf(outfile, "   %4s %20s %12s %12s %s\n", "Iter", "Energy", "Delta E", "Delta D", "Stat");
        fflush(outfile);
 
        bool converged = false;
        for (int iter = 0; iter < options_.get_int("MAXITER") || iter < 2; iter++) {

            // Take an SCF step (including any DIIS)
            std::string status = wfn->step();     

            double dE = wfn->deltaE();  
            double dD = wfn->deltaD();
            double E = wfn->E();

            fprintf(outfile, "   %4d %20.14f %12.5E %12.5E %s\n", iter + 1, E, dE, dD, status.c_str()); fflush(outfile);

            if (fabs(dE) < energy_threshold_ && fabs(dD) < density_threshold_) {
                converged = true;
                break;
            }
        }

        wfn->clear();
    
        if (!converged)
            throw PSIEXCEPTION("OmegaKS: Burn-In iterations did not converge");

        fprintf(outfile, "\n"); fflush(outfile);
    }  

    // Omega procedure 
    fprintf(outfile, "  ==> Omega Optimization Procedure <==\n\n"); 

    bool omega_converged = false;
    for (int omega_iter = 1; omega_iter <= options_.get_int("OMEGA_MAXITER"); omega_iter++) {

        omega_step(omega_iter);

        for (std::map<std::string, boost::shared_ptr<OmegaWavefunction> >::iterator it = wfns_.begin(); it != wfns_.end(); it++) {

            std::string name = (*it).first;
            boost::shared_ptr<OmegaWavefunction> wfn = (*it).second;
          
            // Reset DIIS and whatever 
            wfn->reset();

            fprintf(outfile, "   => %s SCF <=\n\n", name.c_str());
            fprintf(outfile, "   %4s %20s %12s %12s %s\n", "Iter", "Energy", "Delta E", "Delta D", "Stat");
            fflush(outfile);
 
            bool converged = false;
            for (int iter = 0; iter < options_.get_int("MAXITER") || iter < 2; iter++) {

                // Take an SCF step (including any DIIS)
                std::string status = wfn->step();     

                double dE = wfn->deltaE();  
                double dD = wfn->deltaD();
                double E = wfn->E();

                fprintf(outfile, "   %4d %20.14f %12.5E %12.5E %s\n", iter + 1, E, dE, dD, status.c_str()); fflush(outfile);

                if (fabs(dE) < energy_threshold_ && fabs(dD) < density_threshold_) {
                    converged = true;
                    break;
                }
            }

            wfn->clear();
        
            if (!converged)
                throw PSIEXCEPTION("OmegaKS: Burn-In iterations did not converge");

            fprintf(outfile, "\n"); fflush(outfile);
        }  
         
        omega_converged = is_omega_converged();
        if (omega_converged) 
            break;
    }

    if (!omega_converged)
        throw PSIEXCEPTION("OmegaKS: Overall omega optimization did not converge");

    // Final  
    fprintf(outfile, "  ==> Final Convergence Procedure <==\n\n"); 

    for (std::map<std::string, boost::shared_ptr<OmegaWavefunction> >::iterator it = wfns_.begin(); it != wfns_.end(); it++) {

        std::string name = (*it).first;
        boost::shared_ptr<OmegaWavefunction> wfn = (*it).second;
      
        // Reset DIIS and whatever 
        wfn->reset();

        fprintf(outfile, "   => %s SCF <=\n\n", name.c_str());
        fprintf(outfile, "   %4s %20s %12s %12s %s\n", "Iter", "Energy", "Delta E", "Delta D", "Stat");
        fflush(outfile);
 
        bool converged = false;
        for (int iter = 0; iter < options_.get_int("MAXITER") || iter < 2; iter++) {

            // Take an SCF step (including any DIIS)
            std::string status = wfn->step();     

            double dE = wfn->deltaE();  
            double dD = wfn->deltaD();
            double E = wfn->E();

            fprintf(outfile, "   %4d %20.14f %12.5E %12.5E %s\n", iter + 1, E, dE, dD, status.c_str()); fflush(outfile);

            if (fabs(dE) < energy_threshold_ && fabs(dD) < density_threshold_) {
                converged = true;
                break;
            }
        }

        wfn->clear();
  
        fprintf(outfile, "\n   UKS %s Orbital Energies [a.u.]:\n\n", name.c_str()); 
        wfn->print_orbitals();

        fprintf(outfile, "  @UKS %s Final energy: %24.16f\n\n", name.c_str(), wfn->E());
        fflush(outfile);
 
        if (!converged)
            throw PSIEXCEPTION("OmegaKS: Final convergence iterations did not converge");

    } 

    // Printing and finalization
    finalize();
}
boost::shared_ptr<Matrix> OmegaKS::build_S(boost::shared_ptr<BasisSet> primary)
{
    int nso = primary->nbf();
    boost::shared_ptr<Matrix> S(new Matrix("S",nso,nso));

    boost::shared_ptr<IntegralFactory> factory(new IntegralFactory(primary,primary,primary,primary));
    boost::shared_ptr<OneBodyAOInt> Sint(factory->ao_overlap());

    Sint->compute(S);

    return S;
}
boost::shared_ptr<Matrix> OmegaKS::build_X(boost::shared_ptr<BasisSet> primary, double min_S)
{
    int nso = primary->nbf();
    boost::shared_ptr<Matrix> S(new Matrix("S",nso,nso));

    boost::shared_ptr<IntegralFactory> factory(new IntegralFactory(primary,primary,primary,primary));
    boost::shared_ptr<OneBodyAOInt> Sint(factory->ao_overlap());

    Sint->compute(S);

    boost::shared_ptr<Matrix> U(new Matrix("U",nso,nso));
    boost::shared_ptr<Vector> s(new Vector("s",nso));
   
    S->diagonalize(U,s);

    double* sp = s->pointer(); 
    int nmo = 0;
    for (int i = 0; i < nso; i++) {
        if (sp[i] > min_S)
            nmo++;
    }


    boost::shared_ptr<Matrix> X(new Matrix("X", nso,nmo));
    double** Xp = X->pointer();
    double** Up = U->pointer();

    int j = 0;
    for (int i = 0; i < nso; i++) {
        if (sp[i] > min_S) {
            C_DAXPY(nso,pow(sp[i], -1.0/2.0), &Up[0][i], nso, &Xp[0][j], nmo);
            j++;
        }
    }
    return X;
}
boost::shared_ptr<Matrix> OmegaKS::build_H(boost::shared_ptr<BasisSet> primary)
{
    int nso = primary->nbf();
    
    boost::shared_ptr<Matrix> H(new Matrix("H", nso, nso));
    boost::shared_ptr<Matrix> T(new Matrix("T", nso, nso));
    boost::shared_ptr<Matrix> V(new Matrix("V", nso, nso));

    // Integral factory
    boost::shared_ptr<IntegralFactory> integral(new IntegralFactory(primary, primary, primary, primary));
    boost::shared_ptr<OneBodyAOInt>    soT(integral->ao_kinetic());
    boost::shared_ptr<OneBodyAOInt>    soV(integral->ao_potential());

    soT->compute(T);
    soV->compute(V);

    H->copy(T);
    H->add(V);

    return H;
}
OmegaIPKS::OmegaIPKS(Options& options, boost::shared_ptr<PSIO> psio) :
    OmegaKS(options, psio)
{
    common_init();
}
OmegaIPKS::~OmegaIPKS()
{
}
void OmegaIPKS::common_init()
{
    reference_ = Process::environment.reference_wavefunction();
    factory_ = reference_->matrix_factory();

    if (factory_->nirrep() != 1)
        throw PSIEXCEPTION("You want symmetry? Find another coder.");

    molecule_ = reference_->molecule();
    basisset_ = reference_->basisset();

    // If the user doesn't spec a basis name, pick it yourself
    // TODO: Verify that the basis assign does not messs this up
    if (options_.get_str("RI_BASIS_SCF") == "") {
        basisset_->molecule()->set_basis_all_atoms(options_.get_str("BASIS") + "-JKFIT", "RI_BASIS_SCF");
        fprintf(outfile, "  No auxiliary basis selected, defaulting to %s-JKFIT\n\n", options_.get_str("BASIS").c_str()); 
    }

    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    auxiliary_ = BasisSet::construct(parser, basisset_->molecule(), "RI_BASIS_SCF");
    parser.reset();

    if (options_.get_bool("OMEGA_GUESS_INTERPOLATE")) {
        int nso = basisset_->nbf();

        Fa_l_N_ = boost::shared_ptr<Matrix>(new Matrix("F Interpolator", nso, nso));
        Fa_l_M_ = boost::shared_ptr<Matrix>(new Matrix("F Interpolator", nso, nso));
        Fb_l_N_ = boost::shared_ptr<Matrix>(new Matrix("F Interpolator", nso, nso));
        Fb_l_M_ = boost::shared_ptr<Matrix>(new Matrix("F Interpolator", nso, nso));
        Fa_r_N_ = boost::shared_ptr<Matrix>(new Matrix("F Interpolator", nso, nso));
        Fa_r_M_ = boost::shared_ptr<Matrix>(new Matrix("F Interpolator", nso, nso));
        Fb_r_N_ = boost::shared_ptr<Matrix>(new Matrix("F Interpolator", nso, nso));
        Fb_r_M_ = boost::shared_ptr<Matrix>(new Matrix("F Interpolator", nso, nso));

    }

}
void OmegaIPKS::print_header()
{
    fprintf(outfile, "\n");
    fprintf(outfile, "         ---------------------------------------------------------\n");
    fprintf(outfile, "                          Optimized-Omega RC-KS\n");
    fprintf(outfile, "                            IP-Only Algorithm                \n");
    fprintf(outfile, "                               Rob Parrish\n"); 
    fprintf(outfile, "                      %3d Threads, %6ld MiB Core\n", nthread_, memory_ / 1000000L);
    fprintf(outfile, "         ---------------------------------------------------------\n\n");
    
    fprintf(outfile, " ==> Functional <==\n\n");
    fprintf(outfile, "  %s\n\n", functional_->getName().c_str());
    fprintf(outfile, " ==> Geometry <==\n\n");
    molecule_->print();
    fprintf(outfile, "  Nuclear repulsion = %20.15f\n\n", basisset_->molecule()->nuclear_repulsion_energy());
    fprintf(outfile, "  Nalpha = %d\n", reference_->nalphapi()[0]);
    fprintf(outfile, "  Nbeta  = %d\n\n", reference_->nbetapi()[0]);
    fprintf(outfile, "  ==> Primary Basis <==\n\n");
    basisset_->print_by_level(outfile, print_);
    fprintf(outfile, "  ==> Primary Basis <==\n\n");
    auxiliary_->print_by_level(outfile, print_);
    fflush(outfile);

}
void OmegaIPKS::guess_omega() 
{

    if (options_.get_str("OMEGA_GUESS") == "HOMO_SIZE") {

        double alpha = options_.get_double("OMEGA_GUESS_A");
        double beta  = options_.get_double("OMEGA_GUESS_B");

        // HOMO extent
        double Rext = 0.0;
        double x,y,z,xx,yy,zz;
        int nalpha = reference_->nalphapi()[0];
        int nbeta = reference_->nbetapi()[0];
        int nmo = reference_->nmo();
        int nso = basisset_->nbf();

        // Make integral generators
        boost::shared_ptr<IntegralFactory> fact(new IntegralFactory(basisset_,basisset_,basisset_,basisset_));
        boost::shared_ptr<OneBodyAOInt> dipole(fact->ao_dipole());
        boost::shared_ptr<OneBodyAOInt> quadrupole(fact->ao_quadrupole());

        // Find the HOMO
        int HOMO_index;
        double** Cp;

        if (nbeta > 0) {
            double e_a = reference_->epsilon_a()->get(0, nalpha - 1);  
            double e_b = reference_->epsilon_b()->get(0, nbeta - 1);  
            if (e_a > e_b) {
                HOMO_index = nalpha - 1;
                Cp = reference_->Ca()->pointer();
            } else {
                HOMO_index = nbeta - 1;
                Cp = reference_->Cb()->pointer();
            }
        } else {
            HOMO_index = nalpha - 1;
            Cp = reference_->Ca()->pointer();
        }

        // Build a temp for the half-transformed array
        double* T = new double[nso];

        // Get dipole integrals
        std::vector<boost::shared_ptr<Matrix> > dipole_ints;
        dipole_ints.push_back(boost::shared_ptr<Matrix>(new Matrix("Dipole X", nso, nso))); 
        dipole_ints.push_back(boost::shared_ptr<Matrix>(new Matrix("Dipole Y", nso, nso))); 
        dipole_ints.push_back(boost::shared_ptr<Matrix>(new Matrix("Dipole Z", nso, nso))); 
        dipole->compute(dipole_ints);

        if (debug_ > 2) {
            for (int i = 0; i < dipole_ints.size(); i++) {
                dipole_ints[i]->print();
            }
        }

        double** Xp = dipole_ints[0]->pointer();
        double** Yp = dipole_ints[1]->pointer();
        double** Zp = dipole_ints[2]->pointer();

        C_DGEMV('N', nso, nso, 1.0, Xp[0], nso, &Cp[0][HOMO_index], nmo, 0.0, T, 1);
        x = -C_DDOT(nso, T, 1, &Cp[0][HOMO_index], nmo);
    
        C_DGEMV('N', nso, nso, 1.0, Yp[0], nso, &Cp[0][HOMO_index], nmo, 0.0, T, 1);
        y = -C_DDOT(nso, T, 1, &Cp[0][HOMO_index], nmo);
    
        C_DGEMV('N', nso, nso, 1.0, Zp[0], nso, &Cp[0][HOMO_index], nmo, 0.0, T, 1);
        z = -C_DDOT(nso, T, 1, &Cp[0][HOMO_index], nmo);

        dipole_ints.clear();
        
        std::vector<boost::shared_ptr<Matrix> > quadrupole_ints;
        quadrupole_ints.push_back(boost::shared_ptr<Matrix>(new Matrix("Quadrupole XX", nso, nso))); 
        quadrupole_ints.push_back(boost::shared_ptr<Matrix>(new Matrix("Quadrupole XY", nso, nso))); 
        quadrupole_ints.push_back(boost::shared_ptr<Matrix>(new Matrix("Quadrupole XZ", nso, nso))); 
        quadrupole_ints.push_back(boost::shared_ptr<Matrix>(new Matrix("Quadrupole YY", nso, nso))); 
        quadrupole_ints.push_back(boost::shared_ptr<Matrix>(new Matrix("Quadrupole YZ", nso, nso))); 
        quadrupole_ints.push_back(boost::shared_ptr<Matrix>(new Matrix("Quadrupole ZZ", nso, nso))); 
        quadrupole->compute(quadrupole_ints);

        if (debug_ > 2) {
            for (int i = 0; i < quadrupole_ints.size(); i++) {
                quadrupole_ints[i]->print();
            }
        }

        double** XXp = quadrupole_ints[0]->pointer();
        double** YYp = quadrupole_ints[3]->pointer();
        double** ZZp = quadrupole_ints[5]->pointer();

        C_DGEMV('N', nso, nso, 1.0, XXp[0], nso, &Cp[0][HOMO_index], nmo, 0.0, T, 1);
        xx = -C_DDOT(nso, T, 1, &Cp[0][HOMO_index], nmo);

        C_DGEMV('N', nso, nso, 1.0, YYp[0], nso, &Cp[0][HOMO_index], nmo, 0.0, T, 1);
        yy = -C_DDOT(nso, T, 1, &Cp[0][HOMO_index], nmo);

        C_DGEMV('N', nso, nso, 1.0, ZZp[0], nso, &Cp[0][HOMO_index], nmo, 0.0, T, 1);
        zz = -C_DDOT(nso, T, 1, &Cp[0][HOMO_index], nmo);

        quadrupole_ints.clear();
        
        delete[] T;

        double Sx2 = xx - x * x;
        double Sy2 = yy - y * y;
        double Sz2 = zz - z * z;

        Rext = sqrt(Sx2 + Sy2 + Sz2); 

        if (debug_) {
            fprintf(outfile, "  <x> =  %24.14f\n", x);
            fprintf(outfile, "  <y> =  %24.14f\n", y);
            fprintf(outfile, "  <z> =  %24.14f\n", z);
            fprintf(outfile, "  <xx> = %24.14f\n", xx);
            fprintf(outfile, "  <yy> = %24.14f\n", yy);
            fprintf(outfile, "  <zz> = %24.14f\n", zz);
            fprintf(outfile, "  Sx^2 = %24.14f\n", Sx2);
            fprintf(outfile, "  Sy^2 = %24.14f\n", Sy2);
            fprintf(outfile, "  Sz^2 = %24.14f\n", Sz2);
            fprintf(outfile, "  <R> =  %24.14f\n", Rext);
            fprintf(outfile, "  a =    %24.14f\n", alpha);
            fprintf(outfile, "  b =    %24.14f\n\n", beta);
        }

        initial_omega_ = 1.0 / (alpha * Rext + beta);        
   
        fprintf(outfile, "  ==> Omega Guess: HOMO Size <== \n\n");
        fprintf(outfile, "   HOMO <R>:         %12.5E\n", Rext);
        fprintf(outfile, "   Initial Omega^-1: %12.5E\n", 1.0 / initial_omega_);
        fprintf(outfile, "   Initial Omega:    %12.5E\n\n", initial_omega_);
 
    } else {
        initial_omega_ = functional_->getOmega();

        fprintf(outfile, "  ==> Omega Guess: Functional Default <== \n\n");
        fprintf(outfile, "   Initial Omega^-1: %12.5E\n", 1.0 / initial_omega_);
        fprintf(outfile, "   Initial Omega:    %12.5E\n\n", initial_omega_);
    }

    fflush(outfile);
}
void OmegaIPKS::form_DF()
{
    df_ = boost::shared_ptr<OmegaDF>(new OmegaDF(psio_, basisset_, auxiliary_));
    df_->set_omega(initial_omega_);
}
void OmegaIPKS::form_KS()
{
    ks_ = boost::shared_ptr<UKSPotential>(new UKSPotential(functional_, molecule_, basisset_, options_)); 
    functional_->setOmega(initial_omega_);
}
void OmegaIPKS::form_H()
{
    H_ = OmegaKS::build_H(basisset_);

    if (debug_ > 2 ) {
        H_->print(outfile);
    }
}
void OmegaIPKS::form_X()
{
    S_ = OmegaKS::build_S(basisset_);
    X_ = OmegaKS::build_X(basisset_, options_.get_double("S_MIN_EIGENVALUE"));

    if (debug_ > 2)
        X_->print();

    int nmo = X_->colspi()[0];
    int nso = X_->rowspi()[0];

    if (print_)
        fprintf(outfile, "  Canonical Orthogonalization: %d of %d functions retained, %d projected\n\n", nmo, nso, nso - nmo);
}
void OmegaIPKS::finalize() 
{
    fprintf(outfile, "  *** Omega Optimization Profile ***\n\n");
    fprintf(outfile, "  %4s %12s %12s %12s %12s %s\n", "Iter", "Omega", "kIP", "IP", "Delta IP", "Step Type");
    for (int i = 0; i < omega_trace_.size(); i++) {
        fprintf(outfile, "  %4d %12.5E %12.5E %12.5E %12.5E %s\n", i+1, get<0>(omega_trace_[i]),
            get<1>(omega_trace_[i]), get<2>(omega_trace_[i]), get<1>(omega_trace_[i]) -
            get<2>(omega_trace_[i]), get<3>(omega_trace_[i]).c_str());
    }
    fprintf(outfile, "\n");

    Process::environment.globals["OPTIMIZED OMEGA"] = get<0>(omega_trace_[omega_trace_.size()-1]);

    S_.reset();
    X_.reset();
    H_.reset();
}
void OmegaIPKS::populate()
{
    int na_N = reference_->nalphapi()[0];
    int nb_N = reference_->nbetapi()[0];

    int na_M;
    int nb_M;
    if (na_N == nb_N) {
        na_M = na_N;
        nb_M = nb_N - 1;
    } else {
        na_M = na_N - 1;
        nb_M = nb_N;
    }

    boost::shared_ptr<Matrix> Ca = reference_->Ca();
    boost::shared_ptr<Matrix> Cb = reference_->Cb();

    wfns_["N"] = boost::shared_ptr<OmegaWavefunction>(new OmegaWavefunction(
        options_, psio_, basisset_, Ca, na_N, Cb, nb_N, S_, X_, H_, df_, ks_)); 

    wfns_["N-1"] = boost::shared_ptr<OmegaWavefunction>(new OmegaWavefunction(
        options_, psio_, basisset_, Ca, na_M, Cb, nb_M, S_, X_, H_, df_, ks_)); 
}
void OmegaIPKS::omega_step(int iteration)
{
    if (iteration == 1) {
        double old_omega = functional_->getOmega();
        double kIP = wfns_["N"]->koopmansIP();
        double IP = wfns_["N-1"]->E() - wfns_["N"]->E();    
        omega_trace_.push_back(make_tuple(old_omega, kIP, IP, "Initial Guess"));        
        bracketed_ = false;

        if (options_.get_bool("OMEGA_GUESS_INTERPOLATE")) {
            if (kIP - IP > 0.0) {
                Fa_r_N_->copy(wfns_["N"]->Fa());
                Fb_r_N_->copy(wfns_["N"]->Fb());
                Fa_r_M_->copy(wfns_["N-1"]->Fa());
                Fb_r_M_->copy(wfns_["N-1"]->Fb());
            } else {
                Fa_l_N_->copy(wfns_["N"]->Fa());
                Fb_l_N_->copy(wfns_["N"]->Fb());
                Fa_l_M_->copy(wfns_["N-1"]->Fa());
                Fb_l_M_->copy(wfns_["N-1"]->Fb());
            }
        }
    }

    // Bracket first
    if (!bracketed_) {
        int n = omega_trace_.size();
        double old_omega = get<0>(omega_trace_[n-1]);
        double kIP =       get<1>(omega_trace_[n-1]);
        double IP =        get<2>(omega_trace_[n-1]);
        double delta = kIP - IP;
    
        double multiplier = options_.get_double("OMEGA_BRACKET_ALPHA");
        
        // too HF-like
        if (delta > 0) {
            multiplier = 1.0 / multiplier;
        // too KS-like
        } else { 
            multiplier = multiplier;
        }

        double new_omega = old_omega * multiplier;
    
        df_->set_omega(new_omega);
        functional_->setOmega(new_omega);
   
        fprintf(outfile, "  *** Step %4d: Bracketing Step ***\n\n", iteration + 1);
        fprintf(outfile, "    Old Omega:    %15.8E\n", old_omega); 
        fprintf(outfile, "    Old kIP:      %15.8E\n", kIP); 
        fprintf(outfile, "    Old IP:       %15.8E\n", IP); 
        fprintf(outfile, "    Old Delta IP: %15.8E\n", kIP - IP); 
        fprintf(outfile, "    Multipler:    %15.8E\n", multiplier); 
        fprintf(outfile, "    New Omega:    %15.8E\n\n", new_omega); 
 
        return;
    }

    double new_omega = 0.0;
    
    if (options_.get_str("OMEGA_ROOT_ALGORITHM") == "REGULA_FALSI") {
        fprintf(outfile, "  *** Step %4d: Regula Falsi Step ***\n\n", iteration + 1);
        if (fabs(delta_l_ - delta_r_) > 1.0E-14) {
            new_omega = -delta_l_ * (omega_r_ - omega_l_) / (delta_r_ - delta_l_) + omega_l_;
        } else { 
            new_omega = 0.5 * (omega_l_ + omega_r_); 
        }
    } else if (options_.get_str("OMEGA_ROOT_ALGORITHM") == "BISECTION") {
        fprintf(outfile, "  *** Step %4d: Bisection Step ***\n\n", iteration + 1);
        new_omega = 0.5 * (omega_l_ + omega_r_); 
    } else if (options_.get_str("OMEGA_ROOT_ALGORITHM") == "BRENT") {
        fprintf(outfile, "  *** Step %4d: Brent's Method Step ***\n\n", iteration + 1);

        // Assing temps so that larger absolute delta is on the left
        if (fabs(delta_l_) > fabs(delta_r_)) {
            omega_a_ = omega_r_;
            omega_b_ = omega_l_;
            delta_a_ = delta_r_;
            delta_b_ = delta_l_;
        } else { 
            omega_b_ = omega_r_;
            omega_a_ = omega_l_;
            delta_b_ = delta_r_;
            delta_a_ = delta_l_;
        }

        // First iteration
        if (omega_c_ == -1.0) {
            omega_c_ = omega_a_;
            delta_c_ = delta_a_;
        }

        // Quadratic interpolation (if reasonably well-conditioned)
        if (fabs(delta_a_ - delta_b_) > 1.0E-14 &&fabs(delta_a_ - delta_b_) > 1.0E-14 &&
            fabs(delta_a_ - delta_b_) > 1.0E-14) {
            new_omega = omega_a_ * delta_b_ * delta_c_ / ((delta_a_ - delta_b_) * (delta_a_ - delta_c_)) +
                        omega_b_ * delta_a_ * delta_c_ / ((delta_b_ - delta_a_) * (delta_b_ - delta_c_)) +
                        omega_c_ * delta_b_ * delta_a_ / ((delta_c_ - delta_a_) * (delta_c_ - delta_b_));
        } else if (fabs(delta_a_ - delta_b_) > 1.0E-14) {
        // Otherwise linear interpolation
            new_omega = omega_b_ - delta_b_ * (omega_b_ - omega_a_) / (delta_b_ - delta_a_);    
        } else {
        // Otherwise give up
            new_omega = omega_a_;
        }

        // check if we gotta bisect
        bool bisection_required = false;

        // Condition 1: s \in [(3a+b)/4, b]
        if ((3.0 * omega_a_ + omega_b_) / 4.0 < omega_b_) {
            if ((3.0 * omega_a_ + omega_b_) / 4.0 > new_omega ||
                omega_b_ < new_omega)
                bisection_required = true;    
        } else {
            if ((3.0 * omega_a_ + omega_b_) / 4.0 < new_omega ||
                omega_b_ > new_omega)

                bisection_required = true;    
        }

        if (mflag_ && fabs(new_omega - omega_b_) >= 0.5 * fabs(omega_b_ - omega_c_))
            bisection_required = true;
        
        if (!mflag_ && fabs(new_omega - omega_b_) >= 0.5 * fabs(omega_b_ - omega_d_))
            bisection_required = true;
        

        if (bisection_required) {
            new_omega = 0.5 * (omega_a_ + omega_b_); 
            mflag_ = true;
        } else {
            mflag_ = false;
        }

        // Update the old guys
        omega_d_ = omega_c_;
        omega_c_ = omega_b_;
    } else if (options_.get_str("OMEGA_ROOT_ALGORITHM") == "REGULA_FALSI2") {
        fprintf(outfile, "  *** Step %4d: Regula Falsi/Bisection Step ***\n\n", iteration + 1);
        
        // Regula Falsi
        if (left_ < 2 && right_ < 2 && fabs(delta_l_ - delta_r_) > 1.0E-14) {  
            new_omega = omega_l_ - delta_l_ * (omega_r_ - omega_l_) / (delta_r_ - delta_l_);    
        } else {
        // Bisection
            new_omega = 0.5 * (omega_l_ + omega_r_); 
            left_ = right_ = 0;
        }

    } else if (options_.get_str("OMEGA_ROOT_ALGORITHM") == "REGULA_FALSI3") {
        fprintf(outfile, "  *** Step %4d: Regula Falsi 3 Step ***\n\n", iteration + 1);
        
        // Regula Falsi
        if (fabs(delta_l_ - delta_r_) > 1.0E-14) {  
            if (left_ > 1) 
                delta_l_ *= 0.5;
            if (right_ > 1) 
                delta_r_ *= 0.5;
            new_omega = omega_l_ - delta_l_ * (omega_r_ - omega_l_) / (delta_r_ - delta_l_);    
        } else {    
            new_omega = 0.5 * (omega_l_ + omega_r_); 
        }
    }

    // Set the omega
    df_->set_omega(new_omega);
    functional_->setOmega(new_omega);

    // Interpolate the guess
    if (options_.get_bool("OMEGA_GUESS_INTERPOLATE")) {

        boost::shared_ptr<Matrix> Fa(new Matrix("Fa Interpolated", basisset_->nbf(), basisset_->nbf()));
        boost::shared_ptr<Matrix> Fb(new Matrix("Fa Interpolated", basisset_->nbf(), basisset_->nbf()));
        double** Fap = Fa->pointer();
        double** Fbp = Fb->pointer();
        double** Fa_l_Np = Fa_l_N_->pointer();
        double** Fa_l_Mp = Fa_l_M_->pointer();
        double** Fb_l_Np = Fb_l_N_->pointer();
        double** Fb_l_Mp = Fb_l_M_->pointer();
        double** Fa_r_Np = Fa_r_N_->pointer();
        double** Fa_r_Mp = Fa_r_M_->pointer();
        double** Fb_r_Np = Fb_r_N_->pointer();
        double** Fb_r_Mp = Fb_r_M_->pointer();

        double alpha = (new_omega - omega_l_) / (omega_r_ - omega_l_); 
        double beta = 1.0 - alpha;
   
        C_DAXPY(basisset_->nbf() * (ULI) basisset_->nbf(), alpha, Fa_r_Np[0], 1,
            Fap[0], 1);
        C_DAXPY(basisset_->nbf() * (ULI) basisset_->nbf(), beta, Fa_l_Np[0], 1,
            Fap[0], 1);
 
        C_DAXPY(basisset_->nbf() * (ULI) basisset_->nbf(), alpha, Fb_r_Np[0], 1,
            Fbp[0], 1);
        C_DAXPY(basisset_->nbf() * (ULI) basisset_->nbf(), beta, Fb_l_Np[0], 1,
            Fbp[0], 1);
 
        wfns_["N"]->guess(Fa,Fb);

        Fa->zero();
        Fb->zero();
        
        C_DAXPY(basisset_->nbf() * (ULI) basisset_->nbf(), alpha, Fa_r_Mp[0], 1,
            Fap[0], 1);
        C_DAXPY(basisset_->nbf() * (ULI) basisset_->nbf(), beta, Fa_l_Mp[0], 1,
            Fap[0], 1);
 
        C_DAXPY(basisset_->nbf() * (ULI) basisset_->nbf(), alpha, Fb_r_Mp[0], 1,
            Fbp[0], 1);
        C_DAXPY(basisset_->nbf() * (ULI) basisset_->nbf(), beta, Fb_l_Mp[0], 1,
            Fbp[0], 1);

        wfns_["N-1"]->guess(Fa,Fb);
    }

    fprintf(outfile, "    New Omega:    %15.8E\n\n", new_omega); 
}
bool OmegaIPKS::is_omega_converged()
{
    double new_omega = functional_->getOmega();
    double new_kIP = wfns_["N"]->koopmansIP();
    double new_IP = wfns_["N-1"]->E() - wfns_["N"]->E();    
    double new_delta = new_kIP - new_IP; 

    if (options_.get_bool("OMEGA_GUESS_INTERPOLATE")) {
        if (new_delta > 0.0) { 
            Fa_r_N_->copy(wfns_["N"]->Fa());
            Fb_r_N_->copy(wfns_["N"]->Fb());
            Fa_r_M_->copy(wfns_["N-1"]->Fa());
            Fb_r_M_->copy(wfns_["N-1"]->Fb());
        } else {
            Fa_l_N_->copy(wfns_["N"]->Fa());
            Fb_l_N_->copy(wfns_["N"]->Fb());
            Fa_l_M_->copy(wfns_["N-1"]->Fa());
            Fb_l_M_->copy(wfns_["N-1"]->Fb());
        }
    }

    fprintf(outfile, "  *** Step Results ***\n\n");
    fprintf(outfile, "    Omega:    %15.8E\n", new_omega); 
    fprintf(outfile, "    kIP:      %15.8E\n", new_kIP); 
    fprintf(outfile, "    IP:       %15.8E\n", new_IP); 
    fprintf(outfile, "    Delta IP: %15.8E\n\n", new_kIP - new_IP); 
 
    // Has the root been bracketed?
    if (!bracketed_) {

        int n = omega_trace_.size();
        double old_omega = get<0>(omega_trace_[n-1]);
        double old_kIP =   get<1>(omega_trace_[n-1]);
        double old_IP =    get<2>(omega_trace_[n-1]);
        double old_delta = old_kIP - old_IP;

        omega_trace_.push_back(make_tuple(new_omega, new_kIP, new_IP, "Bracket Step"));        

        if (old_delta * new_delta < 0.0) {
            fprintf(outfile, "  *** Omega Bracketed ***\n\n");

            if (old_delta < 0.0) {
                omega_l_ = old_omega;
                delta_l_ = old_delta;
                omega_r_ = new_omega;
                delta_r_ = new_delta;
            } else {
                omega_r_ = old_omega;
                delta_r_ = old_delta;
                omega_l_ = new_omega;
                delta_l_ = new_delta;
            }

            // Brent's method special stuff
            omega_a_ = omega_b_ = omega_c_ = omega_d_ = -1.0;
            omega_a_ = omega_b_ = omega_c_ = -1.0;
            mflag_ = true;

            left_ = 0;
            right_ = 0;

            bracketed_ = true;
        }
        return false;
    }

    if (new_delta < 0.0) {
        omega_l_ = new_omega;
        delta_l_ = new_delta;
        left_ = 0;
        right_++;
    } else {
        omega_r_ = new_omega;
        delta_r_ = new_delta;
        left_++;
        right_ = 0;
    }

    if (options_.get_str("OMEGA_ROOT_ALGORITHM") == "REGULA_FALSI") {
        omega_trace_.push_back(make_tuple(new_omega, new_kIP, new_IP, "Regula-Falsi Step"));        
    } else if (options_.get_str("OMEGA_ROOT_ALGORITHM") == "REGULA_FALSI2") {
        omega_trace_.push_back(make_tuple(new_omega, new_kIP, new_IP, "Regula-Falsi2 Step"));        
    } else if (options_.get_str("OMEGA_ROOT_ALGORITHM") == "REGULA_FALSI3") {
        omega_trace_.push_back(make_tuple(new_omega, new_kIP, new_IP, "Regula-Falsi3 Step"));        
    } else if (options_.get_str("OMEGA_ROOT_ALGORITHM") == "BISECTION") {
        omega_trace_.push_back(make_tuple(new_omega, new_kIP, new_IP, "Bisection Step"));        
    } else if (options_.get_str("OMEGA_ROOT_ALGORITHM") == "BRENT") {
        omega_trace_.push_back(make_tuple(new_omega, new_kIP, new_IP, "Brent's Method Step"));        
    }

    double w_converge = pow(10.0, - (double) options_.get_int("OMEGA_CONVERGE"));
    return (fabs(omega_l_ - omega_r_) < w_converge); 
}    

}} // End Namespaces
