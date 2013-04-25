#include "dftsapt.h"
#include <libmints/mints.h>
#include <libmints/sieve.h>
#include <libfock/jk.h>
#include <libqt/qt.h>
#include <libpsio/psio.hpp>
#include <libpsio/psio.h>
#include <psi4-dec.h>
#include <psifiles.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;
using namespace boost;
using namespace std;

namespace psi {
namespace dftsapt {

DFTSAPT::DFTSAPT()
{
    common_init();
}
DFTSAPT::~DFTSAPT()
{
}
void DFTSAPT::common_init()
{
    print_ = 1;
    debug_ = 0;
    bench_ = 0;
}
boost::shared_ptr<DFTSAPT> DFTSAPT::build(boost::shared_ptr<Wavefunction> d,
                                          boost::shared_ptr<Wavefunction> mA, 
                                          boost::shared_ptr<Wavefunction> mB)
{
    DFTSAPT* sapt = new DFTSAPT();

    Options& options = Process::environment.options; 

    sapt->print_ = options.get_int("PRINT");
    sapt->debug_ = options.get_int("DEBUG");
    sapt->bench_ = options.get_int("BENCH");

    sapt->memory_ = (unsigned long int)(Process::environment.get_memory() * options.get_double("SAPT_MEM_FACTOR") * 0.125);

    sapt->cpks_maxiter_ = options.get_int("MAXITER");
    sapt->cpks_delta_ = options.get_double("D_CONVERGENCE");

    sapt->dimer_     = d->molecule();
    sapt->monomer_A_ = mA->molecule();
    sapt->monomer_B_ = mB->molecule();

    sapt->E_dimer_     = d->reference_energy();
    sapt->E_monomer_A_ = mA->reference_energy();
    sapt->E_monomer_B_ = mB->reference_energy();

    sapt->primary_   = d->basisset();
    sapt->primary_A_ = mA->basisset();
    sapt->primary_B_ = mB->basisset();

    if (sapt->primary_A_->nbf() != sapt->primary_B_->nbf() || sapt->primary_->nbf() != sapt->primary_A_->nbf()) {
        throw PSIEXCEPTION("Monomer-centered bases not allowed in DFT-SAPT");
    }

    sapt->Cocc_A_     = mA->Ca_subset("AO","OCC");
    sapt->Cvir_A_     = mA->Ca_subset("AO","VIR");
    sapt->eps_occ_A_  = mA->epsilon_a_subset("AO","OCC");
    sapt->eps_vir_A_  = mA->epsilon_a_subset("AO","VIR");

    sapt->Caocc_A_    = mA->Ca_subset("AO","ACTIVE_OCC");
    sapt->Cavir_A_    = mA->Ca_subset("AO","ACTIVE_VIR");

    sapt->eps_focc_A_ = mA->epsilon_a_subset("AO","FROZEN_OCC");
    sapt->eps_aocc_A_ = mA->epsilon_a_subset("AO","ACTIVE_OCC");
    sapt->eps_avir_A_ = mA->epsilon_a_subset("AO","ACTIVE_VIR");
    sapt->eps_fvir_A_ = mA->epsilon_a_subset("AO","FROZEN_VIR");

    sapt->Cocc_B_     = mB->Ca_subset("AO","OCC");
    sapt->Cvir_B_     = mB->Ca_subset("AO","VIR");
    sapt->eps_occ_B_  = mB->epsilon_a_subset("AO","OCC");
    sapt->eps_vir_B_  = mB->epsilon_a_subset("AO","VIR");

    sapt->Caocc_B_    = mB->Ca_subset("AO","ACTIVE_OCC");
    sapt->Cavir_B_    = mB->Ca_subset("AO","ACTIVE_VIR");

    sapt->eps_focc_B_ = mB->epsilon_a_subset("AO","FROZEN_OCC");
    sapt->eps_aocc_B_ = mB->epsilon_a_subset("AO","ACTIVE_OCC");
    sapt->eps_avir_B_ = mB->epsilon_a_subset("AO","ACTIVE_VIR");
    sapt->eps_fvir_B_ = mB->epsilon_a_subset("AO","FROZEN_VIR");

    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    sapt->mp2fit_ = BasisSet::construct(parser, sapt->dimer_, "DF_BASIS_SAPT");

    return boost::shared_ptr<DFTSAPT>(sapt);
}
double DFTSAPT::compute_energy()
{
    print_header();

    fock_terms();
    
    mp2_terms();

    print_trailer();
    
    Process::environment.globals["SAPT ENERGY"] = 0.0;

    return 0.0;
}
void DFTSAPT::print_header() const 
{
    fprintf(outfile, "\t --------------------------------------------------------\n");
    fprintf(outfile, "\t                        DF-DFT-SAPT                      \n");
    fprintf(outfile, "\t               Rob Parrish and Ed Hohenstein             \n");
    fprintf(outfile, "\t --------------------------------------------------------\n");
    fprintf(outfile, "\n");

    fprintf(outfile, "  ==> Sizes <==\n");
    fprintf(outfile, "\n");

    fprintf(outfile, "   => Resources <=\n\n");
    
    fprintf(outfile, "    Memory (MB):       %11ld\n", (memory_ *8L) / (1024L * 1024L));
    fprintf(outfile, "\n");

    fprintf(outfile, "   => Orbital Ranges <=\n\n");

    int nmo_A = eps_focc_A_->dim() + eps_aocc_A_->dim() + eps_avir_A_->dim() + eps_fvir_A_->dim();
    int nmo_B = eps_focc_B_->dim() + eps_aocc_B_->dim() + eps_avir_B_->dim() + eps_fvir_B_->dim();

    fprintf(outfile, "    ------------------\n");    
    fprintf(outfile, "    %-6s %5s %5s\n", "Range", "M_A", "M_B");
    fprintf(outfile, "    ------------------\n");    
    fprintf(outfile, "    %-6s %5d %5d\n", "nso", primary_A_->nbf(), primary_B_->nbf());
    fprintf(outfile, "    %-6s %5d %5d\n", "nmo", nmo_A, nmo_B);
    fprintf(outfile, "    %-6s %5d %5d\n", "nocc", eps_aocc_A_->dim() + eps_focc_A_->dim(), eps_aocc_B_->dim() + eps_focc_B_->dim());
    fprintf(outfile, "    %-6s %5d %5d\n", "nvir", eps_avir_A_->dim() + eps_fvir_A_->dim(), eps_avir_B_->dim() + eps_fvir_B_->dim());
    fprintf(outfile, "    %-6s %5d %5d\n", "nfocc", eps_focc_A_->dim(), eps_focc_B_->dim());
    fprintf(outfile, "    %-6s %5d %5d\n", "naocc", eps_aocc_A_->dim(), eps_aocc_B_->dim());
    fprintf(outfile, "    %-6s %5d %5d\n", "navir", eps_avir_A_->dim(), eps_avir_B_->dim());
    fprintf(outfile, "    %-6s %5d %5d\n", "nfvir", eps_fvir_A_->dim(), eps_fvir_B_->dim());
    fprintf(outfile, "    ------------------\n");    
    fprintf(outfile, "\n");

    fprintf(outfile, "   => Primary Basis Set <=\n\n");
    primary_->print_by_level(outfile, print_);

    fflush(outfile);
}
void DFTSAPT::fock_terms()
{
    fprintf(outfile, "  SCF TERMS:\n\n");

    // ==> Setup <== //    

    boost::shared_ptr<JK> jk = JK::build_JK();
   
    // TODO: Account for 2-index overhead in memory
    int nA = Cocc_A_->ncol();
    int nB = Cocc_B_->ncol();
    int nso = Cocc_A_->nrow();
    long int jk_memory = (long int)memory_; 
    jk_memory -= 24 * nso * nso;
    jk_memory -=  4 * nA * nso;
    jk_memory -=  4 * nB * nso;
    if (jk_memory < 0L) {
        throw PSIEXCEPTION("Too little static memory for DFTSAPT::fock_terms");
    }
    jk->set_memory((unsigned long int )jk_memory);

    jk->set_do_J(true);
    jk->set_do_K(true);
    jk->initialize();
    jk->print_header();

    // ==> Generalized Fock Source Terms [Elst/Exch] <== //

    // => Steric Interaction Density Terms (T) <= //

    boost::shared_ptr<Matrix> S = build_S(primary_);
    boost::shared_ptr<Matrix> Sij = build_Sij(S);
    boost::shared_ptr<Matrix> Sij_n = build_Sij_n(Sij);

    std::map<std::string, boost::shared_ptr<Matrix> > Cbar_n = build_Cbar(Sij_n);
    boost::shared_ptr<Matrix> C_T_A_n  = Cbar_n["C_T_A"];
    boost::shared_ptr<Matrix> C_T_B_n  = Cbar_n["C_T_B"];
    boost::shared_ptr<Matrix> C_T_AB_n = Cbar_n["C_T_AB"];
    boost::shared_ptr<Matrix> C_T_BA_n = Cbar_n["C_T_BA"];

    if (debug_ > 2) {
        S->print();
        Sij->print();
        Sij_n->print();
    }

    Sij.reset();
    Sij_n.reset();

    // => Load the JK Object <= //
    
    std::vector<SharedMatrix>& Cl = jk->C_left();
    std::vector<SharedMatrix>& Cr = jk->C_right();
    Cl.clear();
    Cr.clear();

    // J/K[P^A]
    Cl.push_back(Cocc_A_);
    Cr.push_back(Cocc_A_);
    // J/K[T^A, S^\infty]
    Cl.push_back(Cocc_A_);
    Cr.push_back(C_T_A_n);
    // J/K[T^AB, S^\infty]
    Cl.push_back(Cocc_A_);
    Cr.push_back(C_T_AB_n);
    // J/K[P^B]
    Cl.push_back(Cocc_B_);
    Cr.push_back(Cocc_B_);
    
    // => Compute the JK matrices <= //

    jk->compute();

    // => Unload the JK Object <= //

    const std::vector<SharedMatrix>& J = jk->J();
    const std::vector<SharedMatrix>& K = jk->K();

    boost::shared_ptr<Matrix> J_A      = J[0]; 
    boost::shared_ptr<Matrix> J_T_A_n  = J[1]; 
    boost::shared_ptr<Matrix> J_T_AB_n = J[2]; 
    boost::shared_ptr<Matrix> J_B      = J[3]; 

    boost::shared_ptr<Matrix> K_A      = K[0]; 
    boost::shared_ptr<Matrix> K_T_A_n  = K[1]; 
    boost::shared_ptr<Matrix> K_T_AB_n = K[2]; 
    boost::shared_ptr<Matrix> K_B      = K[3]; 

    // => Compute the D matrices <= //
    
    boost::shared_ptr<Matrix> D_A    = build_D(Cocc_A_, Cocc_A_);
    boost::shared_ptr<Matrix> D_B    = build_D(Cocc_B_, Cocc_B_);

    // => Compute the V matrices <= //

    boost::shared_ptr<Matrix> V_A = build_V(primary_A_);
    boost::shared_ptr<Matrix> V_B = build_V(primary_B_);

    // ==> Electrostatic Terms <== //

    double Elst10 = 0.0;

    // Classical physics (watch for cancellation!)
    double Enuc = 0.0;
    Enuc += dimer_->nuclear_repulsion_energy();
    Enuc -= monomer_A_->nuclear_repulsion_energy();
    Enuc -= monomer_B_->nuclear_repulsion_energy();

    std::vector<double> Elst10_terms;
    Elst10_terms.resize(4);
    
    Elst10_terms[0] += 2.0 * D_A->vector_dot(V_B);
    Elst10_terms[1] += 2.0 * D_B->vector_dot(V_A);
    Elst10_terms[2] += 4.0 * D_B->vector_dot(J_A);
    Elst10_terms[3] += 1.0 * Enuc;
    
    for (int k = 0; k < Elst10_terms.size(); k++) {
        Elst10 += Elst10_terms[k];
    }

    if (debug_) {
        for (int k = 0; k < Elst10_terms.size(); k++) {
            fprintf(outfile,"    Elst10 (%1d)          = %18.12lf H\n",k+1,Elst10_terms[k]);
        }
    }

    energies_["Elst10"] = Elst10;

    fprintf(outfile,"    Elst10,r            = %18.12lf H\n",Elst10);
    fflush(outfile);

    // ==> Exchange Terms (S^\infty) <== //
    
    double Exch10_n = 0.0;

    // => Compute the T matrices <= //

    boost::shared_ptr<Matrix> T_A_n  = build_D(Cocc_A_, C_T_A_n);
    boost::shared_ptr<Matrix> T_B_n  = build_D(Cocc_B_, C_T_B_n);
    boost::shared_ptr<Matrix> T_BA_n = build_D(Cocc_B_, C_T_BA_n);
    boost::shared_ptr<Matrix> T_AB_n = build_D(Cocc_A_, C_T_AB_n);

    std::vector<double> Exch10_n_terms;
    Exch10_n_terms.resize(9);

    Exch10_n_terms[0] -= 2.0 * D_A->vector_dot(K_B); 
    
    Exch10_n_terms[1] += 2.0 * T_A_n->vector_dot(V_B);
    Exch10_n_terms[1] += 4.0 * T_A_n->vector_dot(J_B);
    Exch10_n_terms[1] -= 2.0 * T_A_n->vector_dot(K_B);

    Exch10_n_terms[2] += 2.0 * T_B_n->vector_dot(V_A);
    Exch10_n_terms[2] += 4.0 * T_B_n->vector_dot(J_A);
    Exch10_n_terms[2] -= 2.0 * T_B_n->vector_dot(K_A);

    Exch10_n_terms[3] += 2.0 * T_AB_n->vector_dot(V_A); 
    Exch10_n_terms[3] += 4.0 * T_AB_n->vector_dot(J_A);
    Exch10_n_terms[3] -= 2.0 * T_AB_n->vector_dot(K_A);

    Exch10_n_terms[4] += 2.0 * T_AB_n->vector_dot(V_B);
    Exch10_n_terms[4] += 4.0 * T_AB_n->vector_dot(J_B);
    Exch10_n_terms[4] -= 2.0 * T_AB_n->vector_dot(K_B);

    Exch10_n_terms[5] += 4.0 * T_B_n->vector_dot(J_T_AB_n);
    Exch10_n_terms[5] -= 2.0 * T_B_n->vector_dot(K_T_AB_n);

    Exch10_n_terms[6] += 4.0 * T_A_n->vector_dot(J_T_AB_n);
    Exch10_n_terms[6] -= 2.0 * T_A_n->vector_dot(K_T_AB_n);

    Exch10_n_terms[7] += 4.0 * T_B_n->vector_dot(J_T_A_n);
    Exch10_n_terms[7] -= 2.0 * T_B_n->vector_dot(K_T_A_n);

    Exch10_n_terms[8] += 4.0 * T_AB_n->vector_dot(J_T_AB_n);
    Exch10_n_terms[8] -= 2.0 * T_AB_n->vector_dot(K_T_AB_n);

    for (int k = 0; k < Exch10_n_terms.size(); k++) {
        Exch10_n += Exch10_n_terms[k];
    }

    if (debug_) {
        for (int k = 0; k < Exch10_n_terms.size(); k++) {
            fprintf(outfile,"    Exch10 (%1d)          = %18.12lf H\n",k+1,Exch10_n_terms[k]);
        }
    }

    energies_["Exch10"] = Exch10_n;

    fprintf(outfile,"    Exch10              = %18.12lf H\n",Exch10_n);
    fflush(outfile);

    // Clear up some memory
    C_T_A_n.reset();
    C_T_B_n.reset();
    C_T_BA_n.reset();
    C_T_AB_n.reset();
    T_A_n.reset();
    T_B_n.reset();
    T_BA_n.reset();
    T_AB_n.reset();
    J_T_A_n.reset();
    J_T_AB_n.reset();
    K_T_A_n.reset();
    K_T_AB_n.reset();

    // => ESP <= //
        
    // Electrostatic potential of monomer A (AO)
    boost::shared_ptr<Matrix> W_A(V_A->clone());
    W_A->copy(J_A);
    W_A->scale(2.0);
    W_A->add(V_A);

    // Electrostatic potential of monomer B (AO)
    boost::shared_ptr<Matrix> W_B(V_B->clone());
    W_B->copy(J_B);
    W_B->scale(2.0);
    W_B->add(V_B);

    // Electrostatic potential of monomer B (ov space of A)
    boost::shared_ptr<Matrix> w_B = build_w(W_B,Cocc_A_,Cvir_A_);
    // Electrostatic potential of monomer A (ov space of B)
    boost::shared_ptr<Matrix> w_A = build_w(W_A,Cocc_B_,Cvir_B_);

    // ==> Exchange Terms (S^2) <== //
    
    double Exch10_2 = 0.0;

    // => Compute the T matrices <= //

    std::vector<double> Exch10_2_terms;
    Exch10_2_terms.resize(3);

    {
        int vA = Cvir_A_->colspi()[0];
        int vB = Cvir_B_->colspi()[0];

        boost::shared_ptr<Matrix> Sbr(new Matrix("Sbr",nB,vA));
        boost::shared_ptr<Matrix> Sas(new Matrix("Sas",nA,vB));
        boost::shared_ptr<Matrix> Sab(new Matrix("Sab",nA,nB));
        boost::shared_ptr<Matrix> temp(new Matrix("temp",nso,nso));
        double** Sbrp = Sbr->pointer();
        double** Sasp = Sas->pointer();
        double** Sabp = Sab->pointer();
        double** Tp   = temp->pointer();
        double** Sp   = S->pointer();
        double** Cap  = Cocc_A_->pointer();
        double** Cbp  = Cocc_B_->pointer();
        double** Crp  = Cvir_A_->pointer();
        double** Csp  = Cvir_B_->pointer();
        
        C_DGEMM('T','N',nB,nso,nso,1.0,Cbp[0],nB,Sp[0],nso,0.0,Tp[0],nso);
        C_DGEMM('N','N',nB,vA,nso,1.0,Tp[0],nso,Crp[0],vA,0.0,Sbrp[0],vA);

        C_DGEMM('T','N',nA,nso,nso,1.0,Cap[0],nA,Sp[0],nso,0.0,Tp[0],nso);
        C_DGEMM('N','N',nA,vB,nso,1.0,Tp[0],nso,Csp[0],vB,0.0,Sasp[0],vB);

        C_DGEMM('T','N',nA,nso,nso,1.0,Cap[0],nA,Sp[0],nso,0.0,Tp[0],nso);
        C_DGEMM('N','N',nA,nB,nso,1.0,Tp[0],nso,Cbp[0],nB,0.0,Sabp[0],nB);

        boost::shared_ptr<Matrix> Sar2(new Matrix("Sar2",nA,vA));
        boost::shared_ptr<Matrix> Sbs2(new Matrix("Sbs2",nB,vB));
        double** Sar2p = Sar2->pointer();
        double** Sbs2p = Sbs2->pointer();

        C_DGEMM('N','N',nA,vA,nB,1.0,Sabp[0],nB,Sbrp[0],vA,0.0,Sar2p[0],vA);

        C_DGEMM('T','N',nB,vB,nA,1.0,Sabp[0],nB,Sasp[0],vB,0.0,Sbs2p[0],vB);
       
        Exch10_2_terms[0] = -2.0 * w_B->vector_dot(Sar2); 
        Exch10_2_terms[1] = -2.0 * w_A->vector_dot(Sbs2); 

        boost::shared_ptr<Matrix> R(new Matrix("R",nso,nA));
        double** Rp = R->pointer();
        C_DGEMM('N','T',nso,nA,vB,1.0,Csp[0],vB,Sasp[0],vB,0.0,Rp[0],nA);

        Cl.clear();
        Cr.clear();

        Cl.push_back(Cocc_A_);
        Cr.push_back(R);

        jk->compute();

        boost::shared_ptr<Matrix> K2 = K[0];
        boost::shared_ptr<Matrix> Kbr(new Matrix("Kbr",nB,vA));
        double** K2p = K2->pointer();
        double** Kbrp = Kbr->pointer();
    
        C_DGEMM('T','T',nB,nso,nso,1.0,Cbp[0],nB,K2p[0],nso,0.0,Tp[0],nso);
        C_DGEMM('N','N',nB,vA,nso,1.0,Tp[0],nso,Crp[0],vA,0.0,Kbrp[0],vA);

        Exch10_2_terms[2] = -2.0 * Kbr->vector_dot(Sbr);
    }

    for (int k = 0; k < Exch10_2_terms.size(); k++) {
        Exch10_2 += Exch10_2_terms[k];
    }

    if (debug_) {
        for (int k = 0; k < Exch10_2_terms.size(); k++) {
            fprintf(outfile,"    Exch10 (S^2) (%1d)    = %18.12lf H\n",k+1,Exch10_2_terms[k]);
        }
    }

    energies_["Exch10 (S^2)"] = Exch10_2;
    fprintf(outfile,"    Exch10 (S^2)        = %18.12lf H\n",Exch10_2);
    fprintf(outfile, "\n");
    fflush(outfile);

    // Clear up some memory (TODO)

    // ==> Uncorrelated Second-Order Response Terms [Induction] <== //

    // Compute CPHF
    std::pair<boost::shared_ptr<Matrix>, boost::shared_ptr<Matrix> > x_sol = compute_x(jk,w_B,w_A);
    boost::shared_ptr<Matrix> x_A = x_sol.first;
    boost::shared_ptr<Matrix> x_B = x_sol.second;

    // ==> Induction <== //

    double Ind20r_AB = - x_A->vector_dot(w_B);
    double Ind20r_BA = - x_B->vector_dot(w_A);
    double Ind20r = Ind20r_AB + Ind20r_BA;

    energies_["Ind20,r (A<-B)"] = Ind20r_AB;
    energies_["Ind20,r (B->A)"] = Ind20r_BA;
    energies_["Ind20,r"] = Ind20r_AB + Ind20r_BA;

    fprintf(outfile,"    Ind20,r (A<-B)      = %18.12lf H\n",Ind20r_AB);
    fprintf(outfile,"    Ind20,r (A->B)      = %18.12lf H\n",Ind20r_BA);
    fprintf(outfile,"    Ind20,r             = %18.12lf H\n",Ind20r);
    fflush(outfile);

    // ==> Exchange-Induction <== //

    // => Fock-like matrices <= //

    boost::shared_ptr<Matrix> h_A(V_A->clone());
    h_A->copy(J_A);
    h_A->scale(2.0);
    h_A->add(V_A);
    h_A->subtract(K_A);
    J_A.reset();
    V_A.reset();

    boost::shared_ptr<Matrix> h_B(V_B->clone());
    h_B->copy(J_B);
    h_B->scale(2.0);
    h_B->add(V_B);
    h_B->subtract(K_B);
    J_B.reset();
    V_B.reset();

    // => New pertubations <= //
    
    boost::shared_ptr<Matrix> C_O_A = build_C_O(Cocc_A_,S,D_B);
    boost::shared_ptr<Matrix> C_X_A = build_C_X(x_A, Cvir_A_);
    boost::shared_ptr<Matrix> C_X_B = build_C_X(x_B, Cvir_B_);
    
    Cl.clear();
    Cr.clear();

    // J/K[O]
    Cl.push_back(Cocc_A_);
    Cr.push_back(C_O_A);
    // J/K[X_A]
    Cl.push_back(Cocc_A_);
    Cr.push_back(C_X_A);
    // J/K[X_B]
    Cl.push_back(Cocc_B_);
    Cr.push_back(C_X_B);
    
    // => Compute the JK matrices <= //

    jk->compute();

    // => Unload the JK Object <= //

    boost::shared_ptr<Matrix> J_O      = J[0]; 
    boost::shared_ptr<Matrix> K_O      = K[0]; 
    boost::shared_ptr<Matrix> J_X_A    = J[1]; 
    boost::shared_ptr<Matrix> J_X_B    = J[2]; 

    J_O->scale(0.5); 
    K_O->scale(0.5); 
    J_X_A->scale(0.5); 
    J_X_B->scale(0.5); 

    // Clear the JK memory 
    jk.reset();

    // New generalized densities
    boost::shared_ptr<Matrix> O   = build_D(Cocc_A_,C_O_A);
    boost::shared_ptr<Matrix> X_A = build_D(Cocc_A_,C_X_A);
    boost::shared_ptr<Matrix> X_B = build_D(Cocc_B_,C_X_B);

    // Exch-Ind20 (A<-B)
    double ExchInd20_AB = 0.0;
    std::vector<double> ExchInd20_AB_terms;
    ExchInd20_AB_terms.resize(11);
    ExchInd20_AB_terms[0]  += 2.0 * X_A->vector_dot(K_B);
    ExchInd20_AB_terms[1]  += 2.0 * triple(X_A,S,D_B)->vector_dot(h_A);
    ExchInd20_AB_terms[2]  += 2.0 * X_A->vector_dot(J_O);
    ExchInd20_AB_terms[3]  -= 1.0 * X_A->vector_dot(K_O->transpose());
    ExchInd20_AB_terms[4]  += 1.0 * triple(X_B,S,X_A)->vector_dot(h_B);
    ExchInd20_AB_terms[5]  -= 2.0 * triple(D_B,S,triple(X_A,S,D_B))->vector_dot(W_A);
    ExchInd20_AB_terms[6]  -= 2.0 * triple(D_B,S,O)->vector_dot(J_X_A);
    ExchInd20_AB_terms[7]  -= 1.0 * triple(O,S,X_A)->vector_dot(W_B);
    ExchInd20_AB_terms[8]  -= 0.5 * triple(X_A,S,O->transpose())->vector_dot(W_B);
    ExchInd20_AB_terms[9]  -= 1.0 * triple(D_B,S,X_A)->vector_dot(K_O->transpose());
    ExchInd20_AB_terms[10] += 1.0 * triple(X_A,S,D_B)->vector_dot(K_O);
    for (int k = 0; k < ExchInd20_AB_terms.size(); k++) {
        ExchInd20_AB += ExchInd20_AB_terms[k];
    }
    if (debug_) {
        for (int k = 0; k < ExchInd20_AB_terms.size(); k++) {
            fprintf(outfile,"    Exch-Ind20,r (%2d)   = %18.12lf H\n",k+1,ExchInd20_AB_terms[k]);
        }
    }
    fprintf(outfile,"    Exch-Ind20,r (A<-B) = %18.12lf H\n",ExchInd20_AB);
    fflush(outfile);

    K_O->transpose_this();
    O->transpose_this();

    // Exch-Ind20 (B<-A)
    double ExchInd20_BA = 0.0;
    std::vector<double> ExchInd20_BA_terms;
    ExchInd20_BA_terms.resize(11);
    ExchInd20_BA_terms[0]  += 2.0 * X_B->vector_dot(K_A);
    ExchInd20_BA_terms[1]  += 2.0 * triple(X_B,S,D_A)->vector_dot(h_B);
    ExchInd20_BA_terms[2]  += 2.0 * X_B->vector_dot(J_O);
    ExchInd20_BA_terms[3]  -= 1.0 * X_B->vector_dot(K_O->transpose());
    ExchInd20_BA_terms[4]  += 1.0 * triple(X_A,S,X_B)->vector_dot(h_A);
    ExchInd20_BA_terms[5]  -= 2.0 * triple(D_A,S,triple(X_B,S,D_A))->vector_dot(W_B);
    ExchInd20_BA_terms[6]  -= 2.0 * triple(D_A,S,O)->vector_dot(J_X_B);
    ExchInd20_BA_terms[7]  -= 1.0 * triple(O,S,X_B)->vector_dot(W_A);
    ExchInd20_BA_terms[8]  -= 0.5 * triple(X_B,S,O->transpose())->vector_dot(W_A);
    ExchInd20_BA_terms[9]  -= 1.0 * triple(D_A,S,X_B)->vector_dot(K_O->transpose());
    ExchInd20_BA_terms[10] += 1.0 * triple(X_B,S,D_A)->vector_dot(K_O);
    for (int k = 0; k < ExchInd20_BA_terms.size(); k++) {
        ExchInd20_BA += ExchInd20_BA_terms[k];
    }
    if (debug_) {
        for (int k = 0; k < ExchInd20_BA_terms.size(); k++) {
            fprintf(outfile,"    Exch-Ind20,r (%2d)   = %18.12lf H\n",k+1,ExchInd20_BA_terms[k]);
        }
    }
    fprintf(outfile,"    Exch-Ind20,r (B<-A) = %18.12lf H\n",ExchInd20_BA);
    fflush(outfile);

    double ExchInd20 = ExchInd20_AB + ExchInd20_BA;
    fprintf(outfile,"    Exch-Ind20,r        = %18.12lf H\n",ExchInd20);
    fprintf(outfile,"\n");
    fflush(outfile);
}
void DFTSAPT::mp2_terms()
{
    fprintf(outfile, "  PT2 TERMS:\n\n");

    // TODO 
}
void DFTSAPT::print_trailer() const 
{
    fprintf(outfile, "    Oh all the money that e'er I had, I spent it in good company.\n");
    fprintf(outfile, "\n");
}
boost::shared_ptr<Matrix> DFTSAPT::build_S(boost::shared_ptr<BasisSet> basis)
{
    boost::shared_ptr<IntegralFactory> factory(new IntegralFactory(basis));
    boost::shared_ptr<OneBodyAOInt> Sint(factory->ao_overlap());
    boost::shared_ptr<Matrix> S(new Matrix("S (AO)", basis->nbf(), basis->nbf()));
    Sint->compute(S);
    return S;
}
boost::shared_ptr<Matrix> DFTSAPT::build_V(boost::shared_ptr<BasisSet> basis)
{
    boost::shared_ptr<IntegralFactory> factory(new IntegralFactory(basis));
    boost::shared_ptr<OneBodyAOInt> Sint(factory->ao_potential());
    boost::shared_ptr<Matrix> S(new Matrix("V (AO)", basis->nbf(), basis->nbf()));
    Sint->compute(S);
    return S;
}
boost::shared_ptr<Matrix> DFTSAPT::build_Sij(boost::shared_ptr<Matrix> S)
{
    int nso = Cocc_A_->nrow();
    int nocc_A = Cocc_A_->ncol();
    int nocc_B = Cocc_B_->ncol();
    int nocc = nocc_A + nocc_B;

    boost::shared_ptr<Matrix> Sij(new Matrix("Sij (MO)", nocc, nocc));
    boost::shared_ptr<Matrix> T(new Matrix("T", nso, nocc_B));
    
    double** Sp = S->pointer();
    double** Tp = T->pointer();
    double** Sijp = Sij->pointer();
    double** CAp = Cocc_A_->pointer();
    double** CBp = Cocc_B_->pointer();

    C_DGEMM('N','N',nso,nocc_B,nso,1.0,Sp[0],nso,CBp[0],nocc_B,0.0,Tp[0],nocc_B);
    C_DGEMM('T','N',nocc_A,nocc_B,nso,1.0,CAp[0],nocc_A,Tp[0],nocc_B,0.0,&Sijp[0][nocc_A],nocc);

    Sij->copy_upper_to_lower();

    return Sij;
}
boost::shared_ptr<Matrix> DFTSAPT::build_Sij_n(boost::shared_ptr<Matrix> Sij)
{
    int nocc = Sij->nrow();

    boost::shared_ptr<Matrix> Sij2(new Matrix("Sij^inf (MO)", nocc, nocc));

    double** Sijp = Sij->pointer();
    double** Sij2p = Sij2->pointer();

    Sij2->copy(Sij);
    for (int i = 0; i < nocc; i++) {
        Sij2p[i][i] = 1.0;
    }

    int info;
   
    info = C_DPOTRF('L',nocc,Sij2p[0],nocc); 
    if (info) {
        throw PSIEXCEPTION("Sij DPOTRF failed. How far up the steric wall are you?");
    }

    info = C_DPOTRI('L',nocc,Sij2p[0],nocc); 
    if (info) {
        throw PSIEXCEPTION("Sij DPOTRI failed. How far up the steric wall are you?");
    }

    Sij2->copy_upper_to_lower();

    for (int i = 0; i < nocc; i++) {
        Sij2p[i][i] -= 1.0;
    }
   
    return Sij2; 
}
std::map<std::string, boost::shared_ptr<Matrix> > DFTSAPT::build_Cbar(boost::shared_ptr<Matrix> S)
{
    std::map<std::string, boost::shared_ptr<Matrix> > Cbar;

    int nso = Cocc_A_->nrow();
    int nA = Cocc_A_->ncol();
    int nB = Cocc_B_->ncol();
    int no = nA + nB;

    double** Sp = S->pointer();
    double** CAp = Cocc_A_->pointer();
    double** CBp = Cocc_B_->pointer();
    double** Cp;

    Cbar["C_T_A"] = boost::shared_ptr<Matrix>(new Matrix("C_T_A", nso, nA));
    Cp = Cbar["C_T_A"]->pointer();
    C_DGEMM('N','N',nso,nA,nA,1.0,CAp[0],nA,&Sp[0][0],no,0.0,Cp[0],nA);

    Cbar["C_T_B"] = boost::shared_ptr<Matrix>(new Matrix("C_T_B", nso, nB));
    Cp = Cbar["C_T_B"]->pointer();
    C_DGEMM('N','N',nso,nB,nB,1.0,CBp[0],nB,&Sp[nA][nA],no,0.0,Cp[0],nB);

    Cbar["C_T_BA"] = boost::shared_ptr<Matrix>(new Matrix("C_T_BA", nso, nB));   
    Cp = Cbar["C_T_BA"]->pointer();
    C_DGEMM('N','N',nso,nB,nA,1.0,CAp[0],nA,&Sp[0][nA],no,0.0,Cp[0],nB);

    Cbar["C_T_AB"] = boost::shared_ptr<Matrix>(new Matrix("C_T_AB", nso, nA));   
    Cp = Cbar["C_T_AB"]->pointer();
    C_DGEMM('N','N',nso,nA,nB,1.0,CBp[0],nB,&Sp[nA][0],no,0.0,Cp[0],nA);

    return Cbar;
}
boost::shared_ptr<Matrix> DFTSAPT::build_D(boost::shared_ptr<Matrix> L, boost::shared_ptr<Matrix> R)
{
    int nso = L->nrow();
    int nocc = L->ncol();

    boost::shared_ptr<Matrix> D(new Matrix("D", nso, nso));
    
    double** Dp = D->pointer();
    double** Lp = L->pointer();
    double** Rp = R->pointer();

    C_DGEMM('N','T',nso,nso,nocc,1.0,Lp[0],nocc,Rp[0],nocc,0.0,Dp[0],nso);

    return D;
}
boost::shared_ptr<Matrix> DFTSAPT::build_w(boost::shared_ptr<Matrix> W, boost::shared_ptr<Matrix> L, boost::shared_ptr<Matrix> R)
{
    int nso = L->nrow();
    int no = L->ncol();
    int nv = R->ncol();

    boost::shared_ptr<Matrix> w(new Matrix("w", no, nv));
    boost::shared_ptr<Matrix> T(new Matrix("T", no, nso));
    
    double** wp = w->pointer();
    double** Wp = W->pointer();
    double** Lp = L->pointer();
    double** Rp = R->pointer();
    double** Tp = T->pointer();

    C_DGEMM('T','N',no,nso,nso,1.0,Lp[0],no,Wp[0],nso,0.0,Tp[0],nso);
    C_DGEMM('N','N',no,nv,nso,1.0,Tp[0],nso,Rp[0],nv,0.0,wp[0],nv);

    return w;
}
boost::shared_ptr<Matrix> DFTSAPT::build_C_O(boost::shared_ptr<Matrix> C, boost::shared_ptr<Matrix> S, boost::shared_ptr<Matrix> P)
{
    int n = C->nrow();
    int no = C->ncol();

    boost::shared_ptr<Matrix> F(S->clone());
    boost::shared_ptr<Matrix> R(C->clone()); 

    double** Cp = C->pointer();
    double** Sp = S->pointer();
    double** Pp = P->pointer();
    double** Fp = F->pointer();
    double** Rp = R->pointer();

    C_DGEMM('N','N',n,n,n,1.0,Sp[0],n,Pp[0],n,0.0,Fp[0],n);
    C_DGEMM('T','N',n,no,n,1.0,Fp[0],n,Cp[0],no,0.0,Rp[0],no);
    
    return R;
}
boost::shared_ptr<Matrix> DFTSAPT::build_C_X(boost::shared_ptr<Matrix> x, boost::shared_ptr<Matrix> C)
{
    int n = C->nrow();
    int v = C->ncol();
    int o = x->nrow();
    
    boost::shared_ptr<Matrix> R(new Matrix("C_X", n, o));

    double** Cp = C->pointer();
    double** xp = x->pointer();
    double** Rp = R->pointer();

    C_DGEMM('N','T',n,o,v,1.0,Cp[0],v,xp[0],v,0.0,Rp[0],o);

    return R;
}
boost::shared_ptr<Matrix> DFTSAPT::triple(boost::shared_ptr<Matrix> A, boost::shared_ptr<Matrix> B, boost::shared_ptr<Matrix> C)
{
    int n = A->nrow();
    boost::shared_ptr<Matrix> D(A->clone());
    boost::shared_ptr<Matrix> E(A->clone());
    
    double** Ap = A->pointer();    
    double** Bp = B->pointer();    
    double** Cp = C->pointer();    
    double** Dp = D->pointer();    
    double** Ep = E->pointer();    

    C_DGEMM('N','N',n,n,n,1.0,Ap[0],n,Bp[0],n,0.0,Dp[0],n);
    C_DGEMM('N','N',n,n,n,1.0,Dp[0],n,Cp[0],n,0.0,Ep[0],n);

    return E;
}
std::pair<boost::shared_ptr<Matrix>, boost::shared_ptr<Matrix> > DFTSAPT::compute_x(boost::shared_ptr<JK> jk, boost::shared_ptr<Matrix> w_B, boost::shared_ptr<Matrix> w_A)
{
    boost::shared_ptr<CPKS_SAPT> cpks(new CPKS_SAPT);
   
    // Effective constructor 
    cpks->delta_ = cpks_delta_;
    cpks->maxiter_ = cpks_maxiter_;
    cpks->jk_ = jk;

    cpks->w_A_ = w_B; // Reversal of convention
    cpks->Cocc_A_ = Cocc_A_;
    cpks->Cvir_A_ = Cvir_A_;
    cpks->eps_occ_A_ = eps_occ_A_;
    cpks->eps_vir_A_ = eps_vir_A_;

    cpks->w_B_ = w_A; // Reversal of convention
    cpks->Cocc_B_ = Cocc_B_;
    cpks->Cvir_B_ = Cvir_B_;
    cpks->eps_occ_B_ = eps_occ_B_;
    cpks->eps_vir_B_ = eps_vir_B_;

    // Gogo CPKS
    cpks->compute_cpks();

    // Unpack
    std::pair<boost::shared_ptr<Matrix>, boost::shared_ptr<Matrix> > x_sol = make_pair(cpks->x_A_,cpks->x_B_);

    // Scale up for electron pairs
    x_sol.first->scale(2.0);
    x_sol.second->scale(2.0);

    return x_sol;
}

CPKS_SAPT::CPKS_SAPT()
{
}
CPKS_SAPT::~CPKS_SAPT()
{
}
void CPKS_SAPT::compute_cpks()
{
    // Allocate
    x_A_ = boost::shared_ptr<Matrix>(w_A_->clone());
    x_B_ = boost::shared_ptr<Matrix>(w_B_->clone());
    x_A_->zero();
    x_B_->zero();

    boost::shared_ptr<Matrix> r_A(w_A_->clone());
    boost::shared_ptr<Matrix> z_A(w_A_->clone());
    boost::shared_ptr<Matrix> p_A(w_A_->clone());
    boost::shared_ptr<Matrix> r_B(w_B_->clone());
    boost::shared_ptr<Matrix> z_B(w_B_->clone());
    boost::shared_ptr<Matrix> p_B(w_B_->clone());
        
    // Initialize (x_0 = 0)
    r_A->copy(w_A_);
    r_B->copy(w_B_);
    
    preconditioner(r_A,z_A,eps_occ_A_,eps_vir_A_);
    preconditioner(r_B,z_B,eps_occ_B_,eps_vir_B_);
    
    // Uncoupled value
    //fprintf(outfile, "(A<-B): %24.16E\n", -2.0 * z_A->vector_dot(w_A_));
    //fprintf(outfile, "(B<-A): %24.16E\n", -2.0 * z_B->vector_dot(w_B_));

    p_A->copy(z_A);
    p_B->copy(z_B);

    double zr_old_A = z_A->vector_dot(r_A);
    double zr_old_B = z_B->vector_dot(r_B);

    double r2A = 1.0;
    double r2B = 1.0;

    double b2A = sqrt(w_A_->vector_dot(w_A_));
    double b2B = sqrt(w_B_->vector_dot(w_B_));

    fprintf(outfile, "  ==> CPKS Iterations <==\n\n");

    fprintf(outfile, "    Maxiter     = %11d\n", maxiter_);
    fprintf(outfile, "    Convergence = %11.3E\n", delta_);
    fprintf(outfile, "\n");

    fprintf(outfile, "    ------------------------------\n");
    fprintf(outfile, "    %4s %11s  %11s \n", "Iter", "Monomer A", "Monomer B");
    fprintf(outfile, "    ------------------------------\n");
    fflush(outfile);

    int iter;
    for (iter = 0; iter < maxiter_; iter++) {
       
        std::map<std::string, boost::shared_ptr<Matrix> > b;
        if (r2A > delta_) {
            b["A"] = p_A;
        } 
        if (r2B > delta_) {
            b["B"] = p_B;
        } 
    
        std::map<std::string, boost::shared_ptr<Matrix> > s =
            product(b);
    
        if (r2A > delta_) {
            boost::shared_ptr<Matrix> s_A = s["A"];
            double alpha = r_A->vector_dot(z_A) / p_A->vector_dot(s_A);
            if (alpha < 0.0) {
                throw PSIEXCEPTION("Monomer A: A Matrix is not SPD");
            }
            int no = x_A_->nrow();
            int nv = x_A_->ncol();
            double** xp = x_A_->pointer();
            double** rp = r_A->pointer();
            double** pp = p_A->pointer();
            double** sp = s_A->pointer();
            C_DAXPY(no*nv, alpha,pp[0],1,xp[0],1);
            C_DAXPY(no*nv,-alpha,sp[0],1,rp[0],1); 
            r2A = sqrt(C_DDOT(no*nv,rp[0],1,rp[0],1)) / b2A; 
        } 
    
        if (r2B > delta_) {
            boost::shared_ptr<Matrix> s_B = s["B"];
            double alpha = r_B->vector_dot(z_B) / p_B->vector_dot(s_B);
            if (alpha < 0.0) {
                throw PSIEXCEPTION("Monomer B: A Matrix is not SPD");
            }
            int no = x_B_->nrow();
            int nv = x_B_->ncol();
            double** xp = x_B_->pointer();
            double** rp = r_B->pointer();
            double** pp = p_B->pointer();
            double** sp = s_B->pointer();
            C_DAXPY(no*nv, alpha,pp[0],1,xp[0],1);
            C_DAXPY(no*nv,-alpha,sp[0],1,rp[0],1); 
            r2B = sqrt(C_DDOT(no*nv,rp[0],1,rp[0],1)) / b2B; 
        } 

        fprintf(outfile, "    %4d %11.3E%1s %11.3E%1s\n", iter+1,
            r2A, (r2A < delta_ ? "*" : " "), 
            r2B, (r2B < delta_ ? "*" : " ")
            );
        fflush(outfile);

        if (r2A <= delta_ && r2B <= delta_) {
            break;
        }
        
        if (r2A > delta_) {
            preconditioner(r_A,z_A,eps_occ_A_,eps_vir_A_);
            double zr_new = z_A->vector_dot(r_A);
            double beta = zr_new / zr_old_A;
            zr_old_A = zr_new;
            int no = x_A_->nrow();
            int nv = x_A_->ncol();
            double** pp = p_A->pointer();
            double** zp = z_A->pointer();
            C_DSCAL(no*nv,beta,pp[0],1);
            C_DAXPY(no*nv,1.0,zp[0],1,pp[0],1);
        } 
        
        if (r2B > delta_) {
            preconditioner(r_B,z_B,eps_occ_B_,eps_vir_B_);
            double zr_new = z_B->vector_dot(r_B);
            double beta = zr_new / zr_old_B;
            zr_old_B = zr_new;
            int no = x_B_->nrow();
            int nv = x_B_->ncol();
            double** pp = p_B->pointer();
            double** zp = z_B->pointer();
            C_DSCAL(no*nv,beta,pp[0],1);
            C_DAXPY(no*nv,1.0,zp[0],1,pp[0],1);
        } 
    }
    
    fprintf(outfile, "    ------------------------------\n");
    fprintf(outfile, "\n");
    fflush(outfile);

    if (iter == maxiter_) 
        throw PSIEXCEPTION("CPKS did not converge.");
}
void CPKS_SAPT::preconditioner(boost::shared_ptr<Matrix> r,
                               boost::shared_ptr<Matrix> z,
                               boost::shared_ptr<Vector> o,
                               boost::shared_ptr<Vector> v)
{
    int no = o->dim();
    int nv = v->dim();
    
    double** rp = r->pointer();
    double** zp = z->pointer();

    double* op = o->pointer();
    double* vp = v->pointer();

    for (int i = 0; i < no; i++) {
        for (int a = 0; a < nv; a++) {
            zp[i][a] = rp[i][a] / (vp[a] - op[i]);
        }
    }
}
std::map<std::string, boost::shared_ptr<Matrix> > CPKS_SAPT::product(std::map<std::string, boost::shared_ptr<Matrix> > b)
{
    std::map<std::string, boost::shared_ptr<Matrix> > s;

    bool do_A = b.count("A");
    bool do_B = b.count("B");

    std::vector<SharedMatrix>& Cl = jk_->C_left();
    std::vector<SharedMatrix>& Cr = jk_->C_right();
    Cl.clear();
    Cr.clear();

    if (do_A) {
        Cl.push_back(Cocc_A_);    
        int no = b["A"]->nrow();
        int nv = b["A"]->ncol();
        int nso = Cvir_A_->nrow();
        double** Cp = Cvir_A_->pointer();
        double** bp = b["A"]->pointer();
        boost::shared_ptr<Matrix> T(new Matrix("T",nso,no));
        double** Tp = T->pointer();
        C_DGEMM('N','T',nso,no,nv,1.0,Cp[0],nv,bp[0],nv,0.0,Tp[0],no);
        Cr.push_back(T);
    }

    if (do_B) {
        Cl.push_back(Cocc_B_);    
        int no = b["B"]->nrow();
        int nv = b["B"]->ncol();
        int nso = Cvir_B_->nrow();
        double** Cp = Cvir_B_->pointer();
        double** bp = b["B"]->pointer();
        boost::shared_ptr<Matrix> T(new Matrix("T",nso,no));
        double** Tp = T->pointer();
        C_DGEMM('N','T',nso,no,nv,1.0,Cp[0],nv,bp[0],nv,0.0,Tp[0],no);
        Cr.push_back(T);
    }

    jk_->compute();
    
    const std::vector<SharedMatrix>& J = jk_->J();
    const std::vector<SharedMatrix>& K = jk_->K();

    int indA = 0;
    int indB = (do_A ? 1 : 0);

    if (do_A) {
        boost::shared_ptr<Matrix> Jv = J[indA];
        boost::shared_ptr<Matrix> Kv = K[indA];
        Jv->scale(4.0);
        Jv->subtract(Kv);
        Jv->subtract(Kv->transpose());
        
        int no = b["A"]->nrow();
        int nv = b["A"]->ncol();
        int nso = Cvir_A_->nrow();
        boost::shared_ptr<Matrix> T(new Matrix("T", no, nso));
        s["A"] = boost::shared_ptr<Matrix>(new Matrix("S", no, nv));
        double** Cop = Cocc_A_->pointer();
        double** Cvp = Cvir_A_->pointer();
        double** Jp = Jv->pointer();
        double** Tp = T->pointer();
        double** Sp = s["A"]->pointer();
        C_DGEMM('T','N',no,nso,nso,1.0,Cop[0],no,Jp[0],nso,0.0,Tp[0],nso);
        C_DGEMM('N','N',no,nv,nso,1.0,Tp[0],nso,Cvp[0],nv,0.0,Sp[0],nv); 

        double** bp = b["A"]->pointer();
        double* op = eps_occ_A_->pointer();
        double* vp = eps_vir_A_->pointer();
        for (int i = 0; i < no; i++) {
            for (int a = 0; a < nv; a++) {
                Sp[i][a] += bp[i][a] * (vp[a] - op[i]);
            }
        }
    } 

    if (do_B) {
        boost::shared_ptr<Matrix> Jv = J[indB];
        boost::shared_ptr<Matrix> Kv = K[indB];
        Jv->scale(4.0);
        Jv->subtract(Kv);
        Jv->subtract(Kv->transpose());
        
        int no = b["B"]->nrow();
        int nv = b["B"]->ncol();
        int nso = Cvir_B_->nrow();
        boost::shared_ptr<Matrix> T(new Matrix("T", no, nso));
        s["B"] = boost::shared_ptr<Matrix>(new Matrix("S", no, nv));
        double** Cop = Cocc_B_->pointer();
        double** Cvp = Cvir_B_->pointer();
        double** Jp = Jv->pointer();
        double** Tp = T->pointer();
        double** Sp = s["B"]->pointer();
        C_DGEMM('T','N',no,nso,nso,1.0,Cop[0],no,Jp[0],nso,0.0,Tp[0],nso);
        C_DGEMM('N','N',no,nv,nso,1.0,Tp[0],nso,Cvp[0],nv,0.0,Sp[0],nv); 

        double** bp = b["B"]->pointer();
        double* op = eps_occ_B_->pointer();
        double* vp = eps_vir_B_->pointer();
        for (int i = 0; i < no; i++) {
            for (int a = 0; a < nv; a++) {
                Sp[i][a] += bp[i][a] * (vp[a] - op[i]);
            }
        }
    } 

    return s;
}

}}
