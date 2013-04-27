#include "dftsapt.h"
#include <libmints/mints.h>
#include <libmints/sieve.h>
#include <libfock/jk.h>
#include <libthce/thce.h>
#include <libthce/lreri.h>
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

    W_A.reset();
    W_B.reset();

    // ==> Exchange Terms (S^2) <== //
    
    double Exch10_2 = 0.0;

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

    // ==> Uncorrelated Second-Order Response Terms [Induction] <== //

    // Compute CPKS
    std::pair<boost::shared_ptr<Matrix>, boost::shared_ptr<Matrix> > x_sol = compute_x(jk,w_B,w_A);
    boost::shared_ptr<Matrix> x_A = x_sol.first;
    boost::shared_ptr<Matrix> x_B = x_sol.second;

    // Backward in Ed's convention
    x_A->scale(-1.0);
    x_B->scale(-1.0);

    // ==> Induction <== //

    double Ind20r_AB = 2.0 * x_A->vector_dot(w_B);
    double Ind20r_BA = 2.0 * x_B->vector_dot(w_A);
    double Ind20r = Ind20r_AB + Ind20r_BA;

    energies_["Ind20,r (A<-B)"] = Ind20r_AB;
    energies_["Ind20,r (B->A)"] = Ind20r_BA;
    energies_["Ind20,r"] = Ind20r_AB + Ind20r_BA;

    fprintf(outfile,"    Ind20,r (A<-B)      = %18.12lf H\n",Ind20r_AB);
    fprintf(outfile,"    Ind20,r (A->B)      = %18.12lf H\n",Ind20r_BA);
    fprintf(outfile,"    Ind20,r             = %18.12lf H\n",Ind20r);
    fflush(outfile);

    // ==> Exchange-Induction <== //

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
    boost::shared_ptr<Matrix> J_X_A    = J[1]; 
    boost::shared_ptr<Matrix> J_X_B    = J[2]; 

    boost::shared_ptr<Matrix> K_O      = K[0]; 
    boost::shared_ptr<Matrix> K_X_A    = K[1]; 
    boost::shared_ptr<Matrix> K_X_B    = K[2]; 

    // Clear the JK memory 
    jk.reset();
    
    // New generalized densities
    boost::shared_ptr<Matrix> X_A = build_D(Cocc_A_,C_X_A);
    boost::shared_ptr<Matrix> X_B = build_D(Cocc_B_,C_X_B);

    // Don't need these, in AO basis now
    C_O_A.reset(); 
    C_X_A.reset(); 
    C_X_B.reset(); 

    // Exch-Ind20 (A<-B)
    double ExchInd20r_AB = 0.0;
    std::vector<double> ExchInd20r_AB_terms;
    ExchInd20r_AB_terms.resize(18);
    ExchInd20r_AB_terms[0]  -= 2.0 * X_A->vector_dot(K_B);
    ExchInd20r_AB_terms[1]  -= 4.0 * triple(D_B,S,X_A)->vector_dot(J_A);
    ExchInd20r_AB_terms[2]  += 2.0 * triple(D_A,S,D_B)->vector_dot(K_X_A);
    ExchInd20r_AB_terms[3]  -= 4.0 * triple(D_A,S,D_B)->vector_dot(J_X_A);
    ExchInd20r_AB_terms[4]  += 2.0 * triple(D_B,S,X_A)->vector_dot(K_A); 
    ExchInd20r_AB_terms[5]  -= 4.0 * triple(X_A,S,D_B)->vector_dot(J_B);
    ExchInd20r_AB_terms[6]  += 2.0 * triple(X_A,S,D_B)->vector_dot(K_B);
    ExchInd20r_AB_terms[7]  += 4.0 * triple(triple(D_B,S,X_A),S,D_B)->vector_dot(J_A);
    ExchInd20r_AB_terms[8]  += 4.0 * triple(triple(X_A,S,D_B),S,D_A)->vector_dot(J_B);
    ExchInd20r_AB_terms[9]  -= 2.0 * triple(X_A,S,D_B)->vector_dot(K_O);
    ExchInd20r_AB_terms[10] += 4.0 * triple(triple(D_B,S,D_A),S,D_B)->vector_dot(J_X_A);
    ExchInd20r_AB_terms[11] += 4.0 * triple(triple(D_A,S,D_B),S,X_A)->vector_dot(J_B);
    ExchInd20r_AB_terms[12] -= 2.0 * triple(D_B,S,X_A)->vector_dot(K_O->transpose());
    ExchInd20r_AB_terms[13] -= 2.0 * triple(D_B,S,X_A)->vector_dot(V_A); 
    ExchInd20r_AB_terms[14] -= 2.0 * triple(X_A,S,D_B)->vector_dot(V_B);
    ExchInd20r_AB_terms[15] += 2.0 * triple(triple(D_B,S,X_A),S,D_B)->vector_dot(V_A);
    ExchInd20r_AB_terms[16] += 2.0 * triple(triple(X_A,S,D_B),S,D_A)->vector_dot(V_B);
    ExchInd20r_AB_terms[17] += 2.0 * triple(triple(D_A,S,D_B),S,X_A)->vector_dot(V_B);
    for (int k = 0; k < ExchInd20r_AB_terms.size(); k++) {
        ExchInd20r_AB += ExchInd20r_AB_terms[k];
    }
    if (debug_) {
        for (int k = 0; k < ExchInd20r_AB_terms.size(); k++) {
            fprintf(outfile,"    Exch-Ind20,r (%2d)   = %18.12lf H\n",k+1,ExchInd20r_AB_terms[k]);
        }
    }
    fprintf(outfile,"    Exch-Ind20,r (A<-B) = %18.12lf H\n",ExchInd20r_AB);
    fflush(outfile);

    K_O->transpose_this();

    // Exch-Ind20 (B<-A)
    double ExchInd20r_BA = 0.0;
    std::vector<double> ExchInd20r_BA_terms;
    ExchInd20r_BA_terms.resize(18);
    ExchInd20r_BA_terms[0]  -= 2.0 * X_B->vector_dot(K_A);
    ExchInd20r_BA_terms[1]  -= 4.0 * triple(D_A,S,X_B)->vector_dot(J_B);
    ExchInd20r_BA_terms[2]  += 2.0 * triple(D_B,S,D_A)->vector_dot(K_X_B);
    ExchInd20r_BA_terms[3]  -= 4.0 * triple(D_B,S,D_A)->vector_dot(J_X_B);
    ExchInd20r_BA_terms[4]  += 2.0 * triple(D_A,S,X_B)->vector_dot(K_B); 
    ExchInd20r_BA_terms[5]  -= 4.0 * triple(X_B,S,D_A)->vector_dot(J_A);
    ExchInd20r_BA_terms[6]  += 2.0 * triple(X_B,S,D_A)->vector_dot(K_A);
    ExchInd20r_BA_terms[7]  += 4.0 * triple(triple(D_A,S,X_B),S,D_A)->vector_dot(J_B);
    ExchInd20r_BA_terms[8]  += 4.0 * triple(triple(X_B,S,D_A),S,D_B)->vector_dot(J_A);
    ExchInd20r_BA_terms[9]  -= 2.0 * triple(X_B,S,D_A)->vector_dot(K_O);
    ExchInd20r_BA_terms[10] += 4.0 * triple(triple(D_A,S,D_B),S,D_A)->vector_dot(J_X_B);
    ExchInd20r_BA_terms[11] += 4.0 * triple(triple(D_B,S,D_A),S,X_B)->vector_dot(J_A);
    ExchInd20r_BA_terms[12] -= 2.0 * triple(D_A,S,X_B)->vector_dot(K_O->transpose());
    ExchInd20r_BA_terms[13] -= 2.0 * triple(D_A,S,X_B)->vector_dot(V_B); 
    ExchInd20r_BA_terms[14] -= 2.0 * triple(X_B,S,D_A)->vector_dot(V_A);
    ExchInd20r_BA_terms[15] += 2.0 * triple(triple(D_A,S,X_B),S,D_A)->vector_dot(V_B);
    ExchInd20r_BA_terms[16] += 2.0 * triple(triple(X_B,S,D_A),S,D_B)->vector_dot(V_A);
    ExchInd20r_BA_terms[17] += 2.0 * triple(triple(D_B,S,D_A),S,X_B)->vector_dot(V_A);
    for (int k = 0; k < ExchInd20r_BA_terms.size(); k++) {
        ExchInd20r_BA += ExchInd20r_BA_terms[k];
    }
    if (debug_) {
        for (int k = 0; k < ExchInd20r_BA_terms.size(); k++) {
            fprintf(outfile,"    Exch-Ind20,r (%2d)   = %18.12lf H\n",k+1,ExchInd20r_BA_terms[k]);
        }
    }
    fprintf(outfile,"    Exch-Ind20,r (B<-A) = %18.12lf H\n",ExchInd20r_BA);
    fflush(outfile);

    K_O->transpose_this();

    double ExchInd20r = ExchInd20r_AB + ExchInd20r_BA;
    fprintf(outfile,"    Exch-Ind20,r        = %18.12lf H\n",ExchInd20r);
    fprintf(outfile,"\n");
    fflush(outfile);

    energies_["Exch-Ind20,r (A<-B)"] = ExchInd20r_AB;
    energies_["Exch-Ind20,r (B->A)"] = ExchInd20r_BA;
    energies_["Exch-Ind20,r"] = ExchInd20r_AB + ExchInd20r_BA;


    vars_["S"]   = S;
    vars_["D_A"] = D_A;
    vars_["V_A"] = V_A;
    vars_["J_A"] = J_A;
    vars_["K_A"] = K_A;
    vars_["D_B"] = D_B;
    vars_["V_B"] = V_B;
    vars_["J_B"] = J_B;
    vars_["K_B"] = K_B;
    vars_["K_O"] = K_O;
}
void DFTSAPT::mp2_terms()
{
    fprintf(outfile, "  PT2 TERMS:\n\n");

    // => Sizing <= //

    int nn = primary_->nbf();

    int na = Caocc_A_->colspi()[0];
    int nb = Caocc_B_->colspi()[0];
    int nr = Cavir_A_->colspi()[0];
    int ns = Cavir_B_->colspi()[0];
    int nQ = mp2fit_->nbf();
    size_t nrQ = nr * (size_t) nQ;
    size_t nsQ = ns * (size_t) nQ;
    
    int nT = 1;
    #ifdef _OPENMP
        nT = omp_get_max_threads();
    #endif

    // => Stashed Variables <= //

    boost::shared_ptr<Matrix> S   = vars_["S"];
    boost::shared_ptr<Matrix> D_A = vars_["D_A"];
    boost::shared_ptr<Matrix> V_A = vars_["V_A"];
    boost::shared_ptr<Matrix> J_A = vars_["J_A"];
    boost::shared_ptr<Matrix> K_A = vars_["K_A"];
    boost::shared_ptr<Matrix> D_B = vars_["D_B"];
    boost::shared_ptr<Matrix> V_B = vars_["V_B"];
    boost::shared_ptr<Matrix> J_B = vars_["J_B"];
    boost::shared_ptr<Matrix> K_B = vars_["K_B"];
    boost::shared_ptr<Matrix> K_O = vars_["K_O"];

    // => Auxiliary C matrices <= //

    boost::shared_ptr<Matrix> Cr1 = triple(D_B,S,Cavir_A_);
    Cr1->scale(-1.0);
    Cr1->add(Cavir_A_);
    boost::shared_ptr<Matrix> Cs1 = triple(D_A,S,Cavir_B_);
    Cs1->scale(-1.0);
    Cs1->add(Cavir_B_);
    boost::shared_ptr<Matrix> Ca2 = triple(D_B,S,Caocc_A_);
    boost::shared_ptr<Matrix> Cb2 = triple(D_A,S,Caocc_B_);
    boost::shared_ptr<Matrix> Cr3 = triple(D_B,S,Cavir_A_);
    boost::shared_ptr<Matrix> CrX = triple(triple(D_A,S,D_B),S,Cavir_A_);
    Cr3->subtract(CrX);
    Cr3->scale(2.0);
    boost::shared_ptr<Matrix> Cs3 = triple(D_A,S,Cavir_B_);
    boost::shared_ptr<Matrix> CsX = triple(triple(D_B,S,D_A),S,Cavir_B_);
    Cs3->subtract(CsX);
    Cs3->scale(2.0);
    boost::shared_ptr<Matrix> Ca4 = triple(triple(D_A,S,D_B),S,Caocc_A_);
    Ca4->scale(-2.0);    
    boost::shared_ptr<Matrix> Cb4 = triple(triple(D_B,S,D_A),S,Caocc_B_);
    Cb4->scale(-2.0);    

    // => Auxiliary V matrices <= //

    // TODO

    S.reset();
    D_A.reset();
    V_A.reset();
    J_A.reset();
    K_A.reset();
    D_B.reset();
    V_B.reset();
    J_B.reset();
    K_B.reset();
    K_O.reset();

    // => Memory <= //
    
    // TODO: Account for 2-index overhead in memory
    long int mem = memory_;
    mem -= 6L * nn * nn;

    // => Integrals from the THCE <= //

    boost::shared_ptr<DFERI> df = DFERI::build(primary_,mp2fit_,Process::environment.options);
    df->set_memory(mem);
    df->clear();

    std::vector<boost::shared_ptr<Matrix> > Cs;
    Cs.push_back(Caocc_A_);  
    Cs.push_back(Cavir_A_);  
    Cs.push_back(Caocc_B_);  
    Cs.push_back(Cavir_B_);  
    Cs.push_back(Cr1);  
    Cs.push_back(Cs1);  
    Cs.push_back(Ca2);  
    Cs.push_back(Cb2);  
    Cs.push_back(Cr3);  
    Cs.push_back(Cs3);  
    Cs.push_back(Ca4);  
    Cs.push_back(Cb4);  
    boost::shared_ptr<Matrix> Call = Matrix::horzcat(Cs);
    df->set_C(Call);

    int offset = 0;
    df->add_space("a",offset,offset+Caocc_A_->colspi()[0]); offset += Caocc_A_->colspi()[0];
    df->add_space("r",offset,offset+Cavir_A_->colspi()[0]); offset += Cavir_A_->colspi()[0];
    df->add_space("b",offset,offset+Caocc_B_->colspi()[0]); offset += Caocc_B_->colspi()[0];
    df->add_space("s",offset,offset+Cavir_B_->colspi()[0]); offset += Cavir_B_->colspi()[0];
    df->add_space("r1",offset,offset+Cr1->colspi()[0]); offset += Cr1->colspi()[0];
    df->add_space("s1",offset,offset+Cs1->colspi()[0]); offset += Cs1->colspi()[0];
    df->add_space("a2",offset,offset+Ca2->colspi()[0]); offset += Ca2->colspi()[0];
    df->add_space("b2",offset,offset+Cb2->colspi()[0]); offset += Cb2->colspi()[0];
    df->add_space("r3",offset,offset+Cr3->colspi()[0]); offset += Cr3->colspi()[0];
    df->add_space("s3",offset,offset+Cs3->colspi()[0]); offset += Cs3->colspi()[0];
    df->add_space("a4",offset,offset+Ca4->colspi()[0]); offset += Ca4->colspi()[0];
    df->add_space("b4",offset,offset+Cb4->colspi()[0]); offset += Cb4->colspi()[0];

    df->add_pair_space("Aar", "a", "r");
    df->add_pair_space("Abs", "b", "s"); 
    df->add_pair_space("Bas", "a", "s1");
    df->add_pair_space("Bbr", "b", "r1"); 
    df->add_pair_space("Cas", "a2","s");
    df->add_pair_space("Cbr", "b2","r"); 
    df->add_pair_space("Dar", "a", "r3");
    df->add_pair_space("Dbs", "b", "s3"); 
    df->add_pair_space("Ear", "a4","r");
    df->add_pair_space("Ebs", "b4","s"); 

    df->compute();

    std::map<std::string, boost::shared_ptr<Tensor> >& ints = df->ints();

    boost::shared_ptr<Tensor> AarT = ints["Aar"];
    boost::shared_ptr<Tensor> AbsT = ints["Abs"];
    boost::shared_ptr<Tensor> BasT = ints["Bas"];
    boost::shared_ptr<Tensor> BbrT = ints["Bbr"];
    boost::shared_ptr<Tensor> CasT = ints["Cas"];
    boost::shared_ptr<Tensor> CbrT = ints["Cbr"];
    boost::shared_ptr<Tensor> DarT = ints["Dar"];
    boost::shared_ptr<Tensor> DbsT = ints["Dbs"];
    boost::shared_ptr<Tensor> EarT = ints["Ear"];
    boost::shared_ptr<Tensor> EbsT = ints["Ebs"];

    df.reset();

    Cr1.reset();
    Cs1.reset();
    Ca2.reset();
    Cb2.reset();
    Cr3.reset();
    Cs3.reset();
    Ca4.reset();
    Cb4.reset();

    // => Blocking <= //

    long int overhead = 0L;
    overhead += 2L * nT * nr * ns;

    long int rem = mem - overhead;

    if (rem < 0L) {
        throw PSIEXCEPTION("Too little static memory for DFTSAPT::mp2_terms");
    }

    long int cost_a = 2L * nr * nQ + 2L * ns * nQ;
    long int max_a = rem / (2L * cost_a);
    long int max_b = max_a;
    max_a = (max_a > na ? na : max_a);
    max_b = (max_b > nb ? nb : max_b);
    if (max_a < 1L || max_b < 1L) {
        throw PSIEXCEPTION("Too little dynamic memory for DFTSAPT::mp2_terms");
    }

    // => Tensor Slices <= //

    boost::shared_ptr<Matrix> Aar(new Matrix("Aar",max_a*nr,nQ));
    boost::shared_ptr<Matrix> Abs(new Matrix("Abs",max_b*ns,nQ));
    boost::shared_ptr<Matrix> Bas(new Matrix("Bas",max_a*ns,nQ));
    boost::shared_ptr<Matrix> Bbr(new Matrix("Bbr",max_b*nr,nQ));
    boost::shared_ptr<Matrix> Cas(new Matrix("Cas",max_a*ns,nQ));
    boost::shared_ptr<Matrix> Cbr(new Matrix("Cbr",max_b*nr,nQ));
    boost::shared_ptr<Matrix> Dar(new Matrix("Dar",max_a*nr,nQ));
    boost::shared_ptr<Matrix> Dbs(new Matrix("Dbs",max_b*ns,nQ));

    // => Thread Work Arrays <= //

    std::vector<boost::shared_ptr<Matrix> > Trs;
    std::vector<boost::shared_ptr<Matrix> > Vrs;
    for (int t = 0; t < nT; t++) {
        Trs.push_back(boost::shared_ptr<Matrix>(new Matrix("Trs",nr,ns)));
        Vrs.push_back(boost::shared_ptr<Matrix>(new Matrix("Vrs",nr,ns)));
    }

    // => Pointers <= //

    double** Aarp = Aar->pointer();
    double** Absp = Abs->pointer();
    double** Basp = Bas->pointer();
    double** Bbrp = Bbr->pointer();
    double** Casp = Cas->pointer();
    double** Cbrp = Cbr->pointer();
    double** Darp = Dar->pointer();
    double** Dbsp = Dbs->pointer();
    double*  eap  = eps_aocc_A_->pointer();
    double*  ebp  = eps_aocc_B_->pointer();
    double*  erp  = eps_avir_A_->pointer();
    double*  esp  = eps_avir_B_->pointer();

    // => File Pointers <= //
    
    FILE* Aarf = AarT->file_pointer();
    FILE* Absf = AbsT->file_pointer();
    FILE* Basf = BasT->file_pointer();
    FILE* Bbrf = BbrT->file_pointer();
    FILE* Casf = CasT->file_pointer();
    FILE* Cbrf = CbrT->file_pointer();
    FILE* Darf = DarT->file_pointer();
    FILE* Dbsf = DbsT->file_pointer();
    FILE* Earf = EarT->file_pointer();
    FILE* Ebsf = EbsT->file_pointer();

    // => Slice D + E -> D <= //

    fseek(Darf,0L,SEEK_SET);
    fseek(Earf,0L,SEEK_SET);
    for (int astart = 0; astart < na; astart += max_a) {
        int nablock = (astart + max_a >= na ? na - astart : max_a);
        fread(Darp[0],sizeof(double),nablock*nrQ,Darf);
        fread(Aarp[0],sizeof(double),nablock*nrQ,Earf);
        C_DAXPY(nablock*nrQ,1.0,Aarp[0],1,Darp[0],1);
        fseek(Darf,sizeof(double)*astart*nrQ,SEEK_SET);
        fwrite(Darp[0],sizeof(double),nablock*nrQ,Darf); 
    }

    fseek(Dbsf,0L,SEEK_SET);
    fseek(Ebsf,0L,SEEK_SET);
    for (int bstart = 0; bstart < nb; bstart += max_b) {
        int nbblock = (bstart + max_b >= nb ? nb - bstart : max_b);
        fread(Dbsp[0],sizeof(double),nbblock*nsQ,Dbsf);
        fread(Absp[0],sizeof(double),nbblock*nsQ,Ebsf);
        C_DAXPY(nbblock*nsQ,1.0,Absp[0],1,Dbsp[0],1);
        fseek(Dbsf,sizeof(double)*bstart*nsQ,SEEK_SET);
        fwrite(Dbsp[0],sizeof(double),nbblock*nsQ,Dbsf); 
    }

    // => Targets <= //

    double Disp20 = 0.0;
    double ExchDisp20 = 0.0;

    // ==> Master Loop <== //

    fseek(Aarf,0L,SEEK_SET);
    fseek(Basf,0L,SEEK_SET);
    fseek(Casf,0L,SEEK_SET);
    fseek(Darf,0L,SEEK_SET);
    for (int astart = 0; astart < na; astart += max_a) {
        int nablock = (astart + max_a >= na ? na - astart : max_a);

        fread(Aarp[0],sizeof(double),nablock*nrQ,Aarf);
        fread(Basp[0],sizeof(double),nablock*nsQ,Basf);
        fread(Casp[0],sizeof(double),nablock*nsQ,Casf);
        fread(Darp[0],sizeof(double),nablock*nrQ,Darf);

        fseek(Absf,0L,SEEK_SET);
        fseek(Bbrf,0L,SEEK_SET);
        fseek(Cbrf,0L,SEEK_SET);
        fseek(Dbsf,0L,SEEK_SET);
        for (int bstart = 0; bstart < nb; bstart += max_b) {
            int nbblock = (bstart + max_b >= nb ? nb - bstart : max_b);

            fread(Absp[0],sizeof(double),nbblock*nsQ,Absf);
            fread(Bbrp[0],sizeof(double),nbblock*nrQ,Bbrf);
            fread(Cbrp[0],sizeof(double),nbblock*nrQ,Cbrf);
            fread(Dbsp[0],sizeof(double),nbblock*nsQ,Dbsf);

            long int nab = nablock * nbblock;
                
            #pragma omp parallel for schedule(dynamic) reduction(+: Disp20, ExchDisp20)
            for (long int ab = 0L; ab < nab; ab++) {
                int a = ab / nbblock;
                int b = ab % nbblock;

                int thread = 0;
                #ifdef _OPENMP
                    thread = omp_get_thread_num();
                #endif

                double** Trsp = Trs[thread]->pointer();
                double** Vrsp = Vrs[thread]->pointer();

                // => Amplitudes, Disp20 <= //

                C_DGEMM('N','T',nr,ns,nQ,1.0,Aarp[(a + astart)*nr],nQ,Absp[(b + bstart)*ns],nQ,0.0,Vrsp[0],ns);

                for (int r = 0; r < nr; r++) {
                    for (int s = 0; s < ns; s++) {
                        Trsp[r][s] = Vrsp[r][s] / (eap[a + astart] + ebp[b + bstart] - erp[r] - esp[s]);
                        Disp20 += 4.0 * Trsp[r][s] * Vrsp[r][s];
                    }
                }

                // => Q1-Q3 <= //

                C_DGEMM('N','T',nr,ns,nQ,1.0,Bbrp[(b + bstart)*nr],nQ,Basp[(a + astart)*ns],nQ,0.0,Vrsp[0],ns);
                C_DGEMM('N','T',nr,ns,nQ,1.0,Cbrp[(b + bstart)*nr],nQ,Casp[(a + astart)*ns],nQ,1.0,Vrsp[0],ns);
                C_DGEMM('N','T',nr,ns,nQ,1.0,Aarp[(a + astart)*nr],nQ,Dbsp[(b + bstart)*ns],nQ,1.0,Vrsp[0],ns);
                C_DGEMM('N','T',nr,ns,nQ,1.0,Darp[(a + astart)*nr],nQ,Absp[(b + bstart)*ns],nQ,1.0,Vrsp[0],ns);

                for (int r = 0; r < nr; r++) {
                    for (int s = 0; s < ns; s++) {
                        ExchDisp20 -= 2.0 * Trsp[r][s] * Vrsp[r][s];
                    }
                }

                // TODO V,J,K <= //          

            }
        } 
    }

    fprintf(outfile,"    Disp20              = %18.12lf H\n",Disp20);
    fprintf(outfile,"    Exch-Disp20         = %18.12lf H\n",ExchDisp20);
    fprintf(outfile,"\n");
    fflush(outfile);
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
    int Ar = A->nrow();
    int Ac = A->ncol();
    int Br = B->nrow();
    int Bc = B->ncol();
    int Cr = C->nrow();
    int Cc = C->ncol();

    if (Ac != Br) throw PSIEXCEPTION("Nonconforming GEMM");
    if (Bc != Cr) throw PSIEXCEPTION("Nonconforming GEMM");

    boost::shared_ptr<Matrix> D(new Matrix("T",Ar,Bc));
    boost::shared_ptr<Matrix> E(new Matrix("T",Ar,Cc));
    
    double** Ap = A->pointer();    
    double** Bp = B->pointer();    
    double** Cp = C->pointer();    
    double** Dp = D->pointer();    
    double** Ep = E->pointer();    

    C_DGEMM('N','N',Ar,Bc,Ac,1.0,Ap[0],Ac,Bp[0],Bc,0.0,Dp[0],Bc);
    C_DGEMM('N','N',Ar,Cc,Bc,1.0,Dp[0],Bc,Cp[0],Cc,0.0,Ep[0],Cc);

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
