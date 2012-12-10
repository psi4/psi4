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

    sapt->memory_ = (unsigned long int)(Process::environment.get_memory() * options.get_double("SAPT_MEMORY_FACTOR") * 0.125);

    sapt->dimer_     = d->molecule();
    sapt->monomer_A_ = mA->molecule();
    sapt->monomer_B_ = mB->molecule();

    sapt->E_dimer_     = d->reference_energy();
    sapt->E_monomer_A_ = mA->reference_energy();
    sapt->E_monomer_B_ = mB->reference_energy();

    sapt->primary_   = d->basisset();
    sapt->primary_A_ = mA->basisset();
    sapt->primary_B_ = mB->basisset();

    sapt->Cocc_A_     = mA->Ca_subset("AO","OCC");

    sapt->Caocc_A_    = mA->Ca_subset("AO","AOCC");
    sapt->Cavir_A_    = mA->Ca_subset("AO","AVIR");

    sapt->eps_focc_A_ = mA->epsilon_a_subset("AO","FOCC");
    sapt->eps_aocc_A_ = mA->epsilon_a_subset("AO","AOCC");
    sapt->eps_avir_A_ = mA->epsilon_a_subset("AO","AVIR");
    sapt->eps_fvir_A_ = mA->epsilon_a_subset("AO","FVIR");

    sapt->Cocc_B_     = mB->Ca_subset("AO","OCC");

    sapt->Caocc_B_    = mB->Ca_subset("AO","AOCC");
    sapt->Cavir_B_    = mB->Ca_subset("AO","AVIR");

    sapt->eps_focc_B_ = mB->epsilon_a_subset("AO","FOCC");
    sapt->eps_aocc_B_ = mB->epsilon_a_subset("AO","AOCC");
    sapt->eps_avir_B_ = mB->epsilon_a_subset("AO","AVIR");
    sapt->eps_fvir_B_ = mB->epsilon_a_subset("AO","FVIR");

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

    return 0.0;
}
void DFTSAPT::print_header() const 
{
    fprintf(outfile, "\t --------------------------------------------------------\n");
    fprintf(outfile, "\t                        DF-DFT-SAPT                      \n");
    fprintf(outfile, "\t               Rob Parrish and Ed Hohenstein             \n");
    fprintf(outfile, "\t --------------------------------------------------------\n");
    fprintf(outfile, "\n");
}
void DFTSAPT::fock_terms()
{
    // ==> Setup <== //    

    boost::shared_ptr<JK> jk = JK::build_JK();
    // TODO: Account for 2-index overhead in memory
    jk->set_memory(memory_);
    jk->set_do_J(true);
    jk->set_do_K(true);
    jk->initialize();
    jk->print_header();

    // ==> Generalized Fock Source Terms [Elst/Exch] <== //

    // => Steric Interaction Density Terms (T) <= //

    boost::shared_ptr<Matrix> S = build_S(primary_);
    boost::shared_ptr<Matrix> Sij = build_Sij(S);
    boost::shared_ptr<Matrix> Sij_2 = build_Sij_2(Sij);
    boost::shared_ptr<Matrix> Sij_n = build_Sij_n(Sij);

    std::map<std::string, boost::shared_ptr<Matrix> > Cbar_2 = build_Cbar(Sij_2);
    boost::shared_ptr<Matrix> C_T_A_2  = Cbar_2["C_T_A"];
    boost::shared_ptr<Matrix> C_T_B_2  = Cbar_2["C_T_B"];
    boost::shared_ptr<Matrix> C_T_AB_2 = Cbar_2["C_T_AB"];
    boost::shared_ptr<Matrix> C_T_BA_2 = Cbar_2["C_T_BA"];

    std::map<std::string, boost::shared_ptr<Matrix> > Cbar_n = build_Cbar(Sij_n);
    boost::shared_ptr<Matrix> C_T_A_n  = Cbar_n["C_T_A"];
    boost::shared_ptr<Matrix> C_T_B_n  = Cbar_n["C_T_B"];
    boost::shared_ptr<Matrix> C_T_AB_n = Cbar_n["C_T_AB"];
    boost::shared_ptr<Matrix> C_T_BA_n = Cbar_n["C_T_BA"];

    if (debug_ > 2) {
        S->print();
        Sij->print();
        Sij_2->print();
        Sij_n->print();
    }

    Sij.reset();
    Sij_2.reset();
    Sij_n.reset();

    // => Load the JK Object <= //
    
    std::vector<SharedMatrix>& Cl = jk->C_left();
    std::vector<SharedMatrix>& Cr = jk->C_right();
    Cl.clear();
    Cr.clear();

    // J/K[P^A]
    Cl.push_back(Cocc_A_);
    Cr.push_back(Cocc_A_);
    // J/K[T^A, S^2]
    Cl.push_back(Cocc_A_);
    Cr.push_back(C_T_A_2);
    // J/K[T^A, S^\infty]
    Cl.push_back(Cocc_A_);
    Cr.push_back(C_T_A_n);
    // J/K[P^B]
    Cl.push_back(Cocc_B_);
    Cr.push_back(Cocc_B_);
    // J/K[T^BA, S^2]
    Cl.push_back(Cocc_B_);
    Cr.push_back(C_T_BA_2);
    // J/K[T^BA, S^\infty]
    Cl.push_back(Cocc_B_);
    Cr.push_back(C_T_BA_n);
    
    // => Compute the JK matrices <= //

    jk->compute();

    // => Unload the JK Object <= //

    const std::vector<SharedMatrix>& J = jk->J();
    const std::vector<SharedMatrix>& K = jk->K();

    boost::shared_ptr<Matrix> J_A      = J[0]; 
    boost::shared_ptr<Matrix> J_T_A_2  = J[1]; 
    boost::shared_ptr<Matrix> J_T_A_n  = J[2]; 
    boost::shared_ptr<Matrix> J_B      = J[3]; 
    boost::shared_ptr<Matrix> J_T_BA_2 = J[4]; 
    boost::shared_ptr<Matrix> J_T_BA_n = J[5]; 

    boost::shared_ptr<Matrix> K_A      = K[0]; 
    boost::shared_ptr<Matrix> K_T_A_2  = K[1]; 
    boost::shared_ptr<Matrix> K_T_A_n  = K[2]; 
    boost::shared_ptr<Matrix> K_B      = K[3]; 
    boost::shared_ptr<Matrix> K_T_BA_2 = K[4]; 
    boost::shared_ptr<Matrix> K_T_BA_n = K[5]; 

    // => Compute the P matrices <= //
    
    boost::shared_ptr<Matrix> P_A    = build_D(Cocc_A_, Cocc_A_);
    boost::shared_ptr<Matrix> P_B    = build_D(Cocc_B_, Cocc_B_);

    // => Compute the V matrices <= //

    boost::shared_ptr<Matrix> V_A = build_V(primary_A_);
    boost::shared_ptr<Matrix> V_B = build_V(primary_B_);

    // ==> Electrostatic Terms <== //

    double Elst10r = 0.0;

    // Classical physics (watch for cancellation!)
    double Enuc = 0.0;
    Enuc += dimer_->nuclear_repulsion_energy();
    Enuc -= monomer_A_->nuclear_repulsion_energy();
    Enuc -= monomer_B_->nuclear_repulsion_energy();
    
    Elst10r += 1.0 * P_A->vector_dot(V_B);
    Elst10r += 1.0 * P_B->vector_dot(V_A);
    Elst10r += 2.0 * P_B->vector_dot(J_A);
    Elst10r += 1.0 * Enuc;
    
    energies_["Elst10r"] = Elst10r;

    fprintf(outfile,"    Elst10,r            = %18.12lf H\n",Elst10r);
    fflush(outfile);

    // ==> Exchange Terms (S^\infty) <== //
    
    double Exch10_n = 0.0;

    // => Compute the T matrices <= //

    boost::shared_ptr<Matrix> T_A_n  = build_D(Cocc_A_, C_T_A_n);
    boost::shared_ptr<Matrix> T_B_n  = build_D(Cocc_B_, C_T_B_n);
    boost::shared_ptr<Matrix> T_BA_n = build_D(Cocc_B_, C_T_BA_n);
    boost::shared_ptr<Matrix> T_AB_n = build_D(Cocc_A_, C_T_AB_n);

    Exch10_n += 1.0 * P_A->vector_dot(K_B);
    
    Exch10_n += 2.0 * T_A_n->vector_dot(V_B);
    Exch10_n += 4.0 * T_A_n->vector_dot(J_B);
    Exch10_n -= 2.0 * T_A_n->vector_dot(K_B);

    Exch10_n += 2.0 * T_B_n->vector_dot(V_A);
    Exch10_n += 4.0 * T_B_n->vector_dot(J_A);
    Exch10_n -= 2.0 * T_B_n->vector_dot(K_A);

    Exch10_n += 4.0 * T_AB_n->vector_dot(V_A);
    Exch10_n += 8.0 * T_AB_n->vector_dot(J_A);
    Exch10_n -= 4.0 * T_AB_n->vector_dot(K_A);
    Exch10_n += 4.0 * T_AB_n->vector_dot(V_B);
    Exch10_n += 8.0 * T_AB_n->vector_dot(J_B);
    Exch10_n -= 4.0 * T_AB_n->vector_dot(K_B);

    Exch10_n += 8.0 * T_AB_n->vector_dot(J_T_A_n);
    Exch10_n -= 4.0 * T_AB_n->vector_dot(K_T_A_n);
    Exch10_n += 8.0 * T_B_n->vector_dot(J_T_A_n);
    Exch10_n -= 4.0 * T_B_n->vector_dot(K_T_A_n);
    
    Exch10_n += 8.0 * T_BA_n->vector_dot(J_T_BA_n);
    Exch10_n -= 4.0 * T_BA_n->vector_dot(K_T_BA_n);
    Exch10_n += 8.0 * T_B_n->vector_dot(J_T_BA_n);
    Exch10_n -= 4.0 * T_B_n->vector_dot(K_T_BA_n);

    Exch10_n *= -1.0;

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
    J_T_BA_n.reset();
    K_T_A_n.reset();
    K_T_BA_n.reset();

    // ==> Exchange Terms (S^2) <== //
    
    double Exch10_2 = 0.0;

    // => Compute the T matrices <= //

    boost::shared_ptr<Matrix> T_A_2  = build_D(Cocc_A_, C_T_A_2);
    boost::shared_ptr<Matrix> T_B_2  = build_D(Cocc_B_, C_T_B_2);
    boost::shared_ptr<Matrix> T_BA_2 = build_D(Cocc_B_, C_T_BA_2);
    boost::shared_ptr<Matrix> T_AB_2 = build_D(Cocc_A_, C_T_AB_2);

    Exch10_2 += 1.0 * P_A->vector_dot(K_B);
    
    Exch10_2 += 2.0 * T_A_2->vector_dot(V_B);
    Exch10_2 += 4.0 * T_A_2->vector_dot(J_B);
    Exch10_2 -= 2.0 * T_A_2->vector_dot(K_B);

    Exch10_2 += 2.0 * T_B_2->vector_dot(V_A);
    Exch10_2 += 4.0 * T_B_2->vector_dot(J_A);
    Exch10_2 -= 2.0 * T_B_2->vector_dot(K_A);

    Exch10_2 += 4.0 * T_AB_2->vector_dot(V_A);
    Exch10_2 += 8.0 * T_AB_2->vector_dot(J_A);
    Exch10_2 -= 4.0 * T_AB_2->vector_dot(K_A);
    Exch10_2 += 4.0 * T_AB_2->vector_dot(V_B);
    Exch10_2 += 8.0 * T_AB_2->vector_dot(J_B);
    Exch10_2 -= 4.0 * T_AB_2->vector_dot(K_B);

    Exch10_2 += 8.0 * T_AB_2->vector_dot(J_T_A_2);
    Exch10_2 -= 4.0 * T_AB_2->vector_dot(K_T_A_2);
    Exch10_2 += 8.0 * T_B_2->vector_dot(J_T_A_2);
    Exch10_2 -= 4.0 * T_B_2->vector_dot(K_T_A_2);
    
    Exch10_2 += 8.0 * T_BA_2->vector_dot(J_T_BA_2);
    Exch10_2 -= 4.0 * T_BA_2->vector_dot(K_T_BA_2);
    Exch10_2 += 8.0 * T_B_2->vector_dot(J_T_BA_2);
    Exch10_2 -= 4.0 * T_B_2->vector_dot(K_T_BA_2);

    Exch10_2 *= -1.0;

    energies_["Exch10 (S^2)"] = Exch10_2;

    fprintf(outfile,"    Exch10 (S^2)        = %18.12lf H\n",Exch10_2);
    fflush(outfile);

    // Clear up some memory
    C_T_A_2.reset();
    C_T_B_2.reset();
    C_T_BA_2.reset();
    C_T_AB_2.reset();
    T_A_2.reset();
    T_B_2.reset();
    T_BA_2.reset();
    T_AB_2.reset();
    J_T_A_2.reset();
    J_T_BA_2.reset();
    K_T_A_2.reset();
    K_T_BA_2.reset();

    // ==> Uncorrelated Second-Order Response Terms [Induction] <== //

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
    boost::shared_ptr<Matrix> w_B = build_w(W_B,Caocc_A_,Cavir_A_);
    // Electrostatic potential of monomer A (ov space of B)
    boost::shared_ptr<Matrix> w_A = build_w(W_A,Caocc_B_,Cavir_B_);

    // TODO: CPHF

    boost::shared_ptr<Matrix> x_A(w_A->clone());
    boost::shared_ptr<Matrix> x_B(w_B->clone());

    // ==> Induction <== //

    double Ind20r_AB = - x_A->vector_dot(w_B);
    double Ind20r_BA = - x_B->vector_dot(w_A);
    double Ind20r = Ind20r_AB + Ind20r_BA;

    energies_["Ind20,r (A<-B)"] = Ind20r_AB;
    energies_["Ind20,r (B->A)"] = Ind20r_BA;
    energies_["Ind20,r"] = Ind20r_AB + Ind20r_BA;

    fprintf(outfile,"    Ind20 (A<-B)        = %18.12lf H\n",Ind20r_AB);
    fprintf(outfile,"    Ind20 (A->B)        = %18.12lf H\n",Ind20r_BA);
    fprintf(outfile,"    Ind20               = %18.12lf H\n",Ind20r);
    fflush(outfile);

    // ==> Exchange-Induction <== //

    // => AO-basis response <= //    

    // Response of monomber A (AO)
    boost::shared_ptr<Matrix> X_A = build_X(x_A,Caocc_A_,Cavir_A_);
    // Response of monomer B (AO)
    boost::shared_ptr<Matrix> X_B = build_X(x_B,Caocc_B_,Cavir_B_);

    // => O cross densities <= //

    //boost::shared_ptr<Matrix> O_AB = build_triple(P_A,S,P_B);
    //boost::shared_ptr<Matrix> O_BA = build_triple(P_B,S,P_A);

    // TODO Exhange induction
}
void DFTSAPT::mp2_terms()
{
    // TODO 
}
void DFTSAPT::print_trailer() const 
{
    fprintf(outfile, "\t Oh all the money that e'er I had, I spent it in good company.\n");
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
boost::shared_ptr<Matrix> DFTSAPT::build_Sij_2(boost::shared_ptr<Matrix> Sij)
{
    int nocc = Sij->nrow();

    boost::shared_ptr<Matrix> Sij2(new Matrix("Sij2 (MO)", nocc, nocc));

    double** Sijp = Sij->pointer();
    double** Sij2p = Sij2->pointer();

    Sij2->subtract(Sij);
    C_DGEMM('N','N',nocc,nocc,nocc,1.0,Sijp[0],nocc,Sijp[0],nocc,1.0,Sij2p[0],nocc);

    return Sij2;
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
        Sij2p[i][i] = 0.0;
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
    C_DGEMM('N','N',nso,nA,nB,1.0,CBp[0],nB,&Sp[nB][0],no,0.0,Cp[0],nA);

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
boost::shared_ptr<Matrix> DFTSAPT::build_X(boost::shared_ptr<Matrix> x, boost::shared_ptr<Matrix> L, boost::shared_ptr<Matrix> R)
{
    int nso = L->nrow();
    int no = L->ncol();
    int nv = R->ncol();

    boost::shared_ptr<Matrix> X(new Matrix("w", nso, nso));
    boost::shared_ptr<Matrix> T(new Matrix("T", no, nso));
    
    double** xp = x->pointer();
    double** Xp = X->pointer();
    double** Lp = L->pointer();
    double** Rp = R->pointer();
    double** Tp = T->pointer();
    
    C_DGEMM('N','T',no,nso,nv,1.0,xp[0],nv,Rp[0],nv,0.0,Tp[0],nso);
    C_DGEMM('N','N',nso,nso,no,1.0,Lp[0],no,Tp[0],nso,0.0,Xp[0],nso);

    return x;
}

}}
