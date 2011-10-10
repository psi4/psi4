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

OmegaWavefunction::OmegaWavefunction(Options& options,
                                     boost::shared_ptr<PSIO> psio,
                                     boost::shared_ptr<BasisSet> basisset,
                                     boost::shared_ptr<Matrix> Ca, int na,
                                     boost::shared_ptr<Matrix> Cb, int nb,
                                     boost::shared_ptr<Matrix> S, boost::shared_ptr<Matrix> X,
                                     boost::shared_ptr<Matrix> H, 
                                     boost::shared_ptr<OmegaDF> df, boost::shared_ptr<UKSPotential> ks)
 :  options_(options), psio_(psio), basisset_(basisset), nalpha_(na), nbeta_(nb), S_(S), X_(X), H_(H), df_(df), ks_(ks)
{
    nso_ = X_->rowspi()[0];
    nmo_ = X_->colspi()[0];

    Ca_ = boost::shared_ptr<Matrix>(new Matrix("Ca", nso_, nmo_));
    Cb_ = boost::shared_ptr<Matrix>(new Matrix("Cb", nso_, nmo_));

    double** Car = Ca->pointer();
    double** Cbr = Cb->pointer();
    double** Cap = Ca_->pointer();
    double** Cbp = Cb_->pointer();

    int nmo2 = nmo_;
    if (nmo2 > Ca->colspi()[0])
        nmo2 = Ca->colspi()[0];

    for (int i = 0; i < nso_; i++){
        C_DCOPY(nmo2,Car[i],1,Cap[i],1);
        C_DCOPY(nmo2,Cbr[i],1,Cbp[i],1);
    }

    common_init();
}
OmegaWavefunction::~OmegaWavefunction()
{
}
void OmegaWavefunction::common_init()
{
    Enuc_ = basisset_->molecule()->nuclear_repulsion_energy();    

    iteration_ = 0;
    min_diis_vectors_ = options_.get_int("MIN_DIIS_VECTORS");
    max_diis_vectors_ = options_.get_int("MAX_DIIS_VECTORS");
    diis_start_ = options_.get_int("START_DIIS_ITER");
    diis_enabled_ = options_.get_bool("DIIS") && min_diis_vectors_ > 1;

    debug_ = options_.get_int("DEBUG");
    print_ = options_.get_int("PRINT");

    Exc_ = 0.0;
    E_ = 0.0;
    Eold_ = 0.0;
    Drms_ = 0.0;

    epsilon_a_ = boost::shared_ptr<Vector>(new Vector("Epsilon A", nmo_));    
    epsilon_b_ = boost::shared_ptr<Vector>(new Vector("Epsilon B", nmo_));    

    Da_ = boost::shared_ptr<Matrix>(new Matrix("Da", nso_, nso_));    
    Db_ = boost::shared_ptr<Matrix>(new Matrix("Db", nso_, nso_));    
    Dt_ = boost::shared_ptr<Matrix>(new Matrix("Dt", nso_, nso_));    
    Dtold_ = boost::shared_ptr<Matrix>(new Matrix("Dtold", nso_, nso_));    

    J_ = boost::shared_ptr<Matrix>(new Matrix("J", nso_, nso_));    
    Ka_ = boost::shared_ptr<Matrix>(new Matrix("Ka", nso_, nso_));    
    Kb_ = boost::shared_ptr<Matrix>(new Matrix("Kb", nso_, nso_));    
    wKa_ = boost::shared_ptr<Matrix>(new Matrix("wKa", nso_, nso_));    
    wKb_ = boost::shared_ptr<Matrix>(new Matrix("wKb", nso_, nso_));    
    Va_ = boost::shared_ptr<Matrix>(new Matrix("Va", nso_, nso_));    
    Vb_ = boost::shared_ptr<Matrix>(new Matrix("Vb", nso_, nso_));    
    Ga_ = boost::shared_ptr<Matrix>(new Matrix("Ga", nso_, nso_));    
    Gb_ = boost::shared_ptr<Matrix>(new Matrix("Gb", nso_, nso_));    
    Fa_ = boost::shared_ptr<Matrix>(new Matrix("Fa", nso_, nso_));    
    Fb_ = boost::shared_ptr<Matrix>(new Matrix("Fb", nso_, nso_));    

    form_D();

    if (debug_ > 1) {
    
        fprintf(outfile, "  Nalpha = %3d\n", nalpha_);
        fprintf(outfile, "  Nbeta  = %3d\n\n", nbeta_);

        S_->print();
        X_->print();
        H_->print();
        Ca_->print();
        Cb_->print();
        Da_->print();
        Db_->print();
    }
}
std::string OmegaWavefunction::step()
{
    // Save old info
    save_info();

    // Interelectronic Integrals
    form_V();    

    // KS stuff
    form_G();

    // Build KS Matrix
    form_F();   

    // Compute Energy
    compute_E();

    // DIIS
    bool diised = diis();

    // Build C
    form_C(); 
   
    // Build D
    form_D();

    std::string status = (diised ? "DIIS" : "");

    return status;
}
void OmegaWavefunction::guess(boost::shared_ptr<Matrix> Fa, boost::shared_ptr<Matrix> Fb)
{
    Fa_->copy(Fa);
    Fb_->copy(Fb);

    form_C();
    form_D();
}
void OmegaWavefunction::save_info()
{
    Eold_ = E_;
    Dtold_->copy(Dt_);
}
void OmegaWavefunction::form_G()
{
    double a = ks_->functional()->getExactExchange();

    J_ = df_->J(Dt_);
    if (a != 0.0) {
        Ka_ = df_->K(Ca_, nalpha_);
        Kb_ = df_->K(Cb_, nbeta_);
        Ka_->scale(-a);
        Kb_->scale(-a);
    }
    wKa_ = df_->wK(Ca_, nalpha_);
    wKb_ = df_->wK(Cb_, nbeta_);
    wKa_->scale(-(1.0 - a));
    wKb_->scale(-(1.0 - a));
}
void OmegaWavefunction::form_V()
{
    boost::shared_ptr<Dimension> na(new Dimension(1));
    boost::shared_ptr<Dimension> nb(new Dimension(1));

    (*na)[0] = nalpha_;
    (*nb)[0] = nbeta_;

    ks_->buildPotential(Da_, Ca_, na, Db_, Cb_, nb);
    Va_ = ks_->Va_USO();
    Vb_ = ks_->Vb_USO();
    Exc_ = ks_->quadrature_value("FUNCTIONAL");
}
void OmegaWavefunction::form_F()
{
    Ga_->copy(J_);
    Ga_->add(Ka_);
    Ga_->add(wKa_);
    Ga_->add(Va_);

    Gb_->copy(J_);
    Gb_->add(Kb_);
    Gb_->add(wKb_);
    Gb_->add(Vb_);

    Fa_->copy(H_);
    Fa_->add(Ga_);

    Fb_->copy(H_);
    Fb_->add(Gb_);

    if (debug_ > 2) {
        J_->print();
        Ka_->print();
        Kb_->print();
        wKa_->print();
        wKb_->print();
        Va_->print();
        Vb_->print();
        Fa_->print();
        Fb_->print();
    }
}
void OmegaWavefunction::compute_E()
{
    // E_DFT = 2.0 D*H + 2.0 D*J - \alpha D*K + E_xc
    double one_electron_E = Dt_->vector_dot(H_);
    double coulomb_E = Dt_->vector_dot(J_);

    double exchange_E = 0.0;
    exchange_E += Da_->vector_dot(Ka_);
    exchange_E += Db_->vector_dot(Kb_);
    exchange_E += Da_->vector_dot(wKa_);
    exchange_E += Db_->vector_dot(wKb_);

    double Etotal = 0.0;
    Etotal += Enuc_;
    Etotal += one_electron_E;
    Etotal += 0.5 * coulomb_E;
    Etotal += 0.5 * exchange_E;
    Etotal += Exc_;

    double dashD_E = 0.00;
    if (ks_->functional()->isDashD()) {
        dashD_E = ks_->functional()->getDashD()->computeEnergy(basisset_->molecule());
        Etotal += dashD_E;
    }

    if (debug_ > 2) {
        fprintf(outfile, "Nuclear Repulsion Energy = %24.14f\n", Enuc_);
        fprintf(outfile, "One-Electron Energy =      %24.14f\n", one_electron_E);
        fprintf(outfile, "Coulomb Energy =           %24.14f\n", 0.5 * coulomb_E);
        fprintf(outfile, "Hybrid Exchange Energy =   %24.14f\n", 0.5 * exchange_E);
        fprintf(outfile, "XC Functional Energy =     %24.14f\n", Exc_);
        fprintf(outfile, "-D Energy =                %24.14f\n", dashD_E);
    }

    E_ = Etotal;
}
bool OmegaWavefunction::diis()
{
    iteration_++;

    if (diis_enabled_ && iteration_ > 0 && iteration_ >= diis_start_ ) {
        boost::shared_ptr<Matrix> FDSa = form_FDSmSDF(Fa_, Da_);
        boost::shared_ptr<Matrix> FDSb = form_FDSmSDF(Fb_, Db_);
        diis_manager_->add_entry(4, FDSa.get(), FDSb.get(), Fa_.get(), Fb_.get());
    }

    bool diis_performed;
    if (diis_enabled_ == true && iteration_ >= diis_start_ + min_diis_vectors_ - 1) {
        diis_performed = diis_manager_->extrapolate(2, Fa_.get(), Fb_.get());
    } else {
        diis_performed = false;
    }
    return diis_performed;   
}
void OmegaWavefunction::form_C()
{
    boost::shared_ptr<Matrix> T(new Matrix("T1", nso_, nmo_));
    boost::shared_ptr<Matrix> Fp(new Matrix("Fp", nmo_, nmo_));
    boost::shared_ptr<Matrix> Cp(new Matrix("Cp", nmo_, nmo_));

    double** Tp = T->pointer();
    double** Fpp = Fp->pointer();
    double** Cpp = Cp->pointer();
    double** Xp = X_->pointer();
    double** Fap = Fa_->pointer();
    double** Fbp = Fb_->pointer();
    double** Cap = Ca_->pointer();
    double** Cbp = Cb_->pointer();

    C_DGEMM('N','N',nso_,nmo_,nso_,1.0,Fap[0],nso_,Xp[0],nmo_,0.0,Tp[0],nmo_); 
    C_DGEMM('T','N',nmo_,nmo_,nso_,1.0,Xp[0],nmo_,Tp[0],nmo_,0.0,Fpp[0],nmo_); 

    Fp->diagonalize(Cp, epsilon_a_);

    C_DGEMM('N','N',nso_,nmo_,nmo_,1.0,Xp[0],nmo_,Cpp[0],nmo_,0.0,Cap[0],nmo_);

    C_DGEMM('N','N',nso_,nmo_,nso_,1.0,Fbp[0],nso_,Xp[0],nmo_,0.0,Tp[0],nmo_); 
    C_DGEMM('T','N',nmo_,nmo_,nso_,1.0,Xp[0],nmo_,Tp[0],nmo_,0.0,Fpp[0],nmo_); 

    Fp->diagonalize(Cp, epsilon_b_);

    C_DGEMM('N','N',nso_,nmo_,nmo_,1.0,Xp[0],nmo_,Cpp[0],nmo_,0.0,Cbp[0],nmo_);

    if (debug_ > 2) {
        Ca_->print();
        Cb_->print();
        epsilon_a_->print();
        epsilon_b_->print();
    }
}
void OmegaWavefunction::form_D()
{
    double** Dap = Da_->pointer();
    double** Cap = Ca_->pointer();
    double** Dbp = Db_->pointer();
    double** Cbp = Cb_->pointer();
   
    if (nalpha_ > 0) {
        C_DGEMM('N','T',nso_,nso_,nalpha_,1.0,Cap[0],nmo_,Cap[0],nmo_,0.0,Dap[0],nso_); 
    } else {
        Da_->zero();
    }
    if (nbeta_ > 0) {
        C_DGEMM('N','T',nso_,nso_,nbeta_,1.0,Cbp[0],nmo_,Cbp[0],nmo_,0.0,Dbp[0],nso_); 
    } else {
        Db_->zero();
    }

    Dt_->copy(Da_);
    Dt_->add(Db_);

    Dtold_->subtract(Dt_);
    Drms_ = Dtold_->rms();
    Dtold_->add(Dt_);

    if (debug_ > 2) {
        Da_->print();
        Db_->print();
    }
}
void OmegaWavefunction::reset()
{
    iteration_ = 0;
    if (diis_enabled_) {
        if (diis_manager_) diis_manager_.reset();
        diis_manager_ = boost::shared_ptr<DIISManager>(new DIISManager(max_diis_vectors_,
        "HF DIIS VECTORS")); 
    
        boost::shared_ptr<Matrix> FDS = form_FDSmSDF(Fa_, Da_);
        diis_manager_->set_error_vector_size(2, DIISEntry::Matrix, FDS.get(),
                                                DIISEntry::Matrix, FDS.get()); 
        diis_manager_->set_vector_size(2, DIISEntry::Matrix, Fa_.get(),
                                          DIISEntry::Matrix, Fb_.get()); 
    }
}
void OmegaWavefunction::clear()
{
    iteration_ = 0;
    diis_manager_.reset();
}
boost::shared_ptr<Matrix> OmegaWavefunction::form_FDSmSDF(boost::shared_ptr<Matrix> F, boost::shared_ptr<Matrix> D)
{
    boost::shared_ptr<Matrix> FDSmSDF(new Matrix("FDS-SDF", nso_, nso_));
    boost::shared_ptr<Matrix> DS(new Matrix("DS", nso_, nso_));

    DS->gemm(false,false,1.0,D,S_,0.0);
    FDSmSDF->gemm(false,false,1.0,F,DS,0.0);

    boost::shared_ptr<Matrix> SDF(FDSmSDF->transpose());
    FDSmSDF->subtract(SDF);

    DS.reset();
    SDF.reset();

    boost::shared_ptr<Matrix> XP(new Matrix("X'(FDS - SDF)", nmo_, nso_));
    boost::shared_ptr<Matrix> XPX(new Matrix("X'(FDS - SDF)X", nmo_, nmo_));
    XP->gemm(true,false,1.0,X_,FDSmSDF,0.0);
    XPX->gemm(false,false,1.0,XP,X_,0.0);

    //XPX->print();

    return XPX;
}
void OmegaWavefunction::print_orbitals()
{
    fprintf(outfile, "\t%-70s\n\n\t", "Alpha Occupied:");
    int count = 0;
    for (int i = 0; i < nalpha_; i++) {
        fprintf(outfile, "%4d%-4s%11.6f  ", i + 1, "A", epsilon_a_->get(0,i));
        if (count++ % 3 == 2 && count != nalpha_)
            fprintf(outfile, "\n\t");
    }
    fprintf(outfile, "\n\n");

    fprintf(outfile, "\t%-70s\n\n\t", "Alpha Virtual:");
    count = 0;
    for (int i = nalpha_; i < nmo_; i++) {
        fprintf(outfile, "%4d%-4s%11.6f  ", i + 1, "A", epsilon_a_->get(0,i));
        if (count++ % 3 == 2 && count != nmo_ - nalpha_)
            fprintf(outfile, "\n\t");
    }
    fprintf(outfile, "\n\n");

    fprintf(outfile, "\t%-70s\n\n\t", "Beta Occupied:");
    count = 0;
    for (int i = 0; i < nbeta_; i++) {
        fprintf(outfile, "%4d%-4s%11.6f  ", i + 1, "A", epsilon_b_->get(0,i));
        if (count++ % 3 == 2 && count != nbeta_)
            fprintf(outfile, "\n\t");
    }
    fprintf(outfile, "\n\n");

    fprintf(outfile, "\t%-70s\n\n\t", "Beta Virtual:");
    count = 0;
    for (int i = nbeta_; i < nmo_; i++) {
        fprintf(outfile, "%4d%-4s%11.6f  ", i + 1, "A", epsilon_b_->get(0,i));
        if (count++ % 3 == 2 && count != nmo_ - nbeta_)
            fprintf(outfile, "\n\t");
    }
    fprintf(outfile, "\n\n");
}
double OmegaWavefunction::koopmansIP()
{
    double ea = epsilon_a_->get(0, nalpha_ - 1);
    double eb = (nbeta_ > 0 ? epsilon_b_->get(0, nbeta_ - 1) : ea);

    return - (ea > eb ? ea : eb);
}
double OmegaWavefunction::koopmansEA()
{
    double ea = epsilon_a_->get(0, nalpha_);
    double eb = epsilon_b_->get(0, nbeta_);

    return - (ea > eb ? eb : ea);
}

}} // End Namespaces  
