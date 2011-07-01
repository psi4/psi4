#include "3index.h"
#include "libmints/mints.h"
#include <libqt/qt.h>

#include <string>
#include <sstream>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <algorithm>
#include <utility>
#include <ctype.h>

using namespace boost;
using namespace std;
using namespace psi;

namespace psi {

PSTensor::PSTensor(boost::shared_ptr<BasisSet> primary,
                   boost::shared_ptr<Matrix> C,
                   int nocc,
                   int nvir,
                   int naocc,
                   int navir,
                   Options& options,
                   double omega) :
    primary_(primary), C_(C), nocc_(nocc), nvir_(nvir),
    naocc_(naocc), navir_(navir), options_(options), omega_(omega)
{
    common_init();
}
PSTensor::~PSTensor()
{
}
void PSTensor::common_init()
{
    print_ = options_.get_int("PRINT");
    debug_ = options_.get_int("DEBUG");

    use_omega_ = true;
    if (omega_ == -1.0) {
        use_omega_ = false;
    }

    print_header();

    molecule_ = primary_->molecule();

    nfocc_ = nocc_ - naocc_;
    nfvir_ = nvir_ - navir_;

    nso_ = C_->rowspi()[0];
    nmo_ = nmo2_ = C_->colspi()[0];

    Caocc_ = boost::shared_ptr<Matrix>(new Matrix("C active occupied", nso_, naocc_));
    Cavir_ = boost::shared_ptr<Matrix>(new Matrix("C active virtual", nso_, navir_));
    
    double** Cp = C_->pointer();
    double** Cop  = Caocc_->pointer();
    double** Cvp  = Cavir_->pointer();

    for (int m = 0; m < nso_; m++) {
        C_DCOPY(naocc_, &Cp[m][nfocc_],1, Cop[m], 1);
        C_DCOPY(navir_, &Cp[m][nocc_],1, Cvp[m], 1);
    }

    if (debug_) {
        C_->print();
        Caocc_->print();
        Cavir_->print();
    }

    min_S_primary_ = options_.get_double("PS_MIN_S_PRIMARY");
    min_S_dealias_ = options_.get_double("PS_MIN_S_DEALIAS");

    buildDealiasSet();
    buildGrid();
    buildR();
    if (options_.get_str("PS_FITTING_ALGORITHM") == "CONDITIONED") {
        buildQ();
    } else if (options_.get_str("PS_FITTING_ALGORITHM") == "CANONICAL") { 
        buildQ_canonical();
    } else  if (options_.get_str("PS_FITTING_ALGORITHM") == "RENORMALIZED") { 
        buildQ_renormalized();
    } else  if (options_.get_str("PS_FITTING_ALGORITHM") == "QUADRATURE") { 
        buildQ_quadrature();
    }
}
void PSTensor::print_header()
{
    fprintf(outfile,"  ==> PS Tensor (by Rob Parrish) <==\n\n");
    fprintf(outfile,"   => Range-Separation Lowpass Proceedure\n\n");
    if (use_omega_) {
        fprintf(outfile, "    Range separation will be performed with \\omega of %11.4E.\n\n", omega_);
    } else {
        fprintf(outfile, "    No range separation requested.\n\n");
    }
}
void PSTensor::buildDealiasSet() 
{
    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());

    if (print_) {
        fprintf(outfile," => Primary Basis Set <= \n\n");
        primary_->print_by_level(outfile,print_);
        fflush(outfile);
    }

    if (print_) {
        fprintf(outfile," => Dealias Basis Set <= \n\n");
        if (options_.get_str("DEALIAS_BASIS_CC") == "") {

            fprintf(outfile,"  Dealias Basis Automatically Generated\n\n");

            boost::shared_ptr<DealiasBasisSet> d(new DealiasBasisSet(primary_, options_));
            dealias_ = d->dealiasSet();

        } else {
            fprintf(outfile,"  Dealias Basis Read from %s", options_.get_str("DEALIAS_BASIS_CC").c_str()); 
            molecule_->set_basis_all_atoms(options_.get_str("DEALIAS_BASIS_CC"),"DEALIAS_BASIS");
            dealias_ = BasisSet::construct(parser,molecule_,"DEALIAS_BASIS");  
        }
        dealias_->print_by_level(outfile,print_);
    }

    ndealias_ = ndealias2_ = dealias_->nbf();
    naug_ = nmo_ + ndealias_;
    fflush(outfile);
}
void PSTensor::buildGrid()
{
    if (options_.get_str("PS_GRID_FILE") == "") {
        grid_ = boost::shared_ptr<PseudospectralGrid>(new PseudospectralGrid(molecule_,
            primary_, dealias_, options_));
    } else {
        grid_ = boost::shared_ptr<PseudospectralGrid>(new PseudospectralGrid(molecule_,
            primary_, dealias_, options_.get_str("PS_GRID_FILE"), options_));
    }

    grid_->print(outfile, print_);
    fflush(outfile);
    
    naux_ = grid_->npoints();
    w_ = boost::shared_ptr<Vector>(new Vector("Grid Weights", naux_));
    double* wp = w_->pointer();
    
    C_DCOPY(naux_, grid_->w(), 1, wp, 1);
}
void PSTensor::buildR()
{
    // First build the orthonormal dealias set
    form_Spdao(); // <\phi_\mu | \phi_\eta>
    form_Spdmo(); // <\phi_s | \phi_\eta>

    // Now build the collocation matrices
    form_Rpao(); // R_\nu^P
    form_Rdao(); // R_\eta^P
    form_Rpmo(); // R_s^P 
    form_Rdmo(); // R_n^P
    form_Ra(); // R_a^P
}
void PSTensor::form_Spdao()
{
    boost::shared_ptr<IntegralFactory> fact(new IntegralFactory(primary_,dealias_,primary_,primary_));
    boost::shared_ptr<OneBodyAOInt> Sint(fact->ao_overlap());

    Spdao_ = boost::shared_ptr<Matrix>(new Matrix("S (primary x dealias)", nso_, ndealias_));
    Sint->compute(Spdao_);

    if (debug_)
        Spdao_->print();
}
void PSTensor::form_Spdmo()
{
    Spdmo_ = boost::shared_ptr<Matrix>(new Matrix("S (primary' x dealias)", nmo_, ndealias_));
    double** Spd3p = Spdmo_->pointer();
    double** Spdp = Spdao_->pointer();
    double** Xp = C_->pointer();

    C_DGEMM('T','N',nmo_,ndealias_,nso_,1.0,Xp[0],nmo_,Spdp[0],ndealias_,0.0,Spd3p[0],ndealias_);

    if (debug_)
        Spdmo_->print();

    Spdao_.reset();
}
void PSTensor::form_Rpao()
{
    Rpao_ = boost::shared_ptr<Matrix>(new Matrix("R (primary x points)", nso_, naux_));
    double** Rp = Rpao_->pointer();

    boost::shared_ptr<BasisPoints> points(new BasisPoints(primary_, naux_));
    points->setToComputePoints(true);
    double** bpoints = points->getPoints();

    // Compute the basis points
    points->computePoints(grid_->fullGrid());

    // Copy the points in
    for (int i = 0; i < naux_; i++) {
        for (int Q = 0; Q < nso_; Q++)
            Rp[Q][i] = bpoints[i][Q];
    }

    if (debug_)
        Rpao_->print();
}
void PSTensor::form_Rdao()
{
    Rdao_ = boost::shared_ptr<Matrix>(new Matrix("R (dealias x points)", ndealias_, naux_));
    double** Rp = Rdao_->pointer();

    boost::shared_ptr<BasisPoints> points(new BasisPoints(dealias_, naux_));
    points->setToComputePoints(true);
    double** bpoints = points->getPoints();

    // Compute the basis points
    points->computePoints(grid_->fullGrid());

    // Copy the points in
    for (int i = 0; i < naux_; i++) {
        for (int Q = 0; Q < ndealias_; Q++)
            Rp[Q][i] = bpoints[i][Q];
    }

    if (debug_)
        Rdao_->print();
}
void PSTensor::form_Rpmo()
{
    Rpmo_ = boost::shared_ptr<Matrix>(new Matrix("R2 (primary' x points)", nmo_, naux_));
    double** Rp2 = Rpmo_->pointer();
    double** Rp = Rpao_->pointer();
    double** Xp = C_->pointer();

    C_DGEMM('T','N',nmo_,naux_,nso_,1.0,Xp[0],nmo_,Rp[0],naux_,0.0,Rp2[0],naux_);

    if (debug_)
        Rpmo_->print();

    Rmo_ = Rpmo_;
    Rpao_.reset();
}
void PSTensor::form_Rdmo()
{
    Rdmo_ = boost::shared_ptr<Matrix>(new Matrix("R2 (dealias' x points)", ndealias_, naux_));
    Rdmo_->copy(Rdao_);
    double** Rd2p = Rdmo_->pointer();
    double** Rp2p = Rpmo_->pointer();
    double** Cdp = Spdmo_->pointer();

    C_DGEMM('T','N',ndealias_,naux_,nmo_,-1.0,Cdp[0],ndealias_,Rp2p[0],naux_,1.0,Rd2p[0],naux_);

    if (debug_)
        Rdmo_->print();

    Spdmo_.reset();
}
void PSTensor::form_Ra()
{
    Ra_ = boost::shared_ptr<Matrix>(new Matrix("R Augmented (primary' + dealias' x points)", naug_, naux_));
    double** Rap = Ra_->pointer(); 

    double** Rpp = Rpmo_->pointer();
    double** Rdp = Rdmo_->pointer();

    C_DCOPY(nmo_ * naux_, Rpp[0], 1, Rap[0], 1);
    C_DCOPY(ndealias_ * naux_, Rdp[0], 1, Rap[nmo_], 1);

    if (debug_)
        Ra_->print();
}
void PSTensor::buildQ()
{
    fprintf(outfile, "   => Fitting Procedure: Conditioned <=\n\n");

    form_Cpp();
    form_U();

    form_Cpd();    
    form_V();

    form_Cdd();
    form_W();

    form_X();

    if (debug_ > 1) 
        validate_X();
    
    form_Q();
}
void PSTensor::buildQ_canonical()
{
    fprintf(outfile, "   => Fitting Procedure: Canonical <=\n\n");

    nmo2_ = nmo_;
    ndealias2_ = ndealias_;

    boost::shared_ptr<Matrix> C(new Matrix("C" , nmo_ + ndealias_, nmo_ + ndealias_));
    boost::shared_ptr<Matrix> Rw(new Matrix("Rw", nmo_ + ndealias_, naux_));
    double** Cp = C->pointer();
    double** Rp = Ra_->pointer();
    double** Rwp = Rw->pointer();
    double* wp = w_->pointer();

    C_DCOPY((nmo_ + ndealias_) * naux_, Rp[0], 1, Rwp[0], 1);

    for (int Q = 0; Q < naux_; Q++) {
        C_DSCAL(nmo_ + ndealias_, wp[Q], &Rwp[0][Q], naux_);
    } 

    C_DGEMM('N','T',nmo_ + ndealias_,nmo_+ndealias_,naux_,1.0,Rp[0],naux_,Rwp[0],naux_,0.0,Cp[0],nmo_+ndealias_);

    if (debug_)
        C->print();

    C->power(-1.0);     

    if (debug_)
        C->print(outfile, "After Inversion");

    Qmo_ = boost::shared_ptr<Matrix> (new Matrix("Qmo", nmo_, naux_));
    double** Qmop = Qmo_->pointer();
    
    C_DGEMM('N','N',nmo_,naux_,nmo_ + ndealias_,1.0,Cp[0],nmo_ + ndealias_,Rwp[0],naux_,0.0,Qmop[0],naux_);

    if (debug_)
        Qmo_->print();
}
void PSTensor::buildQ_renormalized()
{
    fprintf(outfile, "   => Fitting Procedure: Renormalized <=\n\n");

    nmo2_ = nmo_;
    ndealias2_ = 0;

    form_Cpp();
    double** Cp = Cpp_->pointer();

    Cpp_->power(-1.0);     

    if (debug_)
        Cpp_->print(outfile, "After Inversion");

    boost::shared_ptr<Matrix> Rw(new Matrix("Rw", nmo_, naux_));
    double** Rp = Ra_->pointer();
    double** Rwp = Rw->pointer();
    double* wp = w_->pointer();

    C_DCOPY((nmo_) * naux_, Rp[0], 1, Rwp[0], 1);

    for (int Q = 0; Q < naux_; Q++) {
        C_DSCAL(nmo_, wp[Q], &Rwp[0][Q], naux_);
    } 
    Qmo_ = boost::shared_ptr<Matrix> (new Matrix("Qmo", nmo_, naux_));
    double** Qmop = Qmo_->pointer();
    
    C_DGEMM('N','N',nmo_,naux_,nmo_,1.0,Cp[0],nmo_,Rwp[0],naux_,0.0,Qmop[0],naux_);

    if (debug_)
        Qmo_->print();
}
void PSTensor::buildQ_quadrature()
{
    fprintf(outfile, "   => Fitting Procedure: Quadrature <=\n\n");

    nmo2_ = nmo_;
    ndealias2_ = 0;

    Qmo_ = boost::shared_ptr<Matrix> (new Matrix("Qmo", nmo_, naux_));
    double** Qmop = Qmo_->pointer();
    double** Rp = Ra_->pointer();
    double* wp = w_->pointer();
    C_DCOPY((nmo_) * naux_, Rp[0], 1, Qmop[0], 1);

    for (int Q = 0; Q < naux_; Q++) {
        C_DSCAL(nmo_, wp[Q], &Qmop[0][Q], naux_);
    } 

    if (debug_)
        Qmo_->print();
}
void PSTensor::form_Cpp()
{
    Cpp_ = boost::shared_ptr<Matrix>(new Matrix("Cpp" , nmo_, nmo_));
    boost::shared_ptr<Matrix> Rw(new Matrix("Rw", nmo_, naux_));
    double** Cp = Cpp_->pointer();
    double** Rp = Ra_->pointer();
    double** Rwp = Rw->pointer();
    double* wp = w_->pointer();

    C_DCOPY(nmo_ * naux_, Rp[0], 1, Rwp[0], 1);

    for (int Q = 0; Q < naux_; Q++) {
        C_DSCAL(nmo_, wp[Q], &Rwp[0][Q], naux_);
    } 

    C_DGEMM('N','T',nmo_,nmo_,naux_,1.0,Rp[0],naux_,Rwp[0],naux_,0.0,Cp[0],nmo_);

    if (debug_)
        Cpp_->print();
}
void PSTensor::form_U()
{
    boost::shared_ptr<Matrix> V(new Matrix("Eigvecs", nmo_, nmo_));
    boost::shared_ptr<Vector> c(new Vector("Eigvals", nmo_));
    double** Vp = V->pointer();
    double*  cp = c->pointer();

    Cpp_->diagonalize(V,c);

    if (debug_ > 1)
        V->eivprint(c);
   
    nmo2_ = 0;
    for (int i = 0; i < nmo_; i++) {
        if (cp[i] >= min_S_primary_)
            nmo2_++;
    }

    if (nmo2_ < nmo_) {
        fprintf(outfile, "  WARNING: This grid cannot numerically distinguish primary basis functions.\n"); 
        fprintf(outfile, "           These pseudospectral results may be inaccurate.\n\n");
    }
 
    U_ = boost::shared_ptr<Matrix>(new Matrix("U", nmo_, nmo2_)); 
    double** Up = U_->pointer();
   
    int j = 0; 
    for (int i = 0; i < nmo_; i++) {
        if (cp[i] >= min_S_primary_) {
            C_DAXPY(nmo_, pow(cp[i], -1.0/2.0), &Vp[0][i], nmo_, &Up[0][j], nmo2_); 
            j++;
        }
    }

    if (debug_)
        U_->print();

    Cpp_.reset();
}
void PSTensor::form_Cpd()
{
    boost::shared_ptr<Matrix> Cpdao(new Matrix("Cpdao", nmo_, ndealias_));
    double** Cpdaop = Cpdao->pointer();

    boost::shared_ptr<Matrix> Rw(new Matrix("Rw", nmo_, naux_));
    double** Rp = Ra_->pointer();
    double** Rwp = Rw->pointer();
    double* wp = w_->pointer();

    C_DCOPY(nmo_ * naux_, Rp[0], 1, Rwp[0], 1);

    for (int Q = 0; Q < naux_; Q++) {
        C_DSCAL(nmo_, wp[Q], &Rwp[0][Q], naux_);
    } 

    C_DGEMM('N','T',nmo_,ndealias_,naux_,1.0,Rwp[0],naux_,Rp[nmo_],naux_,0.0,Cpdaop[0],ndealias_);

    Rw.reset();

    if (debug_ > 1)
        Cpdao->print();

    Cpd_ = boost::shared_ptr<Matrix>(new Matrix("Cpd", nmo2_, ndealias_));
    double** Cpdp = Cpd_->pointer();
    double** Xp = U_->pointer();

    C_DGEMM('T','N',nmo2_,ndealias_,nmo_,1.0,Xp[0],nmo2_,Cpdaop[0],ndealias_,0.0,Cpdp[0],ndealias_);


    if (debug_)
        Cpd_->print();

    Cpdao.reset();
}
void PSTensor::form_V()
{
    V_ = boost::shared_ptr<Matrix>(new Matrix("V", nmo2_, ndealias_));
    V_->copy(Cpd_);
    V_->scale(-1.0);
    
    if (debug_) 
        V_->print();
}
void PSTensor::form_Cdd()
{
    Cdd_ = boost::shared_ptr<Matrix>(new Matrix("Cdd",ndealias_,ndealias_));
    double** Cddp = Cdd_->pointer();

    boost::shared_ptr<Matrix> Rw(new Matrix("Rw", ndealias_, naux_));
    double** Rp = Ra_->pointer();
    double** Rwp = Rw->pointer();
    double* wp = w_->pointer();

    C_DCOPY(ndealias_ * naux_, Rp[nmo_], 1, Rwp[0], 1);

    for (int Q = 0; Q < naux_; Q++) {
        C_DSCAL(ndealias_, wp[Q], &Rwp[0][Q], naux_);
    } 

    C_DGEMM('N','T',ndealias_,ndealias_,naux_,1.0,Rwp[0],naux_,Rp[nmo_],naux_,0.0,Cddp[0],ndealias_);
    
    Rw.reset();

    if (debug_ > 1)
        Cdd_->print(outfile, "before orthogonalization");
    
    double** Cpdp = Cpd_->pointer();

    C_DGEMM('T','N',ndealias_,ndealias_,nmo2_,-1.0,Cpdp[0],ndealias_,Cpdp[0],ndealias_,1.0,Cddp[0],ndealias_);

    if (debug_)
        Cdd_->print();
}
void PSTensor::form_W()
{
    boost::shared_ptr<Matrix> V(new Matrix("Eigvecs", ndealias_, ndealias_));
    boost::shared_ptr<Vector> c(new Vector("Eigvals", ndealias_));
    double** Vp = V->pointer();
    double*  cp = c->pointer();

    Cdd_->diagonalize(V,c);

    if (debug_ > 1)
        V->eivprint(c);
   
    ndealias2_ = 0;
    for (int i = 0; i < ndealias_; i++) {
        if (cp[i] >= min_S_dealias_)
            ndealias2_++;
    }

    if (print_) {
        fprintf(outfile, "  %d of %d dealias functions selected, %d projected out.\n\n", ndealias2_, ndealias_,
            ndealias_ - ndealias2_); 
        fflush(outfile);
    }
    

    W_ = boost::shared_ptr<Matrix>(new Matrix("w", ndealias_, ndealias2_)); 
    double** Wp = W_->pointer();
   
    int j = 0; 
    for (int i = 0; i < ndealias_; i++) {
        if (cp[i] >= min_S_dealias_) {
            C_DAXPY(ndealias_, pow(cp[i], -1.0/2.0), &Vp[0][i], ndealias_, &Wp[0][j], ndealias2_); 
            j++;
        }
    }

    if (debug_)
        W_->print();

    Cdd_.reset();
}
void PSTensor::form_X()
{
    X_ = boost::shared_ptr<Matrix>(new Matrix("X", nmo_ + ndealias_, nmo2_ + ndealias2_));
    double** Xp = X_->pointer(); 

    double** Up = U_->pointer();
    double** Vp = V_->pointer();
    double** Wp = W_->pointer();

    for (int m = 0; m < nmo_; m++) {
        C_DCOPY(nmo2_, Up[m], 1, Xp[m], 1);
    }

    boost::shared_ptr<Matrix> T(new Matrix("T",nmo_,ndealias_));   
    double** Tp = T->pointer();

    C_DGEMM('N','N',nmo_,ndealias_,nmo2_,1.0,Up[0],nmo2_,Vp[0],ndealias_,0.0,Tp[0],ndealias_);
    C_DGEMM('N','N',nmo_,ndealias2_,ndealias_,1.0,Tp[0],ndealias_,Wp[0],ndealias2_,0.0,&Xp[0][nmo2_],nmo2_ + ndealias2_);

    for (int d = 0; d < ndealias_; d++) {
        C_DCOPY(ndealias2_, Wp[d], 1, &Xp[nmo2_ + d][nmo_], 1);
    }
    
    if (debug_)
        X_->print();
}
void PSTensor::validate_X()
{
    boost::shared_ptr<Matrix> Rw(new Matrix("Rw", nmo_ + ndealias_, naux_));
    double** Rp = Ra_->pointer();
    double** Rwp = Rw->pointer();
    double* wp = w_->pointer();

    C_DCOPY((nmo_ + ndealias_) * naux_, Rp[0], 1, Rwp[0], 1);

    for (int Q = 0; Q < naux_; Q++) {
        C_DSCAL(nmo_ + ndealias_, wp[Q], &Rwp[0][Q], naux_);
    } 

    boost::shared_ptr<Matrix> C(new Matrix("C", nmo_ + ndealias_, nmo_ + ndealias_));
    double** Cp = C->pointer();

    C_DGEMM('N','T',nmo_ + ndealias_, nmo_ + ndealias_, naux_, 1.0, Rp[0], naux_, Rwp[0], naux_, 
        0.0, Cp[0], nmo_ + ndealias_);

    C->print();

    boost::shared_ptr<Matrix> T(new Matrix("T", nmo_ + ndealias_, nmo2_ + ndealias2_));
    double** Tp = T->pointer();

    double** Xp = X_->pointer();

    C_DGEMM('N','N',nmo_ + ndealias_, nmo2_ + ndealias2_, nmo_ + ndealias_, 1.0, Cp[0], nmo_ + ndealias_,
        Xp[0], nmo2_ + ndealias2_, 0.0, Tp[0], nmo2_ + ndealias2_);

    boost::shared_ptr<Matrix> Cinv(new Matrix("Cinv", nmo2_ + ndealias2_, nmo2_ + ndealias2_));
    double** Cinvp = Cinv->pointer();
    
    C_DGEMM('T','N',nmo2_ + ndealias2_, nmo2_ + ndealias2_, nmo_ + ndealias_, 1.0, Xp[0], nmo2_ + ndealias2_,
        Tp[0], nmo2_ + ndealias2_, 0.0, Cinvp[0], nmo2_ + ndealias2_);

    Cinv->print();
}
void PSTensor::form_Q()
{
    boost::shared_ptr<Matrix> XX(new Matrix("XX^T", nmo_, nmo_ + ndealias_));
    double** XXp = XX->pointer(); 
    double** Xp = X_->pointer();
    
    C_DGEMM('N','T',nmo_,nmo_ + ndealias_,nmo2_ + ndealias2_,1.0,Xp[0],nmo2_ + ndealias2_,
        Xp[0],nmo2_ + ndealias2_,0.0,XXp[0],nmo_ + ndealias_);

    if (debug_ > 1)
        XX->print();

    boost::shared_ptr<Matrix> Rw(new Matrix("Rw", nmo_ + ndealias_, naux_));
    double** Rp = Ra_->pointer();
    double** Rwp = Rw->pointer();
    double* wp = w_->pointer();

    C_DCOPY((nmo_ + ndealias_) * naux_, Rp[0], 1, Rwp[0], 1);

    for (int Q = 0; Q < naux_; Q++) {
        C_DSCAL(nmo_ + ndealias_, wp[Q], &Rwp[0][Q], naux_);
    } 

    Qmo_ = boost::shared_ptr<Matrix> (new Matrix("Qmo", nmo_, naux_));
    double** Qmop = Qmo_->pointer();
    
    C_DGEMM('N','N',nmo_,naux_,nmo_ + ndealias_,1.0,XXp[0],nmo_ + ndealias_,Rwp[0],naux_,0.0,Qmop[0],naux_);

    if (debug_)
        Qmo_->print();
}
boost::shared_ptr<Matrix> PSTensor::Q() 
{
    return Qmo_;    
}
boost::shared_ptr<Matrix> PSTensor::Qocc() 
{
    boost::shared_ptr<Matrix> Q(new Matrix("Qocc", nocc_, naux_));
    double** Qr = Qmo_->pointer();
    double** Qp = Q->pointer();

    C_DCOPY(nocc_ * (ULI) naux_, Qr[0], 1, Qp[0], 1);

    return Q;    
}
boost::shared_ptr<Matrix> PSTensor::Qaocc() 
{
    boost::shared_ptr<Matrix> Q(new Matrix("Qaocc", naocc_, naux_));
    double** Qr = Qmo_->pointer();
    double** Qp = Q->pointer();

    C_DCOPY(naocc_ * (ULI) naux_, Qr[nfocc_], 1, Qp[0], 1);

    return Q;    
}
boost::shared_ptr<Matrix> PSTensor::Qvir() 
{
    boost::shared_ptr<Matrix> Q(new Matrix("Qvir", nvir_, naux_));
    double** Qr = Qmo_->pointer();
    double** Qp = Q->pointer();

    C_DCOPY(nvir_ * (ULI) naux_, Qr[nocc_], 1, Qp[0], 1);

    return Q;    
}
boost::shared_ptr<Matrix> PSTensor::Qavir() 
{
    boost::shared_ptr<Matrix> Q(new Matrix("Qavir", navir_, naux_));
    double** Qr = Qmo_->pointer();
    double** Qp = Q->pointer();

    C_DCOPY(navir_ * (ULI) naux_, Qr[nocc_], 1, Qp[0], 1);

    return Q;    
}
boost::shared_ptr<Matrix> PSTensor::R() 
{
    return Rmo_;    
}
boost::shared_ptr<Matrix> PSTensor::Rocc() 
{
    boost::shared_ptr<Matrix> R(new Matrix("Rocc", nocc_, naux_));
    double** Rr = Rmo_->pointer();
    double** Rp = R->pointer();

    C_DCOPY(nocc_ * (ULI) naux_, Rr[0], 1, Rp[0], 1);

    return R;    
}
boost::shared_ptr<Matrix> PSTensor::Raocc() 
{
    boost::shared_ptr<Matrix> R(new Matrix("Raocc", naocc_, naux_));
    double** Rr = Rmo_->pointer();
    double** Rp = R->pointer();

    C_DCOPY(naocc_ * (ULI) naux_, Rr[nfocc_], 1, Rp[0], 1);

    return R;    
}
boost::shared_ptr<Matrix> PSTensor::Rvir() 
{
    boost::shared_ptr<Matrix> R(new Matrix("Rvir", nvir_, naux_));
    double** Rr = Rmo_->pointer();
    double** Rp = R->pointer();

    C_DCOPY(nvir_ * (ULI) naux_, Rr[nocc_], 1, Rp[0], 1);

    return R;    
}
boost::shared_ptr<Matrix> PSTensor::Ravir() 
{
    boost::shared_ptr<Matrix> R(new Matrix("Ravir", navir_, naux_));
    double** Rr = Rmo_->pointer();
    double** Rp = R->pointer();

    C_DCOPY(navir_ * (ULI) naux_, Rr[nocc_], 1, Rp[0], 1);

    return R;    
}
boost::shared_ptr<Matrix> PSTensor::Aso()
{
    boost::shared_ptr<Matrix> A(new Matrix("Aso",  naux_, nso_ * nso_));
    double** Ap = A->pointer();

    boost::shared_ptr<IntegralFactory> fact(new IntegralFactory(primary_,primary_,primary_,primary_));
    boost::shared_ptr<PseudospectralInt> ints(static_cast<PseudospectralInt*>(fact->ao_pseudospectral()));

    if (use_omega_) {
        ints->set_omega(omega_);
    }

    double* x = grid_->x();
    double* y = grid_->y();
    double* z = grid_->z();

    boost::shared_ptr<Matrix> T(new Matrix("Temp", primary_->nbf(), primary_->nbf()));
    double** Tp = T->pointer();

    for (int P = 0; P < naux_; P++) {
        ints->set_point(x[P], y[P], z[P]);
        T->zero();
        ints->compute(T);

        C_DCOPY(nso_ * nso_, Tp[0], 1, Ap[P], 1);        
    }

    return A;
}
boost::shared_ptr<Matrix> PSTensor::Aoo()
{
    boost::shared_ptr<Matrix> Amn = Aso();
    boost::shared_ptr<Matrix> Ami(new Matrix("Ami", naux_, naocc_ * (ULI) nso_));
    
    double** Amnp = Amn->pointer();
    double** Amip = Ami->pointer();
    double** Cop = Caocc_->pointer();

    C_DGEMM('N','N', naux_ * (ULI) nso_, naocc_, nso_, 1.0, Amnp[0], nso_, Cop[0], naocc_, 
        0.0, Amip[0], naocc_);

    Amn.reset();

    boost::shared_ptr<Matrix> Aia(new Matrix("Aij", naux_, naocc_ * (ULI) naocc_));
    double** Aiap = Aia->pointer();

    for (int Q = 0; Q < naux_; Q++) {
        C_DGEMM('T','N',naocc_,naocc_,nso_,1.0,Amip[Q],naocc_,Cop[0],naocc_, 0.0, Aiap[0], naocc_);
    }    

    return Aia;
}
boost::shared_ptr<Matrix> PSTensor::Aov()
{
    boost::shared_ptr<Matrix> Amn = Aso();
    boost::shared_ptr<Matrix> Ami(new Matrix("Ami", naux_, naocc_ * (ULI) nso_));
    
    double** Amnp = Amn->pointer();
    double** Amip = Ami->pointer();
    double** Cop = Caocc_->pointer();
    double** Cvp = Cavir_->pointer();

    C_DGEMM('N','N', naux_ * (ULI) nso_, naocc_, nso_, 1.0, Amnp[0], nso_, Cop[0], naocc_, 
        0.0, Amip[0], naocc_);

    Amn.reset();

    boost::shared_ptr<Matrix> Aia(new Matrix("Aia", naux_, naocc_ * (ULI) navir_));
    double** Aiap = Aia->pointer();

    for (int Q = 0; Q < naux_; Q++) {
        C_DGEMM('T','N',naocc_,navir_,nso_,1.0,Amip[Q],naocc_,Cvp[0],navir_, 0.0, Aiap[Q], navir_);
    }    

    return Aia;
}
boost::shared_ptr<Matrix> PSTensor::Avv()
{
    boost::shared_ptr<Matrix> Amn = Aso();
    boost::shared_ptr<Matrix> Ami(new Matrix("Ami", naux_, navir_ * (ULI) nso_));
    
    double** Amnp = Amn->pointer();
    double** Amip = Ami->pointer();
    double** Cvp = Cavir_->pointer();

    C_DGEMM('N','N', naux_ * (ULI) nso_, navir_, nso_, 1.0, Amnp[0], nso_, Cvp[0], navir_, 
        0.0, Amip[0], navir_);

    Amn.reset();

    boost::shared_ptr<Matrix> Aia(new Matrix("Aab", naux_, navir_ * (ULI) navir_));
    double** Aiap = Aia->pointer();

    for (int Q = 0; Q < naux_; Q++) {
        C_DGEMM('T','N',navir_,navir_,nso_,1.0,Amip[Q],navir_,Cvp[0],navir_, 0.0, Aiap[Q], navir_);
    }    

    return Aia;
}
boost::shared_ptr<Matrix> PSTensor::Amo()
{
    boost::shared_ptr<Matrix> Amn = Aso();
    boost::shared_ptr<Matrix> Ami(new Matrix("Ami", naux_, nmo_ * (ULI) nso_));
    
    double** Amnp = Amn->pointer();
    double** Amip = Ami->pointer();
    double** Cvp = C_->pointer();

    C_DGEMM('N','N', naux_ * (ULI) nso_, nmo_, nso_, 1.0, Amnp[0], nso_, Cvp[0], nmo_, 
        0.0, Amip[0], nmo_);

    Amn.reset();

    boost::shared_ptr<Matrix> Aia(new Matrix("Amo", naux_, nmo_ * (ULI) nmo_));
    double** Aiap = Aia->pointer();

    for (int Q = 0; Q < naux_; Q++) {
        C_DGEMM('T','N',nmo_,nmo_,nso_,1.0,Amip[Q],nmo_,Cvp[0],nmo_, 0.0, Aiap[Q], nmo_);
    }    

    return Aia;
}
boost::shared_ptr<Matrix> PSTensor::Imo()
{
    boost::shared_ptr<MintsHelper> mints(new MintsHelper());
    return mints->mo_eri(C_,C_);   
}
boost::shared_ptr<Matrix> PSTensor::Ipsmo()
{
    boost::shared_ptr<Matrix> Am = Amo();

    double** Qmop = Qmo_->pointer();
    double** Rmop = Rmo_->pointer();
    double** Amop = Am->pointer();

    boost::shared_ptr<Matrix> Imo(new Matrix("PS MO ERI Tensor", nmo_ * nmo_, nmo_ * nmo_));
    double** Imop = Imo->pointer();

    for (int a = 0; a < nmo_; a++) {
    for (int b = 0; b < nmo_; b++) {
    for (int cd = 0; cd < nmo_ * nmo_; cd++) {
    for (int P = 0; P < naux_; P++) {
        
        Imop[a * nmo_ + b][cd] += Qmop[a][P] * Rmop[b][P] * Amop[P][cd];
    
    }}}}    

    return Imo; 
}

}
