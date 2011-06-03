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
                   Options& options) :
    primary_(primary), C_(C), nocc_(nocc), nvir_(nvir),
    naocc_(naocc), navir_(navir), options_(options)
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

    print_header();

    molecule_ = primary_->molecule();

    nfocc_ = nocc_ - naocc_;
    nfvir_ = nvir_ - navir_;

    nso_ = C_->rowspi()[0];
    nmo_ = C_->colspi()[0];

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

    min_S_dealias_ = options_.get_double("PS_MIN_S_DEALIAS");

    buildDealiasSet();
    buildGrid();
    buildR();
    buildQ();
}
void PSTensor::print_header()
{
    fprintf(outfile,"  ==> PS Tensor (by Rob Parrish) <==\n\n");
}
void PSTensor::buildDealiasSet() 
{
    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());

    fprintf(outfile," => Primary Basis Set <= \n\n");
    primary_->print_by_level(outfile,print_);
    fflush(outfile);

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
    ndealias_ = ndealias2_ = dealias_->nbf();
    naug_ = naug2_ = nso_ + ndealias_;
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
    form_Spp();
    form_Spd();
    form_Sdd();

    form_Spd3();
    form_Cdp();
    form_Sdd4();
    form_Xdd();

    // Check if it worked
    if (debug_) {
        form_Sa();
        form_Sa3();
        form_Sa4();
        form_Sa2();
    }

    // Now build the collocation matrices
    form_Rp();
    form_Rd();
    form_Rpmo();
    form_Rdmo();
    form_Ramo();
}
void PSTensor::buildQ()
{
    if (options_.get_str("PS_INVERSION") == "CHOLESKY") {
        form_Q_cholesky();
    } else if (options_.get_str("PS_INVERSION") == "SVD") {
        form_Q_SVD();
    } else if (options_.get_str("PS_INVERSION") == "QR") {
        form_Q_QR();
    }
}
void PSTensor::form_Spp()
{
    boost::shared_ptr<IntegralFactory> fact(new IntegralFactory(primary_,primary_,primary_,primary_));
    boost::shared_ptr<OneBodyAOInt> Sint(fact->ao_overlap());

    Spp_ = boost::shared_ptr<Matrix>(new Matrix("S (primary x primary)", nso_, nso_));
    Sint->compute(Spp_);

    if (debug_)
        Spp_->print();
}
void PSTensor::form_Spd()
{
    boost::shared_ptr<IntegralFactory> fact(new IntegralFactory(primary_,dealias_,primary_,primary_));
    boost::shared_ptr<OneBodyAOInt> Sint(fact->ao_overlap());

    Spd_ = boost::shared_ptr<Matrix>(new Matrix("S (primary x dealias)", nso_, ndealias_));
    Sint->compute(Spd_);

    if (debug_)
        Spd_->print();
}
void PSTensor::form_Sdd()
{
    boost::shared_ptr<IntegralFactory> fact(new IntegralFactory(dealias_,dealias_,primary_,primary_));
    boost::shared_ptr<OneBodyAOInt> Sint(fact->ao_overlap());

    Sdd_ = boost::shared_ptr<Matrix>(new Matrix("S (dealias x dealias)", ndealias_, ndealias_));
    Sint->compute(Sdd_);

    if (debug_)
        Sdd_->print();
}
void PSTensor::form_Spd3()
{
    Spd3_ = boost::shared_ptr<Matrix>(new Matrix("S (primary' x dealias)", nmo_, ndealias_));
    double** Spd3p = Spd3_->pointer();
    double** Spdp = Spd_->pointer();
    double** Xp = C_->pointer();

    C_DGEMM('T','N',nmo_,ndealias_,nso_,1.0,Xp[0],nmo_,Spdp[0],ndealias_,0.0,Spd3p[0],ndealias_);

    if (debug_)
        Spd3_->print();
}
void PSTensor::form_Cdp()
{
    Cdp_ = boost::shared_ptr<Matrix>(new Matrix("Orthogonalization coefficients (dealias x primary')", ndealias_, nmo_));
    double** Cp = Cdp_->pointer();
    double** Sp = Spd3_->pointer();

    for (int i = 0; i < ndealias_; i++)
        C_DCOPY(nmo_,&Sp[0][i],ndealias_,Cp[i],1);

    Cdp_->scale(-1.0);

    if (debug_)
        Cdp_->print();
}
void PSTensor::form_Sdd4()
{
    Sdd4_ = boost::shared_ptr<Matrix>(new Matrix("S Separated (dealias x dealias)", ndealias_, ndealias_));
    double** Sdd4p = Sdd4_->pointer();
    double** Spdp = Spd3_->pointer();
    double** Cp = Cdp_->pointer();

    Sdd4_->copy(Sdd_);

    C_DGEMM('T','T',ndealias_,ndealias_,nmo_,1.0,Spdp[0],ndealias_,Cp[0],nmo_,1.0,Sdd4p[0],ndealias_);
    C_DGEMM('N','N',ndealias_,ndealias_,nmo_,1.0,Cp[0],nmo_,Spdp[0],ndealias_,1.0,Sdd4p[0],ndealias_);
    C_DGEMM('N','T',ndealias_,ndealias_,nmo_,1.0,Cp[0],nmo_,Cp[0],nmo_,1.0,Sdd4p[0],ndealias_);

    if (debug_)
        Sdd4_->print();
}
void PSTensor::form_Xdd()
{
    boost::shared_ptr<Matrix> St(new Matrix("Temporary S", ndealias_, ndealias_));        
    boost::shared_ptr<Matrix> Xt(new Matrix("Temporary X", ndealias_, ndealias_));        
    boost::shared_ptr<Vector> st(new Vector("s", ndealias_));

    double** Stp = St->pointer();
    double** Xtp = Xt->pointer();
    double* stp = st->pointer();

    St->copy(Sdd4_);

    St->diagonalize(Xt, st);

    if (debug_) 
        Xt->eivprint(st);

    ndealias2_ = 0;
    for (int i = 0; i < ndealias_; i++)
    {
        if (stp[i] > min_S_dealias_)
            ndealias2_++;
    }
    naug2_ = nmo_ + ndealias2_;

    Xdd_ = boost::shared_ptr<Matrix>(new Matrix("X Matrix (dealias x dealias')", ndealias_, ndealias2_)); 
    double** Xp = Xdd_->pointer();    
 
    int m = 0;
    for (int i = 0; i < ndealias_; i++) {
        if (stp[i] > min_S_dealias_) {
            C_DCOPY(ndealias_, &Xtp[0][i], ndealias_, &Xp[0][m], ndealias2_);
            C_DSCAL(ndealias_, pow(stp[i], -1.0 / 2.0), &Xp[0][m], ndealias2_);
            m++;
        }
    }

    if (debug_)
        Xdd_->print();
}
void PSTensor::form_Sa()
{
    Sa_ = boost::shared_ptr<Matrix>(new Matrix("S Augmented, Raw (primary + dealias x primary + dealias)", naug_, naug_));
    double** Sap = Sa_->pointer();
    double** Sppp = Spp_->pointer();
    double** Spdp = Spd_->pointer();
    double** Sddp = Sdd_->pointer();

    for (int m = 0; m < nso_; m++) {
        C_DCOPY(nso_, Sppp[m], 1, Sap[m], 1);
    } 

    for (int m = 0; m < nso_; m++) {
        C_DCOPY(ndealias_, Spdp[m], 1, &Sap[m][nso_], 1);
    } 
    
    for (int m = 0; m < nso_; m++) {
        C_DCOPY(ndealias_, Spdp[m], 1, &Sap[nso_][m], naug_);
    }

    for (int a = 0; a < ndealias_; a++) {
        C_DCOPY(ndealias_, Sddp[a], 1, &Sap[nso_ + a][nso_], 1);
    } 
   
    if (debug_)
        Sa_->print();
}
void PSTensor::form_Sa3()
{
    Sa3_ = boost::shared_ptr<Matrix>(new Matrix("S3 Augmented, Raw (primary' + dealias x primary' + dealias)", nmo_ + ndealias_, nmo_ + ndealias_));

    double** Sap = Sa3_->pointer();

    double** Sppp = Spp_->pointer();
    double** Xp   = C_->pointer();
    double** Spdp = Spd_->pointer();
    double** Sddp = Sdd_->pointer();

    boost::shared_ptr<Matrix> T(new Matrix("Temp",nmo_,nso_));
    double** Tp = T->pointer();

    C_DGEMM('T','N',nmo_,nso_,nso_,1.0,Xp[0],nmo_,Sppp[0],nso_,0.0,Tp[0],nso_);
    C_DGEMM('N','N',nmo_,nmo_,nso_,1.0,Tp[0],nso_,Xp[0],nmo_,0.0,Sap[0],nmo_ + ndealias_);

    C_DGEMM('T','N',nmo_,ndealias_,nso_,1.0,Xp[0],nmo_,Spdp[0],ndealias_,0.0,&Sap[0][nmo_], nmo_ + ndealias_);
    C_DGEMM('T','N',ndealias_,nmo_,nso_,1.0,Spdp[0],ndealias_,Xp[0],nmo_,0.0,&Sap[nmo_][0], nmo_ + ndealias_);
    
    for (int a = 0; a < ndealias_; a++) {
        C_DCOPY(ndealias_, Sddp[a], 1, &Sap[nmo_ + a][nmo_], 1);
    } 
   
    if (debug_)
        Sa3_->print();
}
void PSTensor::form_Sa4()
{
    Sa4_ = boost::shared_ptr<Matrix>(new Matrix("S4 Augmented, Raw (primary' + dealias x primary' + dealias)", nmo_ + ndealias_, nmo_ + ndealias_));
    Sa4_->copy(Sa3_);

    double** Sap = Sa4_->pointer();
    double** Sppp = Spp_->pointer();
    double** Spdp = Spd3_->pointer();
    double** Sddp = Sdd_->pointer();

    double** Cp = Cdp_->pointer();

    C_DGEMM('N','T',nmo_,ndealias_,nmo_,1.0,Sap[0], nmo_ + ndealias_, Cp[0], nmo_,1.0,&Sap[0][nmo_], nmo_ + ndealias_);
    C_DGEMM('N','N',ndealias_,nmo_,nmo_,1.0,Cp[0],nmo_,Sap[0],nmo_ + ndealias_,1.0,&Sap[nmo_][0], nmo_ + ndealias_);

    C_DGEMM('T','T',ndealias_,ndealias_,nmo_,1.0,Spdp[0],ndealias_,Cp[0],nmo_,1.0,&Sap[nmo_][nmo_], nmo_ + ndealias_); 
    C_DGEMM('N','N',ndealias_,ndealias_,nmo_,1.0,Cp[0],nmo_,Spdp[0],ndealias_,1.0,&Sap[nmo_][nmo_], nmo_ + ndealias_); 
    C_DGEMM('N','T',ndealias_,ndealias_,nmo_,1.0,Cp[0],nmo_,Cp[0],nmo_,1.0,&Sap[nmo_][nmo_], nmo_ + ndealias_); 
 
    if (debug_)
        Sa4_->print();
}
void PSTensor::form_Sa2()
{
    Sa2_ = boost::shared_ptr<Matrix>(new Matrix("S2 Augmented, Finished (primary' + dealias' x primary' + dealias')", naug2_, naug2_));

    double** Sap = Sa2_->pointer();

    double** Sppp = Sa3_->pointer();
    double** Sddp = Sdd4_->pointer();

    for (int i = 0; i < nmo_; i++)
        C_DCOPY(nmo_,Sppp[i],1,Sap[i],1);

    boost::shared_ptr<Matrix> T(new Matrix("Temp", ndealias2_, ndealias_));
    double** Tp = T->pointer();

    double** Xp = Xdd_->pointer();

    C_DGEMM('T','N',ndealias2_,ndealias_,ndealias_,1.0,Xp[0],ndealias2_,Sddp[0],ndealias_,0.0,Tp[0],ndealias_);
    C_DGEMM('N','N',ndealias2_,ndealias2_,ndealias_,1.0,Tp[0],ndealias_,Xp[0],ndealias2_,0.0,&Sap[nmo_][nmo_],naug2_);

    if (debug_)
        Sa2_->print();
}
void PSTensor::form_Rp()
{
    Rp_ = boost::shared_ptr<Matrix>(new Matrix("R (primary x points)", nso_, naux_));
    double** Rp = Rp_->pointer();

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
        Rp_->print();
}
void PSTensor::form_Rd()
{
    Rd_ = boost::shared_ptr<Matrix>(new Matrix("R (dealias x points)", ndealias_, naux_));
    double** Rp = Rd_->pointer();

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
        Rd_->print();
}
void PSTensor::form_Rpmo()
{
    Rpmo_ = boost::shared_ptr<Matrix>(new Matrix("R2 (primary' x points)", nmo_, naux_));
    double** Rp2 = Rpmo_->pointer();
    double** Rp = Rp_->pointer();
    double** Xp = C_->pointer();

    C_DGEMM('T','N',nmo_,naux_,nso_,1.0,Xp[0],nmo_,Rp[0],naux_,0.0,Rp2[0],naux_);

    if (debug_)
        Rpmo_->print();

    Rmo_ = Rpmo_;
}
void PSTensor::form_Rdmo()
{
    Rdmo_ = boost::shared_ptr<Matrix>(new Matrix("R2 (dealias' x points)", ndealias2_, naux_));
    double** Rd2p = Rdmo_->pointer();
    double** Rp2p = Rpmo_->pointer();
    double** Rdp = Rd_->pointer();

    double** Xdp = Xdd_->pointer();
    double** Cdp = Cdp_->pointer();

    C_DGEMM('T','N',ndealias2_,naux_,ndealias_,1.0,Xdp[0],ndealias2_,Rdp[0],naux_,0.0,Rd2p[0],naux_);

    boost::shared_ptr<Matrix> T(new Matrix("Temp",ndealias_,naux_));
    double** Tp = T->pointer();

    C_DGEMM('N','N',ndealias_,naux_,nmo_,1.0,Cdp[0],nmo_,Rp2p[0],naux_,0.0,Tp[0],naux_);
    C_DGEMM('T','N',ndealias2_,naux_,ndealias_,1.0,Xdp[0],ndealias2_,Tp[0],naux_,1.0,Rd2p[0],naux_);

    if (debug_)
        Rdmo_->print();
}
void PSTensor::form_Ramo()
{
    Ra_ = boost::shared_ptr<Matrix>(new Matrix("R Augmented (primary' + dealias' x points)", naug2_, naux_));
    double** Rap = Ra_->pointer(); 

    double** Rpp = Rpmo_->pointer();
    double** Rdp = Rdmo_->pointer();

    C_DCOPY(nmo_ * naux_, Rpp[0], 1, Rap[0], 1);
    C_DCOPY(ndealias2_ * naux_, Rdp[0], 1, Rap[nmo_], 1);

    if (debug_)
        Ra_->print();
}
void PSTensor::form_Q_cholesky()
{
    boost::shared_ptr<Matrix> T(new Matrix("Rw (points x  naug')", naux_, naug2_));
    double** Rp = Ra_->pointer();
    double** Tp = T->pointer();
    double*  wp = w_->pointer();
   
    for (int Q = 0; Q < naux_; Q++) {
        C_DCOPY(naug2_, &Rp[0][Q], naux_, Tp[Q], 1);
        C_DSCAL(naug2_, wp[Q], Tp[Q], 1);
    } 

    boost::shared_ptr<Matrix> C(new Matrix("C (naug' x naug')", naug2_, naug2_));
    double** Cp = C->pointer(); 

    C_DGEMM('N','N',naug2_, naug2_, naux_, 1.0, Rp[0], naux_, Tp[0], naug2_, 0.0, Cp[0], naug2_);

    if (debug_) 
        C->print();

    int info = C_DPOTRF('L',naug2_,Cp[0],naug2_);

    if (info != 0)
        throw PSIEXCEPTION("PSTensor::form_Q_cholesky: C_DPOTRF failed.");

    info = C_DPOTRS('L',naug2_,naux_,Cp[0],naug2_,Tp[0],naug2_);

    if (info != 0)
        throw PSIEXCEPTION("PSTensor::form_Q_cholesky: C_DPOTRS failed.");

    Qmo_ = boost::shared_ptr<Matrix>(new Matrix("Q2 (primary' x points", nmo_, naux_));
    double** Qp = Qmo_->pointer();
         
    for (int i = 0; i < nmo_; i++) {
        C_DCOPY(naux_, &Tp[0][i], naug2_, Qp[i], 1);
    }
    
    if (debug_)
        Qmo_->print();
}
void PSTensor::form_Q_SVD()
{
    throw FeatureNotImplemented("PSTensor","form_Q_SVD",__FILE__,__LINE__);
}
void PSTensor::form_Q_QR()
{
    throw FeatureNotImplemented("PSTensor","form_Q_QR",__FILE__,__LINE__);
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
