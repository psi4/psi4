#include "3index.h"
#include <psi4-dec.h>
#include <libmints/mints.h>
#include <libfock/cubature.h>
#include <libfock/points.h>
#include <libqt/qt.h>

#include <string>
#include <sstream>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <algorithm>
#include <utility>
#include <ctype.h>

//MKL Header
#include <psiconfig.h>
#if HAVE_MKL
#include <mkl.h>
#endif

//OpenMP Header
#if _OPENMP
#include <omp.h>
#endif

#include <boost/shared_ptr.hpp>

using namespace boost;
using namespace std;
using namespace psi;

namespace psi {

PseudoTrial::PseudoTrial() : 
    options_(Process::environment.options)
{
    common_init();
}

PseudoTrial::~PseudoTrial()
{
}

void PseudoTrial::common_init()
{
    print_header();
    fflush(outfile);

    debug_ = options_.get_int("DEBUG");
    print_ = options_.get_int("PRINT");   
    min_S_primary_ = options_.get_double("PS_MIN_S_PRIMARY");
    min_S_dealias_ = options_.get_double("PS_MIN_S_DEALIAS");

    form_molecule();
    fflush(outfile);

    form_bases();
    form_grid();
    fflush(outfile);

    form_Spp();
    form_Spd();
    form_Sdd();

    form_Xpp();

    if (do_dealias_) {

        form_Spd3();
        form_Cdp();
        form_Sdd4();
        form_Xdd();

        // Verify
        form_Sa();
        form_Sa3();
        form_Sa4();
        form_Sa2();
    }

    form_Rp();
    form_Rd();
    form_Rp2();
    form_Rd2();
    form_Ra();

    form_P();
    form_SX();

    form_Q();
    form_A();

    form_Ips();
    form_I();
    verify();
}

void PseudoTrial::print_header()
{
    fprintf(outfile,"\t\t--------------------------------------------------\n");
    fprintf(outfile,"\t\t                                                  \n");
    fprintf(outfile,"\t\t      PseudoTrial: A Pseudospectral Sandbox       \n");
    fprintf(outfile,"\t\t                  Rob Parrish                     \n");
    fprintf(outfile,"\t\t                  21 May 2011                     \n");
    fprintf(outfile,"\t\t                                                  \n");
    fprintf(outfile,"\t\t--------------------------------------------------\n\n");
    fflush(outfile);
}

void PseudoTrial::form_molecule()
{
    molecule_ = Process::environment.molecule(); 
    fprintf(outfile," => Molecule <= \n\n");
    molecule_->print();
}

void PseudoTrial::form_bases()
{
    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());

    // Primary
    molecule_->set_basis_all_atoms(options_.get_str("BASIS"),"BASIS");
    primary_ = BasisSet::construct(parser,molecule_,"BASIS");  
    nso_ = primary_->nbf();   
 
    fprintf(outfile," => Primary Basis Set <= \n\n");
    primary_->print_by_level(outfile,print_);

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
    do_dealias_ = true;
    ndealias_ = dealias_->nbf();
    naug_ = nso_ + ndealias_;

    dealias_->print_by_level(outfile,print_);
}

void PseudoTrial::form_grid()
{

    if (options_.get_str("PS_GRID_FILE") != "") {
        grid_ = boost::shared_ptr<PseudospectralGrid>(new PseudospectralGrid(molecule_, primary_, options_.get_str("PS_GRID_FILE"), options_));
    } else {
        grid_ = boost::shared_ptr<PseudospectralGrid>(new PseudospectralGrid(molecule_, primary_, options_));
    }

    grid_->print();

    naux_ = grid_->npoints();

    double* w = grid_->w();

    w_ = boost::shared_ptr<Vector> (new Vector("Grid Weights", naux_));
    double* wp = w_->pointer();

    for (int Q = 0; Q < naux_; Q++)
        wp[Q] = w[Q];
}

void PseudoTrial::form_Spp()
{
    boost::shared_ptr<IntegralFactory> fact(new IntegralFactory(primary_,primary_,primary_,primary_));
    boost::shared_ptr<OneBodyAOInt> Sint(fact->ao_overlap());

    Spp_ = SharedMatrix(new Matrix("S (primary x primary)", nso_, nso_));
    Sint->compute(Spp_);

    if (debug_)
        Spp_->print();
}

void PseudoTrial::form_Spd()
{
    if (!do_dealias_) return;

    boost::shared_ptr<IntegralFactory> fact(new IntegralFactory(primary_,dealias_,primary_,primary_));
    boost::shared_ptr<OneBodyAOInt> Sint(fact->ao_overlap());

    Spd_ = SharedMatrix(new Matrix("S (primary x dealias)", nso_, ndealias_));
    Sint->compute(Spd_);

    if (debug_)
        Spd_->print();
}

void PseudoTrial::form_Sdd()
{
    if (!do_dealias_) return;

    boost::shared_ptr<IntegralFactory> fact(new IntegralFactory(dealias_,dealias_,primary_,primary_));
    boost::shared_ptr<OneBodyAOInt> Sint(fact->ao_overlap());

    Sdd_ = SharedMatrix(new Matrix("S (dealias x dealias)", ndealias_, ndealias_));
    Sint->compute(Sdd_);

    if (debug_)
        Sdd_->print();
}

void PseudoTrial::form_Spd3()
{
    Spd3_ = SharedMatrix(new Matrix("S (primary' x dealias)", nmo_, ndealias_));
    double** Spd3p = Spd3_->pointer();
    double** Spdp = Spd_->pointer();
    double** Xp = Xpp_->pointer();

    C_DGEMM('T','N',nmo_,ndealias_,nso_,1.0,Xp[0],nmo_,Spdp[0],ndealias_,0.0,Spd3p[0],ndealias_);

    if (debug_)
        Spd3_->print();
}

void PseudoTrial::form_Sdd4()
{
    Sdd4_ = SharedMatrix(new Matrix("S Separated (dealias x dealias)", ndealias_, ndealias_));
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

void PseudoTrial::form_Xpp()
{
    SharedMatrix St(new Matrix("Temporary S", nso_, nso_));        
    SharedMatrix Xt(new Matrix("Temporary X", nso_, nso_));        
    boost::shared_ptr<Vector> st(new Vector("s", nso_));

    double** Stp = St->pointer();
    double** Xtp = Xt->pointer();
    double* stp = st->pointer();

    St->copy(Spp_);

    St->diagonalize(Xt, st);

    if (debug_) 
        Xt->eivprint(st);

    nmo_ = 0;
    for (int i = 0; i < nso_; i++)
    {
        if (stp[i] > min_S_primary_)
            nmo_++;
    }

    Xpp_ = SharedMatrix(new Matrix("X Matrix (primary x primary')", nso_, nmo_)); 
    double** Xp = Xpp_->pointer();    
 
    int m = 0;
    for (int i = 0; i < nso_; i++) {
        if (stp[i] > min_S_primary_) {
            C_DCOPY(nso_, &Xtp[0][i], nso_, &Xp[0][m], nmo_);
            C_DSCAL(nso_, pow(stp[i], -1.0 / 2.0), &Xp[0][m], nmo_);
            m++;
        }
    }

    if (debug_)
        Xpp_->print();

    ndealias2_ = 0;
    naug2_ = nmo_;
}

void PseudoTrial::form_Cdp()
{
    Cdp_ = SharedMatrix(new Matrix("Orthogonalization coefficients (dealias x primary')", ndealias_, nmo_));
    double** Cp = Cdp_->pointer();
    double** Sp = Spd3_->pointer();

    for (int i = 0; i < ndealias_; i++)
        C_DCOPY(nmo_,&Sp[0][i],ndealias_,Cp[i],1);

    Cdp_->scale(-1.0);

    if (debug_)
        Cdp_->print();
}

void PseudoTrial::form_Xdd()
{
    if (!do_dealias_) {
        ndealias2_ = 0;
        naug2_ = nmo_;
        return;
    }

    SharedMatrix St(new Matrix("Temporary S", ndealias_, ndealias_));        
    SharedMatrix Xt(new Matrix("Temporary X", ndealias_, ndealias_));        
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

    Xdd_ = SharedMatrix(new Matrix("X Matrix (dealias x dealias')", ndealias_, ndealias2_)); 
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

void PseudoTrial::form_Sa()
{
    Sa_ = SharedMatrix(new Matrix("S Augmented, Raw (primary + dealias x primary + dealias)", naug_, naug_));
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

void PseudoTrial::form_Sa3()
{
    Sa3_ = SharedMatrix(new Matrix("S3 Augmented, Raw (primary' + dealias x primary' + dealias)", nmo_ + ndealias_, nmo_ + ndealias_));

    double** Sap = Sa3_->pointer();

    double** Sppp = Spp_->pointer();
    double** Xp   = Xpp_->pointer();
    double** Spdp = Spd_->pointer();
    double** Sddp = Sdd_->pointer();

    SharedMatrix T(new Matrix("Temp",nmo_,nso_));
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

void PseudoTrial::form_Sa4()
{
    Sa4_ = SharedMatrix(new Matrix("S4 Augmented, Raw (primary' + dealias x primary' + dealias)", nmo_ + ndealias_, nmo_ + ndealias_));
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

void PseudoTrial::form_Sa2()
{
    Sa2_ = SharedMatrix(new Matrix("S2 Augmented, Finished (primary' + dealias' x primary' + dealias')", naug2_, naug2_));

    double** Sap = Sa2_->pointer();

    double** Sppp = Sa3_->pointer();
    double** Sddp = Sdd4_->pointer();

    for (int i = 0; i < nmo_; i++)
        C_DCOPY(nmo_,Sppp[i],1,Sap[i],1);

    SharedMatrix T(new Matrix("Temp", ndealias2_, ndealias_));
    double** Tp = T->pointer();

    double** Xp = Xdd_->pointer();

    C_DGEMM('T','N',ndealias2_,ndealias_,ndealias_,1.0,Xp[0],ndealias2_,Sddp[0],ndealias_,0.0,Tp[0],ndealias_);
    C_DGEMM('N','N',ndealias2_,ndealias2_,ndealias_,1.0,Tp[0],ndealias_,Xp[0],ndealias2_,0.0,&Sap[nmo_][nmo_],naug2_);

    if (debug_)
        Sa2_->print();
}

void PseudoTrial::form_Rp()
{
    Rp_ = SharedMatrix(new Matrix("R (primary x points)", nso_, naux_));
    double** Rp = Rp_->pointer();

    #if 0

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

    #endif

    if (debug_)
        Rp_->print();

    R_ = Rp_;
}

void PseudoTrial::form_Rd()
{
    if (!do_dealias_) {
        Rd_ = Rp_;
        return;
    }

    Rd_ = SharedMatrix(new Matrix("R (dealias x points)", ndealias_, naux_));
    double** Rp = Rd_->pointer();

    #if 0

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
    
    #endif

    if (debug_)
        Rd_->print();
}

void PseudoTrial::form_Rp2()
{
    Rp2_ = SharedMatrix(new Matrix("R2 (primary' x points)", nmo_, naux_));
    double** Rp2 = Rp2_->pointer();
    double** Rp = Rp_->pointer();
    double** Xp = Xpp_->pointer();

    C_DGEMM('T','N',nmo_,naux_,nso_,1.0,Xp[0],nmo_,Rp[0],naux_,0.0,Rp2[0],naux_);

    if (debug_)
        Rp2_->print();

    R_ = Rp_;
}

void PseudoTrial::form_Rd2()
{
    if (!do_dealias_) {
        Rd2_ = Rp2_;
        return;
    }

    Rd2_ = SharedMatrix(new Matrix("R2 (dealias' x points)", ndealias2_, naux_));
    double** Rd2p = Rd2_->pointer();
    double** Rp2p = Rp2_->pointer();
    double** Rdp = Rd_->pointer();

    double** Xdp = Xdd_->pointer();
    double** Cdp = Cdp_->pointer();

    C_DGEMM('T','N',ndealias2_,naux_,ndealias_,1.0,Xdp[0],ndealias2_,Rdp[0],naux_,0.0,Rd2p[0],naux_);

    SharedMatrix T(new Matrix("Temp",ndealias_,naux_));
    double** Tp = T->pointer();

    C_DGEMM('N','N',ndealias_,naux_,nmo_,1.0,Cdp[0],nmo_,Rp2p[0],naux_,0.0,Tp[0],naux_);
    C_DGEMM('T','N',ndealias2_,naux_,ndealias_,1.0,Xdp[0],ndealias2_,Tp[0],naux_,1.0,Rd2p[0],naux_);

    if (debug_)
        Rd2_->print();
}

void PseudoTrial::form_Ra()
{
    if (!do_dealias_) {
        Ra_ = Rp2_;
        return;
    }

    Ra_ = SharedMatrix(new Matrix("R Augmented (primary' + dealias' x points)", naug2_, naux_));
    double** Rap = Ra_->pointer(); 

    double** Rpp = Rp2_->pointer();
    double** Rdp = Rd2_->pointer();

    C_DCOPY(nmo_ * naux_, Rpp[0], 1, Rap[0], 1);
    C_DCOPY(ndealias2_ * naux_, Rdp[0], 1, Rap[nmo_], 1);

    if (debug_)
        Ra_->print();
}

void PseudoTrial::form_Q()
{
    C_ = SharedMatrix(new Matrix("C Matrix (primary' + dealias' x primary' + dealias'", naug2_, naug2_));
    Cinv_ = SharedMatrix(new Matrix("C^-1 Matrix (primary' + dealias' x primary' + dealias'", naug2_, naug2_));
    Qfull_ = SharedMatrix(new Matrix("Full Q Matrix (primary' + dealias' x points", naug2_, naux_));
    Qmo_ = SharedMatrix(new Matrix("Q Matrix (primary' x points)", nmo_, naux_));
    Q_ = SharedMatrix(new Matrix("Q Matrix (primary x points)", nso_, naux_));

    double** Cp = C_->pointer();
    double** Cinvp = Cinv_->pointer();
    double** Qfullp = Qfull_->pointer();
    double** Qmop = Qmo_->pointer();
    double** Qp = Q_->pointer();
    double** Pp = P_->pointer();
    double** SXp = SX_->pointer();
    double** Rp = Ra_->pointer();
    double* w = w_->pointer();

    SharedMatrix Rt(new Matrix("Shared R matrix for scaling",naug2_, naux_));
    Rt->copy(Ra_);
    double** Rtp = Rt->pointer();

    if (debug_)
        w_->print();

    for (int Q = 0; Q < naux_; Q++)
        C_DSCAL(naug2_, w[Q], &Rtp[0][Q], naux_);

    C_DGEMM('N','T',naug2_,naug2_,naux_,1.0,Rtp[0],naux_,Rp[0],naux_,0.0,Cp[0],naug2_);

    if (debug_)
        C_->print();

    Cinv_->copy(C_);

    C_DPOTRF('L',naug2_,Cinvp[0],naug2_);
    C_DPOTRI('L',naug2_,Cinvp[0],naug2_);
    Cinv_->copy_upper_to_lower();
    
    if (debug_)
        Cinv_->print();

    C_DGEMM('N','N',naug2_,naux_,naug2_,1.0,Cinvp[0],naug2_,Rtp[0],naux_,0.0,Qfullp[0],naux_);

    if (debug_)
        Qfull_->print(); 

    C_DGEMM('N','N',nmo_,naux_,naug2_,1.0,Pp[0],naug2_,Qfullp[0],naux_,0.0,Qmop[0],naux_);

    if (debug_)
        Qmo_->print();

    C_DGEMM('N','N',nso_,naux_,nmo_,1.0,SXp[0],nmo_,Qmop[0],naux_,0.0,Qp[0],naux_);

    if (debug_)
        Q_->print();
}

void PseudoTrial::form_P()
{
    P_ = SharedMatrix(new Matrix("Projector Matrix (primary' x primary' + dealias')", nmo_, naug2_));
    double** Pp = P_->pointer();

    // First try: [1 0]
    for (int i = 0; i < nmo_; i++)
        Pp[i][i] = 1.0;

    // Next try: [S 0]
    //double** Sp = Spp_->pointer();
    //for (int i = 0; i < nso_; i++)
    //    C_DCOPY(nso_,Sp[i],1,Pp[i],1);

    if (debug_)
        P_->print();
}

void PseudoTrial::form_SX()
{
    SX_ = SharedMatrix(new Matrix("SX (primary x primary')", nso_, nmo_));

    double** SXp = SX_->pointer();
    double** Sp = Spp_->pointer();
    double** Xp = Xpp_->pointer();
    
    C_DGEMM('N','N',nso_,nmo_,nso_,1.0,Sp[0],nso_,Xp[0],nmo_,0.0,SXp[0],nmo_);

    if (debug_)
        SX_->pointer();
}   
    
void PseudoTrial::form_A()
{
    A_ = SharedMatrix(new Matrix("A (primary-primary x points)", nso_ * nso_, naux_));
    double** Ap = A_->pointer();

    boost::shared_ptr<IntegralFactory> fact(new IntegralFactory(primary_,primary_,primary_,primary_));
    boost::shared_ptr<PseudospectralInt> ints(static_cast<PseudospectralInt*>(fact->ao_pseudospectral()));

    double* x = grid_->x();
    double* y = grid_->y();
    double* z = grid_->z();

    SharedMatrix T(new Matrix("Temp", primary_->nbf(), primary_->nbf()));
    double** Tp = T->pointer();

    for (int P = 0; P < naux_; P++) {
        ints->set_point(x[P], y[P], z[P]);
        T->zero();
        ints->compute(T);

        C_DCOPY(nso_ * nso_, Tp[0], 1, &Ap[0][P], naux_);        
    }

    if (debug_)
        A_->print();
}

void PseudoTrial::form_I()
{
    boost::shared_ptr<MintsHelper> mints(new MintsHelper());
    I_ = mints->ao_eri();
    I_->print();
}

void PseudoTrial::form_Ips()
{
    Ips_ = SharedMatrix(new Matrix("PS AO ERI Tensor", nso_ * nso_, nso_ * nso_));
    double** Ip = Ips_->pointer();

    T_ = SharedMatrix(new Matrix("QR product", nso_ * nso_, naux_));
    double** Tp = T_->pointer();

    double** Qp = Q_->pointer();
    double** Rp = R_->pointer();
    double** Ap = A_->pointer();

    for (int m = 0, mn = 0; m < nso_; m++) {
        for (int n = 0; n < nso_; n++, mn++) {
            for (int Q = 0; Q < naux_; Q++) {
                Tp[mn][Q] = Qp[m][Q] * Rp[n][Q];
            }
        }
    }

    if (debug_)
        T_->print();

    C_DGEMM('N','T',nso_*nso_,nso_*nso_,naux_,1.0,Tp[0],naux_,Ap[0],naux_,0.0,Ip[0],nso_*nso_); 
   
    Ips_->print(); 
}

void PseudoTrial::verify()
{
    SharedMatrix E(new Matrix("Error in AO TEI tensor", nso_ * nso_ , nso_ * nso_));
    double** Ep = E->pointer();
    double** Ip = I_->pointer();
    double** Ipsp = Ips_->pointer();

    C_DCOPY(nso_*nso_*nso_*nso_,Ipsp[0],1,Ep[0],1);
    C_DAXPY(nso_*nso_*nso_*nso_,-1.0,Ip[0],1,Ep[0],1);

    E->print();
}

SharedMatrix PseudoTrial::getI() const { return I_; }    
SharedMatrix PseudoTrial::getIPS() const { return Ips_; }    

SharedMatrix PseudoTrial::getQ() const { return Q_; }    
SharedMatrix PseudoTrial::getR() const { return R_; }    
SharedMatrix PseudoTrial::getA() const { return A_; }    

}
