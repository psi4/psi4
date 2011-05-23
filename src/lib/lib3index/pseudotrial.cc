#include "3index.h"
#include <psi4-dec.h>
#include <libmints/mints.h>
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
#ifdef HAVE_MKL
#include <mkl.h>
#endif

//OpenMP Header
#ifdef _OPENMP
#include <omp.h>
#endif

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
    min_S_ = options_.get_double("S_MIN_EIGENVALUE");

    form_molecule();
    fflush(outfile);

    form_bases();
    form_grid();
    fflush(outfile);

    form_Spp();
    form_Spd();
    form_Sdd();
    form_Sa();

    form_Xpp();

    if (do_dealias_) {
        form_Cdp();
        form_Xdd();
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

    if (options_.get_str("DEALIAS_BASIS_CC") == "") {
        do_dealias_ = false;
        ndealias_ = 0;
        naug_ = nso_;
        fprintf(outfile," => Dealias Basis Set <= \n\n");
        fprintf(outfile, "  No dealiasing basis provided.\n\n");
    } else {
        do_dealias_ = true;
        // Dealias 
        molecule_->set_basis_all_atoms(options_.get_str("DEALIAS_BASIS_CC"),"DEALIAS_BASIS");
        dealias_ = BasisSet::construct(parser,molecule_,"DEALIAS_BASIS");  
        ndealias_ = dealias_->nbf();
        naug_ = nso_ + ndealias_;

        fprintf(outfile," => Dealias Basis Set <= \n\n");
        dealias_->print_by_level(outfile,print_);
    }
}

void PseudoTrial::form_grid()
{
    grid_ = shared_ptr<PseudoGrid>(new PseudoGrid(molecule_, options_.get_str("PS_GRID_FILE")));
    grid_->parse(options_.get_str("PS_GRID_FILE"));
    naux_ = grid_->getBlock()->getTruePoints();

    fprintf(outfile," => Pseudospectral Grid [a.u.] <= \n\n");

    double* x = grid_->getBlock()->getX();
    double* y = grid_->getBlock()->getY();
    double* z = grid_->getBlock()->getZ();
    double* w = grid_->getBlock()->getWeights();

    fprintf(outfile," %6s %16s %16s %16s %16s\n","N","x","y","z","w");
    for (int Q = 0; Q < naux_; Q++) {
        fprintf(outfile," %6d %16.10f %16.10f %16.10f %16.10f\n",Q+1,x[Q],y[Q],z[Q],w[Q]);
    }
    fprintf(outfile,"\n"); 

    w_ = shared_ptr<Vector> (new Vector("Grid Weights", naux_));
    double* wp = w_->pointer();

    for (int Q = 0; Q < naux_; Q++)
        wp[Q] = w[Q];
}

void PseudoTrial::form_Spp()
{
    shared_ptr<IntegralFactory> fact(new IntegralFactory(primary_,primary_,primary_,primary_));
    shared_ptr<OneBodyAOInt> Sint(fact->ao_overlap());

    Spp_ = shared_ptr<Matrix>(new Matrix("S (primary x primary)", nso_, nso_));
    Sint->compute(Spp_);

    if (debug_)
        Spp_->print();
}

void PseudoTrial::form_Spd()
{
    if (!do_dealias_) return;

    shared_ptr<IntegralFactory> fact(new IntegralFactory(primary_,dealias_,primary_,primary_));
    shared_ptr<OneBodyAOInt> Sint(fact->ao_overlap());

    Spd_ = shared_ptr<Matrix>(new Matrix("S (primary x dealias)", nso_, ndealias_));
    Sint->compute(Spd_);

    if (debug_)
        Spd_->print();
}

void PseudoTrial::form_Sdd()
{
    if (!do_dealias_) return;

    shared_ptr<IntegralFactory> fact(new IntegralFactory(dealias_,dealias_,primary_,primary_));
    shared_ptr<OneBodyAOInt> Sint(fact->ao_overlap());

    Sdd_ = shared_ptr<Matrix>(new Matrix("S (dealias x dealias)", ndealias_, ndealias_));
    Sint->compute(Sdd_);

    if (debug_)
        Sdd_->print();
}

void PseudoTrial::form_Sa()
{
    Sa_ = shared_ptr<Matrix>(new Matrix("S Augmented, Raw (primary + dealias x primary + dealias)", naug_, naug_));
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

void PseudoTrial::form_Xpp()
{
    
}

void PseudoTrial::form_Cdp()
{
    Cdp_ = shared_ptr<Matrix>(new Matrix("Orthogonalization coefficients (dealias x primary)", ndealias_, nso_));
    double** Cp = Cdp_->pointer();
    double** Spdp = Spd_->pointer();
    for (int i = 0; i < ndealias_; i++)
        C_DCOPY(nso_, &Spdp[0][i], ndealias_, Cp[i], 1);

    if (debug_)
        Cdp_->print(outfile, " Before solution");

    shared_ptr<Matrix> St(new Matrix("S Primary Temp", nso_, nso_));
    St->copy(Spp_);
    double** Stp = St->pointer();

    C_DPOTRF('L', nso_, Stp[0], nso_);
    St->zero_lower();   

    if (debug_)
        St->print(outfile, " after Cholesky");

    C_DPOTRS('L',nso_,ndealias_,Stp[0],nso_,Cp[0],nso_);
    C_DSCAL(nso_*ndealias_, -1.0, Cp[0], 1);

    if (debug_)
        Cdp_->print();
}

void PseudoTrial::form_Xdd()
{
}

void PseudoTrial::form_Rp()
{
    Rp_ = shared_ptr<Matrix>(new Matrix("R (primary x points)", nso_, naux_));
    double** Rp = Rp_->pointer();

    boost::shared_ptr<BasisPoints> points(new BasisPoints(primary_, naux_));
    points->setToComputePoints(true);
    double** bpoints = points->getPoints();

    // Compute the basis points
    points->computePoints(grid_->getBlock());

    // Copy the points in
    for (int i = 0; i < naux_; i++) {
        for (int Q = 0; Q < nso_; Q++)
            Rp[Q][i] = bpoints[i][Q];
    }

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

    Rd_ = shared_ptr<Matrix>(new Matrix("R (dealias x points)", ndealias_, naux_));
    double** Rp = Rd_->pointer();

    boost::shared_ptr<BasisPoints> points(new BasisPoints(dealias_, naux_));
    points->setToComputePoints(true);
    double** bpoints = points->getPoints();

    // Compute the basis points
    points->computePoints(grid_->getBlock());

    // Copy the points in
    for (int i = 0; i < naux_; i++) {
        for (int Q = 0; Q < ndealias_; Q++)
            Rp[Q][i] = bpoints[i][Q];
    }

    if (debug_)
        Rd_->print();
}

void PseudoTrial::form_Rp2()
{
}

void PseudoTrial::form_Rd2()
{
}

void PseudoTrial::form_Ra()
{
    if (!do_dealias_) {
        Ra_ = Rp2_;
        return;
    }

    Ra_ = shared_ptr<Matrix>(new Matrix("R Augmented (primary + dealias x points)", naug_, naux_));
    double** Rap = Ra_->pointer(); 
    double** Rpp = Rp_->pointer();
    double** Rdp = Rd_->pointer();

    double** Cp = Cdp_->pointer();

    C_DCOPY(nso_ * naux_, Rpp[0], 1, Rap[0], 1);
    C_DCOPY(ndealias_ * naux_, Rdp[0], 1, Rap[nso_], 1);

    C_DGEMM('N','N',ndealias_,naux_,nso_,1.0,Cp[0],nso_,Rap[0],naux_,1.0,Rap[nso_],naux_);

    if (debug_)
        Ra_->print();
}

void PseudoTrial::form_Sa2()
{
    Sa2_ = shared_ptr<Matrix>(new Matrix("S Augmented (primary + dealias x primary + dealias)", naug_, naug_));
    double** Sap = Sa2_->pointer();
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
        Sa_->print(outfile, "Before Orthogonalization");

    double** Cp = Cdp_->pointer();

    // pd block
    C_DGEMM('N','N',ndealias_,nso_,nso_,1.0,Cp[0],nso_,Sap[0],naug_,1.0,Sap[nso_],naug_); 
    
    for (int m = 0; m < nso_; m++) {
        C_DCOPY(ndealias_, &Sap[nso_][m], naug_, &Sap[m][nso_], 1);
    }

    // dd block
    C_DGEMM('N','N',ndealias_,ndealias_,nso_,1.0,Cp[0],nso_,Spdp[0],ndealias_,1.0,&Sap[nso_][nso_],naug_);   
    C_DGEMM('T','T',ndealias_,ndealias_,nso_,1.0,Spdp[0],ndealias_,Cp[0],nso_,1.0,&Sap[nso_][nso_],naug_);   

    shared_ptr<Matrix> T(new Matrix("Temp dp matrix", ndealias_, nso_));
    double** Tp = T->pointer();

    C_DGEMM('N','N',ndealias_,nso_,nso_,1.0,Cp[0],nso_,Sppp[0],nso_,0.0,Tp[0],nso_);
    C_DGEMM('N','T',ndealias_,ndealias_,nso_,1.0,Tp[0],nso_,Cp[0],nso_,1.0,&Sap[nso_][nso_],naug_);

    if (debug_)
        Sa_->print(); 
}

void PseudoTrial::form_Q()
{
    C_ = shared_ptr<Matrix>(new Matrix("C Matrix (primary + dealias x primary + dealias", naug_, naug_));
    Cinv_ = shared_ptr<Matrix>(new Matrix("C^-1 Matrix (primary + dealias x primary + dealias", naug_, naug_));
    Qfull_ = shared_ptr<Matrix>(new Matrix("Full Q Matrix (primary + dealias x points", naug_, naux_));
    Q_ = shared_ptr<Matrix>(new Matrix("Q Matrix (primary x points)", nso_, naux_));
    double** Cp = C_->pointer();
    double** Cinvp = Cinv_->pointer();
    double** Qfullp = Qfull_->pointer();
    double** Qp = Q_->pointer();
    double** Pp = P_->pointer();
    double** Rp = Ra_->pointer();
    double* w = w_->pointer();

    shared_ptr<Matrix> Rt(new Matrix("Shared R matrix for scaling",naug_, naux_));
    Rt->copy(Ra_);
    double** Rtp = Rt->pointer();

    if (debug_)
        w_->print();

    for (int Q = 0; Q < naux_; Q++)
        C_DSCAL(naug_, w[Q], &Rtp[0][Q], naux_);

    C_DGEMM('N','T',naug_,naug_,naux_,1.0,Rtp[0],naux_,Rp[0],naux_,0.0,Cp[0],naug_);

    if (debug_)
        C_->print();

    Cinv_->copy(C_);

    C_DPOTRF('L',naug_,Cinvp[0],naug_);
    C_DPOTRI('L',naug_,Cinvp[0],naug_);
    Cinv_->copy_upper_to_lower();
    
    if (debug_)
        Cinv_->print();

    C_DGEMM('N','N',naug_,naux_,naug_,1.0,Cinvp[0],naug_,Rtp[0],naux_,0.0,Qfullp[0],naux_);

    if (debug_)
        Qfull_->print(); 

    C_DGEMM('N','N',nso_,naux_,naug_,1.0,Pp[0],naug_,Qfullp[0],naux_,0.0,Qp[0],naux_);

    if (debug_)
        Q_->print();
}

void PseudoTrial::form_P()
{
    P_ = shared_ptr<Matrix>(new Matrix("Projector Matrix (primary x primary + dealias)", nmo_, naug2_));
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
}
    
void PseudoTrial::form_A()
{
    A_ = shared_ptr<Matrix>(new Matrix("A (primary-primary x points)", nso_ * nso_, naux_));
    double** Ap = A_->pointer();

    shared_ptr<IntegralFactory> fact(new IntegralFactory(primary_,primary_,primary_,primary_));
    boost::shared_ptr<PseudospectralInt> ints(static_cast<PseudospectralInt*>(fact->ao_pseudospectral()));

    boost::shared_ptr<GridBlock> block = grid_->getBlock();
    double* x = block->getX();
    double* y = block->getY();
    double* z = block->getZ();

    boost::shared_ptr<Matrix> T(new Matrix("Temp", primary_->nbf(), primary_->nbf()));
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
    shared_ptr<MintsHelper> mints(new MintsHelper());
    I_ = mints->ao_eri();
    I_->print();
}

void PseudoTrial::form_Ips()
{
    Ips_ = shared_ptr<Matrix>(new Matrix("PS AO ERI Tensor", nso_ * nso_, nso_ * nso_));
    double** Ip = Ips_->pointer();

    T_ = shared_ptr<Matrix>(new Matrix("QR product", nso_ * nso_, naux_));
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
    shared_ptr<Matrix> E(new Matrix("Error in AO TEI tensor", nso_ * nso_ , nso_ * nso_));
    double** Ep = E->pointer();
    double** Ip = I_->pointer();
    double** Ipsp = Ips_->pointer();

    C_DCOPY(nso_*nso_*nso_*nso_,Ipsp[0],1,Ep[0],1);
    C_DAXPY(nso_*nso_*nso_*nso_,-1.0,Ip[0],1,Ep[0],1);

    E->print();
}

}
