#include "3index.h"
#include <libmints/mints.h>
#include <libqt/qt.h>
#include <boost/regex.hpp>
#include <boost/xpressive/xpressive.hpp>
#include <boost/xpressive/regex_actions.hpp>
#include <boost/algorithm/string.hpp>

#include <string>
#include <sstream>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <algorithm>
#include <utility>
#include <ctype.h>

using namespace std;
using namespace psi;

namespace psi {

Pseudospectral::Pseudospectral(shared_ptr<PSIO> psio, shared_ptr<BasisSet> primary, shared_ptr<BasisSet> dealias, shared_ptr<PseudoGrid> g) :
    psio_(psio), primary_(primary), dealias_(dealias), grid_(g), npoints_(g->getBlock()->getTruePoints())
{
} 
Pseudospectral::~Pseudospectral()
{
}
shared_ptr<Matrix> Pseudospectral::form_X()
{
    shared_ptr<Matrix> X(new Matrix("Primary Collocation Metric (Dirac)", primary_->nbf(), npoints_));
    double** Xp = X->pointer();
    
    shared_ptr<BasisPoints> points(new BasisPoints(primary_, npoints_));    
    points->setToComputePoints(true);
    double** bpoints = points->getPoints();

    // Compute the basis points
    points->computePoints(grid_->getBlock());
 
    // Copy the points in
    for (int i = 0; i < npoints_; i++) {
        for (int Q = 0; Q < primary_->nbf(); Q++)
            Xp[Q][i] = bpoints[i][Q]; 
    }

    return X;
}
shared_ptr<Matrix> Pseudospectral::form_X_dealias()
{
    shared_ptr<Matrix> X(new Matrix("Dealias Collocation Metric (Dirac)", dealias_->nbf(), npoints_));
    double** Xp = X->pointer();
    
    shared_ptr<BasisPoints> points(new BasisPoints(dealias_, npoints_));    
    points->setToComputePoints(true);
    double** bpoints = points->getPoints();

    // Compute the basis points
    points->computePoints(grid_->getBlock());
 
    // Copy the points in
    for (int i = 0; i < npoints_; i++) {
        for (int Q = 0; Q < dealias_->nbf(); Q++)
            Xp[Q][i] = bpoints[i][Q]; 
    }

    return X;
}
shared_ptr<Matrix> Pseudospectral::form_S()
{
    shared_ptr<IntegralFactory> fact(new IntegralFactory(primary_, primary_, primary_, primary_));
    shared_ptr<OneBodyAOInt> o2(fact->ao_overlap());
    shared_ptr<Matrix> S(new Matrix("Primary Overlap Matrix", primary_->nbf(), primary_->nbf()));
    o2->compute(S);

    return S;
}
shared_ptr<Matrix> Pseudospectral::form_S_dealias()
{
    shared_ptr<IntegralFactory> fact(new IntegralFactory(dealias_, primary_, primary_, primary_));
    shared_ptr<OneBodyAOInt> o2(fact->ao_overlap());
    shared_ptr<Matrix> S(new Matrix("Dealias Overlap Matrix", dealias_->nbf(), primary_->nbf()));
    o2->compute(S);

    return S;
}
shared_ptr<Matrix> Pseudospectral::form_Q()
{
    int nbf = primary_->nbf();
    int ndf = dealias_->nbf();
    int ntotal = nbf + ndf;

    shared_ptr<Matrix> Q(new Matrix("Q (Least-Squares)", nbf, npoints_));
    double** Qp = Q->pointer();     
    
    shared_ptr<Matrix> X = form_X();
    shared_ptr<Matrix> Xd = form_X_dealias();
    shared_ptr<Matrix> S = form_S();
    shared_ptr<Matrix> Sd = form_S_dealias();

    X->print();
    Xd->print();
    S->print();
    Sd->print();

    double** Xp = X->pointer();
    double** Xdp = Xd->pointer();
    double** Sp = S->pointer();
    double** Sdp = Sd->pointer();

    // Copy S for later use (final projection)
    shared_ptr<Matrix> S2(new Matrix("S (Copy)", nbf, nbf));
    double** S2p = S->pointer();
    memcpy(static_cast<void*> (S2p[0]), static_cast<void*> (Sp[0]), nbf*nbf*sizeof(double)); 

    // Orthogonalize the primary and dealias bases
    C_DPOTRF('L', nbf, Sp[0], nbf);
    S->print();
    C_DPOTRS('L', nbf, ndf, Sp[0], nbf, Sdp[0], nbf);
    Sd->print();    

    // Build the Collocation matrix
    shared_ptr<Matrix> R(new Matrix("R (Full Collocation Matrix)", ntotal, npoints_));
    double** Rp = R->pointer();
    double** Rdp = &Rp[nbf];
    memcpy(static_cast<void*> (Rp[0]), static_cast<void*> (Xp[0]), nbf*npoints_*sizeof(double)); 
    memcpy(static_cast<void*> (Rp[nbf]), static_cast<void*> (Xdp[0]), ndf*npoints_*sizeof(double)); 
    R->print();

    // Remove the overlap from the primary basis on the dealias collocation partition  
    C_DGEMM('N','N', ndf, npoints_, nbf, -1.0, Sdp[0], nbf, Xp[0], npoints_, 1.0, Rdp[0], npoints_); 
    R->print();

    // Roll the square root of the weights into R (one square root w will chill until the end)
    double* w = new double[npoints_];
    double* wp = grid_->getBlock()->getWeights();
    for (int P = 0; P < npoints_; P++)
        w[P] = sqrt(wp[P]);

    for (int P = 0; P < npoints_; P++) 
        C_DSCAL(ntotal, w[P], &Rp[0][P], npoints_);
    R->print();

    // Form (X'X)X'\sqrt(w) via QR Decomposition
    double* tau = new double[ntotal];

    // First, find out how much workspace to provide
    double work_size;
    C_DGEQRF(npoints_,ntotal,Rp[0],npoints_,tau,&work_size, -1);  

    // Now, do the QR decomposition
    int lwork = (int)work_size;
    double *work = new double[lwork];
    C_DGEQRF(npoints_,ntotal,Rp[0],npoints_,tau,work, lwork);  
    R->set_name("Q (of QR Decomposition");
    delete[] work;

    R->print();

    // Put R in the upper triangle where it belongs
    shared_ptr<Matrix> r(new Matrix("R (of QR decomposition)", ntotal, ntotal));
    double** rp = r->pointer();
    for (int i = 0; i < ntotal; i++)
        for (int j = i; j < ntotal; j++) {
            rp[i][j] = Rp[j][i]; 
        }
 
    r->print();
    
    // First, find out how much workspace to provide
    C_DORGQR(npoints_,ntotal,ntotal,Rp[0],npoints_,tau,&work_size,-1); 

    // Now, form Q
    lwork = (int)work_size;
    work = new double[lwork];
    C_DORGQR(npoints_,ntotal,ntotal,Rp[0],npoints_,tau,work,lwork); 
    delete[] work;
    delete[] tau;

    R->print();
    
    // Scale Q' by sqrt w
    for (int P = 0; P < npoints_; P++) 
        C_DSCAL(ntotal, w[P], &Rp[0][P], npoints_);
    delete[] w;

    R->print();
    
    // Backsolve R^-1 Q' sqrt w
    C_DTRSM('L','U','N','N', ntotal, npoints_, 1.0, rp[0], ntotal, Rp[0], npoints_);

    R->print();

    // Do P R^-1 Q' sqrt w, using S (this bit is confusing, Friesner never gets it straight)
    C_DGEMM('N','N', nbf, npoints_, nbf, 1.0, S2p[0], nbf, Rp[0], npoints_, 0.0, Qp[0], npoints_);

    Q->print();

    return Q;
}
shared_ptr<Matrix> Pseudospectral::form_A()
{
    shared_ptr<IntegralFactory> fact(new IntegralFactory(primary_, primary_, primary_, primary_));
    shared_ptr<Matrix> A(new Matrix("A Integrals", npoints_, primary_->nbf()*primary_->nbf()));
    double** Ap = A->pointer();
    shared_ptr<Matrix> T(new Matrix("Temp", primary_->nbf(), primary_->nbf()));
    double** Tp = T->pointer();
    
    shared_ptr<PseudospectralInt> ints(static_cast<PseudospectralInt*>(fact->ao_pseudospectral()));
    const double* buffer = ints->buffer();

    shared_ptr<GridBlock> block = grid_->getBlock();
    double* x = block->getX();
    double* y = block->getY();
    double* z = block->getZ();

    for (int P = 0; P < npoints_; P++) {
        ints->set_point(x[P], y[P], z[P]);
        T->zero();
        ints->compute(T);
    
        fprintf(outfile, " Amn integrals for P = %d\n", P);
        T->print(); 
    
        memcpy(static_cast<void*>(Ap[P]), static_cast<void*>(Tp[0]), primary_->nbf()*primary_->nbf()*sizeof(double));
    }
        
    return A;
}
void Pseudospectral::form_A_disk()
{
    shared_ptr<IntegralFactory> fact(new IntegralFactory(primary_, primary_, primary_, primary_));
    shared_ptr<Matrix> T(new Matrix("Temp", primary_->nbf(), primary_->nbf()));
    double** Tp = T->pointer();
    
    shared_ptr<PseudospectralInt> ints(static_cast<PseudospectralInt*>(fact->ao_pseudospectral()));
    const double* buffer = ints->buffer();

    shared_ptr<GridBlock> block = grid_->getBlock();
    double* x = block->getX();
    double* y = block->getY();
    double* z = block->getZ();

    psio_address address = PSIO_ZERO;
    for (int P = 0; P < npoints_; P++) {
        ints->set_point(x[P], y[P], z[P]);
        ints->compute(T);
        psio_->write(PSIF_PS_TENSOR, "AO PS Integrals", (char*) Tp[0],primary_->nbf()*primary_->nbf()*sizeof(double), address, &address); 
    }
}
shared_ptr<Matrix> Pseudospectral::form_I()
{
    shared_ptr<Matrix> X = form_X();
    shared_ptr<Matrix> Q = form_Q();
    shared_ptr<Matrix> A = form_A();

    X->print();
    Q->print();
    A->print();

    int nbf = primary_->nbf();

    shared_ptr<Matrix> I(new Matrix("PS ERIs", nbf*nbf, nbf*nbf));

    double** Xp = X->pointer();
    double** Ap = A->pointer();
    double** Qp = Q->pointer();
    double** Ip = I->pointer();

    shared_ptr<Matrix> temp(new Matrix("Temp", npoints_, nbf*nbf));
    double** Tp = temp->pointer();

    for (int P = 0; P < npoints_; P++) {
        for (int m = 0; m < nbf; m++) {
            for (int n = 0; n < nbf; n++) {
                Tp[P][m*nbf + n] = Qp[m][P] * Xp[n][P];
            }
        }
    } 

    C_DGEMM('T', 'N', nbf*nbf, nbf*nbf, npoints_, 1.0, Tp[0], nbf*nbf, Ap[0], nbf*nbf, 0.0, Ip[0], nbf*nbf);

    return I;
}



}
