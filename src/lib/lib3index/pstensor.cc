#include "3index.h"
#include <libmints/mints.h>
#include <libmints/cubature.h>
#include <libqt/qt.h>

#include <string>
#include <sstream>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <algorithm>
#include <utility>
#include <ctype.h>
#include <cmath>

using namespace boost;
using namespace std;
using namespace psi;

namespace psi {

PSTensorII::PSTensorII(boost::shared_ptr<BasisSet> primary,
                   SharedMatrix C,
                   int nocc,
                   int nvir,
                   int naocc,
                   int navir,
                   ULI memory,
                   Options& options) :
    primary_(primary), C_(C), nocc_(nocc), nvir_(nvir),
    naocc_(naocc), navir_(navir), memory_(memory), options_(options)
{
    common_init();
}
PSTensorII::~PSTensorII()
{
}
void PSTensorII::common_init()
{
    print_ = options_.get_int("PRINT");
    debug_ = options_.get_int("DEBUG");

    do_omega_ = options_.get_bool("PS_USE_OMEGA");
    omega_ = options_.get_double("PS_OMEGA");

    min_S_dealias_ = options_.get_double("PS_MIN_S_DEALIAS");

    alpha_ = options_.get_double("PS_ALPHA");

    if (options_.get_str("PS_FITTING_ALGORITHM") == "QUADRATURE") {
        do_renormalize_ = false;
        do_dealias_ = false;
    } else if (options_.get_str("PS_FITTING_ALGORITHM") == "RENORMALIZED") {
        do_renormalize_ = true;
        do_dealias_ = false;
    } else if (options_.get_str("PS_FITTING_ALGORITHM") == "DEALIASED") {
        do_renormalize_ = true;
        do_dealias_ = true;
    }

    molecule_ = primary_->molecule();

    nfocc_ = nocc_ - naocc_;
    nfvir_ = nvir_ - navir_;

    nso_ = C_->rowspi()[0];
    nmo_ = C_->colspi()[0];
    ndso_ = ndmo_ = 0;

    Caocc_ = SharedMatrix(new Matrix("C active occupied", nso_, naocc_));
    Cavir_ = SharedMatrix(new Matrix("C active virtual", nso_, navir_));

    double** Cp = C_->pointer();
    double** Cop  = Caocc_->pointer();
    double** Cvp  = Cavir_->pointer();

    for (int m = 0; m < nso_; m++) {
        C_DCOPY(naocc_, &Cp[m][nfocc_],1, Cop[m], 1);
        C_DCOPY(navir_, &Cp[m][nocc_],1, Cvp[m], 1);
    }

    buildGrid();

    print_header();

    if (options_.get_str("PS_FITTING_ALGORITHM") == "DEALIASED") {
        buildDealiased();
    } else  if (options_.get_str("PS_FITTING_ALGORITHM") == "RENORMALIZED") {
        buildRenormalized();
    } else  if (options_.get_str("PS_FITTING_ALGORITHM") == "QUADRATURE") {
        buildQuadrature();
    }
}
void PSTensorII::buildQuadrature()
{
    form_Rpao();
    form_Rpmo();

    form_Q_quadrature();
}
void PSTensorII::buildRenormalized()
{
    form_Rpao();
    form_Rpmo();

    form_Q_renormalized();
}
void PSTensorII::buildDealiased()
{
    buildDealiasSet();

    form_Spdao();
    form_Spdmo();
    form_Sddao();
    form_Sddoo();
    form_Cdd();

    form_Rpao();
    form_Rpmo();
    form_Rdao();
    form_Rdmo();

    form_Q_dealiased();
}
void PSTensorII::print_header()
{
    fprintf(outfile,"  ==> PS Tensor II (by Rob Parrish) <==\n\n");

    fprintf(outfile, "    %s fitting algorithm will be applied.\n", options_.get_str("PS_FITTING_ALGORITHM").c_str());

    if (do_omega_) {
        fprintf(outfile, "    Range separation will be performed with \\omega of %11.4E.\n\n", omega_);
    } else {
        fprintf(outfile, "    No range separation requested.\n\n");
    }

    if (print_) {
        molecule_->print();
        primary_->print_by_level(outfile, print_);
    }

    if (debug_ > 1) {
        C_->print();
        Caocc_->print();
        Cavir_->print();
    }

    if (print_) {
        grid_->print(outfile, print_);
        fflush(outfile);
    }
}
void PSTensorII::buildGrid()
{
    if (options_.get_str("PS_GRID_FILE") == "") {
        grid_ = boost::shared_ptr<PseudospectralGrid>(new PseudospectralGrid(molecule_,
            primary_, options_));
    } else {
        grid_ = boost::shared_ptr<PseudospectralGrid>(new PseudospectralGrid(molecule_,
            primary_, options_.get_str("PS_GRID_FILE"), options_));
    }

    naux_ = grid_->npoints();
    w_ = boost::shared_ptr<Vector>(new Vector("Grid Weights", naux_));
    double* wp = w_->pointer();

    C_DCOPY(naux_, grid_->w(), 1, wp, 1);
}
void PSTensorII::buildDealiasSet()
{
    if (print_)
        fprintf(outfile," => Dealias Basis Set <= \n\n");

    if (options_.get_str("DEALIAS_BASIS_CC") == "") {
        if (print_)
            fprintf(outfile,"  Dealias Basis Automatically Generated.\n\n");

        boost::shared_ptr<DealiasBasisSet> d(new DealiasBasisSet(primary_, options_));
        dealias_ = d->dealiasSet();
    } else {
        if (print_)
            fprintf(outfile,"  Dealias Basis Read from %s.\n\n", options_.get_str("DEALIAS_BASIS_CC").c_str());

        boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
        molecule_->set_basis_all_atoms(options_.get_str("DEALIAS_BASIS_CC"),"DEALIAS_BASIS");
        dealias_ = BasisSet::construct(parser,molecule_,"DEALIAS_BASIS");
    }

    if (print_) {
        dealias_->print_by_level(outfile,print_);
        fflush(outfile);
    }

    ndso_ = ndmo_ = dealias_->nbf();
}
void PSTensorII::form_Q_quadrature()
{
    Rmo_ = Rpmo_;
    Rmo_->set_name("Rmo (primary x points)");

    Qmo_ = SharedMatrix(new Matrix("Qmo (primary x points)", nmo_, naux_));

    double** Qp = Qmo_->pointer();
    double** Rp = Rmo_->pointer();
    double* wp = w_->pointer();

    for (int Q = 0; Q < naux_ ; Q++) {
        C_DAXPY(nmo_, wp[Q], &Rp[0][Q], naux_, &Qp[0][Q], naux_);
    }

    if (debug_ > 1) {
        Qmo_->print();
        Rmo_->print();
    }

    if (debug_ > 2) {
        SharedMatrix S(new Matrix("Sbar", nmo_, nmo_));
        double** Sp = S->pointer();
        double** Qp = Qmo_->pointer();
        double** Rp = Rmo_->pointer();
        C_DGEMM('N','T',nmo_,nmo_,naux_,1.0,Qp[0],naux_,Rp[0],naux_,0.0,Sp[0],nmo_);
        S->print();
    }
}
void PSTensorII::form_Q_renormalized()
{
    SharedMatrix Sbar(new Matrix("Sbar", nmo_, nmo_));
    SharedMatrix Rw(new Matrix("Rw", nmo_, naux_));

    double** Rwp = Rw->pointer();
    double** Rp = Rpmo_->pointer();
    double* wp = w_->pointer();

    //  Rw
    for (int Q = 0; Q < naux_ ; Q++) {
        C_DAXPY(nmo_, wp[Q], &Rp[0][Q], naux_, &Rwp[0][Q], naux_);
    }

    double** Sp = Sbar->pointer();

    //  Sbar
    C_DGEMM('N','T',nmo_,nmo_,naux_,1.0,Rwp[0],naux_,Rp[0],naux_,0.0,Sp[0],nmo_);

    // Sbar ^ -(alpha) Rw, Sbar ^ -(1-alpha) R
    if (alpha_ == 1.0) {
        int info = C_DPOTRF('L',nmo_,Sp[0],nmo_);
        if (info) {
            throw PSIEXCEPTION("Sbar Cholesky Factorization not SPD");
        }
        info = C_DPOTRI('L',nmo_,Sp[0],nmo_);
        if (info) {
            throw PSIEXCEPTION("Sbar Cholesky Inverse not SPD");
        }
        Sbar->copy_upper_to_lower();

        Qmo_ = SharedMatrix(new Matrix("Qmo (primary x points)", nmo_, naux_));
        double** Qp = Qmo_->pointer();

        C_DGEMM('N','N',nmo_,naux_,nmo_,1.0,Sp[0],nmo_,Rwp[0],naux_,0.0,Qp[0],naux_);

        Rmo_ = Rpmo_;
        Rmo_->set_name("Rmo (primary x points)");
    } else {
        SharedMatrix Sbarl(new Matrix("Sbarl", nmo_, nmo_));
        Sbarl->copy(Sbar);
        double** Slp = Sbarl->pointer();

        Sbar->power(-alpha_);
        Sbarl->power(-(1.0-alpha_));

        Qmo_ = SharedMatrix(new Matrix("Qmo (primary x points)", nmo_, naux_));
        Rmo_ = SharedMatrix(new Matrix("Rmo (primary x points)", nmo_, naux_));
        double** Rmop = Rmo_->pointer();
        double** Qmop = Qmo_->pointer();

        C_DGEMM('N','N',nmo_,naux_,nmo_,1.0,Sp[0],nmo_,Rwp[0],naux_,0.0,Qmop[0],naux_);
        C_DGEMM('N','N',nmo_,naux_,nmo_,1.0,Slp[0],nmo_,Rp[0],naux_,0.0,Rmop[0],naux_);
    }

    Rpmo_.reset();
}
void PSTensorII::form_Q_dealiased()
{
    SharedMatrix Sbar(new Matrix("Sbar", nmo_ + ndmo_, nmo_ + ndmo_));
    SharedMatrix Rw(new Matrix("Rw", nmo_ + ndmo_, naux_));

    double** Rwp = Rw->pointer();
    double** Rp = Rpmo_->pointer();
    double** Rdp = Rdmo_->pointer();
    double* wp = w_->pointer();

    //  Rw
    for (int Q = 0; Q < naux_ ; Q++) {
        C_DAXPY(nmo_, wp[Q], &Rp[0][Q], naux_, &Rwp[0][Q], naux_);
        C_DAXPY(ndmo_, wp[Q], &Rdp[0][Q], naux_, &Rwp[nmo_][Q], naux_);
    }

    double** Sp = Sbar->pointer();

    //  Sbar
    C_DGEMM('N','T',nmo_ + ndmo_,nmo_,naux_,1.0,Rwp[0],naux_,Rp[0],naux_,0.0,Sp[0],ndmo_ + nmo_);
    C_DGEMM('N','T',nmo_ + ndmo_,ndmo_,naux_,1.0,Rwp[0],naux_,Rdp[0],naux_,0.0,&Sp[0][nmo_],ndmo_ + nmo_);

    // Sbar ^ -(alpha) Rw, Sbar ^ -(1-alpha) R
    if (alpha_ == 1.0) {
        int info = C_DPOTRF('L',ndmo_ + nmo_,Sp[0],ndmo_ + nmo_);
        if (info) {
            throw PSIEXCEPTION("Sbar Cholesky Factorization not SPD");
        }
        info = C_DPOTRI('L',ndmo_ + nmo_,Sp[0],ndmo_ + nmo_);
        if (info) {
            throw PSIEXCEPTION("Sbar Cholesky Inverse not SPD");
        }
        Sbar->copy_upper_to_lower();

        Qmo_ = SharedMatrix(new Matrix("Qmo (primary x points)", nmo_, naux_));
        double** Qp = Qmo_->pointer();

        C_DGEMM('N','N',nmo_,naux_,nmo_ + ndmo_,1.0,Sp[0],nmo_ + ndmo_,Rwp[0],naux_,0.0,Qp[0],naux_);

        Rmo_ = Rpmo_;
        Rmo_->set_name("Rmo (primary x points)");
    } else {
        SharedMatrix Sbarl(new Matrix("Sbarl", ndmo_ + nmo_, ndmo_ + nmo_));
        Sbarl->copy(Sbar);
        double** Slp = Sbarl->pointer();

        Sbar->power(-alpha_);
        Sbarl->power(-(1.0-alpha_));

        Qmo_ = SharedMatrix(new Matrix("Qmo (primary x points)", nmo_, naux_));
        Rmo_ = SharedMatrix(new Matrix("Rmo (primary x points)", nmo_, naux_));
        double** Rmop = Rmo_->pointer();
        double** Qmop = Qmo_->pointer();

        C_DGEMM('N','N',nmo_,naux_,ndmo_ + nmo_,1.0,Sp[0],ndmo_ + nmo_,Rwp[0],naux_,0.0,Qmop[0],naux_);
        C_DGEMM('N','N',nmo_,naux_,nmo_,1.0,Slp[0],ndmo_ + nmo_,Rp[0],naux_,0.0,Rmop[0],naux_);
        C_DGEMM('N','N',nmo_,naux_,ndmo_,1.0,&Slp[0][nmo_],ndmo_ + nmo_,Rdp[0],naux_,1.0,Rmop[0],naux_);
    }

    Rpmo_.reset();
    Rdmo_.reset();
}
SharedMatrix PSTensorII::O()
{
    if (!do_omega_) {
        throw PSIEXCEPTION("O Double-Pseudospectral Operator Requested, but Range-Separation is not applied.");
    }

    if (O_.get())
        return O_;

    O_ = SharedMatrix(new Matrix("O (naux x naux)", naux_, naux_));
    double** Op = O_->pointer();
    double* xp = grid_->x();
    double* yp = grid_->y();
    double* zp = grid_->z();

    #ifndef HAVE_FUNC_ERF
        throw PSIEXCEPTION("erf not implemented. Get a C99 compiler.");
    #else
    for (int Q = 0; Q < naux_; Q++) {
        for (int P = Q; P < naux_; P++) {
            if (P == Q) {
                // M_2_SQRTPI is defined in math.h
                Op[P][P] = M_2_SQRTPI * omega_;
            } else {
                double R = 1.0 / sqrt((xp[P] - xp[Q]) * (xp[P] - xp[Q]) +
                                      (yp[P] - yp[Q]) * (yp[P] - yp[Q]) +
                                      (zp[P] - zp[Q]) * (zp[P] - zp[Q]));
                // erf is defined in math.h if C99 or better
                Op[P][Q] = Op[Q][P] = erf(omega_ * R) / R;
            }
        }
    }
    #endif

    return O_;
}
SharedMatrix PSTensorII::Q()
{
    return Qmo_;
}
SharedMatrix PSTensorII::Qocc()
{
    SharedMatrix Q(new Matrix("Qocc", nocc_, naux_));
    double** Qr = Qmo_->pointer();
    double** Qp = Q->pointer();

    C_DCOPY(nocc_ * (ULI) naux_, Qr[0], 1, Qp[0], 1);

    return Q;
}
SharedMatrix PSTensorII::Qaocc()
{
    SharedMatrix Q(new Matrix("Qaocc", naocc_, naux_));
    double** Qr = Qmo_->pointer();
    double** Qp = Q->pointer();

    C_DCOPY(naocc_ * (ULI) naux_, Qr[nfocc_], 1, Qp[0], 1);

    return Q;
}
SharedMatrix PSTensorII::Qvir()
{
    SharedMatrix Q(new Matrix("Qvir", nvir_, naux_));
    double** Qr = Qmo_->pointer();
    double** Qp = Q->pointer();

    C_DCOPY(nvir_ * (ULI) naux_, Qr[nocc_], 1, Qp[0], 1);

    return Q;
}
SharedMatrix PSTensorII::Qavir()
{
    SharedMatrix Q(new Matrix("Qavir", navir_, naux_));
    double** Qr = Qmo_->pointer();
    double** Qp = Q->pointer();

    C_DCOPY(navir_ * (ULI) naux_, Qr[nocc_], 1, Qp[0], 1);

    return Q;
}
SharedMatrix PSTensorII::R()
{
    return Rmo_;
}
SharedMatrix PSTensorII::Rocc()
{
    SharedMatrix R(new Matrix("Rocc", nocc_, naux_));
    double** Rr = Rmo_->pointer();
    double** Rp = R->pointer();

    C_DCOPY(nocc_ * (ULI) naux_, Rr[0], 1, Rp[0], 1);

    return R;
}
SharedMatrix PSTensorII::Raocc()
{
    SharedMatrix R(new Matrix("Raocc", naocc_, naux_));
    double** Rr = Rmo_->pointer();
    double** Rp = R->pointer();

    C_DCOPY(naocc_ * (ULI) naux_, Rr[nfocc_], 1, Rp[0], 1);

    return R;
}
SharedMatrix PSTensorII::Rvir()
{
    SharedMatrix R(new Matrix("Rvir", nvir_, naux_));
    double** Rr = Rmo_->pointer();
    double** Rp = R->pointer();

    C_DCOPY(nvir_ * (ULI) naux_, Rr[nocc_], 1, Rp[0], 1);

    return R;
}
SharedMatrix PSTensorII::Ravir()
{
    SharedMatrix R(new Matrix("Ravir", navir_, naux_));
    double** Rr = Rmo_->pointer();
    double** Rp = R->pointer();

    C_DCOPY(navir_ * (ULI) naux_, Rr[nocc_], 1, Rp[0], 1);

    return R;
}
SharedMatrix PSTensorII::Aso()
{
    SharedMatrix A(new Matrix("Aso",  naux_, nso_ * nso_));
    double** Ap = A->pointer();

    boost::shared_ptr<IntegralFactory> fact(new IntegralFactory(primary_,primary_,primary_,primary_));
    boost::shared_ptr<PseudospectralInt> ints(static_cast<PseudospectralInt*>(fact->ao_pseudospectral()));

    if (do_omega_) {
        ints->set_omega(omega_);
    }

    double* x = grid_->x();
    double* y = grid_->y();
    double* z = grid_->z();

    SharedMatrix T(new Matrix("Temp", primary_->nbf(), primary_->nbf()));
    double** Tp = T->pointer();

    for (int P = 0; P < naux_; P++) {
        ints->set_point(x[P], y[P], z[P]);
        T->zero();
        ints->compute(T);

        C_DCOPY(nso_ * nso_, Tp[0], 1, Ap[P], 1);
    }

    return A;
}
SharedMatrix PSTensorII::Aoo()
{
    SharedMatrix Amn = Aso();
    SharedMatrix Ami(new Matrix("Ami", naux_, naocc_ * (ULI) nso_));

    double** Amnp = Amn->pointer();
    double** Amip = Ami->pointer();
    double** Cop = Caocc_->pointer();

    C_DGEMM('N','N', naux_ * (ULI) nso_, naocc_, nso_, 1.0, Amnp[0], nso_, Cop[0], naocc_,
        0.0, Amip[0], naocc_);

    Amn.reset();

    SharedMatrix Aia(new Matrix("Aij", naux_, naocc_ * (ULI) naocc_));
    double** Aiap = Aia->pointer();

    for (int Q = 0; Q < naux_; Q++) {
        C_DGEMM('T','N',naocc_,naocc_,nso_,1.0,Amip[Q],naocc_,Cop[0],naocc_, 0.0, Aiap[0], naocc_);
    }

    return Aia;
}
SharedMatrix PSTensorII::Aov()
{
    SharedMatrix Amn = Aso();
    SharedMatrix Ami(new Matrix("Ami", naux_, naocc_ * (ULI) nso_));

    double** Amnp = Amn->pointer();
    double** Amip = Ami->pointer();
    double** Cop = Caocc_->pointer();
    double** Cvp = Cavir_->pointer();

    C_DGEMM('N','N', naux_ * (ULI) nso_, naocc_, nso_, 1.0, Amnp[0], nso_, Cop[0], naocc_,
        0.0, Amip[0], naocc_);

    Amn.reset();

    SharedMatrix Aia(new Matrix("Aia", naux_, naocc_ * (ULI) navir_));
    double** Aiap = Aia->pointer();

    for (int Q = 0; Q < naux_; Q++) {
        C_DGEMM('T','N',naocc_,navir_,nso_,1.0,Amip[Q],naocc_,Cvp[0],navir_, 0.0, Aiap[Q], navir_);
    }

    return Aia;
}
SharedMatrix PSTensorII::Avv()
{
    SharedMatrix Amn = Aso();
    SharedMatrix Ami(new Matrix("Ami", naux_, navir_ * (ULI) nso_));

    double** Amnp = Amn->pointer();
    double** Amip = Ami->pointer();
    double** Cvp = Cavir_->pointer();

    C_DGEMM('N','N', naux_ * (ULI) nso_, navir_, nso_, 1.0, Amnp[0], nso_, Cvp[0], navir_,
        0.0, Amip[0], navir_);

    Amn.reset();

    SharedMatrix Aia(new Matrix("Aab", naux_, navir_ * (ULI) navir_));
    double** Aiap = Aia->pointer();

    for (int Q = 0; Q < naux_; Q++) {
        C_DGEMM('T','N',navir_,navir_,nso_,1.0,Amip[Q],navir_,Cvp[0],navir_, 0.0, Aiap[Q], navir_);
    }

    return Aia;
}
SharedMatrix PSTensorII::Amo()
{
    SharedMatrix Amn = Aso();
    SharedMatrix Ami(new Matrix("Ami", naux_, nmo_ * (ULI) nso_));

    double** Amnp = Amn->pointer();
    double** Amip = Ami->pointer();
    double** Cvp = C_->pointer();

    C_DGEMM('N','N', naux_ * (ULI) nso_, nmo_, nso_, 1.0, Amnp[0], nso_, Cvp[0], nmo_,
        0.0, Amip[0], nmo_);

    Amn.reset();

    SharedMatrix Aia(new Matrix("Amo", naux_, nmo_ * (ULI) nmo_));
    double** Aiap = Aia->pointer();

    for (int Q = 0; Q < naux_; Q++) {
        C_DGEMM('T','N',nmo_,nmo_,nso_,1.0,Amip[Q],nmo_,Cvp[0],nmo_, 0.0, Aiap[Q], nmo_);
    }

    return Aia;
}
SharedMatrix PSTensorII::Imo()
{
    boost::shared_ptr<MintsHelper> mints(new MintsHelper());
    return mints->mo_eri(C_,C_);
}
SharedMatrix PSTensorII::Ipsmo()
{
    SharedMatrix Am = Amo();

    double** Qmop = Qmo_->pointer();
    double** Rmop = Rmo_->pointer();
    double** Amop = Am->pointer();

    SharedMatrix QR(new Matrix("QR", nmo_ * nmo_, naux_));
    double** QRp = QR->pointer();

    for (int a = 0; a < nmo_; a++) {
    for (int b = 0; b < nmo_; b++) {
    for (int Q = 0; Q < naux_; Q++) {
        QRp[a * nmo_ + b][Q] = Qmop[a][Q] * Rmop[b][Q];
    }}}

    SharedMatrix Imo(new Matrix("PS MO ERI Tensor", nmo_ * nmo_, nmo_ * nmo_));
    double** Imop = Imo->pointer();

    C_DGEMM('N','N',nmo_ * nmo_, nmo_ * nmo_, naux_, 1.0, QRp[0], naux_, Amop[0], nmo_ * nmo_, 0.0, Imop[0], nmo_ * nmo_);

    return Imo;
}
SharedMatrix PSTensorII::Idpsmo()
{
    O();

    double** Qmop = Qmo_->pointer();
    double** Rmop = Rmo_->pointer();
    double** Op = O_->pointer();

    SharedMatrix QR(new Matrix("QR", nmo_ * nmo_, naux_));
    double** QRp = QR->pointer();

    for (int a = 0; a < nmo_; a++) {
    for (int b = 0; b < nmo_; b++) {
    for (int Q = 0; Q < naux_; Q++) {
        QRp[a * nmo_ + b][Q] = Qmop[a][Q] * Rmop[b][Q];
    }}}

    SharedMatrix QRO(new Matrix("QRO", nmo_ * nmo_, naux_));
    double** QROp = QRO->pointer();

    C_DGEMM('N','N',nmo_ * nmo_, naux_, naux_, 1.0, QRp[0], naux_, Op[0], naux_, 0.0, QROp[0], naux_);

    O_.reset();

    SharedMatrix Imo(new Matrix("Double PS MO ERI Tensor", nmo_ * nmo_, nmo_ * nmo_));
    double** Imop = Imo->pointer();

    C_DGEMM('N','T',nmo_ * nmo_, nmo_ * nmo_, naux_, 1.0, QROp[0], naux_, QRp[0], naux_, 0.0, Imop[0], nmo_ * nmo_);

    return Imo;
}
void PSTensorII::form_Rpao()
{
    Rpao_ = SharedMatrix(new Matrix("R (primary x points)", nso_, naux_));
    double** Rp = Rpao_->pointer();

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

    if (debug_ > 1)
        Rpao_->print();
}
void PSTensorII::form_Rdao()
{
    Rdao_ = SharedMatrix(new Matrix("R (dealias x points)", ndso_, naux_));
    double** Rp = Rdao_->pointer();

    #if 0

    boost::shared_ptr<BasisPoints> points(new BasisPoints(dealias_, naux_));
    points->setToComputePoints(true);
    double** bpoints = points->getPoints();

    // Compute the basis points
    points->computePoints(grid_->fullGrid());

    // Copy the points in
    for (int i = 0; i < naux_; i++) {
        for (int Q = 0; Q < ndso_; Q++)
            Rp[Q][i] = bpoints[i][Q];
    }

    #endif

    if (debug_ > 1)
        Rdao_->print();
}
void PSTensorII::form_Rpmo()
{
    Rpmo_ = SharedMatrix(new Matrix("R2 (primary' x points)", nmo_, naux_));
    double** Rp2 = Rpmo_->pointer();
    double** Rp = Rpao_->pointer();
    double** Xp = C_->pointer();

    C_DGEMM('T','N',nmo_,naux_,nso_,1.0,Xp[0],nmo_,Rp[0],naux_,0.0,Rp2[0],naux_);

    if (debug_ > 1)
        Rpmo_->print();

    Rpao_.reset();
}
void PSTensorII::form_Rdmo()
{
    SharedMatrix Rdso(new Matrix("R2 (dealias x points)", ndso_, naux_));
    Rdso->copy(Rdao_);
    double** Rd4p = Rdso->pointer();
    double** Rp2p = Rpmo_->pointer();
    double** Cdp = Spdmo_->pointer();
    double** Cdd = Cdd_->pointer();

    C_DGEMM('T','N',ndso_,naux_,nmo_,-1.0,Cdp[0],ndso_,Rp2p[0],naux_,1.0,Rd4p[0],naux_);

    Rdao_.reset();
    Spdmo_.reset();

    Rdmo_ = SharedMatrix(new Matrix("R2 (dealias' x points)", ndmo_, naux_));
    double** Rd3p = Rdmo_->pointer();

    C_DGEMM('T','N',ndmo_,naux_,ndso_,1.0,Cdd[0],ndmo_,Rd4p[0],naux_,0.0,Rd3p[0],naux_);

    Cdd_.reset();

    if (debug_ > 1)
        Rdmo_->print();
}
void PSTensorII::form_Spdao()
{
    boost::shared_ptr<IntegralFactory> fact(new IntegralFactory(primary_,dealias_,primary_,primary_));
    boost::shared_ptr<OneBodyAOInt> Sint(fact->ao_overlap());

    Spdao_ = SharedMatrix(new Matrix("S (primary x dealias)", nso_, ndso_));
    Sint->compute(Spdao_);

    if (debug_ > 1)
        Spdao_->print();
}
void PSTensorII::form_Spdmo()
{
    Spdmo_ = SharedMatrix(new Matrix("S (primary' x dealias)", nmo_, ndso_));
    double** Spd3p = Spdmo_->pointer();
    double** Spdp = Spdao_->pointer();
    double** Xp = C_->pointer();

    C_DGEMM('T','N',nmo_,ndso_,nso_,1.0,Xp[0],nmo_,Spdp[0],ndso_,0.0,Spd3p[0],ndso_);

    if (debug_ > 1)
        Spdmo_->print();

    Spdao_.reset();
}
void PSTensorII::form_Sddao()
{
    boost::shared_ptr<IntegralFactory> fact(new IntegralFactory(dealias_,dealias_,primary_,primary_));
    boost::shared_ptr<OneBodyAOInt> Sint(fact->ao_overlap());

    Sddao_ = SharedMatrix(new Matrix("S (dealias x dealias)", ndso_, ndso_));
    Sint->compute(Sddao_);

    if (debug_ > 1)
        Sddao_->print();
}
void PSTensorII::form_Sddoo()
{
    Sddoo_ = SharedMatrix(new Matrix("S (dealias x dealias)", ndso_, ndso_));
    Sddoo_->copy(Sddao_);
    Sddao_.reset();

    double** Sp = Sddoo_->pointer();
    double** Cpdp = Spdmo_->pointer();

    C_DGEMM('T','N',ndso_,ndso_,nmo_,-1.0,Cpdp[0],ndso_,Cpdp[0],ndso_,1.0,Sp[0],ndso_);
}
void PSTensorII::form_Cdd()
{
    SharedMatrix V(new Matrix("Eigvecs", ndso_, ndso_));
    boost::shared_ptr<Vector> c(new Vector("Eigvals", ndso_));
    double** Vp = V->pointer();
    double*  cp = c->pointer();

    Sddoo_->diagonalize(V,c);
    Sddoo_.reset();

    if (debug_ > 2)
        V->eivprint(c);

    ndmo_ = 0;
    for (int i = 0; i < ndso_; i++) {
        if (cp[i] >= min_S_dealias_)
            ndmo_++;
    }

    if (print_) {
        fprintf(outfile, "  %d of %d dealias functions selected, %d projected out.\n\n", ndmo_, ndso_,
            ndso_ - ndmo_);
        fflush(outfile);
    }


    Cdd_ = SharedMatrix(new Matrix("Cdd", ndso_, ndmo_));
    double** Wp = Cdd_->pointer();

    int j = 0;
    for (int i = 0; i < ndso_; i++) {
        if (cp[i] >= min_S_dealias_) {
            C_DAXPY(ndso_, pow(cp[i], -1.0/2.0), &Vp[0][i], ndso_, &Wp[0][j], ndmo_);
            j++;
        }
    }

    if (debug_ > 1)
        Cdd_->print();
}
// Older PSTensor
PSTensor::PSTensor(boost::shared_ptr<BasisSet> primary,
                   SharedMatrix C,
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

    Caocc_ = SharedMatrix(new Matrix("C active occupied", nso_, naocc_));
    Cavir_ = SharedMatrix(new Matrix("C active virtual", nso_, navir_));

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
            primary_, options_));
    } else {
        grid_ = boost::shared_ptr<PseudospectralGrid>(new PseudospectralGrid(molecule_,
            primary_, options_.get_str("PS_GRID_FILE"), options_));
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

    Spdao_ = SharedMatrix(new Matrix("S (primary x dealias)", nso_, ndealias_));
    Sint->compute(Spdao_);

    if (debug_)
        Spdao_->print();
}
void PSTensor::form_Spdmo()
{
    Spdmo_ = SharedMatrix(new Matrix("S (primary' x dealias)", nmo_, ndealias_));
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
    Rpao_ = SharedMatrix(new Matrix("R (primary x points)", nso_, naux_));
    double** Rp = Rpao_->pointer();

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
        Rpao_->print();
}
void PSTensor::form_Rdao()
{
    Rdao_ = SharedMatrix(new Matrix("R (dealias x points)", ndealias_, naux_));
    double** Rp = Rdao_->pointer();

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
        Rdao_->print();
}
void PSTensor::form_Rpmo()
{
    Rpmo_ = SharedMatrix(new Matrix("R2 (primary' x points)", nmo_, naux_));
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
    Rdmo_ = SharedMatrix(new Matrix("R2 (dealias' x points)", ndealias_, naux_));
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
    Ra_ = SharedMatrix(new Matrix("R Augmented (primary' + dealias' x points)", naug_, naux_));
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

    SharedMatrix C(new Matrix("C" , nmo_ + ndealias_, nmo_ + ndealias_));
    SharedMatrix Rw(new Matrix("Rw", nmo_ + ndealias_, naux_));
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

    Qmo_ = SharedMatrix (new Matrix("Qmo", nmo_, naux_));
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

    SharedMatrix Rw(new Matrix("Rw", nmo_, naux_));
    double** Rp = Ra_->pointer();
    double** Rwp = Rw->pointer();
    double* wp = w_->pointer();

    C_DCOPY((nmo_) * naux_, Rp[0], 1, Rwp[0], 1);

    for (int Q = 0; Q < naux_; Q++) {
        C_DSCAL(nmo_, wp[Q], &Rwp[0][Q], naux_);
    }
    Qmo_ = SharedMatrix (new Matrix("Qmo", nmo_, naux_));
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

    Qmo_ = SharedMatrix (new Matrix("Qmo", nmo_, naux_));
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
    Cpp_ = SharedMatrix(new Matrix("Cpp" , nmo_, nmo_));
    SharedMatrix Rw(new Matrix("Rw", nmo_, naux_));
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
    SharedMatrix V(new Matrix("Eigvecs", nmo_, nmo_));
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

    U_ = SharedMatrix(new Matrix("U", nmo_, nmo2_));
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
    SharedMatrix Cpdao(new Matrix("Cpdao", nmo_, ndealias_));
    double** Cpdaop = Cpdao->pointer();

    SharedMatrix Rw(new Matrix("Rw", nmo_, naux_));
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

    Cpd_ = SharedMatrix(new Matrix("Cpd", nmo2_, ndealias_));
    double** Cpdp = Cpd_->pointer();
    double** Xp = U_->pointer();

    C_DGEMM('T','N',nmo2_,ndealias_,nmo_,1.0,Xp[0],nmo2_,Cpdaop[0],ndealias_,0.0,Cpdp[0],ndealias_);


    if (debug_)
        Cpd_->print();

    Cpdao.reset();
}
void PSTensor::form_V()
{
    V_ = SharedMatrix(new Matrix("V", nmo2_, ndealias_));
    V_->copy(Cpd_);
    V_->scale(-1.0);

    if (debug_)
        V_->print();
}
void PSTensor::form_Cdd()
{
    Cdd_ = SharedMatrix(new Matrix("Cdd",ndealias_,ndealias_));
    double** Cddp = Cdd_->pointer();

    SharedMatrix Rw(new Matrix("Rw", ndealias_, naux_));
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
    SharedMatrix V(new Matrix("Eigvecs", ndealias_, ndealias_));
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


    W_ = SharedMatrix(new Matrix("w", ndealias_, ndealias2_));
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
    X_ = SharedMatrix(new Matrix("X", nmo_ + ndealias_, nmo2_ + ndealias2_));
    double** Xp = X_->pointer();

    double** Up = U_->pointer();
    double** Vp = V_->pointer();
    double** Wp = W_->pointer();

    for (int m = 0; m < nmo_; m++) {
        C_DCOPY(nmo2_, Up[m], 1, Xp[m], 1);
    }

    SharedMatrix T(new Matrix("T",nmo_,ndealias_));
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
    SharedMatrix Rw(new Matrix("Rw", nmo_ + ndealias_, naux_));
    double** Rp = Ra_->pointer();
    double** Rwp = Rw->pointer();
    double* wp = w_->pointer();

    C_DCOPY((nmo_ + ndealias_) * naux_, Rp[0], 1, Rwp[0], 1);

    for (int Q = 0; Q < naux_; Q++) {
        C_DSCAL(nmo_ + ndealias_, wp[Q], &Rwp[0][Q], naux_);
    }

    SharedMatrix C(new Matrix("C", nmo_ + ndealias_, nmo_ + ndealias_));
    double** Cp = C->pointer();

    C_DGEMM('N','T',nmo_ + ndealias_, nmo_ + ndealias_, naux_, 1.0, Rp[0], naux_, Rwp[0], naux_,
        0.0, Cp[0], nmo_ + ndealias_);

    C->print();

    SharedMatrix T(new Matrix("T", nmo_ + ndealias_, nmo2_ + ndealias2_));
    double** Tp = T->pointer();

    double** Xp = X_->pointer();

    C_DGEMM('N','N',nmo_ + ndealias_, nmo2_ + ndealias2_, nmo_ + ndealias_, 1.0, Cp[0], nmo_ + ndealias_,
        Xp[0], nmo2_ + ndealias2_, 0.0, Tp[0], nmo2_ + ndealias2_);

    SharedMatrix Cinv(new Matrix("Cinv", nmo2_ + ndealias2_, nmo2_ + ndealias2_));
    double** Cinvp = Cinv->pointer();

    C_DGEMM('T','N',nmo2_ + ndealias2_, nmo2_ + ndealias2_, nmo_ + ndealias_, 1.0, Xp[0], nmo2_ + ndealias2_,
        Tp[0], nmo2_ + ndealias2_, 0.0, Cinvp[0], nmo2_ + ndealias2_);

    Cinv->print();
}
void PSTensor::form_Q()
{
    SharedMatrix XX(new Matrix("XX^T", nmo_, nmo_ + ndealias_));
    double** XXp = XX->pointer();
    double** Xp = X_->pointer();

    C_DGEMM('N','T',nmo_,nmo_ + ndealias_,nmo2_ + ndealias2_,1.0,Xp[0],nmo2_ + ndealias2_,
        Xp[0],nmo2_ + ndealias2_,0.0,XXp[0],nmo_ + ndealias_);

    if (debug_ > 1)
        XX->print();

    SharedMatrix Rw(new Matrix("Rw", nmo_ + ndealias_, naux_));
    double** Rp = Ra_->pointer();
    double** Rwp = Rw->pointer();
    double* wp = w_->pointer();

    C_DCOPY((nmo_ + ndealias_) * naux_, Rp[0], 1, Rwp[0], 1);

    for (int Q = 0; Q < naux_; Q++) {
        C_DSCAL(nmo_ + ndealias_, wp[Q], &Rwp[0][Q], naux_);
    }

    Qmo_ = SharedMatrix (new Matrix("Qmo", nmo_, naux_));
    double** Qmop = Qmo_->pointer();

    C_DGEMM('N','N',nmo_,naux_,nmo_ + ndealias_,1.0,XXp[0],nmo_ + ndealias_,Rwp[0],naux_,0.0,Qmop[0],naux_);

    if (debug_)
        Qmo_->print();
}
SharedMatrix PSTensor::Q()
{
    return Qmo_;
}
SharedMatrix PSTensor::Qocc()
{
    SharedMatrix Q(new Matrix("Qocc", nocc_, naux_));
    double** Qr = Qmo_->pointer();
    double** Qp = Q->pointer();

    C_DCOPY(nocc_ * (ULI) naux_, Qr[0], 1, Qp[0], 1);

    return Q;
}
SharedMatrix PSTensor::Qaocc()
{
    SharedMatrix Q(new Matrix("Qaocc", naocc_, naux_));
    double** Qr = Qmo_->pointer();
    double** Qp = Q->pointer();

    C_DCOPY(naocc_ * (ULI) naux_, Qr[nfocc_], 1, Qp[0], 1);

    return Q;
}
SharedMatrix PSTensor::Qvir()
{
    SharedMatrix Q(new Matrix("Qvir", nvir_, naux_));
    double** Qr = Qmo_->pointer();
    double** Qp = Q->pointer();

    C_DCOPY(nvir_ * (ULI) naux_, Qr[nocc_], 1, Qp[0], 1);

    return Q;
}
SharedMatrix PSTensor::Qavir()
{
    SharedMatrix Q(new Matrix("Qavir", navir_, naux_));
    double** Qr = Qmo_->pointer();
    double** Qp = Q->pointer();

    C_DCOPY(navir_ * (ULI) naux_, Qr[nocc_], 1, Qp[0], 1);

    return Q;
}
SharedMatrix PSTensor::R()
{
    return Rmo_;
}
SharedMatrix PSTensor::Rocc()
{
    SharedMatrix R(new Matrix("Rocc", nocc_, naux_));
    double** Rr = Rmo_->pointer();
    double** Rp = R->pointer();

    C_DCOPY(nocc_ * (ULI) naux_, Rr[0], 1, Rp[0], 1);

    return R;
}
SharedMatrix PSTensor::Raocc()
{
    SharedMatrix R(new Matrix("Raocc", naocc_, naux_));
    double** Rr = Rmo_->pointer();
    double** Rp = R->pointer();

    C_DCOPY(naocc_ * (ULI) naux_, Rr[nfocc_], 1, Rp[0], 1);

    return R;
}
SharedMatrix PSTensor::Rvir()
{
    SharedMatrix R(new Matrix("Rvir", nvir_, naux_));
    double** Rr = Rmo_->pointer();
    double** Rp = R->pointer();

    C_DCOPY(nvir_ * (ULI) naux_, Rr[nocc_], 1, Rp[0], 1);

    return R;
}
SharedMatrix PSTensor::Ravir()
{
    SharedMatrix R(new Matrix("Ravir", navir_, naux_));
    double** Rr = Rmo_->pointer();
    double** Rp = R->pointer();

    C_DCOPY(navir_ * (ULI) naux_, Rr[nocc_], 1, Rp[0], 1);

    return R;
}
SharedMatrix PSTensor::Aso()
{
    SharedMatrix A(new Matrix("Aso",  naux_, nso_ * nso_));
    double** Ap = A->pointer();

    boost::shared_ptr<IntegralFactory> fact(new IntegralFactory(primary_,primary_,primary_,primary_));
    boost::shared_ptr<PseudospectralInt> ints(static_cast<PseudospectralInt*>(fact->ao_pseudospectral()));

    if (use_omega_) {
        ints->set_omega(omega_);
    }

    double* x = grid_->x();
    double* y = grid_->y();
    double* z = grid_->z();

    SharedMatrix T(new Matrix("Temp", primary_->nbf(), primary_->nbf()));
    double** Tp = T->pointer();

    for (int P = 0; P < naux_; P++) {
        ints->set_point(x[P], y[P], z[P]);
        T->zero();
        ints->compute(T);

        C_DCOPY(nso_ * nso_, Tp[0], 1, Ap[P], 1);
    }

    return A;
}
SharedMatrix PSTensor::Aoo()
{
    SharedMatrix Amn = Aso();
    SharedMatrix Ami(new Matrix("Ami", naux_, naocc_ * (ULI) nso_));

    double** Amnp = Amn->pointer();
    double** Amip = Ami->pointer();
    double** Cop = Caocc_->pointer();

    C_DGEMM('N','N', naux_ * (ULI) nso_, naocc_, nso_, 1.0, Amnp[0], nso_, Cop[0], naocc_,
        0.0, Amip[0], naocc_);

    Amn.reset();

    SharedMatrix Aia(new Matrix("Aij", naux_, naocc_ * (ULI) naocc_));
    double** Aiap = Aia->pointer();

    for (int Q = 0; Q < naux_; Q++) {
        C_DGEMM('T','N',naocc_,naocc_,nso_,1.0,Amip[Q],naocc_,Cop[0],naocc_, 0.0, Aiap[0], naocc_);
    }

    return Aia;
}
SharedMatrix PSTensor::Aov()
{
    SharedMatrix Amn = Aso();
    SharedMatrix Ami(new Matrix("Ami", naux_, naocc_ * (ULI) nso_));

    double** Amnp = Amn->pointer();
    double** Amip = Ami->pointer();
    double** Cop = Caocc_->pointer();
    double** Cvp = Cavir_->pointer();

    C_DGEMM('N','N', naux_ * (ULI) nso_, naocc_, nso_, 1.0, Amnp[0], nso_, Cop[0], naocc_,
        0.0, Amip[0], naocc_);

    Amn.reset();

    SharedMatrix Aia(new Matrix("Aia", naux_, naocc_ * (ULI) navir_));
    double** Aiap = Aia->pointer();

    for (int Q = 0; Q < naux_; Q++) {
        C_DGEMM('T','N',naocc_,navir_,nso_,1.0,Amip[Q],naocc_,Cvp[0],navir_, 0.0, Aiap[Q], navir_);
    }

    return Aia;
}
SharedMatrix PSTensor::Avv()
{
    SharedMatrix Amn = Aso();
    SharedMatrix Ami(new Matrix("Ami", naux_, navir_ * (ULI) nso_));

    double** Amnp = Amn->pointer();
    double** Amip = Ami->pointer();
    double** Cvp = Cavir_->pointer();

    C_DGEMM('N','N', naux_ * (ULI) nso_, navir_, nso_, 1.0, Amnp[0], nso_, Cvp[0], navir_,
        0.0, Amip[0], navir_);

    Amn.reset();

    SharedMatrix Aia(new Matrix("Aab", naux_, navir_ * (ULI) navir_));
    double** Aiap = Aia->pointer();

    for (int Q = 0; Q < naux_; Q++) {
        C_DGEMM('T','N',navir_,navir_,nso_,1.0,Amip[Q],navir_,Cvp[0],navir_, 0.0, Aiap[Q], navir_);
    }

    return Aia;
}
SharedMatrix PSTensor::Amo()
{
    SharedMatrix Amn = Aso();
    SharedMatrix Ami(new Matrix("Ami", naux_, nmo_ * (ULI) nso_));

    double** Amnp = Amn->pointer();
    double** Amip = Ami->pointer();
    double** Cvp = C_->pointer();

    C_DGEMM('N','N', naux_ * (ULI) nso_, nmo_, nso_, 1.0, Amnp[0], nso_, Cvp[0], nmo_,
        0.0, Amip[0], nmo_);

    Amn.reset();

    SharedMatrix Aia(new Matrix("Amo", naux_, nmo_ * (ULI) nmo_));
    double** Aiap = Aia->pointer();

    for (int Q = 0; Q < naux_; Q++) {
        C_DGEMM('T','N',nmo_,nmo_,nso_,1.0,Amip[Q],nmo_,Cvp[0],nmo_, 0.0, Aiap[Q], nmo_);
    }

    return Aia;
}
SharedMatrix PSTensor::Imo()
{
    boost::shared_ptr<MintsHelper> mints(new MintsHelper());
    return mints->mo_eri(C_,C_);
}
SharedMatrix PSTensor::Ipsmo()
{
    SharedMatrix Am = Amo();

    double** Qmop = Qmo_->pointer();
    double** Rmop = Rmo_->pointer();
    double** Amop = Am->pointer();

    SharedMatrix Imo(new Matrix("PS MO ERI Tensor", nmo_ * nmo_, nmo_ * nmo_));
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
