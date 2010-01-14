#include <stdexcept>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.hpp>
#include <psifiles.h>

#include <libmints/symmetry.h>

using namespace psi;

SymmetryOperation::SymmetryOperation()
{
    zero();
}

SymmetryOperation::SymmetryOperation(const SymmetryOperation &so)
{
    for (int i=0; i<3; ++i) {
        for (int j=0; j<3; ++j) {
            d_[i][j] = so.d_[i][j];
        }
    }
}

SymmetryOperation::~SymmetryOperation()
{
}

SymmetryOperation SymmetryOperation::operate(const SymmetryOperation& r) const
{
    SymmetryOperation ret;
    for (int i=0; i<3; ++i) {
        for (int j=0; j<3; ++j) {
            double t=0;
            for (int k=0; k<3; ++k) {
                t += r.d_[i][k] * d_[k][j];
            }
            ret. d_[i][j] = t;
        }
    }
    return ret;
}

SymmetryOperation SymmetryOperation::transform(const SymmetryOperation& r) const
{
    int i,j,k;
    SymmetryOperation ret,foo;

    // foo = r * d
    for (i=0; i < 3; i++) {
        for (j=0; j < 3; j++) {
            double t=0;
            for (k=0; k < 3; k++)
                t += r.d_[i][k] * d_[k][j];
            foo.d_[i][j] = t;
        }
    }

    // ret = (r*d)*r~ = foo*r~
    for (i=0; i < 3; i++) {
        for (j=0; j < 3; j++) {
            double t=0;
            for (k=0; k < 3; k++)
                t += foo.d_[i][k]*r.d_[j][k];
            ret.d_[i][j]=t;
        }
    }

    return ret;
}

void SymmetryOperation::rotation(int n)
{
    double theta = (n) ? 2.0 * M_PI / n : 2.0 * M_PI;
    rotation(theta);
}

void SymmetryOperation::rotation(double theta)
{
    zero();

    double ctheta = cos(theta);
    double stheta = sin(theta);

    d_[0][0] = ctheta;
    d_[0][1] = stheta;
    d_[1][0] = -stheta;
    d_[1][1] = ctheta;
    d_[2][2] = 1.0;
}

void SymmetryOperation::transpose()
{
    for (int i=1; i<3; ++i) {
        for (int j=0; j<i; ++j) {
            double tmp = d_[i][j];
            d_[i][j] = d_[j][i];
            d_[j][i] = tmp;
        }
    }
}

//
// SymRep
//
SymRep::SymRep(int i) : n_(i)
{
    zero();
}

SymRep::SymRep(const SymmetryOperation& so) : n_(3)
{
    memset(d_, 0, sizeof(double)*25);
    for (int i=0; i<3; ++i)
        for (int j=0; j<3; ++j)
            d_[i][j] = so[i][j];
}

SymRep::~SymRep()
{
    n_=0;
}

SymRep::operator SymmetryOperation() const
{
    if (n_ != 3)
        throw PSIEXCEPTION("Trying to cast to SymmetryOperation when n != 3");

    SymmetryOperation so;

    for (int i=0; i<3; ++i)
        for (int j=0; j<3; ++j)
            so[i][j] = d_[i][j];

    return so;
}

SymRep SymRep::operate(const SymRep& r) const
{
    if (r.n_ != n_)
        throw PSIEXCEPTION("SymRep::operate(): dimensions don't match");
    
    SymRep ret(n_);

    for (int i=0; i<n_; ++i) {
        for (int j=0; j<n_; ++j) {
            double t=0;
            for (int k=0; k<n_; ++k) {
                t += r[i][k] * d_[k][j];
            }
            ret[i][j] = t;
        }
    }

    return ret;
}

SymRep SymRep::transform(const SymRep& r) const
{
    int i, j, k;

    if (r.n_ != n_)
        throw PSIEXCEPTION("SymRep::transform(): dimensions don't match");

    SymRep ret(n_), foo(n_);

    // foo = r * d_
    for (i=0; i<n_; ++i) {
        for (j=0; j<n_; ++j) {
            double t=0;
            for (k=0; k<n_; ++k) {
                t += r[i][k] * d_[k][j];
            }
            foo[i][j] = t;
        }
    }

    // ret = (r*d)*r~ = foo*r~
    for (i=0; i<n_; ++i) {
        for (j=0; j<n_; ++j) {
            double t=0;
            for (k=0; k<n_; ++k) {
                t += foo[i][k] * r[j][k];
            }
            ret[i][j] = t;
        }
    }

    return ret;
}

void SymRep::sigma_h()
{
    unit();

    if (n_== 3) {
        d_[2][2] = -1.0;
    }
    else if (n_ == 5) {
        d_[3][3] = d_[4][4] = -1.0;
    }
}

void SymRep::sigma_xz()
{
    unit();

    if (n_==2 || n_==3 || n_==4) {
        d_[1][1] = -1.0;
        if (n_==4)
            d_[2][2] = -1.0;
    }
    else if (n_==5) {
        d_[2][2] = d_[4][4] = -1.0;
    }
}

void SymRep::sigma_yz()
{
    unit();

    if (n_==2 || n_==3 || n_==4) {
        d_[0][0] = -1.0;
        if (n_==4)
            d_[3][3] = -1.0;
    }
    else if (n_==5) {
        d_[2][2] = d_[3][3] = -1.0;
    }
}

void SymRep::rotation(int nt)
{
    double theta = (nt) ? 2.0*M_PI/nt : 2.0*M_PI;
    rotation(theta);
}

void SymRep::rotation(double theta)
{
    zero();

    double ctheta = cos(theta);
    double stheta = sin(theta);
    double c2theta = cos(2*theta);
    double s2theta = sin(2*theta);

    switch (n_) {
    case 1:
        d_[0][0] = 1.0;
        break;

    case 3:
        d_[0][0] = ctheta;
        d_[0][1] = stheta;
        d_[1][0] = -stheta;
        d_[1][1] = ctheta;
        d_[2][2] = 1.0;
        break;

    case 4:
    case 2:
        d_[0][0] = ctheta;
        d_[0][1] = stheta;
        d_[1][0] = -stheta;
        d_[1][1] = ctheta;

        // this is okay since d_ is hardwired
        d_[2][2] = c2theta;
        d_[2][3] = -s2theta;
        d_[3][2] = s2theta;
        d_[3][3] = c2theta;
        break;

    case 5:
        d_[0][0] = 1.0;

        d_[1][1] = c2theta;
        d_[1][2] = s2theta;
        d_[2][1] = -s2theta;
        d_[2][2] = c2theta;

        d_[3][3] = ctheta;
        d_[3][4] = -stheta;
        d_[4][3] = stheta;
        d_[4][4] = ctheta;
        break;

    default:
        throw PSIEXCEPTION("SymRep::n_ > 5");
    }
}

void SymRep::c2_x()
{
    i();

    if (n_ == 2 || n_ == 3 || n_ == 4) {
        d_[0][0] = 1.0;
        if (n_ == 4)
            d_[3][3] = 1.0;
    }
    else if (n_ == 5) {
        d_[0][0] = d_[1][1] = d_[4][4] = 1.0;
    }
}

void SymRep::c2_y()
{
    i();

    if (n_ == 2 || n_ == 3 || n_ == 4) {
        d_[1][1] = 1.0;
        if (n_ == 4)
            d_[2][2] = 1.0;
    }
    else if (n_ == 5) {
        d_[0][0] = d_[1][1] = d_[3][3] = 1.0;
    }
}

//
// IrreducibleRepresentation
//

IrreducibleRepresentation::IrreducibleRepresentation() :
  g(0), degen(0), nrot_(0), ntrans_(0), complex_(0), symb(0), rep(0), csymb(0)
{
}

IrreducibleRepresentation::IrreducibleRepresentation(
  int order, int d, const char *lab, const char *clab) :
  g_(0), degen_(0), nrot_(0), ntrans_(0), complex_(0), symb_(0), rep_(0), csymb_(0)
{
    init(order,d,lab,clab);
}


IrreducibleRepresentation::IrreducibleRepresentation(
  const IrreducibleRepresentation& ir) :
  g_(0), degen_(0), nrot_(0), ntrans_(0), complex_(0), symb_(0), rep_(0), csymb_(0)
{
    *this = ir;
}

IrreducibleRepresentation::~IrreducibleRepresentation()
{
    init();
}

IrreducibleRepresentation&
IrreducibleRepresentation::operator=(const IrreducibleRepresentation& ir)
{
    init(ir.g_,ir.degen_,ir.symb_,ir.csymb_);

    nrot_ = ir.nrot_;
    ntrans_ = ir.ntrans_;
    complex_ = ir.complex_;

    for (int i=0; i < g_; i++)
        rep_[i]=ir.rep_[i];

    return *this;
}

void
IrreducibleRepresentation::init(int order, int d, const char *lab,
                                const char *clab)
{
    g_=order;
    degen_=d;
    ntrans_=nrot_=complex_=0;

    delete[] symb_;
    symb_ = new_string(lab);

    delete[] csymb_;
    if (clab) csymb_ = new_string(clab);
    else csymb_ = 0;

    if (rep_) {
        delete[] rep_;
        rep_=0;
    }

    if (g_) {
        rep_ = new SymRep[g_];
        for (int i=0; i < g_; i++)
            rep[i].set_dim(d);
    }
}

