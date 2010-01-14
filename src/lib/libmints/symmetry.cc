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
		t += r.d[i][k] * d[k][j];
	    foo.d[i][j] = t;
	}
    }

    // ret = (r*d)*r~ = foo*r~
    for (i=0; i < 3; i++) {
	for (j=0; j < 3; j++) {
	    double t=0;
	    for (k=0; k < 3; k++)
		t += foo.d[i][k]*r.d[j][k];
	    ret.d[i][j]=t;
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

