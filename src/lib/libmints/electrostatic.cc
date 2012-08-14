#include <libciomr/libciomr.h>

#include "mints.h"

#include <physconst.h>

#define MAX(a, b) ((a) > (b) ? (a) : (b))

using namespace boost;
using namespace psi;

// Initialize potential_recur_ to +1 basis set angular momentum
ElectrostaticInt::ElectrostaticInt(std::vector<SphericalTransform>& st, boost::shared_ptr<BasisSet> bs1, boost::shared_ptr<BasisSet> bs2, int deriv) :
    PotentialInt(st, bs1, bs2, deriv)
{
}

ElectrostaticInt::~ElectrostaticInt()
{

}

void ElectrostaticInt::compute(SharedMatrix &result, const Vector3& C)
{
    // Do not worry about zeroing out result
    int ns1 = bs1_->nshell();
    int ns2 = bs2_->nshell();

    int i_offset=0;
    double *location;

    // Leave as this full double for loop. We could be computing nonsymmetric integrals
    for (int i=0; i<ns1; ++i) {
        int ni = force_cartesian_ ? bs1_->shell(i).ncartesian() : bs1_->shell(i).nfunction();
        int j_offset=0;
        for (int j=0; j<ns2; ++j) {
            int nj = force_cartesian_ ? bs2_->shell(j).ncartesian() : bs2_->shell(j).nfunction();

            // Compute the shell (automatically transforms to pure am in needed)
            compute_shell(i, j, C);

            // For each integral that we got put in its contribution
            location = buffer_;
            for (int p=0; p<ni; ++p) {
                for (int q=0; q<nj; ++q) {
                    result->add(0, i_offset+p, j_offset+q, *location);
                    location++;
                }
            }

            j_offset += nj;
        }
        i_offset += ni;
    }
}

void ElectrostaticInt::compute_shell(int sh1, int sh2, const Vector3& C)
{
    const GaussianShell& s1 = bs1_->shell(sh1);
    const GaussianShell& s2 = bs2_->shell(sh2);

    // Call the child's compute_pair method, results better be in buffer_.
    compute_pair(s1, s2, C);

    // Normalize for angular momentum
    normalize_am(s1, s2, nchunk_);
    if(!force_cartesian_){
        // Pure angular momentum (6d->5d, ...) transformation
        pure_transform(s1, s2, nchunk_);
    }
}

// The engine only supports segmented basis sets
void ElectrostaticInt::compute_pair(const GaussianShell& s1, const GaussianShell& s2, const Vector3& C)
{
    int ao12;
    int am1 = s1.am();
    int am2 = s2.am();
    int nprim1 = s1.nprimitive();
    int nprim2 = s2.nprimitive();
    double A[3], B[3];
    A[0] = s1.center()[0];
    A[1] = s1.center()[1];
    A[2] = s1.center()[2];
    B[0] = s2.center()[0];
    B[1] = s2.center()[1];
    B[2] = s2.center()[2];

    int izm = 1;
    int iym = am1 + 1;
    int ixm = iym * iym;
    int jzm = 1;
    int jym = am2 + 1;
    int jxm = jym * jym;

    // compute intermediates
    double AB2 = 0.0;
    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
    AB2 += (A[2] - B[2]) * (A[2] - B[2]);

    memset(buffer_, 0, s1.ncartesian() * s2.ncartesian() * sizeof(double));

    double ***vi = potential_recur_->vi();

    for (int p1=0; p1<nprim1; ++p1) {
        double a1 = s1.exp(p1);
        double c1 = s1.coef(p1);
        for (int p2=0; p2<nprim2; ++p2) {
            double a2 = s2.exp(p2);
            double c2 = s2.coef(p2);
            double gamma = a1 + a2;
            double oog = 1.0/gamma;

            double PA[3], PB[3];
            double P[3];

            P[0] = (a1*A[0] + a2*B[0])*oog;
            P[1] = (a1*A[1] + a2*B[1])*oog;
            P[2] = (a1*A[2] + a2*B[2])*oog;
            PA[0] = P[0] - A[0];
            PA[1] = P[1] - A[1];
            PA[2] = P[2] - A[2];
            PB[0] = P[0] - B[0];
            PB[1] = P[1] - B[1];
            PB[2] = P[2] - B[2];

            double over_pf = exp(-a1*a2*AB2*oog) * sqrt(M_PI*oog) * M_PI * oog * c1 * c2;

            // Loop over atoms of basis set 1 (only works if bs1_ and bs2_ are on the same
            // molecule)
            double PC[3];

            PC[0] = P[0] - C[0];
            PC[1] = P[1] - C[1];
            PC[2] = P[2] - C[2];

            // Do recursion
            potential_recur_->compute(PA, PB, PC, gamma, am1, am2);

            ao12 = 0;
            for(int ii = 0; ii <= am1; ii++) {
                int l1 = am1 - ii;
                for(int jj = 0; jj <= ii; jj++) {
                    int m1 = ii - jj;
                    int n1 = jj;
                    /*--- create all am components of sj ---*/
                    for(int kk = 0; kk <= am2; kk++) {
                        int l2 = am2 - kk;
                        for(int ll = 0; ll <= kk; ll++) {
                            int m2 = kk - ll;
                            int n2 = ll;

                            // Compute location in the recursion
                            int iind = l1 * ixm + m1 * iym + n1 * izm;
                            int jind = l2 * jxm + m2 * jym + n2 * jzm;

                            buffer_[ao12++] += -vi[iind][jind][0] * over_pf;
                        }
                    }
                }
            }
        }
    }
}

SharedVector ElectrostaticInt::nuclear_contribution(boost::shared_ptr<Molecule> mol)
{
    boost::shared_ptr<Vector> sret(new Vector(mol->natom()));
    double *ret = sret->pointer();

    int natom = mol->natom();
    for(int k=0;k<natom;k++) {
        Vector3 kgeom = mol->xyz(k);
        for(int i=0;i<natom;i++) {
            if (i != k) {
                Vector3 igeom = mol->xyz(i);

                double x = kgeom[0] - igeom[0];
                double y = kgeom[1] - igeom[1];
                double z = kgeom[2] - igeom[2];
                double r2 = x*x+y*y+z*z;
                double r = std::sqrt(r2);
                ret[k]  += mol->Z(i)/r;
            }
        }
    }

    return sret;
}
