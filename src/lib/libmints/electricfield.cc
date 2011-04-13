#include <stdexcept>
#include <vector>
#include <libciomr/libciomr.h>

#include "mints.h"
#include <physconst.h>

#define MAX(a, b) ((a) > (b) ? (a) : (b))

using namespace boost;
using namespace psi;
using namespace std;

ElectricFieldInt::ElectricFieldInt(vector<SphericalTransform>& spherical_transforms, boost::shared_ptr<BasisSet> bs1, boost::shared_ptr<BasisSet> bs2, int nderiv) :
    OneBodyAOInt(spherical_transforms, bs1, bs2, nderiv), efield_recur_(bs1->max_am(), bs2->max_am()), natom_(bs1->molecule()->natom())
{
    int maxam1 = bs1_->max_am();
    int maxam2 = bs2_->max_am();

    int maxnao1 = INT_NCART(maxam1);
    int maxnao2 = INT_NCART(maxam2);

    buffer_ = new double[3*natom_*maxnao1*maxnao2];

    set_chunks(3*natom_);
}

ElectricFieldInt::~ElectricFieldInt()
{
    delete[] buffer_;
}

void ElectricFieldInt::compute_pair(const boost::shared_ptr<GaussianShell>& s1,
                                    const boost::shared_ptr<GaussianShell>& s2)
{
    int ao12;
    int am1 = s1->am();
    int am2 = s2->am();
    int nprim1 = s1->nprimitive();
    int nprim2 = s2->nprimitive();
    double A[3], B[3];
    A[0] = s1->center()[0];
    A[1] = s1->center()[1];
    A[2] = s1->center()[2];
    B[0] = s2->center()[0];
    B[1] = s2->center()[1];
    B[2] = s2->center()[2];

    int izm = 1;
    int iym = am1 + 1;
    int ixm = iym * iym;
    int jzm = 1;
    int jym = am2 + 1;
    int jxm = jym * jym;

    // Not sure if these are needed.
    int size =  INT_NCART(am1) * INT_NCART(am2);
    int ydisp = size;
    int zdisp = ydisp + size;

    // compute intermediates
    double AB2 = 0.0;
    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
    AB2 += (A[2] - B[2]) * (A[2] - B[2]);
    double AB = sqrt(AB2);

    memset(buffer_, 0, 3 * natom_ * s1->ncartesian() * s2->ncartesian() * sizeof(double));

    double ***ex = efield_recur_.ex();
    double ***ey = efield_recur_.ey();
    double ***ez = efield_recur_.ez();

//    double ***ex = efield_recur_.vx();
//    double ***ey = efield_recur_.vy();
//    double ***ez = efield_recur_.vz();

    for (int p1=0; p1<nprim1; ++p1) {
        double a1 = s1->exp(p1);
        double c1 = s1->coef(p1);
        for (int p2=0; p2<nprim2; ++p2) {
            double a2 = s2->exp(p2);
            double c2 = s2->coef(p2);
            double gamma = a1 + a2;
            double oog = 1.0 / gamma;

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

            // Loop over atoms of basis set 1 (only works if bs1_ and bs2_ are on the same molecule
            for (int atom=0; atom<bs1_->molecule()->natom(); ++atom) {
                int atomoffset = atom * 3 * size;

                double PC[3];
                double Z = (double)bs1_->molecule()->Z(atom);
                Vector3 C = bs1_->molecule()->xyz(atom);

                PC[0] = P[0] - C[0];
                PC[1] = P[1] - C[1];
                PC[2] = P[2] - C[2];

                // Get recursive
                efield_recur_.compute(PA, PB, PC, gamma, am1, am2);

                // Gather contributions.
                ao12 = 0;
                for (int ii = 0; ii <= am1; ++ii) {
                    int l1 = am1 - ii;
                    for (int jj = 0; jj <= ii; ++jj) {
                        int m1 = ii - jj;
                        int n1 = jj;

                        for (int kk = 0; kk <= am2; ++kk) {
                            int l2 = am2 - kk;
                            for (int ll = 0; ll <= kk; ++ll) {
                                int m2 = kk - ll;
                                int n2 = ll;

                                // Compute location in the recursion
                                int iind = l1 * ixm + m1 * iym + n1 * izm;
                                int jind = l2 * jxm + m2 * jym + n2 * jzm;

                                buffer_[atomoffset+ao12]       += ex[iind][jind][0] * over_pf;
                                buffer_[atomoffset+ao12+ydisp] += ey[iind][jind][0] * over_pf;
                                buffer_[atomoffset+ao12+zdisp] += ez[iind][jind][0] * over_pf;

//                                fprintf(outfile, "final: over_pf = %10.8lf AIX[%d][%d][0] = %10.8lf\tAIY = %10.8lf\tAIZ = %10.8lf\n", over_pf, iind, jind, ex[iind][jind][0], ey[iind][jind][0], ez[iind][jind][0]);

                                ao12++;
                            }
                        }
                    }
                }
            }
        }
    }
}
