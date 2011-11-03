#include <libciomr/libciomr.h>

#include "mints.h"

#include <physconst.h>

#define MAX(a, b) ((a) > (b) ? (a) : (b))

using namespace boost;
using namespace psi;

// Initialize potential_recur_ to +1 basis set angular momentum
PotentialInt::PotentialInt(std::vector<SphericalTransform>& st, boost::shared_ptr<BasisSet> bs1, boost::shared_ptr<BasisSet> bs2, int deriv) :
    OneBodyAOInt(st, bs1, bs2, deriv), potential_recur_(bs1->max_am()+1, bs2->max_am()+1),
    potential_deriv_recur_(bs1->max_am()+2, bs2->max_am()+2)
{
    if (deriv == 0)
        potential_recur_ = ObaraSaikaTwoCenterVIRecursion(bs1->max_am()+1, bs2->max_am()+1);
    else if (deriv == 1)
        potential_recur_ = ObaraSaikaTwoCenterVIDerivRecursion(bs1->max_am()+2, bs2->max_am()+2);
    else if (deriv == 2)
        potential_recur_ = ObaraSaikaTwoCenterVIDeriv2Recursion(bs1->max_am()+3, bs2->max_am()+3);
    else
        throw PSIEXCEPTION("PotentialInt: deriv > 2 is not supported.");

    int maxam1 = bs1_->max_am();
    int maxam2 = bs2_->max_am();

    int maxnao1 = INT_NCART(maxam1);
    int maxnao2 = INT_NCART(maxam2);

    if (deriv == 1) {
        // We set chunk count for normalize_am and pure_transform
        // We can't use the trick of using less memory that I implemented in overlap & kinetic
        // since potential integral derivatives also have a contribution to center c...which is
        // over all atoms.
        set_chunks(3*natom_);

        maxnao1 *= 3*natom_;
        maxnao2 *= 3*natom_;
    }

    buffer_ = new double[maxnao1*maxnao2];

    // Setup the initial field of partial charges
    Zxyz_ = SharedMatrix (new Matrix("Partial Charge Field (Z,x,y,z)", bs1_->molecule()->natom(), 4));
    double** Zxyzp = Zxyz_->pointer();

    for (int A = 0; A < bs1_->molecule()->natom(); A++) {
        Zxyzp[A][0] = (double) bs1_->molecule()->Z(A);
        Zxyzp[A][1] = bs1_->molecule()->x(A);
        Zxyzp[A][2] = bs1_->molecule()->y(A);
        Zxyzp[A][3] = bs1_->molecule()->z(A);
    }
}

PotentialInt::~PotentialInt()
{
    delete[] buffer_;
    delete potential_recur_;
}

// The engine only supports segmented basis sets
void PotentialInt::compute_pair(const boost::shared_ptr<GaussianShell>& s1,
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

    // compute intermediates
    double AB2 = 0.0;
    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
    AB2 += (A[2] - B[2]) * (A[2] - B[2]);

    memset(buffer_, 0, s1->ncartesian() * s2->ncartesian() * sizeof(double));

    double ***vi = potential_recur_->vi();

    double** Zxyzp = Zxyz_->pointer();
    int ncharge = Zxyz_->rowspi()[0];

    for (int p1=0; p1<nprim1; ++p1) {
        double a1 = s1->exp(p1);
        double c1 = s1->coef(p1);
        for (int p2=0; p2<nprim2; ++p2) {
            double a2 = s2->exp(p2);
            double c2 = s2->coef(p2);
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
            for (int atom=0; atom<ncharge; ++atom) {
                double PC[3];

                double Z = Zxyzp[atom][0];

                PC[0] = P[0] - Zxyzp[atom][1];
                PC[1] = P[1] - Zxyzp[atom][2];
                PC[2] = P[2] - Zxyzp[atom][3];

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

                                buffer_[ao12++] += -vi[iind][jind][0] * over_pf * Z;

//                                fprintf(outfile, "ao12=%d, vi[%d][%d][0] = %20.14f, over_pf = %20.14f, Z = %f\n", ao12-1, iind, jind, vi[iind][jind][0], over_pf, Z);
                            }
                        }
                    }
                }
            }
        }
    }
}

void PotentialInt::compute_pair_deriv1(const boost::shared_ptr<GaussianShell>& s1, const boost::shared_ptr<GaussianShell>& s2)
{
    int ao12;
    const int am1 = s1->am();
    const int am2 = s2->am();
    const int nprim1 = s1->nprimitive();
    const int nprim2 = s2->nprimitive();
    const int ncenteri = s1->ncenter();
    const int ncenterj = s2->ncenter();

    double A[3], B[3];
    A[0] = s1->center()[0];
    A[1] = s1->center()[1];
    A[2] = s1->center()[2];
    B[0] = s2->center()[0];
    B[1] = s2->center()[1];
    B[2] = s2->center()[2];

    // size of the length of a perturbation
    const size_t size = s1->ncartesian() * s2->ncartesian();
    const int center_i = ncenteri * 3 * size;
    const int center_j = ncenterj * 3 * size;

    const int izm1 = 1;
    const int iym1 = am1 + 1 + 1;  // extra 1 for derivative
    const int ixm1 = iym1 * iym1;
    const int jzm1 = 1;
    const int jym1 = am2 + 1 + 1;  // extra 1 for derivative
    const int jxm1 = jym1 * jym1;

    // compute intermediates
    double AB2 = 0.0;
    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
    AB2 += (A[2] - B[2]) * (A[2] - B[2]);

    memset(buffer_, 0, 3 * natom_ * size * sizeof(double));

    double ***vi = potential_recur_->vi();
    double ***vx = potential_recur_->vx();
    double ***vy = potential_recur_->vy();
    double ***vz = potential_recur_->vz();

    double** Zxyzp = Zxyz_->pointer();
    int ncharge = Zxyz_->rowspi()[0];

    for (int p1=0; p1<nprim1; ++p1) {
        double a1 = s1->exp(p1);
        double c1 = s1->coef(p1);
        for (int p2=0; p2<nprim2; ++p2) {
            double a2 = s2->exp(p2);
            double c2 = s2->coef(p2);
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
            for (int atom=0; atom<ncharge; ++atom) {
                double PC[3];

                double Z = Zxyzp[atom][0];

                PC[0] = P[0] - Zxyzp[atom][1];
                PC[1] = P[1] - Zxyzp[atom][2];
                PC[2] = P[2] - Zxyzp[atom][3];

                // Do recursion
                potential_recur_->compute(PA, PB, PC, gamma, am1+1, am2+1);

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
                                int iind = l1 * ixm1 + m1 * iym1 + n1 * izm1;
                                int jind = l2 * jxm1 + m2 * jym1 + n2 * jzm1;

                                const double pfac = over_pf * Z;

                                // x
                                double temp = 2.0*a1*vi[iind+ixm1][jind][0];
                                if (l1)
                                    temp -= l1*vi[iind-ixm1][jind][0];
                                buffer_[center_i+(0*size)+ao12] -= temp * pfac;
                                // printf("ix temp = %f ", temp);

                                temp = 2.0*a2*vi[iind][jind+jxm1][0];
                                if (l2)
                                    temp -= l2*vi[iind][jind-jxm1][0];
                                buffer_[center_j+(0*size)+ao12] -= temp * pfac;
                                // printf("jx temp = %f ", temp);

                                buffer_[3*size*atom+ao12] -= vx[iind][jind][0] * pfac;

                                // y
                                temp = 2.0*a1*vi[iind+iym1][jind][0];
                                if (m1)
                                    temp -= m1*vi[iind-iym1][jind][0];
                                buffer_[center_i+(1*size)+ao12] -= temp * pfac;
                                // printf("iy temp = %f ", temp);

                                temp = 2.0*a2*vi[iind][jind+jym1][0];
                                if (m2)
                                    temp -= m2*vi[iind][jind-jym1][0];
                                buffer_[center_j+(1*size)+ao12] -= temp * pfac;
                                // printf("jy temp = %f ", temp);

                                buffer_[3*size*atom+size+ao12] -= vy[iind][jind][0] * pfac;

                                // z
                                temp = 2.0*a1*vi[iind+izm1][jind][0];
                                if (n1)
                                    temp -= n1*vi[iind-izm1][jind][0];
                                buffer_[center_i+(2*size)+ao12] -= temp * pfac;
                                // printf("iz temp = %f ", temp);

                                temp = 2.0*a2*vi[iind][jind+jzm1][0];
                                if (n2)
                                    temp -= n2*vi[iind][jind-jzm1][0];
                                buffer_[center_j+(2*size)+ao12] -= temp * pfac;
                                // printf("jz temp = %f \n", temp);

                                buffer_[3*size*atom+2*size+ao12] -= vz[iind][jind][0] * pfac;

                                ao12++;
                            }
                        }
                    }
                }
            }
        }
    }
}

void PotentialInt::compute_pair_deriv2(const boost::shared_ptr<GaussianShell>& s1, const boost::shared_ptr<GaussianShell>& s2)
{
    int ao12;
    const int am1 = s1->am();
    const int am2 = s2->am();
    const int nprim1 = s1->nprimitive();
    const int nprim2 = s2->nprimitive();
    const int ncenteri = s1->ncenter();
    const int ncenterj = s2->ncenter();

    double A[3], B[3];
    A[0] = s1->center()[0];
    A[1] = s1->center()[1];
    A[2] = s1->center()[2];
    B[0] = s2->center()[0];
    B[1] = s2->center()[1];
    B[2] = s2->center()[2];

    // size of the length of a perturbation
    const size_t size = s1->ncartesian() * s2->ncartesian();
    const int center_i = ncenteri * 3 * size;
    const int center_j = ncenterj * 3 * size;

    const int izm1 = 1;
    const int iym1 = am1 + 1 + 1;  // extra 1 for derivative
    const int ixm1 = iym1 * iym1;
    const int jzm1 = 1;
    const int jym1 = am2 + 1 + 1;  // extra 1 for derivative
    const int jxm1 = jym1 * jym1;

    // compute intermediates
    double AB2 = 0.0;
    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
    AB2 += (A[2] - B[2]) * (A[2] - B[2]);

    memset(buffer_, 0, 3 * natom_ * size * sizeof(double));

    double ***vi = potential_recur_->vi();
    double ***vx = potential_recur_->vx();
    double ***vy = potential_recur_->vy();
    double ***vz = potential_recur_->vz();
    double ***vxx = potential_recur_->vxx();
    double ***vxy = potential_recur_->vxx();
    double ***vxz = potential_recur_->vxx();
    double ***vyy = potential_recur_->vxx();
    double ***vyz = potential_recur_->vxx();
    double ***vzz = potential_recur_->vxx();

    double** Zxyzp = Zxyz_->pointer();
    int ncharge = Zxyz_->rowspi()[0];

    for (int p1=0; p1<nprim1; ++p1) {
        double a1 = s1->exp(p1);
        double c1 = s1->coef(p1);
        for (int p2=0; p2<nprim2; ++p2) {
            double a2 = s2->exp(p2);
            double c2 = s2->coef(p2);
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
            for (int atom=0; atom<ncharge; ++atom) {
                double PC[3];

                double Z = Zxyzp[atom][0];

                PC[0] = P[0] - Zxyzp[atom][1];
                PC[1] = P[1] - Zxyzp[atom][2];
                PC[2] = P[2] - Zxyzp[atom][3];

                // Do recursion
                potential_recur_->compute(PA, PB, PC, gamma, am1+1, am2+1);

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
                                int iind = l1 * ixm1 + m1 * iym1 + n1 * izm1;
                                int jind = l2 * jxm1 + m2 * jym1 + n2 * jzm1;

                                const double pfac = over_pf * Z;

                                // x
                                double temp = 2.0*a1*vi[iind+ixm1][jind][0];
                                if (l1)
                                    temp -= l1*vi[iind-ixm1][jind][0];
                                buffer_[center_i+(0*size)+ao12] -= temp * pfac;
                                // printf("ix temp = %f ", temp);

                                temp = 2.0*a2*vi[iind][jind+jxm1][0];
                                if (l2)
                                    temp -= l2*vi[iind][jind-jxm1][0];
                                buffer_[center_j+(0*size)+ao12] -= temp * pfac;
                                // printf("jx temp = %f ", temp);

                                buffer_[3*size*atom+ao12] -= vx[iind][jind][0] * pfac;

                                // y
                                temp = 2.0*a1*vi[iind+iym1][jind][0];
                                if (m1)
                                    temp -= m1*vi[iind-iym1][jind][0];
                                buffer_[center_i+(1*size)+ao12] -= temp * pfac;
                                // printf("iy temp = %f ", temp);

                                temp = 2.0*a2*vi[iind][jind+jym1][0];
                                if (m2)
                                    temp -= m2*vi[iind][jind-jym1][0];
                                buffer_[center_j+(1*size)+ao12] -= temp * pfac;
                                // printf("jy temp = %f ", temp);

                                buffer_[3*size*atom+size+ao12] -= vy[iind][jind][0] * pfac;

                                // z
                                temp = 2.0*a1*vi[iind+izm1][jind][0];
                                if (n1)
                                    temp -= n1*vi[iind-izm1][jind][0];
                                buffer_[center_i+(2*size)+ao12] -= temp * pfac;
                                // printf("iz temp = %f ", temp);

                                temp = 2.0*a2*vi[iind][jind+jzm1][0];
                                if (n2)
                                    temp -= n2*vi[iind][jind-jzm1][0];
                                buffer_[center_j+(2*size)+ao12] -= temp * pfac;
                                // printf("jz temp = %f \n", temp);

                                buffer_[3*size*atom+2*size+ao12] -= vz[iind][jind][0] * pfac;

                                ao12++;
                            }
                        }
                    }
                }
            }
        }
    }
}

void PotentialInt::compute_deriv1(std::vector<boost::shared_ptr<SimpleMatrix> > &result)
{
    if (deriv_ < 1)
        throw SanityCheckError("PotentialInt::compute_deriv1(result): integral object not created to handle derivatives.", __FILE__, __LINE__);

    // Do not worry about zeroing out result
    int ns1 = bs1_->nshell();
    int ns2 = bs2_->nshell();
    int result_size = result.size();
    int i_offset = 0;
    double *location = 0;

    // Check the length of result, must be 3*natom_
    if (result.size() != 3*natom_)
        throw SanityCheckError("PotentialInt::compute_derv1(result): result must be 3 * natom in length.", __FILE__, __LINE__);

    for (int i=0; i<ns1; ++i) {
        int ni = force_cartesian_ ? bs1_->shell(i)->ncartesian() : bs1_->shell(i)->nfunction();
        int j_offset=0;
        for (int j=0; j<ns2; ++j) {
            int nj = force_cartesian_ ? bs2_->shell(j)->ncartesian() : bs2_->shell(j)->nfunction();

            // Compute the shell
            compute_shell_deriv1(i, j);

            // For each integral that we got put in its contribution
            location = buffer_;
            for (int r=0; r<result_size; ++r) {
                for (int p=0; p<ni; ++p) {
                    for (int q=0; q<nj; ++q) {
                        result[r]->add(i_offset+p, j_offset+q, *location);
                        location++;
                    }
                }
            }
            j_offset += nj;
        }
        i_offset += ni;
    }
}

void PotentialInt::compute_deriv1(std::vector<SharedMatrix > &result)
{
    if (deriv_ < 1)
        throw SanityCheckError("PotentialInt::compute_deriv1(result): integral object not created to handle derivatives.", __FILE__, __LINE__);

    // Do not worry about zeroing out result
    int ns1 = bs1_->nshell();
    int ns2 = bs2_->nshell();
    int result_size = result.size();
    int i_offset = 0;
    double *location = 0;

    // Check the length of result, must be 3*natom_
    if (result.size() != 3*natom_)
        throw SanityCheckError("PotentialInt::compute_derv1(result): result must be 3 * natom in length.", __FILE__, __LINE__);

    for (int i=0; i<ns1; ++i) {
        int ni = force_cartesian_ ? bs1_->shell(i)->ncartesian() : bs1_->shell(i)->nfunction();
        int j_offset=0;
        for (int j=0; j<ns2; ++j) {
            int nj = force_cartesian_ ? bs2_->shell(j)->ncartesian() : bs2_->shell(j)->nfunction();

            // Compute the shell
            compute_shell_deriv1(i, j);

            // For each integral that we got put in its contribution
            location = buffer_;
            for (int r=0; r<result_size; ++r) {
                for (int p=0; p<ni; ++p) {
                    for (int q=0; q<nj; ++q) {
                        result[r]->add(0, i_offset+p, j_offset+q, *location);
                        location++;
                    }
                }
            }
            j_offset += nj;
        }
        i_offset += ni;
    }
}

PotentialSOInt::PotentialSOInt(const boost::shared_ptr<OneBodyAOInt> &aoint, const boost::shared_ptr<IntegralFactory> &fact)
    : OneBodySOInt(aoint, fact)
{
    natom_ = ob_->basis1()->molecule()->natom();
}

PotentialSOInt::PotentialSOInt(const boost::shared_ptr<OneBodyAOInt> &aoint, const IntegralFactory *fact)
    : OneBodySOInt(aoint, fact)
{
    natom_ = ob_->basis1()->molecule()->natom();
}

void PotentialSOInt::compute_deriv1(std::vector<SharedMatrix > result,
                                    const CdSalcList &cdsalcs)
{
    // Do not worry about zeroing out result.

    // Do some checks:
    if (deriv_ < 1)
        throw SanityCheckError("OneBodySOInt::compute_deriv1: integral object not created to handle derivatives.", __FILE__, __LINE__);

    if (result.size() != cdsalcs.ncd())
        throw SanityCheckError("OneBodySOInt::compute_deriv1: result vector size does not match SALC size.", __FILE__, __LINE__);

    int ns1 = b1_->nshell();
    int ns2 = b2_->nshell();
    const double *aobuf = ob_->buffer();

    // Loop over unique SO shells.
    for (int ish=0; ish<ns1; ++ish) {
        const SOTransform& t1 = b1_->sotrans(ish);
        int nso1 = b1_->nfunction(ish);
        int nao1 = b1_->naofunction(ish);

        for (int jsh=0; jsh<ns2; ++jsh) {
            const SOTransform& t2= b2_->sotrans(jsh);
            int nso2 = b2_->nfunction(jsh);
            int nao2 = b2_->naofunction(jsh);

            int nao12 = nao1 * nao2;
            int nso12 = nso1 * nso2;

            // loop through the AO shells that make up this SO shell
            // by the end of these 4 for loops we will have our final integral in buffer_
            for (int i=0; i<t1.naoshell; ++i) {
                const SOTransformShell &s1 = t1.aoshell[i];
                int center_i = ob_->basis1()->shell(s1.aoshell)->ncenter();

                for (int j=0; j<t2.naoshell; ++j) {
                    const SOTransformShell &s2 = t2.aoshell[j];
                    int center_j = ob_->basis2()->shell(s2.aoshell)->ncenter();

                    // If we're working on the same atomic center, don't even bother with the derivative
                    // Does this still hold for potentials? nope
//                    if (center_i == center_j)
//                        continue;

                    ob_->compute_shell_deriv1(s1.aoshell, s2.aoshell);

                    // handle SO transform
                    for (int itr=0; itr<s1.nfunc; ++itr) {
                        const SOTransformFunction &ifunc = s1.func[itr];
                        // SO transform coefficient
                        double icoef = ifunc.coef;
                        // AO function offset in a linear array
                        int iaofunc  = ifunc.aofunc;
                        // SO function offset in a linear array
                        int isofunc  = b1_->function_offset_within_shell(ish, ifunc.irrep) + ifunc.sofunc;
                        // AO function offset in a linear array
                        int iaooff   = iaofunc;
                        // SO function offset in a lienar array
                        int isooff   = isofunc;
                        // Relative position of the SO function within its irrep
                        int irel     = b1_->function_within_irrep(ish, isofunc);
                        int iirrep   = ifunc.irrep;

                        for (int jtr=0; jtr<s2.nfunc; ++jtr) {
                            const SOTransformFunction &jfunc = s2.func[jtr];
                            double jcoef = jfunc.coef * icoef;
                            int jaofunc  = jfunc.aofunc;
                            int jsofunc  = b2_->function_offset_within_shell(jsh, jfunc.irrep) + jfunc.sofunc;
                            int jaooff   = iaooff*nao2 + jaofunc;
                            int jsooff   = isooff*nso2 + jsofunc;
                            int jrel     = b2_->function_within_irrep(jsh, jsofunc);
                            int jirrep   = jfunc.irrep;

                            // Need to loop over the cdsalcs over ALL atoms
                            // Potential integral derivatives include contribution
                            // to a third atom.

                            // third atom loop (actually goes over all atoms)
                            for (int a=0; a<natom_; ++a) {
                                const CdSalcWRTAtom& cdsalc1 = cdsalcs.atom_salc(a);
                                int offset = jaooff + 3*a*nao12;

                                double jcoef_aobuf = jcoef * aobuf[offset+(0*nao12)];
                                for (int nx=0; nx<cdsalc1.nx(); ++nx) {
                                    const CdSalcWRTAtom::Component element = cdsalc1.x(nx);
                                    double temp = jcoef_aobuf * element.coef;
                                    if ((iirrep ^ jirrep) == element.irrep && fabs(temp) > 1.0e-10) {
                                        result[element.salc]->add(iirrep, irel, jrel, temp);
                                    }
                                }

                                jcoef_aobuf = jcoef * aobuf[offset+(1*nao12)];
                                for (int ny=0; ny<cdsalc1.ny(); ++ny) {
                                    const CdSalcWRTAtom::Component element = cdsalc1.y(ny);
                                    double temp = jcoef_aobuf * element.coef;
                                    if ((iirrep ^ jirrep) == element.irrep && fabs(temp) > 1.0e-10) {
                                        result[element.salc]->add(iirrep, irel, jrel, temp);
                                    }
                                }

                                jcoef_aobuf = jcoef * aobuf[offset+(2*nao12)];
                                for (int nz=0; nz<cdsalc1.nz(); ++nz) {
                                    const CdSalcWRTAtom::Component element = cdsalc1.z(nz);
                                    double temp = jcoef_aobuf * element.coef;
                                    if ((iirrep ^ jirrep) == element.irrep && fabs(temp) > 1.0e-10) {
                                        result[element.salc]->add(iirrep, irel, jrel, temp);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
