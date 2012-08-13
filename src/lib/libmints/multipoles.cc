#include "mints.h"

using namespace psi;
using namespace boost;

MultipoleInt::MultipoleInt(std::vector<SphericalTransform>& spherical_transforms, boost::shared_ptr<BasisSet> bs1, boost::shared_ptr<BasisSet> bs2, int order, int nderiv) :
    order_(order),
    OneBodyAOInt(spherical_transforms, bs1, bs2, nderiv), mi_recur_(bs1->max_am()+2, bs2->max_am()+2, order)
{
    int maxam1 = bs1_->max_am();
    int maxam2 = bs2_->max_am();

    int maxnao1 = (maxam1+1)*(maxam1+2)/2;
    int maxnao2 = (maxam2+1)*(maxam2+2)/2;

    // The number of multipole components to compute.  N.B. we don't compute the 0th one
    int n_mult = (order_+1)*(order_+2)*(order_+3)/6 - 1;

    // Increase buffer size to handle x, y, and z components
    if (deriv_ == 0) {
        buffer_ = new double[n_mult*maxnao1*maxnao2];
        set_chunks(n_mult);
    }
    else {
        throw PSIEXCEPTION("Derivatives are NYI for arbitrary-order multipoles");
    }
}

MultipoleInt::~MultipoleInt()
{
    delete[] buffer_;
}

SharedVector MultipoleInt::nuclear_contribution(boost::shared_ptr<Molecule> mol, int order)
{
    int ntot = (order+1)*(order+2)*(order+3)/6 - 1;
    boost::shared_ptr<Vector> sret(new Vector(ntot));
    double *ret = sret->pointer();

    int address = 0;
    for(int l = 1; l <= order; ++l){
        for(int ii = 0; ii <= l; ii++) {
            int lx = l - ii;
            for(int lz = 0; lz <= ii; lz++) {
                int ly = ii - lz;
                for(int atom = 0; atom < mol->natom(); ++atom) {
                    Vector3 geom = mol->xyz(atom);
                    ret[address] += mol->Z(atom)*pow(geom[0], lx)*pow(geom[1], ly)*pow(geom[2], lz);
                }
                ++address;
            }
        }
    }

    return sret;
}

inline uint64_t binomial(int n, int c1)
{
    uint64_t num = 1;
    uint64_t den = 1;
    int c2 = n - c1;
    int i;
    for (i=c2+1; i<=n; i++) {
        num *= i;
    }
    for (i=2; i<=c1; i++) {
        den *= i;
    }
    return num/den;
}


// The engine only supports segmented basis sets
void MultipoleInt::compute_pair(const GaussianShell& s1, const GaussianShell& s2)
{
    int am1 = s1.am();
    int am2 = s2.am();
    int nprim1 = s1.nprimitive();
    int nprim2 = s2.nprimitive();

    // The number of bf components in each shell pair
    int stride = INT_NCART(am1) * INT_NCART(am2);

    memset(buffer_, 0, nchunk_ * INT_NCART(am1) * INT_NCART(am2) * sizeof(double));

    // Buffers to hold the {x,y,z} moments
    double *Xpowers = new double[order_ + 1];
    double *Ypowers = new double[order_ + 1];
    double *Zpowers = new double[order_ + 1];

    double A[3], B[3];
    A[0] = s1.center()[0];
    A[1] = s1.center()[1];
    A[2] = s1.center()[2];
    B[0] = s2.center()[0];
    B[1] = s2.center()[1];
    B[2] = s2.center()[2];

    // compute intermediates
    double AB2 = 0.0;
    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
    AB2 += (A[2] - B[2]) * (A[2] - B[2]);

    double ***x = mi_recur_.x();
    double ***y = mi_recur_.y();
    double ***z = mi_recur_.z();


    for (int p1=0; p1<nprim1; ++p1) {
        double a1 = s1.exp(p1);
        double c1 = s1.coef(p1);
        for (int p2=0; p2<nprim2; ++p2) {
            double a2 = s2.exp(p2);
            double c2 = s2.coef(p2);
            double gamma = a1 + a2;
            double oog = 1.0/gamma;

            double PA[3], PB[3], PC[3];
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
            double PCx = (P[0] - origin_[0]);
            double PCy = (P[1] - origin_[1]);
            double PCz = (P[2] - origin_[2]);

            double over_pf = exp(-a1*a2*AB2*oog) * sqrt(M_PI*oog) * M_PI * oog * c1 * c2;

            // Do recursion
            mi_recur_.compute(PA, PB, gamma, am1+2, am2+2);
            int bf = 0;

            // Basis function A.M. components on center 1
            for(int ii = 0; ii <= am1; ii++) {
                int lx1 = am1 - ii;
                for(int lz1 = 0; lz1 <= ii; lz1++) {
                    int ly1 = ii - lz1;
                    // Basis function A.M. components on center 2
                    for(int kk = 0; kk <= am2; kk++) {
                        int lx2 = am2 - kk;
                        for(int lz2 = 0; lz2 <= kk; lz2++) {
                            int ly2 = kk - lz2;

                            int chunk = 0;

                            /*
                             * Define X_C = X_P + X_PC, then expand X_C in powers,
                             * up to the requested order; these will be combined to
                             * form the total moment operators below
                             */
                            // Xpowers[0] = (l) * X_C^0 * X_PC^0
                            //              (0)
                            Xpowers[0] = x[lx1][lx2][0];
                            Ypowers[0] = y[ly1][ly2][0];
                            Zpowers[0] = z[lz1][lz2][0];
                            for(int l = 1; l <= order_; ++l){
                                double px = PCx;
                                double py = PCy;
                                double pz = PCz;
                                // Xpowers[l] = (l) * X_C^0 * X_PC^0  +  (l) * X_C^0 * X_PC^l
                                //              (0)                      (l)
                                Xpowers[l] = x[lx1][lx2][l] + pow(px, l) * x[lx1][lx2][0];
                                Ypowers[l] = y[ly1][ly2][l] + pow(py, l) * y[ly1][ly2][0];
                                Zpowers[l] = z[lz1][lz2][l] + pow(pz, l) * z[lz1][lz2][0];
                                for(int i = 1; i < l; ++i){
                                    double coef = (double) binomial(l, i);
                                    // Xpowers[l] = (l) * X_C^(l-i) * X_PC^i
                                    //              (i)
                                    Xpowers[l] += coef * px * x[lx1][lx2][l-i];
                                    Ypowers[l] += coef * py * y[ly1][ly2][l-i];
                                    Zpowers[l] += coef * pz * z[lz1][lz2][l-i];
                                    px *= PCx;
                                    py *= PCy;
                                    pz *= PCz;
                                }

                                /*
                                 * Loop over the multipole components, and combine the moments computed above
                                 */
                                for(int ii = 0; ii <= l; ii++) {
                                    int lx = l - ii;
                                    for(int lz = 0; lz <= ii; lz++) {
                                        int ly = ii - lz;
                                        double val = Xpowers[lx] * Ypowers[ly] * Zpowers[lz] * over_pf;
                                        buffer_[chunk*stride + bf] -= val;
                                        ++chunk;
                                    }
                                }
                            } //End loop over l
                            bf++;
                        }
                    } //End loop over shell 2 A.M. components
                }
            } //End loop over shell 1 A.M. components
        } //End loop over shell 2 primitives
    } //End loop over shell 1 primitives
    delete [] Xpowers;
    delete [] Ypowers;
    delete [] Zpowers;

}

