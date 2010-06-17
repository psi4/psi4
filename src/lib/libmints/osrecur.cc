#include <stdexcept>
#include <libciomr/libciomr.h>
#include <libmints/integral.h>
#include <libmints/wavefunction.h>   // for df
#include <libmints/osrecur.h>
#include <exception.h>

using namespace psi;

double ***init_box(int a, int b, int c)
{
    int i,j;
    double ***box;

    box = (double ***) malloc(sizeof(double **)*a);
    for(i=0;i<a;i++)
        box[i] = (double **) malloc(sizeof(double *)*b);
    for(i=0;i<a;i++)
        for(j=0;j<b;j++) {
            box[i][j] = (double *) malloc(sizeof(double)*c);
            memset((void *)box[i][j], '\0', sizeof(double)*c);
        }

    return box;
}

void free_box(double ***box, int a, int b)
{
    int i,j;

    for(i=0;i<a;i++)
        for(j=0;j<b;j++)
            free(box[i][j]);

    for(i=0;i<a;i++)
        free(box[i]);

    free(box);
}

ObaraSaikaTwoCenterMIRecursion::ObaraSaikaTwoCenterMIRecursion(int max_am1, int max_am2, int max_m):
    max_am1_(max_am1), max_am2_(max_am2), max_m_(max_m)
{
    if (max_am1 < 0)
        throw SanityCheckError("ObaraSaikaTwoCenterMIRecursion -- max_am1 must be nonnegative", __FILE__, __LINE__);
    if (max_am2 < 0)
        throw SanityCheckError("ERROR: ObaraSaikaTwoCenterMIRecursion -- max_am2 must be nonnegative", __FILE__, __LINE__);
    if (max_m > 3)
        throw SanityCheckError("ERROR: ObaraSaikaTwoCenterMIRecursion -- max_m must be nonnegative and less than 4", __FILE__, __LINE__);

    x_ = init_box(max_am1+3, max_am2+3, max_m+1);
    y_ = init_box(max_am1+3, max_am2+3, max_m+1);
    z_ = init_box(max_am1+3, max_am2+3, max_m+1);
}

ObaraSaikaTwoCenterMIRecursion::~ObaraSaikaTwoCenterMIRecursion()
{
    free_box(x_, max_am1_+3, max_am2_+3);
    free_box(y_, max_am1_+3, max_am2_+3);
    free_box(z_, max_am1_+3, max_am2_+3);
}

void ObaraSaikaTwoCenterMIRecursion::compute(double PA[3], double PB[3], double gamma, int am1, int am2)
{
    if (am1 < 0 || am1 > max_am1_)
        throw SanityCheckError("ERROR: ObaraSaikaTwoCenterMIRecursion::compute -- am1 out of bounds", __FILE__, __LINE__);
    if (am2 < 0 || am2 > max_am2_)
        throw SanityCheckError("ERROR: ObaraSaikaTwoCenterMIRecursion::compute -- am2 out of bounds", __FILE__, __LINE__);

    int i, j, k;
    double oog = 1.0 / (2.0 * gamma);

    if (max_m_) {
        x_[0][0][2] = y_[0][0][2] = z_[0][0][2] = oog;
    }

    // Upward recursion in j for i=0
    for (j=0; j<am1; ++j) {
        for (k=0; k<=max_m_; ++k) {
            x_[0][j+1][k] = PB[0] * x_[0][j][k];
            y_[0][j+1][k] = PB[1] * y_[0][j][k];
            z_[0][j+1][k] = PB[2] * z_[0][j][k];

            if (j > 0) {
                x_[0][j+1][k] += j * oog * x_[0][j-1][k];
                y_[0][j+1][k] += j * oog * y_[0][j-1][k];
                z_[0][j+1][k] += j * oog * z_[0][j-1][k];
            }

            if (k > 0) {
                x_[0][j+1][k] += k * oog * x_[0][j][k-1];
                y_[0][j+1][k] += k * oog * y_[0][j][k-1];
                z_[0][j+1][k] += k * oog * z_[0][j][k-1];
            }
        }
    }

    // Upward recursion in i for all j's
    for (i=0; i<am1; ++i) {
        for (j=0; j<=am2; ++j) {
            x_[i+1][j][k] = PA[0] * x_[i][j][k];
            y_[i+1][j][k] = PA[1] * y_[i][j][k];
            z_[i+1][j][k] = PA[2] * z_[i][j][k];

            if (i > 0) {
                x_[i+1][j][k] += i * oog * x_[i-1][j][k];
                y_[i+1][j][k] += i * oog * y_[i-1][j][k];
                z_[i+1][j][k] += i * oog * z_[i-1][j][k];
            }

            if (j > 0) {
                x_[i+1][j][k] += j * oog * x_[i][j-1][k];
                y_[i+1][j][k] += j * oog * y_[i][j-1][k];
                z_[i+1][j][k] += j * oog * z_[i][j-1][k];
            }

            if (k > 0) {
                x_[i+1][j][k] += k * oog * x_[i][j][k-1];
                y_[i+1][j][k] += k * oog * y_[i][j][k-1];
                z_[i+1][j][k] += k * oog * z_[i][j][k-1];
            }
        }
    }
}

ObaraSaikaTwoCenterRecursion::ObaraSaikaTwoCenterRecursion(int max_am1, int max_am2):
    max_am1_(max_am1), max_am2_(max_am2)
{
    if (max_am1 < 0)
        throw SanityCheckError("ERROR: ObaraSaikaTwoCenterRecursion -- max_am1 must be nonnegative", __FILE__, __LINE__);
    if (max_am2 < 0)
        throw SanityCheckError("ERROR: ObaraSaikaTwoCenterRecursion -- max_am2 must be nonnegative", __FILE__, __LINE__);

    int am1 = max_am1;
    int am2 = max_am2;
    if (max_am1_ == 0)
        am1 = 1;
    if (max_am2_ == 0)
        am2 = 1;

    x_ = block_matrix(am1+1, am2+1);
    y_ = block_matrix(am1+1, am2+1);
    z_ = block_matrix(am1+1, am2+1);
}

ObaraSaikaTwoCenterRecursion::~ObaraSaikaTwoCenterRecursion()
{
    free_block(x_);
    free_block(y_);
    free_block(z_);
}

void ObaraSaikaTwoCenterRecursion::compute(double PA[3], double PB[3], double gamma, int am1, int am2)
{
    if (am1 < 0 || am1 > max_am1_)
        throw SanityCheckError("ERROR: ObaraSaikaTwoCenterRecursion::compute -- am1 out of bounds", __FILE__, __LINE__);
    if (am2 < 0 || am2 > max_am2_)
        throw SanityCheckError("ERROR: ObaraSaikaTwoCenterRecursion::compute -- am2 out of bounds", __FILE__, __LINE__);

    int i,j;
    double pp = 1/(2*gamma);
    int lmaxi = am1;
    int lmaxj = am2;

    // Try zeroing out the matrices to fix dipoles
    memset(x_[0], 0, sizeof(double) * (max_am1_+1) * (max_am2_+1));
    memset(y_[0], 0, sizeof(double) * (max_am1_+1) * (max_am2_+1));
    memset(z_[0], 0, sizeof(double) * (max_am1_+1) * (max_am2_+1));

    x_[0][0] = y_[0][0] = z_[0][0] = 1.0;

    /* Upward recursion in j for i=0 */

    x_[0][1] = PB[0];
    y_[0][1] = PB[1];
    z_[0][1] = PB[2];

    for(j=1;j<lmaxj;j++) {
        x_[0][j+1] = PB[0]*x_[0][j];
        y_[0][j+1] = PB[1]*y_[0][j];
        z_[0][j+1] = PB[2]*z_[0][j];
        x_[0][j+1] += j*pp*x_[0][j-1];
        y_[0][j+1] += j*pp*y_[0][j-1];
        z_[0][j+1] += j*pp*z_[0][j-1];
    }

  /* Upward recursion in i for all j's */
    if (lmaxi > 0) {
        x_[1][0] = PA[0];
        y_[1][0] = PA[1];
        z_[1][0] = PA[2];
        for(j=1;j<=lmaxj;j++) {
            x_[1][j] = PA[0]*x_[0][j];
            y_[1][j] = PA[1]*y_[0][j];
            z_[1][j] = PA[2]*z_[0][j];
            x_[1][j] += j*pp*x_[0][j-1];
            y_[1][j] += j*pp*y_[0][j-1];
            z_[1][j] += j*pp*z_[0][j-1];
        }
        for(i=1;i<lmaxi;i++) {
            x_[i+1][0] = PA[0]*x_[i][0];
            y_[i+1][0] = PA[1]*y_[i][0];
            z_[i+1][0] = PA[2]*z_[i][0];
            x_[i+1][0] += i*pp*x_[i-1][0];
            y_[i+1][0] += i*pp*y_[i-1][0];
            z_[i+1][0] += i*pp*z_[i-1][0];
            for(j=1;j<=lmaxj;j++) {
                x_[i+1][j] = PA[0]*x_[i][j];
                y_[i+1][j] = PA[1]*y_[i][j];
                z_[i+1][j] = PA[2]*z_[i][j];
                x_[i+1][j] += i*pp*x_[i-1][j];
                y_[i+1][j] += i*pp*y_[i-1][j];
                z_[i+1][j] += i*pp*z_[i-1][j];
                x_[i+1][j] += j*pp*x_[i][j-1];
                y_[i+1][j] += j*pp*y_[i][j-1];
                z_[i+1][j] += j*pp*z_[i][j-1];
            }
        }
    }
}

ObaraSaikaTwoCenterVIRecursion::ObaraSaikaTwoCenterVIRecursion(int max_am1, int max_am2):
    max_am1_(max_am1), max_am2_(max_am2)
{
    if (max_am1 < 0)
        throw SanityCheckError("ERROR: ObaraSaikaTwoCenterVIRecursion -- max_am1 must be nonnegative", __FILE__, __LINE__);
    if (max_am2 < 0)
        throw SanityCheckError("ERROR: ObaraSaikaTwoCenterVIRecursion -- max_am2 must be nonnegative", __FILE__, __LINE__);

    size_ = max_am1 > max_am2 ? max_am1 : max_am2;
    size_ += 1;
    size_ = (size_-1)*size_*(size_+1)+1;
    vi_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
}

ObaraSaikaTwoCenterVIRecursion::~ObaraSaikaTwoCenterVIRecursion()
{
    free_box(vi_, size_, size_);
}

#define EPS 1.0e-17

void ObaraSaikaTwoCenterVIRecursion::calculate_f(double *F, int n, double t)
{
  int i, m, k;
  int m2;
  double t2;
  double num;
  double sum;
  double term1, term2;
  static double K = 1.0/M_2_SQRTPI;
  double et;


  if (t>20.0){
    t2 = 2*t;
    et = exp(-t);
    t = sqrt(t);
    F[0] = K*erf(t)/t;
    for(m=0; m<=n-1; m++){
      F[m+1] = ((2*m + 1)*F[m] - et)/(t2);
      }
    }
  else {
    et = exp(-t);
    t2 = 2*t;
    m2 = 2*n;
    num = df[m2];
    i=0;
    sum = 1.0/(m2+1);
    do{
      i++;
      num = num*t2;
      term1 = num/df[m2+2*i+2];
      sum += term1;
      } while (fabs(term1) > EPS && i < MAX_FAC);
    F[n] = sum*et;
    for(m=n-1;m>=0;m--){
      F[m] = (t2*F[m+1] + et)/(2*m+1);
      }
    }
}

void ObaraSaikaTwoCenterVIRecursion::compute(double PA[3], double PB[3], double PC[3], double zeta, int am1, int am2)
{
    int a, b, m;
    int azm = 1;
    int aym = am1 + 1;
    int axm = aym * aym;
    int bzm = 1;
    int bym = am2 + 1;
    int bxm = bym * bym;
    int ax, ay, az, bx, by, bz;
    int aind, bind;
    double ooz = 1.0/(2.0 * zeta);
    int mmax = max_am1_ + max_am2_;

    // Prefactor from A20
    double tmp = sqrt(zeta) * M_2_SQRTPI;
    // U from A21
    double u = zeta * (PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2]);
    double *F = new double[mmax+1];

    // Form Fm(U) from A20
    calculate_f(F, mmax, u);

    // Perform recursion in m for (a|A(0)|s) using A20
    for (m=0; m<=mmax; ++m) {
        vi_[0][0][m] = tmp * F[m];
    }

    // Perform recursion in b with a=0
    //  subset of A19
    for (b=1; b<=am2; ++b) {
        for (bx=0; bx<=b; ++bx) {
            for (by=0; by<=b-bx; ++by) {
                bz = b-bx-by;

                // Compute the index into VI for bx,by,bz
                bind = bx*bxm + by*bym + bz*bzm;

                // Compute each x, y, z contribution
                if (bz > 0) {
                    for (m=0; m<=mmax-b; ++m) {
                        vi_[0][bind][m] = PB[2] * vi_[0][bind-bzm][m] - PC[2] * vi_[0][bind-bzm][m+1];
                    }
                    if (bz > 1) {
                        for (m=0; m<=mmax-b; ++m) {
                            vi_[0][bind][m] += ooz * (bz-1) * (vi_[0][bind-2*bzm][m] - vi_[0][bind-2*bzm][m+1]);
                        }
                    }
                }
                else if (by > 0) {
                    for (m=0; m<=mmax-b; ++m) {
                        vi_[0][bind][m] = PB[1] * vi_[0][bind-bym][m] - PC[1] * vi_[0][bind-bym][m+1];
                    }
                    if (by > 1) {
                        for (m=0; m<=mmax-b; ++m) {
                            vi_[0][bind][m] += ooz * (by-1) * (vi_[0][bind-2*bym][m] - vi_[0][bind-2*bym][m+1]);
                        }
                    }
                }
                else if (bx > 0) {
                    for (m=0; m<=mmax-b; ++m) {
                        vi_[0][bind][m] = PB[0] * vi_[0][bind-bxm][m] - PC[0] * vi_[0][bind-bxm][m+1];
                    }
                    if (bx > 1) {
                        for (m=0; m<=mmax-b; ++m) {
                            vi_[0][bind][m] += ooz * (bx-1) * (vi_[0][bind-2*bxm][m] - vi_[0][bind-2*bxm][m+1]);
                        }
                    }
                }
            }
        }
    }

    // Perform upward recursion in a with all b's
    for (b=0; b<=am2; b++) {
        for (bx=0; bx<=b; bx++) {
            for (by=0; by<=b-bx;by++) {
                bz = b-bx-by;
                bind = bx*bxm + by*bym + bz*bzm;

                for (a=1; a<=am1; a++) {
                    // This next for loop was for (ax=0; ax<=b; ax++)
                    // this could explain why dx2 was not being computed.
                    // change for for(ax=0; ax<a; ax++) on 2005-09-15 4:11pm
                    for (ax=0; ax<=a; ax++) {
                        for (ay=0; ay<=a-ax; ay++) {
                            az = a-ax-ay;
                            aind = ax*axm + ay*aym + az*azm;

                            if (az > 0) {
                                for (m=0; m<=mmax-a-b; m++) {
                                    vi_[aind][bind][m] = PA[2] * vi_[aind-azm][bind][m] - PC[2] * vi_[aind-azm][bind][m+1];
                                }

                                if (az > 1) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        vi_[aind][bind][m] += ooz * (az-1) * (vi_[aind-2*azm][bind][m] - vi_[aind-2*azm][bind][m+1]);
                                    }
                                }
                                if (bz > 0) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        vi_[aind][bind][m] += ooz * bz * (vi_[aind-azm][bind-bzm][m] - vi_[aind-azm][bind-bzm][m+1]);
                                    }
                                }
                            }
                            else if (ay > 0) {
                                for (m=0; m<=mmax-a-b; m++) {
                                    vi_[aind][bind][m] = PA[1] * vi_[aind-aym][bind][m] - PC[1] * vi_[aind-aym][bind][m+1];
                                }
                                if (ay > 1) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        vi_[aind][bind][m] += ooz * (ay-1) * (vi_[aind-2*aym][bind][m] - vi_[aind-2*aym][bind][m+1]);
                                    }
                                }
                                if (by > 0) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        vi_[aind][bind][m] += ooz * by * (vi_[aind-aym][bind-bym][m] - vi_[aind-aym][bind-bym][m+1]);
                                    }
                                }
                            }
                            else if (ax > 0) {
                                for (m=0; m<=mmax-a-b; m++) {
                                    vi_[aind][bind][m] = PA[0] * vi_[aind-axm][bind][m] - PC[0] * vi_[aind-axm][bind][m+1];
                                }

                                if (ax > 1) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        vi_[aind][bind][m] += ooz * (ax-1) * (vi_[aind-2*axm][bind][m] - vi_[aind-2*axm][bind][m+1]);
                                    }
                                }
                                if (bx > 0) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        vi_[aind][bind][m] += ooz * bx * (vi_[aind-axm][bind-bxm][m] - vi_[aind-axm][bind-bxm][m+1]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    delete[] F;
}

ObaraSaikaTwoCenterVIDerivRecursion::ObaraSaikaTwoCenterVIDerivRecursion(int max_am1, int max_am2)
    : ObaraSaikaTwoCenterVIRecursion(max_am1+1, max_am2+1)
{
    vx_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    vy_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    vz_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
}

ObaraSaikaTwoCenterVIDerivRecursion::~ObaraSaikaTwoCenterVIDerivRecursion()
{
    free_box(vx_, size_, size_);
    free_box(vy_, size_, size_);
    free_box(vz_, size_, size_);
}

void ObaraSaikaTwoCenterVIDerivRecursion::compute(double PA[3], double PB[3], double PC[3], double zeta, int am1, int am2)
{
    int a, b, m;
    int azm = 1;
    int aym = am1 + 1;
    int axm = aym * aym;
    int bzm = 1;
    int bym = am2 + 1;
    int bxm = bym * bym;
    int ax, ay, az, bx, by, bz;
    int aind, bind;
    double ooz = 1.0/(2.0 * zeta);
    int mmax = am1 + am2;

    // Prefactor from A20
    double tmp = sqrt(zeta) * M_2_SQRTPI;
    // U from A21
    double u = zeta * (PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2]);
    double *F = new double[mmax+1];

    // Zero out F
    memset(F, 0, sizeof(double) * (mmax+1));

    // Form Fm(U) from A20
    calculate_f(F, mmax, u);

    // Perform recursion in m for (a|A(0)|s) using A20
    for (m=0; m<=mmax; ++m) {
        vi_[0][0][m] = tmp * F[m];
    }
    for (m=0; m<=mmax-1; ++m) {
        vx_[0][0][m] = 2.0*zeta*PC[0]*vi_[0][0][m+1];
        vy_[0][0][m] = 2.0*zeta*PC[1]*vi_[0][0][m+1];
        vz_[0][0][m] = 2.0*zeta*PC[2]*vi_[0][0][m+1];
    }

    // Perform recursion in b with a=0
    //  subset of A19
    for (b=1; b<=am2; ++b) {
        for (bx=0; bx<=b; ++bx) {
            for (by=0; by<=b-bx; ++by) {
                bz = b-bx-by;

                // Compute the index into VI for bx,by,bz
                bind = bx*bxm + by*bym + bz*bzm;

                // Compute each x, y, z contribution
                if (bz > 0) {
                    for (m=0; m<=mmax-b; ++m) {
                        vi_[0][bind][m] = PB[2] * vi_[0][bind-bzm][m] - PC[2] * vi_[0][bind-bzm][m+1];
                    }
                    for (m=0; m<=mmax-b-1; ++m) {
                        vx_[0][bind][m] = PB[2] * vx_[0][bind-bzm][m] - PC[2] * vx_[0][bind-bzm][m+1];
                        vy_[0][bind][m] = PB[2] * vy_[0][bind-bzm][m] - PC[2] * vy_[0][bind-bzm][m+1];
                        vz_[0][bind][m] = PB[2] * vz_[0][bind-bzm][m] - PC[2] * vz_[0][bind-bzm][m+1] + vi_[0][bind-bzm][m+1];
                    }
                    if (bz > 1) {
                        for (m=0; m<=mmax-b; ++m) {
                            vi_[0][bind][m] += ooz * (bz-1) * (vi_[0][bind-2*bzm][m] - vi_[0][bind-2*bzm][m+1]);
                        }
                        for (m=0; m<=mmax-b-1; ++m) {
                            vx_[0][bind][m] += ooz * (bz-1) * (vx_[0][bind-2*bzm][m] - vx_[0][bind-2*bzm][m+1]);
                            vy_[0][bind][m] += ooz * (bz-1) * (vy_[0][bind-2*bzm][m] - vy_[0][bind-2*bzm][m+1]);
                            vz_[0][bind][m] += ooz * (bz-1) * (vz_[0][bind-2*bzm][m] - vz_[0][bind-2*bzm][m+1]);
                        }
                    }
                }
                else if (by > 0) {
                    for (m=0; m<=mmax-b; ++m) {
                        vi_[0][bind][m] = PB[1] * vi_[0][bind-bym][m] - PC[1] * vi_[0][bind-bym][m+1];
                    }
                    for (m=0; m<=mmax-b-1; ++m) {
                        vx_[0][bind][m] = PB[1] * vx_[0][bind-bym][m] - PC[1] * vx_[0][bind-bym][m+1];
                        vy_[0][bind][m] = PB[1] * vy_[0][bind-bym][m] - PC[1] * vy_[0][bind-bym][m+1] + vi_[0][bind-bym][m+1];
                        vz_[0][bind][m] = PB[1] * vz_[0][bind-bym][m] - PC[1] * vz_[0][bind-bym][m+1];
                    }
                    if (by > 1) {
                        for (m=0; m<=mmax-b; ++m) {
                            vi_[0][bind][m] += ooz * (by-1) * (vi_[0][bind-2*bym][m] - vi_[0][bind-2*bym][m+1]);
                        }
                        for (m=0; m<=mmax-b-1; ++m) {
                            vx_[0][bind][m] += ooz * (by-1) * (vx_[0][bind-2*bym][m] - vx_[0][bind-2*bym][m+1]);
                            vy_[0][bind][m] += ooz * (by-1) * (vy_[0][bind-2*bym][m] - vy_[0][bind-2*bym][m+1]);
                            vz_[0][bind][m] += ooz * (by-1) * (vz_[0][bind-2*bym][m] - vz_[0][bind-2*bym][m+1]);
                        }
                    }
                }
                else if (bx > 0) {
                    for (m=0; m<=mmax-b; ++m) {
                        vi_[0][bind][m] = PB[0] * vi_[0][bind-bxm][m] - PC[0] * vi_[0][bind-bxm][m+1];
                    }
                    for (m=0; m<=mmax-b-1; ++m) {
                        vx_[0][bind][m] = PB[0] * vx_[0][bind-bxm][m] - PC[0] * vx_[0][bind-bxm][m+1] + vi_[0][bind-bxm][m+1];
                        vy_[0][bind][m] = PB[0] * vy_[0][bind-bxm][m] - PC[0] * vy_[0][bind-bxm][m+1];
                        vz_[0][bind][m] = PB[0] * vz_[0][bind-bxm][m] - PC[0] * vz_[0][bind-bxm][m+1];
                    }
                    if (bx > 1) {
                        for (m=0; m<=mmax-b; ++m) {
                            vi_[0][bind][m] += ooz * (bx-1) * (vi_[0][bind-2*bxm][m] - vi_[0][bind-2*bxm][m+1]);
                        }
                        for (m=0; m<=mmax-b-1; ++m) {
                            vx_[0][bind][m] += ooz * (bx-1) * (vx_[0][bind-2*bxm][m] - vx_[0][bind-2*bxm][m+1]);
                            vy_[0][bind][m] += ooz * (bx-1) * (vy_[0][bind-2*bxm][m] - vy_[0][bind-2*bxm][m+1]);
                            vz_[0][bind][m] += ooz * (bx-1) * (vz_[0][bind-2*bxm][m] - vz_[0][bind-2*bxm][m+1]);
                        }
                    }
                }
            }
        }
    }

    // Perform upward recursion in a with all b's
    for (b=0; b<=am2; b++) {
        for (bx=0; bx<=b; bx++) {
            for (by=0; by<=b-bx;by++) {
                bz = b-bx-by;
                bind = bx*bxm + by*bym + bz*bzm;

                for (a=1; a<=am1; a++) {
                    // This next for loop was for (ax=0; ax<=b; ax++)
                    // this could explain why dx2 was not being computed.
                    // change for for(ax=0; ax<a; ax++) on 2005-09-15 4:11pm
                    for (ax=0; ax<=a; ax++) {
                        for (ay=0; ay<=a-ax; ay++) {
                            az = a-ax-ay;
                            aind = ax*axm + ay*aym + az*azm;

                            if (az > 0) {
                                for (m=0; m<=mmax-a-b; m++) {
                                    vi_[aind][bind][m] = PA[2] * vi_[aind-azm][bind][m] - PC[2] * vi_[aind-azm][bind][m+1];
                                }
                                for (m=0; m<=mmax-a-b-1; ++m) {
                                    vx_[aind][bind][m] = PA[2] * vx_[aind-azm][bind][m] - PC[2] * vx_[aind-azm][bind][m+1];
                                    vy_[aind][bind][m] = PA[2] * vy_[aind-azm][bind][m] - PC[2] * vy_[aind-azm][bind][m+1];
                                    vz_[aind][bind][m] = PA[2] * vz_[aind-azm][bind][m] - PC[2] * vz_[aind-azm][bind][m+1] + vi_[aind-azm][bind][m+1];
                                }

                                if (az > 1) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        vi_[aind][bind][m] += ooz * (az-1) * (vi_[aind-2*azm][bind][m] - vi_[aind-2*azm][bind][m+1]);
                                    }
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        vx_[aind][bind][m] += ooz * (az-1) * (vx_[aind-2*azm][bind][m] - vx_[aind-2*azm][bind][m+1]);
                                        vy_[aind][bind][m] += ooz * (az-1) * (vy_[aind-2*azm][bind][m] - vy_[aind-2*azm][bind][m+1]);
                                        vz_[aind][bind][m] += ooz * (az-1) * (vz_[aind-2*azm][bind][m] - vz_[aind-2*azm][bind][m+1]);
                                    }
                                }
                                if (bz > 0) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        vi_[aind][bind][m] += ooz * bz * (vi_[aind-azm][bind-bzm][m] - vi_[aind-azm][bind-bzm][m+1]);
                                    }
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        vx_[aind][bind][m] += ooz * bz * (vx_[aind-azm][bind-bzm][m] - vx_[aind-azm][bind-bzm][m+1]);
                                        vy_[aind][bind][m] += ooz * bz * (vy_[aind-azm][bind-bzm][m] - vy_[aind-azm][bind-bzm][m+1]);
                                        vz_[aind][bind][m] += ooz * bz * (vz_[aind-azm][bind-bzm][m] - vz_[aind-azm][bind-bzm][m+1]);
                                    }
                                }
                            }
                            else if (ay > 0) {
                                for (m=0; m<=mmax-a-b; m++) {
                                    vi_[aind][bind][m] = PA[1] * vi_[aind-aym][bind][m] - PC[1] * vi_[aind-aym][bind][m+1];
                                }
                                for (m=0; m<=mmax-a-b-1; ++m) {
                                    vx_[aind][bind][m] = PA[1] * vx_[aind-aym][bind][m] - PC[1] * vx_[aind-aym][bind][m+1];
                                    vy_[aind][bind][m] = PA[1] * vy_[aind-aym][bind][m] - PC[1] * vy_[aind-aym][bind][m+1] + vi_[aind-aym][bind][m+1];
                                    vz_[aind][bind][m] = PA[1] * vz_[aind-aym][bind][m] - PC[1] * vz_[aind-aym][bind][m+1];
                                }
                                if (ay > 1) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        vi_[aind][bind][m] += ooz * (ay-1) * (vi_[aind-2*aym][bind][m] - vi_[aind-2*aym][bind][m+1]);
                                    }
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        vx_[aind][bind][m] += ooz * (ay-1) * (vx_[aind-2*aym][bind][m] - vx_[aind-2*aym][bind][m+1]);
                                        vy_[aind][bind][m] += ooz * (ay-1) * (vy_[aind-2*aym][bind][m] - vy_[aind-2*aym][bind][m+1]);
                                        vz_[aind][bind][m] += ooz * (ay-1) * (vz_[aind-2*aym][bind][m] - vz_[aind-2*aym][bind][m+1]);
                                    }
                                }
                                if (by > 0) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        vi_[aind][bind][m] += ooz * by * (vi_[aind-aym][bind-bym][m] - vi_[aind-aym][bind-bym][m+1]);
                                    }
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        vx_[aind][bind][m] += ooz * by * (vx_[aind-aym][bind-bym][m] - vx_[aind-aym][bind-bym][m+1]);
                                        vy_[aind][bind][m] += ooz * by * (vy_[aind-aym][bind-bym][m] - vy_[aind-aym][bind-bym][m+1]);
                                        vz_[aind][bind][m] += ooz * by * (vz_[aind-aym][bind-bym][m] - vz_[aind-aym][bind-bym][m+1]);
                                    }
                                }
                            }
                            else if (ax > 0) {
                                for (m=0; m<=mmax-a-b; m++) {
                                    vi_[aind][bind][m] = PA[0] * vi_[aind-axm][bind][m] - PC[0] * vi_[aind-axm][bind][m+1];
                                }
                                for (m=0; m<=mmax-a-b-1; ++m) {
                                    vx_[aind][bind][m] = PA[0] * vx_[aind-axm][bind][m] - PC[0] * vx_[aind-axm][bind][m+1] + vi_[aind-axm][bind][m+1];
                                    vy_[aind][bind][m] = PA[0] * vy_[aind-axm][bind][m] - PC[0] * vy_[aind-axm][bind][m+1];
                                    vz_[aind][bind][m] = PA[0] * vz_[aind-axm][bind][m] - PC[0] * vz_[aind-axm][bind][m+1];
                                }

                                if (ax > 1) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        vi_[aind][bind][m] += ooz * (ax-1) * (vi_[aind-2*axm][bind][m] - vi_[aind-2*axm][bind][m+1]);
                                    }
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        vx_[aind][bind][m] += ooz * (ax-1) * (vx_[aind-2*axm][bind][m] - vx_[aind-2*axm][bind][m+1]);
                                        vy_[aind][bind][m] += ooz * (ax-1) * (vy_[aind-2*axm][bind][m] - vy_[aind-2*axm][bind][m+1]);
                                        vz_[aind][bind][m] += ooz * (ax-1) * (vz_[aind-2*axm][bind][m] - vz_[aind-2*axm][bind][m+1]);
                                    }
                                }
                                if (bx > 0) {
                                    for (m=0; m<= mmax-a-b; m++) {
                                        vi_[aind][bind][m] += ooz * bx * (vi_[aind-axm][bind-bxm][m] - vi_[aind-axm][bind-bxm][m+1]);
                                    }
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        vx_[aind][bind][m] += ooz * bx * (vx_[aind-axm][bind-bxm][m] - vx_[aind-axm][bind-bxm][m+1]);
                                        vy_[aind][bind][m] += ooz * bx * (vy_[aind-axm][bind-bxm][m] - vy_[aind-axm][bind-bxm][m+1]);
                                        vz_[aind][bind][m] += ooz * bx * (vz_[aind-axm][bind-bxm][m] - vz_[aind-axm][bind-bxm][m+1]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    delete[] F;
}

ObaraSaikaTwoCenterElectricField::ObaraSaikaTwoCenterElectricField(int max_am1, int max_am2)
    : ObaraSaikaTwoCenterVIRecursion(max_am1+1, max_am2+1)
{
    ex_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    ey_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    ez_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
}

ObaraSaikaTwoCenterElectricField::~ObaraSaikaTwoCenterElectricField()
{
    free_box(ex_, size_, size_);
    free_box(ey_, size_, size_);
    free_box(ez_, size_, size_);
}

void ObaraSaikaTwoCenterElectricField::compute(double PA[3], double PB[3], double PC[3], double zeta, int am1, int am2)
{
    int a, b, m;
    int azm = 1;
    int aym = am1 + 1;
    int axm = aym * aym;
    int bzm = 1;
    int bym = am2 + 1;
    int bxm = bym * bym;
    int ax, ay, az, bx, by, bz;
    int aind, bind;
    double ooz = 1.0/(2.0 * zeta);
    int mmax = am1 + am2;

    // Call our super class to compute potential integrals
    ObaraSaikaTwoCenterVIRecursion::compute(PA, PB, PC, zeta, am1+1, am2+1);

    // Compute starting integrals using A25
    fprintf(outfile, "setup:\n");
    for (m=0; m<=mmax; ++m) {
        ex_[0][0][m] = 2.0 * zeta * PC[0] * vi_[0][0][m+1];
        ey_[0][0][m] = 2.0 * zeta * PC[1] * vi_[0][0][m+1];
        ez_[0][0][m] = 2.0 * zeta * PC[2] * vi_[0][0][m+1];

        fprintf(outfile, "1 ex_[0][0][%d] = %lf\tey_[0][0][%d] = %lf\tez_[0][0][%d] = %lf\n", m, ex_[0][0][m], m, ey_[0][0][m], m, ez_[0][0][m]);
    }

    // Perform recursion in b with a=0
    //    subset of A24
    fprintf(outfile, "upward in j with i=0:\n");
    for (b=1; b<=am2; ++b) {
        for (bx=0; bx<=b; ++bx) {
            for (by=0; by<=b-bx; ++by) {
                bz = b-bx-by;

                // Compute the index into ex_, ey_, ez_ for bx, by, bz
                bind = bx*bxm + by*bym + bz*bzm;

                // Compute each x, y, z contribution
                if (bz > 0) {
                    for (m=0; m<=am2-b; ++m) {
                        ex_[0][bind][m] = PB[2] * ex_[0][bind-bzm][m] - PC[2] * ex_[0][bind-bzm][m+1];
                        ey_[0][bind][m] = PB[2] * ey_[0][bind-bzm][m] - PC[2] * ey_[0][bind-bzm][m+1];
                        ez_[0][bind][m] = PB[2] * ez_[0][bind-bzm][m] - PC[2] * ez_[0][bind-bzm][m+1] + vi_[0][bind-bzm][m+1];

                        fprintf(outfile, "2 ex_[0][%d][%d] = %lf\tey_[0][%d][%d] = %lf\tez_[0][%d][%d] = %lf\n", bind, m, ex_[0][bind][m], bind, m, ey_[0][bind][m], bind, m, ez_[0][bind][m]);
                    }
                    if (bz > 1) {
                        for (m=0; m<=am2-b; ++m) {
                            ex_[0][bind][m] += ooz * (bz-1) * (ex_[0][bind-2*bzm][m] - ex_[0][bind-2*bzm][m+1]);
                            ey_[0][bind][m] += ooz * (bz-1) * (ey_[0][bind-2*bzm][m] - ey_[0][bind-2*bzm][m+1]);
                            ez_[0][bind][m] += ooz * (bz-1) * (ez_[0][bind-2*bzm][m] - ez_[0][bind-2*bzm][m+1]);

                            fprintf(outfile, "3 ex_[0][%d][%d] = %lf\tey_[0][%d][%d] = %lf\tez_[0][%d][%d] = %lf\n", bind, m, ex_[0][bind][m], bind, m, ey_[0][bind][m], bind, m, ez_[0][bind][m]);
                        }
                    }
                }
                else if (by > 0) {
                    for (m=0; m<=am2-b; ++m) {
                        ex_[0][bind][m] = PB[1] * ex_[0][bind-bym][m] - PC[1] * ex_[0][bind-bym][m+1];
                        ey_[0][bind][m] = PB[1] * ey_[0][bind-bym][m] - PC[1] * ey_[0][bind-bym][m+1] + vi_[0][bind-bym][m+1];
                        ez_[0][bind][m] = PB[1] * ez_[0][bind-bym][m] - PC[1] * ez_[0][bind-bym][m+1];

                        fprintf(outfile, "4 ex_[0][%d][%d] = %lf\tey_[0][%d][%d] = %lf\tez_[0][%d][%d] = %lf\n", bind, m, ex_[0][bind][m], bind, m, ey_[0][bind][m], bind, m, ez_[0][bind][m]);
                    }
                    if (by > 1) {
                        for (m=0; m<=am2-b; ++m) {
                            ex_[0][bind][m] += ooz * (by-1) * (ex_[0][bind-2*bym][m] - ex_[0][bind-2*bym][m+1]);
                            ey_[0][bind][m] += ooz * (by-1) * (ey_[0][bind-2*bym][m] - ey_[0][bind-2*bym][m+1]);
                            ez_[0][bind][m] += ooz * (by-1) * (ez_[0][bind-2*bym][m] - ez_[0][bind-2*bym][m+1]);

                            fprintf(outfile, "5 ex_[0][%d][%d] = %lf\tey_[0][%d][%d] = %lf\tez_[0][%d][%d] = %lf\n", bind, m, ex_[0][bind][m], bind, m, ey_[0][bind][m], bind, m, ez_[0][bind][m]);
                        }
                    }
                }
                else if (bx > 0) {
                    for (m=0; m<=am2-b; ++m) {
                        ex_[0][bind][m] = PB[0] * ex_[0][bind-bxm][m] - PC[0] * ex_[0][bind-bxm][m+1] + vi_[0][bind-bxm][m+1];
                        ey_[0][bind][m] = PB[0] * ey_[0][bind-bxm][m] - PC[0] * ey_[0][bind-bxm][m+1];
                        ez_[0][bind][m] = PB[0] * ez_[0][bind-bxm][m] - PC[0] * ez_[0][bind-bxm][m+1];

                        fprintf(outfile, "6 ex_[0][%d][%d] = %lf\tey_[0][%d][%d] = %lf\tez_[0][%d][%d] = %lf\n", bind, m, ex_[0][bind][m], bind, m, ey_[0][bind][m], bind, m, ez_[0][bind][m]);
                    }
                    if (bx > 1) {
                        for (m=0; m<=am2-b; ++m) {
                            ex_[0][bind][m] += ooz * (bx-1) * (ex_[0][bind-2*bxm][m] - ex_[0][bind-2*bxm][m+1]);
                            ey_[0][bind][m] += ooz * (bx-1) * (ey_[0][bind-2*bxm][m] - ey_[0][bind-2*bxm][m+1]);
                            ez_[0][bind][m] += ooz * (bx-1) * (ez_[0][bind-2*bxm][m] - ez_[0][bind-2*bxm][m+1]);

                            fprintf(outfile, "7 ex_[0][%d][%d] = %lf\tey_[0][%d][%d] = %lf\tez_[0][%d][%d] = %lf\n", bind, m, ex_[0][bind][m], bind, m, ey_[0][bind][m], bind, m, ez_[0][bind][m]);
                        }
                    }
                }
            }
        }
    }

    // Perform upward recursion in a with all b's
    for (b=0; b<=am2; ++b) {
        for (bx=0; bx<=b; bx++) {
            for (by=0; by<=b-bx;by++) {
                bz = b-bx-by;
                bind = bx*bxm + by*bym + bz*bzm;

                for (a=1; a<=am1; a++) {
                    for (ax=0; ax<=a; ax++) {
                        for (ay=0; ay<=a-ax; ay++) {
                            az = a-ax-ay;
                            aind = ax*axm + ay*aym + az*azm;

                            if (az > 0) {
                                for (m=0; m<=mmax-a-b; ++m) {
//                                    AIX[iind][jind][m] = paz*AIX[iind-izm][jind][m] - pcz*AIX[iind-izm][jind][m+1];
//                                    AIY[iind][jind][m] = paz*AIY[iind-izm][jind][m] - pcz*AIY[iind-izm][jind][m+1];
//                                    AIZ[iind][jind][m] = paz*AIZ[iind-izm][jind][m] - pcz*AIZ[iind-izm][jind][m+1] + AI0[iind-izm][jind][m+1];

                                    ex_[aind][bind][m] = PA[2] * ex_[aind-azm][bind][m] - PC[2] * ex_[aind-azm][bind][m+1];
                                    ey_[aind][bind][m] = PA[2] * ey_[aind-azm][bind][m] - PC[2] * ey_[aind-azm][bind][m+1];
                                    ez_[aind][bind][m] = PA[2] * ez_[aind-azm][bind][m] - PC[2] * ez_[aind-azm][bind][m+1] + vi_[aind-azm][bind][m+1];

                                    fprintf(outfile, "8 ex_[%d][%d][%d] = %lf\tey_[%d][%d][%d] = %lf\tez_[%d][%d][%d] = %lf\tvi_ = %lf %d %d %d\n", aind, bind, m, ex_[aind][bind][m], aind, bind, m, ey_[aind][bind][m], aind, bind, m, ez_[aind][bind][m], vi_[aind-azm][bind][m+1], aind-azm, bind, m+1);
                                }
                                if (az > 1) {
                                    for (m=0; m<=mmax-a-b; ++m) {
                                        ex_[aind][bind][m] += ooz * (az-1) * (ex_[aind-2*azm][bind][m] - ex_[aind-2*azm][bind][m+1]);
                                        ey_[aind][bind][m] += ooz * (az-1) * (ey_[aind-2*azm][bind][m] - ey_[aind-2*azm][bind][m+1]);
                                        ez_[aind][bind][m] += ooz * (az-1) * (ez_[aind-2*azm][bind][m] - ez_[aind-2*azm][bind][m+1]);
                                        fprintf(outfile, "9 ex_[%d][%d][%d] = %lf\tey_[%d][%d][%d] = %lf\tez_[%d][%d][%d] = %lf\n", aind, bind, m, ex_[aind][bind][m], aind, bind, m, ey_[aind][bind][m], aind, bind, m, ez_[aind][bind][m]);
                                    }
                                }
                                if (bz > 0) {
                                    for (m=0; m<=mmax-a-b; ++m) {
                                        ex_[aind][bind][m] += ooz * bz * (ex_[aind-azm][bind-bzm][m] - ex_[aind-azm][bind-bzm][m+1]);
                                        ey_[aind][bind][m] += ooz * bz * (ey_[aind-azm][bind-bzm][m] - ey_[aind-azm][bind-bzm][m+1]);
                                        ez_[aind][bind][m] += ooz * bz * (ez_[aind-azm][bind-bzm][m] - ez_[aind-azm][bind-bzm][m+1]);
                                        fprintf(outfile, "10 ex_[%d][%d][%d] = %lf\tey_[%d][%d][%d] = %lf\tez_[%d][%d][%d] = %lf\n", aind, bind, m, ex_[aind][bind][m], aind, bind, m, ey_[aind][bind][m], aind, bind, m, ez_[aind][bind][m]);
                                    }
                                }
                            }
                            else if (ay > 0) {
                                for (m=0; m<=mmax-a-b; ++m) {
                                    ex_[aind][bind][m] = PA[1] * ex_[aind-aym][bind][m] - PC[1] * ex_[aind-aym][bind][m+1];
                                    ey_[aind][bind][m] = PA[1] * ey_[aind-aym][bind][m] - PC[1] * ey_[aind-aym][bind][m+1] + vi_[aind-aym][bind][m+1];
                                    ez_[aind][bind][m] = PA[1] * ez_[aind-aym][bind][m] - PC[1] * ez_[aind-aym][bind][m+1];
                                    fprintf(outfile, "11 ex_[%d][%d][%d] = %lf\tey_[%d][%d][%d] = %lf\tez_[%d][%d][%d] = %lf\n", aind, bind, m, ex_[aind][bind][m], aind, bind, m, ey_[aind][bind][m], aind, bind, m, ez_[aind][bind][m]);
                                }
                                if (ay > 1) {
                                    for (m=0; m<=mmax-a-b; ++m) {
                                        ex_[aind][bind][m] += ooz * (ay-1) * (ex_[aind-2*aym][bind][m] - ex_[aind-2*aym][bind][m+1]);
                                        ey_[aind][bind][m] += ooz * (ay-1) * (ey_[aind-2*aym][bind][m] - ey_[aind-2*aym][bind][m+1]);
                                        ez_[aind][bind][m] += ooz * (ay-1) * (ez_[aind-2*aym][bind][m] - ez_[aind-2*aym][bind][m+1]);
                                        fprintf(outfile, "12 ex_[%d][%d][%d] = %lf\tey_[%d][%d][%d] = %lf\tez_[%d][%d][%d] = %lf\n", aind, bind, m, ex_[aind][bind][m], aind, bind, m, ey_[aind][bind][m], aind, bind, m, ez_[aind][bind][m]);
                                    }
                                }
                                if (by > 0) {
                                    for (m=0; m<=mmax-a-b; ++m) {
                                        ex_[aind][bind][m] += ooz * by * (ex_[aind-aym][bind-bym][m] - ex_[aind-aym][bind-bym][m+1]);
                                        ey_[aind][bind][m] += ooz * by * (ey_[aind-aym][bind-bym][m] - ey_[aind-aym][bind-bym][m+1]);
                                        ez_[aind][bind][m] += ooz * by * (ez_[aind-aym][bind-bym][m] - ez_[aind-aym][bind-bym][m+1]);
                                        fprintf(outfile, "13 ex_[%d][%d][%d] = %lf\tey_[%d][%d][%d] = %lf\tez_[%d][%d][%d] = %lf\n", aind, bind, m, ex_[aind][bind][m], aind, bind, m, ey_[aind][bind][m], aind, bind, m, ez_[aind][bind][m]);
                                    }
                                }
                            }
                            else if (ax > 0) {
                                for (m=0; m<=mmax-a-b; ++m) {
                                    ex_[aind][bind][m] = PA[0] * ex_[aind-axm][bind][m] - PC[0] * ex_[aind-axm][bind][m+1] + vi_[aind-axm][bind][m+1];
                                    ey_[aind][bind][m] = PA[0] * ey_[aind-axm][bind][m] - PC[0] * ey_[aind-axm][bind][m+1];
                                    ez_[aind][bind][m] = PA[0] * ez_[aind-axm][bind][m] - PC[0] * ez_[aind-axm][bind][m+1];
                                    fprintf(outfile, "14 ex_[%d][%d][%d] = %lf\tey_[%d][%d][%d] = %lf\tez_[%d][%d][%d] = %lf\n", aind, bind, m, ex_[aind][bind][m], aind, bind, m, ey_[aind][bind][m], aind, bind, m, ez_[aind][bind][m]);
                                }
                                if (ax > 1) {
                                    for (m=0; m<=mmax-a-b; ++m) {
                                        ex_[aind][bind][m] += ooz * (ax-1) * (ex_[aind-2*axm][bind][m] - ex_[aind-2*axm][bind][m+1]);
                                        ey_[aind][bind][m] += ooz * (ax-1) * (ey_[aind-2*axm][bind][m] - ey_[aind-2*axm][bind][m+1]);
                                        ez_[aind][bind][m] += ooz * (ax-1) * (ez_[aind-2*axm][bind][m] - ez_[aind-2*axm][bind][m+1]);
                                        fprintf(outfile, "15 ex_[%d][%d][%d] = %lf\tey_[%d][%d][%d] = %lf\tez_[%d][%d][%d] = %lf\n", aind, bind, m, ex_[aind][bind][m], aind, bind, m, ey_[aind][bind][m], aind, bind, m, ez_[aind][bind][m]);
                                    }
                                }
                                if (bx > 0) {
                                    for (m=0; m<=mmax-a-b; ++m) {
                                        ex_[aind][bind][m] += ooz * bx * (ex_[aind-axm][bind-bxm][m] - ex_[aind-axm][bind-bxm][m+1]);
                                        ey_[aind][bind][m] += ooz * bx * (ey_[aind-axm][bind-bxm][m] - ey_[aind-axm][bind-bxm][m+1]);
                                        ez_[aind][bind][m] += ooz * bx * (ez_[aind-axm][bind-bxm][m] - ez_[aind-axm][bind-bxm][m+1]);
                                        fprintf(outfile, "16 ex_[%d][%d][%d] = %lf\tey_[%d][%d][%d] = %lf\tez_[%d][%d][%d] = %lf\n", aind, bind, m, ex_[aind][bind][m], aind, bind, m, ey_[aind][bind][m], aind, bind, m, ez_[aind][bind][m]);
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

ObaraSaikaTwoCenterElectricFieldGradient::ObaraSaikaTwoCenterElectricFieldGradient(int max_am1, int max_am2)
    : ObaraSaikaTwoCenterElectricField(max_am1+1, max_am2+1)
{
    // Check memory allocation for 3rd parameter. I think it should be max_am1_ + max_am2_ - 2
    exx_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    eyy_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    ezz_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    exy_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    exz_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
    eyz_ = init_box(size_, size_, max_am1_ + max_am2_ + 1);
}

ObaraSaikaTwoCenterElectricFieldGradient::~ObaraSaikaTwoCenterElectricFieldGradient()
{
    free_box(exx_, size_, size_);
    free_box(eyy_, size_, size_);
    free_box(ezz_, size_, size_);
    free_box(exy_, size_, size_);
    free_box(exz_, size_, size_);
    free_box(eyz_, size_, size_);
}

void ObaraSaikaTwoCenterElectricFieldGradient::compute(double PA[3], double PB[3], double PC[3], double zeta, int am1, int am2)
{
    int a, b, m;
    int azm = 1;
    int aym = am1 + 1;
    int axm = aym * aym;
    int bzm = 1;
    int bym = am2 + 1;
    int bxm = bym * bym;
    int ax, ay, az, bx, by, bz;
    int aind, bind;
    double ooz = 1.0/(2.0 * zeta);
    int mmax = am1 + am2;

    // Call our super class to compute electric field and potential integrals
    ObaraSaikaTwoCenterElectricField::compute(PA, PB, PC, zeta, am1+1, am2+1);

    // Compute starting integrals using A26
    for (m=0; m<=mmax-2; ++m) {
        exx_[0][0][m] = 4.0 * zeta * zeta * PC[0] * PC[0] * vi_[0][0][m+2] - 2.0 * zeta * vi_[0][0][m+1];
        eyy_[0][0][m] = 4.0 * zeta * zeta * PC[1] * PC[1] * vi_[0][0][m+2] - 2.0 * zeta * vi_[0][0][m+1];
        ezz_[0][0][m] = 4.0 * zeta * zeta * PC[2] * PC[2] * vi_[0][0][m+2] - 2.0 * zeta * vi_[0][0][m+1];
        exy_[0][0][m] = 4.0 * zeta * zeta * PC[0] * PC[1] * vi_[0][0][m+2];
        exz_[0][0][m] = 4.0 * zeta * zeta * PC[0] * PC[2] * vi_[0][0][m+2];
        eyz_[0][0][m] = 4.0 * zeta * zeta * PC[1] * PC[2] * vi_[0][0][m+2];
    }

    // Perform recursion in b with a=0
    //    subset of A24
    for (b=1; b<=am2; ++b) {
        for (bx=0; bx<=b; ++bx) {
            for (by=0; by<=b-bx; ++by) {
                bz = b-bx-by;

                // Compute the index into exx_, eyy_, ezz_, exy_, exz_, eyz
                bind = bx*bxm + by*bym + bz*bzm;

                // Compute each x, y, z contribution
                if (bz > 0) {
                    for (m=0; m<=mmax-b-2; ++m) {
                        exx_[0][bind][m] = PB[2] * exx_[0][bind-bzm][m] - PC[2] * exx_[0][bind-bzm][m+1];
                        eyy_[0][bind][m] = PB[2] * eyy_[0][bind-bzm][m] - PC[2] * eyy_[0][bind-bzm][m+1];
                        ezz_[0][bind][m] = PB[2] * ezz_[0][bind-bzm][m] - PC[2] * ezz_[0][bind-bzm][m+1] + 2.0 * ez_[0][bind-bzm][m+1];
                        exy_[0][bind][m] = PB[2] * exy_[0][bind-bzm][m] - PC[2] * exy_[0][bind-bzm][m+1];
                        exz_[0][bind][m] = PB[2] * exz_[0][bind-bzm][m] - PC[2] * exz_[0][bind-bzm][m+1] + ex_[0][bind-bzm][m+1];
                        eyz_[0][bind][m] = PB[2] * eyz_[0][bind-bzm][m] - PC[2] * eyz_[0][bind-bzm][m+1] + ey_[0][bind-bzm][m+1];
                    }
                    if (bz > 1) {
                        for (m=0; m<=mmax-b-2; ++m) {
                            exx_[0][bind][m] += ooz * (bz-1) * (exx_[0][bind-2*bzm][m] - exx_[0][bind-2*bzm][m+1]);
                            eyy_[0][bind][m] += ooz * (bz-1) * (eyy_[0][bind-2*bzm][m] - eyy_[0][bind-2*bzm][m+1]);
                            ezz_[0][bind][m] += ooz * (bz-1) * (ezz_[0][bind-2*bzm][m] - ezz_[0][bind-2*bzm][m+1]);
                            exy_[0][bind][m] += ooz * (bz-1) * (exy_[0][bind-2*bzm][m] - exy_[0][bind-2*bzm][m+1]);
                            exz_[0][bind][m] += ooz * (bz-1) * (exz_[0][bind-2*bzm][m] - exz_[0][bind-2*bzm][m+1]);
                            eyz_[0][bind][m] += ooz * (bz-1) * (eyz_[0][bind-2*bzm][m] - eyz_[0][bind-2*bzm][m+1]);
                        }
                    }
                }
                else if (by > 0) {
                    for (m=0; m<=mmax-b-2; ++m) {
                        exx_[0][bind][m] = PB[1] * exx_[0][bind-bym][m] - PC[1] * exx_[0][bind-bym][m+1];
                        eyy_[0][bind][m] = PB[1] * eyy_[0][bind-bym][m] - PC[1] * eyy_[0][bind-bym][m+1] + 2.0 * ey_[0][bind-bym][m+1];
                        ezz_[0][bind][m] = PB[1] * ezz_[0][bind-bym][m] - PC[1] * ezz_[0][bind-bym][m+1];
                        exy_[0][bind][m] = PB[1] * exy_[0][bind-bym][m] - PC[1] * exy_[0][bind-bym][m+1] - ex_[0][bind-bym][m+1];
                        exz_[0][bind][m] = PB[1] * exz_[0][bind-bym][m] - PC[1] * exz_[0][bind-bym][m+1];
                        eyz_[0][bind][m] = PB[1] * eyz_[0][bind-bym][m] - PC[1] * eyz_[0][bind-bym][m+1] - ez_[0][bind-bym][m+1];
                    }
                    if (by > 1) {
                        for (m=0; m<=mmax-b-2; ++m) {
                            exx_[0][bind][m] += ooz * (by-1) * (exx_[0][bind-2*bym][m] - exx_[0][bind-2*bym][m+1]);
                            eyy_[0][bind][m] += ooz * (by-1) * (eyy_[0][bind-2*bym][m] - eyy_[0][bind-2*bym][m+1]);
                            ezz_[0][bind][m] += ooz * (by-1) * (ezz_[0][bind-2*bym][m] - ezz_[0][bind-2*bym][m+1]);
                            exy_[0][bind][m] += ooz * (by-1) * (exy_[0][bind-2*bym][m] - exy_[0][bind-2*bym][m+1]);
                            exz_[0][bind][m] += ooz * (by-1) * (exz_[0][bind-2*bym][m] - exz_[0][bind-2*bym][m+1]);
                            eyz_[0][bind][m] += ooz * (by-1) * (eyz_[0][bind-2*bym][m] - eyz_[0][bind-2*bym][m+1]);
                        }
                    }
                }
                else if (bx > 0) {
                    for (m=0; m<=mmax-b-2; ++m) {
                        exx_[0][bind][m] = PB[0] * exx_[0][bind-bxm][m] - PC[0] * exx_[0][bind-bxm][m+1] + 2.0 * ex_[0][bind-bxm][m+1];
                        eyy_[0][bind][m] = PB[0] * eyy_[0][bind-bxm][m] - PC[0] * eyy_[0][bind-bxm][m+1];
                        ezz_[0][bind][m] = PB[0] * ezz_[0][bind-bxm][m] - PC[0] * ezz_[0][bind-bxm][m+1];
                        exy_[0][bind][m] = PB[0] * exy_[0][bind-bxm][m] - PC[0] * exy_[0][bind-bxm][m+1] - ey_[0][bind-bxm][m+1];
                        exz_[0][bind][m] = PB[0] * exz_[0][bind-bxm][m] - PC[0] * exz_[0][bind-bxm][m+1] - ez_[0][bind-bxm][m+1];
                        eyz_[0][bind][m] = PB[0] * eyz_[0][bind-bxm][m] - PC[0] * eyz_[0][bind-bxm][m+1];
                    }
                    if (bx > 1) {
                        for (m=0; m<=mmax-b-2; ++m) {
                            exx_[0][bind][m] += ooz * (bx-1) * (exx_[0][bind-2*bxm][m] - exx_[0][bind-2*bxm][m+1]);
                            eyy_[0][bind][m] += ooz * (bx-1) * (eyy_[0][bind-2*bxm][m] - eyy_[0][bind-2*bxm][m+1]);
                            ezz_[0][bind][m] += ooz * (bx-1) * (ezz_[0][bind-2*bxm][m] - ezz_[0][bind-2*bxm][m+1]);
                            exy_[0][bind][m] += ooz * (bx-1) * (exy_[0][bind-2*bxm][m] - exy_[0][bind-2*bxm][m+1]);
                            exz_[0][bind][m] += ooz * (bx-1) * (exz_[0][bind-2*bxm][m] - exz_[0][bind-2*bxm][m+1]);
                            eyz_[0][bind][m] += ooz * (bx-1) * (eyz_[0][bind-2*bxm][m] - eyz_[0][bind-2*bxm][m+1]);
                        }
                    }
                }
            }
        }
    }

    // Perform upward recursion in a with all b's
    for (b=0; b<=am2; ++b) {
        for (bx=0; bx<=b; bx++) {
            for (by=0; by<=b-bx;by++) {
                bz = b-bx-by;
                bind = bx*bxm + by*bym + bz*bzm;

                for (a=1; a<=am1; a++) {
                    for (ax=0; ax<=a; ax++) {
                        for (ay=0; ay<=a-ax; ay++) {
                            az = a-ax-ay;
                            aind = ax*axm + ay*aym + az*azm;

                            if (az > 0) {
                                for (m=0; m<=mmax-a-b-2; ++m) {
                                    exx_[aind][bind][m] = PA[2] * exx_[aind-azm][bind][m] - PC[2] * exx_[aind-azm][bind][m+1];
                                    eyy_[aind][bind][m] = PA[2] * eyy_[aind-azm][bind][m] - PC[2] * eyy_[aind-azm][bind][m+1];
                                    ezz_[aind][bind][m] = PA[2] * ezz_[aind-azm][bind][m] - PC[2] * ezz_[aind-azm][bind][m+1] - 2.0 * ez_[aind-azm][bind][m+1];
                                    exy_[aind][bind][m] = PA[2] * exy_[aind-azm][bind][m] - PC[2] * exy_[aind-azm][bind][m+1];
                                    exz_[aind][bind][m] = PA[2] * exz_[aind-azm][bind][m] - PC[2] * exz_[aind-azm][bind][m+1] + ex_[aind-azm][bind][m+1];
                                    eyz_[aind][bind][m] = PA[2] * eyz_[aind-azm][bind][m] - PC[2] * eyz_[aind-azm][bind][m+1] + ey_[aind-azm][bind][m+1];
                                }
                                if (az > 1) {
                                    for (m=0; m<=mmax-a-b-2; ++m) {
                                        exx_[aind][bind][m] += ooz * (az-1) * (exx_[aind-2*azm][bind][m] - exx_[aind-2*azm][bind][m+1]);
                                        eyy_[aind][bind][m] += ooz * (az-1) * (eyy_[aind-2*azm][bind][m] - eyy_[aind-2*azm][bind][m+1]);
                                        ezz_[aind][bind][m] += ooz * (az-1) * (ezz_[aind-2*azm][bind][m] - ezz_[aind-2*azm][bind][m+1]);
                                        exy_[aind][bind][m] += ooz * (az-1) * (exy_[aind-2*azm][bind][m] - exy_[aind-2*azm][bind][m+1]);
                                        exz_[aind][bind][m] += ooz * (az-1) * (exz_[aind-2*azm][bind][m] - exz_[aind-2*azm][bind][m+1]);
                                        eyz_[aind][bind][m] += ooz * (az-1) * (eyz_[aind-2*azm][bind][m] - eyz_[aind-2*azm][bind][m+1]);
                                    }
                                }
                                if (bz > 0) {
                                    for (m=0; m<=mmax-a-b-2; ++m) {
                                        exx_[aind][bind][m] += ooz * bz * (exx_[aind-azm][bind-bzm][m] - exx_[aind-azm][bind-bzm][m+1]);
                                        eyy_[aind][bind][m] += ooz * bz * (eyy_[aind-azm][bind-bzm][m] - eyy_[aind-azm][bind-bzm][m+1]);
                                        ezz_[aind][bind][m] += ooz * bz * (ezz_[aind-azm][bind-bzm][m] - ezz_[aind-azm][bind-bzm][m+1]);
                                        exy_[aind][bind][m] += ooz * bz * (exy_[aind-azm][bind-bzm][m] - exy_[aind-azm][bind-bzm][m+1]);
                                        exz_[aind][bind][m] += ooz * bz * (exz_[aind-azm][bind-bzm][m] - exz_[aind-azm][bind-bzm][m+1]);
                                        eyz_[aind][bind][m] += ooz * bz * (eyz_[aind-azm][bind-bzm][m] - eyz_[aind-azm][bind-bzm][m+1]);
                                    }
                                }
                            }
                            else if (ay > 0) {
                                for (m=0; m<=mmax-a-b-2; ++m) {
                                    exx_[aind][bind][m] = PA[1] * exx_[aind-aym][bind][m] - PC[1] * exx_[aind-aym][bind][m+1];
                                    eyy_[aind][bind][m] = PA[1] * eyy_[aind-aym][bind][m] - PC[1] * eyy_[aind-aym][bind][m+1] - 2.0 * ey_[aind-aym][bind][m+1];
                                    ezz_[aind][bind][m] = PA[1] * ezz_[aind-aym][bind][m] - PC[1] * ezz_[aind-aym][bind][m+1];
                                    exy_[aind][bind][m] = PA[1] * exy_[aind-aym][bind][m] - PC[1] * exy_[aind-aym][bind][m+1] + ex_[aind-aym][bind][m+1];
                                    exz_[aind][bind][m] = PA[1] * exz_[aind-aym][bind][m] - PC[1] * exz_[aind-aym][bind][m+1];
                                    eyz_[aind][bind][m] = PA[1] * eyz_[aind-aym][bind][m] - PC[1] * eyz_[aind-aym][bind][m+1] + ez_[aind-aym][bind][m+1];
                                }
                                if (ay > 1) {
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        exx_[aind][bind][m] += ooz * (ay-1) * (exx_[aind-2*aym][bind][m] - exx_[aind-2*aym][bind][m+1]);
                                        eyy_[aind][bind][m] += ooz * (ay-1) * (eyy_[aind-2*aym][bind][m] - eyy_[aind-2*aym][bind][m+1]);
                                        ezz_[aind][bind][m] += ooz * (ay-1) * (ezz_[aind-2*aym][bind][m] - ezz_[aind-2*aym][bind][m+1]);
                                        exy_[aind][bind][m] += ooz * (ay-1) * (exy_[aind-2*aym][bind][m] - exy_[aind-2*aym][bind][m+1]);
                                        exz_[aind][bind][m] += ooz * (ay-1) * (exz_[aind-2*aym][bind][m] - exz_[aind-2*aym][bind][m+1]);
                                        eyz_[aind][bind][m] += ooz * (ay-1) * (eyz_[aind-2*aym][bind][m] - eyz_[aind-2*aym][bind][m+1]);
                                    }
                                }
                                if (by > 0) {
                                    for (m=0; m<=mmax-a-b-2; ++m) {
                                        exx_[aind][bind][m] += ooz * by * (exx_[aind-aym][bind-bym][m] - exx_[aind-aym][bind-bym][m+1]);
                                        eyy_[aind][bind][m] += ooz * by * (eyy_[aind-aym][bind-bym][m] - eyy_[aind-aym][bind-bym][m+1]);
                                        ezz_[aind][bind][m] += ooz * by * (ezz_[aind-aym][bind-bym][m] - ezz_[aind-aym][bind-bym][m+1]);
                                        exy_[aind][bind][m] += ooz * by * (exy_[aind-aym][bind-bym][m] - exy_[aind-aym][bind-bym][m+1]);
                                        exz_[aind][bind][m] += ooz * by * (exz_[aind-aym][bind-bym][m] - exz_[aind-aym][bind-bym][m+1]);
                                        eyz_[aind][bind][m] += ooz * by * (eyz_[aind-aym][bind-bym][m] - eyz_[aind-aym][bind-bym][m+1]);
                                    }
                                }
                            }
                            else if (ax > 0) {
                                for (m=0; m<=mmax-a-b-2; ++m) {
                                    exx_[aind][bind][m] = PA[0] * exx_[aind-axm][bind][m] - PC[0] * exx_[aind-axm][bind][m+1] - 2.0 * ex_[aind-axm][bind][m+1];
                                    eyy_[aind][bind][m] = PA[0] * eyy_[aind-axm][bind][m] - PC[0] * eyy_[aind-axm][bind][m+1];
                                    ezz_[aind][bind][m] = PA[0] * ezz_[aind-axm][bind][m] - PC[0] * ezz_[aind-axm][bind][m+1];
                                    exy_[aind][bind][m] = PA[0] * exy_[aind-axm][bind][m] - PC[0] * exy_[aind-axm][bind][m+1] + ey_[aind-axm][bind][m+1];
                                    exz_[aind][bind][m] = PA[0] * exz_[aind-axm][bind][m] - PC[0] * exz_[aind-axm][bind][m+1] + ez_[aind-axm][bind][m+1];
                                    eyz_[aind][bind][m] = PA[0] * eyz_[aind-axm][bind][m] - PC[0] * eyz_[aind-axm][bind][m+1];
                                }
                                if (ax > 1) {
                                    for (m=0; m<=mmax-a-b-1; ++m) {
                                        exx_[aind][bind][m] += ooz * (ax-1) * (exx_[aind-2*axm][bind][m] - exx_[aind-2*axm][bind][m+1]);
                                        eyy_[aind][bind][m] += ooz * (ax-1) * (eyy_[aind-2*axm][bind][m] - eyy_[aind-2*axm][bind][m+1]);
                                        ezz_[aind][bind][m] += ooz * (ax-1) * (ezz_[aind-2*axm][bind][m] - ezz_[aind-2*axm][bind][m+1]);
                                        exy_[aind][bind][m] += ooz * (ax-1) * (exy_[aind-2*axm][bind][m] - exy_[aind-2*axm][bind][m+1]);
                                        exz_[aind][bind][m] += ooz * (ax-1) * (exz_[aind-2*axm][bind][m] - exz_[aind-2*axm][bind][m+1]);
                                        eyz_[aind][bind][m] += ooz * (ax-1) * (eyz_[aind-2*axm][bind][m] - eyz_[aind-2*axm][bind][m+1]);
                                    }
                                }
                                if (bx > 0) {
                                    for (m=0; m<=mmax-a-b-2; ++m) {
                                        exx_[aind][bind][m] += ooz * bx * (exx_[aind-axm][bind-bxm][m] - exx_[aind-axm][bind-bxm][m+1]);
                                        eyy_[aind][bind][m] += ooz * bx * (eyy_[aind-axm][bind-bxm][m] - eyy_[aind-axm][bind-bxm][m+1]);
                                        ezz_[aind][bind][m] += ooz * bx * (ezz_[aind-axm][bind-bxm][m] - ezz_[aind-axm][bind-bxm][m+1]);
                                        exy_[aind][bind][m] += ooz * bx * (exy_[aind-axm][bind-bxm][m] - exy_[aind-axm][bind-bxm][m+1]);
                                        exz_[aind][bind][m] += ooz * bx * (exz_[aind-axm][bind-bxm][m] - exz_[aind-axm][bind-bxm][m+1]);
                                        eyz_[aind][bind][m] += ooz * bx * (eyz_[aind-axm][bind-bxm][m] - eyz_[aind-axm][bind-bxm][m+1]);
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

