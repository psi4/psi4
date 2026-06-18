#include "psi4/pragma.h"

#include "psi4/libfmm/multipoles_helper.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/vector3.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include <functional>
#include <memory>
#include <tuple>
#include <vector>
#include <cmath>
#include <unordered_map>

namespace psi {

MultipoleRotationFactory::MultipoleRotationFactory(Vector3 R_a, Vector3 R_b, int lmax) {
    lmax_ = lmax;
    Vector3 R_ab = R_a - R_b;
    double R = R_ab.norm();

    Vector3 z_axis = R_ab / R;
    Vector3 ca(z_axis);

    if (R_a[1] == R_b[1] && R_a[2] == R_b[2]) {
        ca[1] += 1.0;
    } else {
        ca[0] += 1.0;
    }

    double this_dot = ca.dot(z_axis);
    ca -= (z_axis * this_dot);

    Vector3 x_axis = ca / ca.norm();
    Vector3 y_axis = z_axis.cross(x_axis);

    Uz_ = std::make_shared<Matrix>("Uz Matrix", 3, 3);

    for (int i = 0; i < 3; i++) {
        Uz_->set(0, i, x_axis[i]);
        Uz_->set(1, i, y_axis[i]);
        Uz_->set(2, i, z_axis[i]);
    }

    for (int l = 0; l <= lmax_; l++) {
        D_cache_.push_back(nullptr);
    }

}

double MultipoleRotationFactory::U(int l, int m, int M) {
    return P(0, l, m, M);
}

double MultipoleRotationFactory::V(int l, int m, int M) {
    double val;
    if (m == 0) {
        val = P(1, l, 1, M) + P(-1, l, -1, M);
    } else if (m > 0) {
        if (m == 1) {
            val = std::sqrt(2.0) * P(1, l, m-1, M);
        } else {
            val = P(1, l, m-1, M) - P(-1, l, -m+1, M);
        }
    } else {
        if (m == -1) {
            val = std::sqrt(2.0) * P(-1, l, -m-1, M);
        } else {
            val = P(1, l, m+1, M) + P(-1, l, -m-1, M);
        }
    }
    return val;
}

double MultipoleRotationFactory::W(int l, int m, int M) {
    double val;
    if (m == 0) {
        val = 0.0;
    } else if (m > 0) {
        val = P(1, l, m+1, M) + P(-1, l, -m-1, M);
    } else {
        val = P(1, l, m-1, M) - P(-1, l, -m+1, M);
    }
    return val;
}

double MultipoleRotationFactory::P(int i, int l, int mu, int M) {
    int I = m_addr(i);
    SharedMatrix D1 = get_D(1);
    SharedMatrix Dl1 = get_D(l-1);
    if (std::abs(M) < l) {
        return D1->get(I, 0) * Dl1->get(m_addr(mu), m_addr(M));
    } else if (M == l) {
        return D1->get(I, 1) * Dl1->get(m_addr(mu), m_addr(M-1)) - D1->get(I, 2) * Dl1->get(m_addr(mu), m_addr(-M+1));
    } else {
        return D1->get(I, 1) * Dl1->get(m_addr(mu), m_addr(M+1)) + D1->get(I, 2) * Dl1->get(m_addr(mu), m_addr(-M-1));
    }
}

SharedMatrix MultipoleRotationFactory::get_D(int l) {
    if (l > lmax_) {
        throw PsiException("Input l is larger than the set lmax for the Rotation Matrix", __FILE__, __LINE__);
    }

    if (D_cache_[l]) {
        return D_cache_[l];
    }

    SharedMatrix Drot;

    if (l == 0) {
        Drot = std::make_shared<Matrix>("D Rotation Matrix", 1, 1);
        Drot->set(0, 0, 1.0);
    } else if (l == 1) {
        Drot = std::make_shared<Matrix>("D Rotation Matrix", 3, 3);
        std::vector<int> permute {2, 0, 1};
        for (int i = 0; i < 3; i++) {
            int ip = permute[i];
            for (int j = 0; j < 3; j++) {
                int jp = permute[j];
                Drot->set(i, j, Uz_->get(ip, jp));
            }
        }
    } else {
        Drot = std::make_shared<Matrix>("D Rotation Matrix", 2*l+1, 2*l+1);

        for (int m1 = -l; m1 <= l; m1++) {
            int k1 = m_addr(m1);
            for (int m2 = -l; m2 <= l; m2++) {
                int k2 = m_addr(m2);
                double Uterm = u(l, m1, m2);
                if (Uterm != 0.0) Uterm *= U(l, m1, m2);
                double Vterm = v(l, m1, m2);
                if (Vterm != 0.0) Vterm *= V(l, m1, m2);
                double Wterm = w(l, m1, m2);
                if (Wterm != 0.0) Wterm *= W(l, m1, m2);
                Drot->set(k1, k2, Uterm + Vterm + Wterm);
            }
        }
    }

    D_cache_[l] = Drot;
    return D_cache_[l];
}

HarmonicCoefficients::HarmonicCoefficients(int lmax, SolidHarmonicsType type) {
    lmax_ = lmax;
    type_ = type;

    Rc_.resize(lmax_+1);
    Rs_.resize(lmax_+1);
    mpole_terms_.resize(lmax_+1);

    for (int l = 0; l <= lmax_; l++) {
        int nc = ncart(l);
        Rc_[l].resize(l+1);
        Rs_[l].resize(l+1);
        mpole_terms_[l].resize(2*l+1);
    }

    if (type_ == Regular) compute_terms_regular();
    if (type_ == Irregular) compute_terms_irregular();
}

/// TODO: Implement Helgaker 9.13.85 - 9.13.89
void HarmonicCoefficients::compute_terms_irregular() {
    throw FeatureNotImplemented("libfmm", "RealSolidHarmonics::compute_terms_irregular()", __FILE__, __LINE__);
}

/// Helgaker 9.13.78 - 9.13.82
void HarmonicCoefficients::compute_terms_regular() {

    for (int l = 0; l <= lmax_; l++) {
        int ncl = ncart(l);

        if (l == 0) {
            Rc_[0][0][0] = 1.0;
            Rs_[0][0][0] = 0.0;
        } else {
            // m < l-1 terms
            for (int m = 0; m < l-1; m++) {
                double denom = (l+m)*(l-m);

                // Rc_[l-1][m] contributions to Rc_[l][m]
                for (const std::pair<int, double>& rpair : Rc_[l-1][m]) {
                    int ncl1 = ncart(l-1);
                    int ind = rpair.first;
                    double coef = rpair.second;
                    if (std::abs(coef) < 1.0e-16) continue;
                    int a = ind / (ncl1*ncl1);
                    int bc = ind % (ncl1*ncl1);
                    int b = bc / ncl1;
                    int c = bc % ncl1;

                    int newind = a * ncl * ncl + b * ncl + (c+1);

                    Rc_[l][m][newind] += coef * (2*l-1) / denom;

                }

                // Rc_[l-2][m] contributions to Rc_[l][m]
                for (const std::pair<int, double>& rpair : Rc_[l-2][m]) {
                    int ncl2 = ncart(l-2);
                    int ind = rpair.first;
                    double coef = rpair.second;
                    if (std::abs(coef) < 1.0e-16) continue;
                    int a = ind / (ncl2*ncl2);
                    int bc = ind % (ncl2*ncl2);
                    int b = bc / ncl2;
                    int c = bc % ncl2;

                    int newind1 = (a+2) * ncl * ncl + b * ncl + c;
                    Rc_[l][m][newind1] += -coef / denom;
                    int newind2 = a * ncl * ncl + (b+2) * ncl + c;
                    Rc_[l][m][newind2] += -coef / denom;
                    int newind3 = a * ncl * ncl + b * ncl + (c+2);
                    Rc_[l][m][newind3] += -coef / denom;

                }

                // Rs_[l-1][m] contributions to Rs_[l][m]
                for (const std::pair<int, double>& rpair : Rs_[l-1][m]) {
                    int ncl1 = ncart(l-1);
                    int ind = rpair.first;
                    double coef = rpair.second;
                    if (std::abs(coef) < 1.0e-16) continue;
                    int a = ind / (ncl1*ncl1);
                    int bc = ind % (ncl1*ncl1);
                    int b = bc / ncl1;
                    int c = bc % ncl1;

                    int newind = a * ncl * ncl + b * ncl + (c+1);

                    Rs_[l][m][newind] += coef * (2*l-1) / denom;

                }

                // Rs_[l-2][m] contributions to Rs_[l][m]
                for (const std::pair<int, double>& rpair : Rs_[l-2][m]) {
                    int ncl2 = ncart(l-2);
                    int ind = rpair.first;
                    double coef = rpair.second;
                    if (std::abs(coef) < 1.0e-16) continue;
                    int a = ind / (ncl2*ncl2);
                    int bc = ind % (ncl2*ncl2);
                    int b = bc / ncl2;
                    int c = bc % ncl2;

                    int newind1 = (a+2) * ncl * ncl + b * ncl + c;
                    Rs_[l][m][newind1] += -coef / denom;
                    int newind2 = a * ncl * ncl + (b+2) * ncl + c;
                    Rs_[l][m][newind2] += -coef / denom;
                    int newind3 = a * ncl * ncl + b * ncl + (c+2);
                    Rs_[l][m][newind3] += -coef / denom;

                }
            }

            // => m = l-1 <= //

            // Rc[l][l-1]
            for (const std::pair<int, double>& rpair : Rc_[l-1][l-1]) {
                int ncl1 = ncart(l-1);
                int ind = rpair.first;
                double coef = rpair.second;
                if (std::abs(coef) < 1.0e-16) continue;
                int a = ind / (ncl1*ncl1);
                int bc = ind % (ncl1*ncl1);
                int b = bc / ncl1;
                int c = bc % ncl1;

                int newind = a * ncl * ncl + b * ncl + (c+1);

                Rc_[l][l-1][newind] += coef;
            }

            // Rs[l][l-1]
            for (const std::pair<int, double>& rpair : Rs_[l-1][l-1]) {
                int ncl1 = ncart(l-1);
                int ind = rpair.first;
                double coef = rpair.second;
                if (std::abs(coef) < 1.0e-16) continue;
                int a = ind / (ncl1*ncl1);
                int bc = ind % (ncl1*ncl1);
                int b = bc / ncl1;
                int c = bc % ncl1;

                int newind = a * ncl * ncl + b * ncl + (c+1);

                Rs_[l][l-1][newind] += coef;
            }

            // => m = l <= //

            // Rc[l-1][l-1] contribution to Rc[l][l]
            for (const std::pair<int, double>& rpair : Rc_[l-1][l-1]) {
                int ncl1 = ncart(l-1);
                int ind = rpair.first;
                double coef = rpair.second;
                if (std::abs(coef) < 1.0e-16) continue;
                int a = ind / (ncl1*ncl1);
                int bc = ind % (ncl1*ncl1);
                int b = bc / ncl1;
                int c = bc % ncl1;

                int newind = (a+1) * ncl * ncl + b * ncl + c;

                Rc_[l][l][newind] += -coef/(2*l);

            }

            // Rs[l-1][l-1] contribution to Rc[l][l]
            for (const std::pair<int, double>& rpair : Rs_[l-1][l-1]) {
                int ncl1 = ncart(l-1);
                int ind = rpair.first;
                double coef = rpair.second;
                if (std::abs(coef) < 1.0e-16) continue;
                int a = ind / (ncl1*ncl1);
                int bc = ind % (ncl1*ncl1);
                int b = bc / ncl1;
                int c = bc % ncl1;

                int newind = a * ncl * ncl + (b+1) * ncl + c;

                Rc_[l][l][newind] += coef/(2*l);

            }

            // Rc[l-1][l-1] contribution to Rs[l][l]
            for (const std::pair<int, double>& rpair : Rc_[l-1][l-1]) {
                int ncl1 = ncart(l-1);
                int ind = rpair.first;
                double coef = rpair.second;
                if (std::abs(coef) < 1.0e-16) continue;
                int a = ind / (ncl1*ncl1);
                int bc = ind % (ncl1*ncl1);
                int b = bc / ncl1;
                int c = bc % ncl1;

                int newind = a * ncl * ncl + (b+1) * ncl + c;

                Rs_[l][l][newind] += -coef/(2*l);

            }

            // Rs[l-1][l-1] contribution to Rs[l][l]
            for (const std::pair<int, double>& rpair : Rs_[l-1][l-1]) {
                int ncl1 = ncart(l-1);
                int ind = rpair.first;
                double coef = rpair.second;
                if (std::abs(coef) < 1.0e-16) continue;
                int a = ind / (ncl1*ncl1);
                int bc = ind % (ncl1*ncl1);
                int b = bc / ncl1;
                int c = bc % ncl1;

                int newind = (a+1) * ncl * ncl + b * ncl + c;

                Rs_[l][l][newind] += -coef/(2*l);

            }
        }

        // => Renormalization <= //
        
        // Convert FROM Helgaker Convention (Equ. 9.13.14) TO Stone Convention (Equ B.1.3)
        for (int m = -l; m <= l; m++) {
            // m is signed address
            // mu is unsigned address
            int mu = m_addr(m);
            double prefactor = 1.0;

            if ((mu == 0) || (mu % 2 == 1)) {
                if (mu == 0) {
                    prefactor = factorial(l);
                } else {
                    prefactor = std::pow(-1.0, (double) m) * std::sqrt(2.0 * factorial(l-m) * factorial(l+m));
                }
                for (const std::pair<int, double>& rpair : Rc_[l][m]) {
                    int ind = rpair.first;
                    double coef = rpair.second;
                    if (std::abs(coef) < 1.0e-16) continue;
                    mpole_terms_[l][mu][ind] = prefactor * coef;
                }

            } else {
                prefactor = std::pow(-1.0, (double) m) * std::sqrt(2.0 * factorial(l-m) * factorial(l+m));
                for (const std::pair<int, double>& rpair : Rs_[l][-m]) {
                    int ind = rpair.first;
                    double coef = rpair.second;
                    if (std::abs(coef) < 1.0e-16) continue;
                    mpole_terms_[l][mu][ind] = prefactor * coef;
                }
            }
        }
    }
}

RealSolidHarmonics::RealSolidHarmonics(int lmax, Vector3 center, SolidHarmonicsType type) {
    lmax_ = lmax;
    center_ = center;
    type_ = type;

    // Set all multipoles to zero
    Ylm_.resize(lmax_+1);

    for (int l = 0; l <= lmax_; l++) {
        Ylm_[l].resize(2*l+1, 0.0);
    }
    
}

std::shared_ptr<RealSolidHarmonics> RealSolidHarmonics::copy() {
    std::shared_ptr<RealSolidHarmonics> new_harm = std::make_shared<RealSolidHarmonics>(lmax_, center_, type_);
    for (int l = 0; l <= lmax_; l++) {
        for (int mu = 0; mu < 2*l+1; mu++) {
            new_harm->Ylm_[l][mu] = Ylm_[l][mu];
        }
    }
    return new_harm;
}

void RealSolidHarmonics::add(const RealSolidHarmonics& rsh) {
    for (int l = 0; l <= lmax_; l++) {
        for (int mu = 0; mu < 2*l+1; mu++) {
#pragma omp atomic
            Ylm_[l][mu] += rsh.Ylm_[l][mu];
        }
    }
}

void RealSolidHarmonics::add(const std::shared_ptr<RealSolidHarmonics>& rsh) {
    this->add(*rsh);
}

double RealSolidHarmonics::dot(const RealSolidHarmonics& rsh) {
    double result;
    for (int l = 0; l <= lmax_; l++) {
        for (int mu = 0; mu < 2*l+1; mu++) {
            result += Ylm_[l][mu] * rsh.Ylm_[l][mu];
        }
    }
    return result;
}

double RealSolidHarmonics::dot(const std::shared_ptr<RealSolidHarmonics>& rsh) {
    return this->dot(*rsh);
}

void RealSolidHarmonics::scale(double val) {
    for (int l = 0; l <= lmax_; l++) {
        for (int mu = 0; mu < 2*l+1; mu++) {
            Ylm_[l][mu] *= val;
        }
    }
}

std::shared_ptr<RealSolidHarmonics> RealSolidHarmonics::translate(const Vector3& new_center) {
    if (type_ == Regular) return translate_regular(new_center);
    if (type_ == Irregular) return translate_irregular(new_center);
}

std::shared_ptr<RealSolidHarmonics> RealSolidHarmonics::translate_irregular(Vector3 new_center) {
    auto trans_harmonics = std::make_shared<RealSolidHarmonics>(lmax_, new_center, Irregular);
    auto rotation_factory = std::make_shared<MultipoleRotationFactory>(center_, new_center, lmax_);

    Vector3 R_ab = new_center - center_;
    double R = R_ab.norm();

    std::vector<std::vector<double>> rot_mpoles(lmax_+1);
    std::vector<std::vector<double>> trans_rot_mpoles(lmax_+1);

    // Rotate Multipoles to direction of translation
    for (int l = 0; l <= lmax_; l++) {
        rot_mpoles[l].resize(2*l+1, 0.0);
        int dim = 2*l+1;
        SharedMatrix Dmat = rotation_factory->get_D(l);
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                rot_mpoles[l][i] += Dmat->get(i, j) * Ylm_[l][j];
            }
        }
    }

    // Translate Rotated Multipoles
    for (int l = 0; l <= lmax_; l++) {
        trans_rot_mpoles[l].resize(2*l+1, 0.0);
        for (int j = l; j <= lmax_; j++) {
            for (int m = -l; m <= l; m++) {
                int mu = m_addr(m);
                double coef = std::sqrt((double) choose(j+m,l+m)*choose(j-m,l-m));

                trans_rot_mpoles[l][mu] += coef * std::pow(R, j-l) * rot_mpoles[j][mu];
            }
        }
    }

    // Backrotation of Multipoles
    for (int l = 0; l <= lmax_; l++) {
        SharedMatrix Dmat = rotation_factory->get_D(l);
        int dim = 2*l+1;
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                trans_harmonics->Ylm_[l][i] += Dmat->get(j, i) * trans_rot_mpoles[l][j];
            }
        }
    }

    return trans_harmonics;
}

std::shared_ptr<RealSolidHarmonics> RealSolidHarmonics::translate_regular(Vector3 new_center) {
    auto trans_harmonics = std::make_shared<RealSolidHarmonics>(lmax_, new_center, Regular);
    auto rotation_factory = std::make_shared<MultipoleRotationFactory>(center_, new_center, lmax_);

    Vector3 R_ab = new_center - center_;
    double R = R_ab.norm();

    std::vector<std::vector<double>> rot_mpoles(lmax_+1);
    std::vector<std::vector<double>> trans_rot_mpoles(lmax_+1);

    for (int l = 0; l <= lmax_; l++) {
        rot_mpoles[l].resize(2*l+1, 0.0);
        trans_rot_mpoles[l].resize(2*l+1, 0.0);
    }

    // Rotate Multipoles to direction of translation
    for (int l = 0; l <= lmax_; l++) {
        int dim = 2*l+1;
        SharedMatrix Dmat = rotation_factory->get_D(l);
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                rot_mpoles[l][i] += Dmat->get(i, j) * Ylm_[l][j];
            }
        }
    }

    // Translate Rotated Multipoles
    for (int l = 0; l <= lmax_; l++) {
        for (int j = 0; j <= l; j++) {
            for (int m = -j; m <= j; m++) {
                int mu = m_addr(m);
                double coef = std::sqrt((double) choose(l+m,j+m)*choose(l-m,j-m));
                trans_rot_mpoles[l][mu] += coef * std::pow(-R, l-j) * rot_mpoles[j][mu];
            }
        }
    }

    // Backrotation of Multipoles
    for (int l = 0; l <= lmax_; l++) {
        SharedMatrix Dmat = rotation_factory->get_D(l);
        int dim = 2*l+1;
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                trans_harmonics->Ylm_[l][i] += Dmat->get(j, i) * trans_rot_mpoles[l][j];
            }
        }
    }

    return trans_harmonics;
}

// A helper method to compute the interaction tensor between aligned multipoles after rotation
SharedVector RealSolidHarmonics::build_T_spherical(int la, int lb, double R) {
    int lmin = std::min(la, lb);
    SharedVector Tvec = std::make_shared<Vector>(2*lmin+1);
    double denom = std::pow(R, (double) la+lb+1);

    for (int m = -lmin; m <= lmin; m++) {
        int mu = m_addr(m);
        double Tval = std::pow(-1.0, (double) (lb-m)) * std::sqrt((double) choose(la+lb, la+m) * choose(la+lb, la-m)) / denom;
        Tvec->set(mu, Tval);
    }

    return Tvec;
}

std::shared_ptr<RealSolidHarmonics> RealSolidHarmonics::far_field_vector(const Vector3& far_center) {
    
    Vector3 R_ab = far_center - center_;
    double R = R_ab.norm();

    auto Vff = std::make_shared<RealSolidHarmonics>(lmax_, far_center, Irregular);
    auto rotation_factory = std::make_shared<MultipoleRotationFactory>(far_center, center_, lmax_);

    for (int l = 0; l <= lmax_; l++) {
        for (int j = 0; j <= lmax_; j++) {
            SharedVector Tvec = build_T_spherical(l, j, R);
            int nterms = 2*std::min(l,j)+1;

            std::vector<double> rotated_mpole(2*j+1, 0.0);
            SharedMatrix Dmat = rotation_factory->get_D(j);
            for (int u = 0; u < 2*j+1; u++) {
                for (int v = 0; v < 2*j+1; v++) {
                    rotated_mpole[u] += Dmat->get(u, v) * Ylm_[j][v];
                }
            }

            std::vector<double> temp(nterms, 0.0);
            for (int u = 0; u < nterms; u++) {
                temp[u] = Tvec->get(u) * rotated_mpole[u];
            }

            SharedMatrix Dl = rotation_factory->get_D(l);
            for (int r = 0; r < 2*l+1; r++) {
                for (int s = 0; s < nterms; s++) {
                    Vff->Ylm_[l][r] += Dl->get(s, r) * temp[s];
                }
            }

        }
    }

    return Vff;
}

} // namespace psi
