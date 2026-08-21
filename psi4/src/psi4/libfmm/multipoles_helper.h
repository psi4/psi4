/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef PSI4_LIBFMM_MULTIPOLES_HELPER_H
#define PSI4_LIBFMM_MULTIPOLES_HELPER_H

#include "psi4/pragma.h"

#include "psi4/libmints/vector.h"
#include "psi4/libmints/vector3.h"
#include "psi4/libmints/matrix.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <unordered_map>
#include <vector>

namespace psi {

enum SolidHarmonicsType { Regular, Irregular };

// Indexing conventions shared with libmints/solidharmonics.cc.
static inline int npure(int l) { return 2 * l + 1; }
static inline int icart(int a, int b, int c) { return (((((a + b + c + 1) << 1) - a) * (a + 1)) >> 1) - b - 1; }
static inline int ipure(int, int m) { return m < 0 ? 2 * -m : (m == 0 ? 0 : 2 * m - 1); }

static inline int ncart(int l) { return (l + 1) * (l + 2) / 2; }

static inline double choose(int n, int r) {
    if (r < 0 || r > n) {
        return 0.0;
    }
    int small = std::min(n, n - r);
    double nCr = 1.0;
    for (int t = 0; t < small; t++) {
        nCr *= n;
        nCr /= (t + 1);
        n -= 1;
    }
    return nCr;
}

static inline int m_addr(int m) {
    /*- Return the unsigned (array) address of m -*/
    if (m <= 0) {
        // 0, 1s, 2s, 3s, ...
        return 2 * (-m);
    } else {
        // 1c, 2c, 3c, ...
        return 2 * m - 1;
    }
}

static inline double cfmm_factorial(const double n) {
    return std::tgamma(n + 1.0);
}

class PSI_API MultipoleRotationFactory {

    protected:
      // New Z axis in rotated frame of reference
      SharedMatrix Uz_;

      // Maximal Angular Momentum
      int lmax_;

      // Cached Rotation Matrices in a vector of Matrices
      std::vector<SharedMatrix> D_cache_;

      double U(int l, int m, int M);
      double V(int l, int m, int M);
      double W(int l, int m, int M);
      double P(int i, int l, int mu, int M);

    inline double u(int l, int m, int M) {
        if (std::abs(M) < l) {
            return std::sqrt((double)(l+m)*(l-m) /((l+M)*(l-M)));
        } else {
            return std::sqrt((double)(l+m)*(l-m)/(2*l*(2*l-1)));
        }
    }

    inline double v(int l, int m, int M) {
        double dm0 = 0.0;
        if (m == 0) dm0 = 1.0;
        if (std::abs(M) < l) {
            return 0.5 * (1.0 - 2.0*dm0) * std::sqrt((1.0+dm0)*(l+std::abs(m)-1)*(l+std::abs(m)) / ((l+M)*(l-M)));
        } else {
            return 0.5 * (1.0 - 2.0*dm0) * std::sqrt((1.0+dm0)*(l+std::abs(m)-1)*(l+std::abs(m)) / ((2*l)*(2*l-1)));
        }
    }

    inline double w(int l, int m, int M) {
        double dm0 = 0.0;
        if (m == 0) dm0 = 1.0;
        if (std::abs(M) < l) {
            return 0.5 * (dm0 - 1) * std::sqrt((double)(l-std::abs(m)-1)*(l-std::abs(m)) / ((l+M)*(l-M)));
        } else {
            return 0.5 * (dm0 - 1) * std::sqrt((double)(l-std::abs(m)-1)*(l-std::abs(m)) / ((2*l)*(2*l-1)));
        }
    }

    public:
      // Constructor
      MultipoleRotationFactory(Vector3 source_center, Vector3 target_center, int lmax);

      SharedMatrix get_D(int l);

}; // End MultipoleRotationFactory

class PSI_API HarmonicCoefficients {
    protected:
      // Ylm[l][m] = sum(coeff * x^a * y^b * z^c). Each map stores a flattened
      // Cartesian exponent index -> coefficient and uses Stone's normalization.
      std::vector<std::vector<std::unordered_map<int, double>>> mpole_terms_;
      // Cosine-like real regular-solid-harmonic coefficients.
      std::vector<std::vector<std::unordered_map<int, double>>> Rc_;
      // Sine-like real regular-solid-harmonic coefficients.
      std::vector<std::vector<std::unordered_map<int, double>>> Rs_;
      // Maximum angular momentum
      int lmax_;
      // Regular or Irregular?
      SolidHarmonicsType type_;

      // Compute regular-solid-harmonic terms.
      void compute_terms_regular();
      // Compute irregular-solid-harmonic terms.
      void compute_terms_irregular();

    public:
      // Constructor
      HarmonicCoefficients(int lmax, SolidHarmonicsType type);
      // Returns a reference to the terms
      const std::unordered_map<int, double>& get_terms(int l, int mu) const { return mpole_terms_[l][mu]; }

};

class PSI_API RealSolidHarmonics {

    protected:
      // Values of the Real Solid Harmonics, normalized according to Stone's convention
      std::vector<std::vector<double>> Ylm_;
      // Maximum angular momentum
      int lmax_;
      // Regular or Irregular?
      SolidHarmonicsType type_;
      // Center of the Harmonics
      Vector3 center_;

      // Return a translated copy of a regular expansion.
      std::shared_ptr<RealSolidHarmonics> translate_regular(Vector3 new_center);
      // Return a translated copy of an irregular expansion.
      std::shared_ptr<RealSolidHarmonics> translate_irregular(Vector3 new_center);

    public:
      // Constructor
      RealSolidHarmonics(int lmax, Vector3 center, SolidHarmonicsType type);

      // Returns a copy of a RealSolidHarmonics object
      std::shared_ptr<RealSolidHarmonics> copy();
      // Adds two harmonics together
      void add(const RealSolidHarmonics& other);
      void add(const std::shared_ptr<RealSolidHarmonics>& other);
      // Contract two multipole expansions element by element.
      double dot(const RealSolidHarmonics& other);
      double dot(const std::shared_ptr<RealSolidHarmonics>& other);
      // Scale the harmonics by a constant
      void scale(double factor);
      // Adds to a specific harmonic term
      void add(int l, int mu, double val) { Ylm_[l][mu] += val; }
      // Get a specific multipole term
      double get(int l, int mu) { return Ylm_[l][mu]; }

      // Returns a reference of Ylm, to be computed by something else
      std::vector<std::vector<double>>& get_multipoles() { return Ylm_; }

      // Get an "internuclear" interaction tensor between two points separated by a distance R
      SharedVector build_T_spherical(int la, int lb, double distance);

      // Translate the solid harmonics
      std::shared_ptr<RealSolidHarmonics> translate(const Vector3& new_center);
      // Calculate the far-field potential generated by this multipole expansion.
      std::shared_ptr<RealSolidHarmonics> far_field_vector(const Vector3& far_center);

}; // End RealSolidHarmonics class

} // namespace psi

#endif  // PSI4_LIBFMM_MULTIPOLES_HELPER_H
