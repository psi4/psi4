/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
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

#ifndef _psi_src_lib_libmints_osrecur_h
#define _psi_src_lib_libmints_osrecur_h

#include "psi4/pragma.h"

namespace psi {

/*! \ingroup MINTS
 *  \class ObaraSaikaTwoCenterRecursion
 *  \brief Generic Obara and Saika recursion object.
 */
class ObaraSaikaTwoCenterRecursion {
    int max_am1_;
    int max_am2_;

    double **x_;
    double **y_;
    double **z_;

    // No default constructor
    ObaraSaikaTwoCenterRecursion();
    // No assignment operator
    ObaraSaikaTwoCenterRecursion &operator=(const ObaraSaikaTwoCenterRecursion &);

   public:
    /// Constructor, max_am1 and max_am2 are the max angular momentum on center 1 and 2.
    /// Needed to allocate enough memory.
    ObaraSaikaTwoCenterRecursion(int max_am1, int max_am2);
    ~ObaraSaikaTwoCenterRecursion();

    /// Returns the x recursion matrix.
    double **x() const { return x_; }
    /// Returns the y recursion matrix.
    double **y() const { return y_; }
    /// Returns the z recursion matrix.
    double **z() const { return z_; }

    /// Computes the recursion matrices for the data provided.
    void compute(double PA[3], double PB[3], double gamma, int am1, int am2);
};

/*! \ingroup MINTS
 *  \class ObaraSaikaTwoCenterMIRecursion
 *  \brief Obara and Saika recursion object for moment integrals. Currently not used by DipoleInt, hopefully soon.
 *  THIS CLASS HAS NOT BEEN TESTED!!!
 */
class ObaraSaikaTwoCenterMIRecursion {
    int max_am1_;
    int max_am2_;
    int max_m_;

    double ***x_;
    double ***y_;
    double ***z_;

    // No default constructor
    ObaraSaikaTwoCenterMIRecursion();
    // No assignment operator
    ObaraSaikaTwoCenterMIRecursion &operator=(const ObaraSaikaTwoCenterMIRecursion &);

   public:
    ObaraSaikaTwoCenterMIRecursion(int max_am1, int max_am2, int max_m);
    ~ObaraSaikaTwoCenterMIRecursion();

    double ***x() const { return x_; }
    double ***y() const { return y_; }
    double ***z() const { return z_; }
    void compute(double PA[3], double PB[3], double gamma, int am1, int am2);
};

/*! \ingroup MINTS
 *  \class ObaraSaikaTwoCenterVIRecursion
 *  \brief Obara and Saika recursion object for potential integrals.
 */
class ObaraSaikaTwoCenterVIRecursion {
   protected:
    int max_am1_;
    int max_am2_;
    int size_;

    double ***vi_;

    // Forms Fm(U) from A20 (OS 1986)
    void calculate_f(double *F, int n, double t);

   private:
    // No default constructor
    ObaraSaikaTwoCenterVIRecursion();
    // No assignment operator
    ObaraSaikaTwoCenterVIRecursion &operator=(const ObaraSaikaTwoCenterVIRecursion &);

   public:
    /// Constructor, max_am1 and max_am2 are the max angular momentum on center 1 and 2.
    /// Needed to allocate enough memory.
    ObaraSaikaTwoCenterVIRecursion(int max_am1, int max_am2);
    virtual ~ObaraSaikaTwoCenterVIRecursion();

    /// Returns the potential integral 3D matrix
    double ***vi() const { return vi_; }

    virtual double ***vx() const { return nullptr; }
    virtual double ***vy() const { return nullptr; }
    virtual double ***vz() const { return nullptr; }

    virtual double ***vxx() const { return nullptr; }
    virtual double ***vxy() const { return nullptr; }
    virtual double ***vxz() const { return nullptr; }
    virtual double ***vyy() const { return nullptr; }
    virtual double ***vyz() const { return nullptr; }
    virtual double ***vzz() const { return nullptr; }

    /// Computes the potential integral 3D matrix using the data provided.
    virtual void compute(double PA[3], double PB[3], double PC[3], double zeta, int am1, int am2);
    /// Computes the Ewald potential integral with modified zeta -> zetam 3D matrix using the data provided.
    virtual void compute_erf(double PA[3], double PB[3], double PC[3], double zeta, int am1, int am2, double zetam);
};

/*! \ingroup MINTS
 *  \class ObaraSaikaTwoCenterVIDerivRecursion
 *  \brief Obara and Saika recursion object for computing potential derivatives.
 */
class ObaraSaikaTwoCenterVIDerivRecursion : public ObaraSaikaTwoCenterVIRecursion {
   protected:
    double ***vx_;
    double ***vy_;
    double ***vz_;

   private:
    // No default constructor();
    ObaraSaikaTwoCenterVIDerivRecursion();
    // No assignment operator
    ObaraSaikaTwoCenterVIDerivRecursion &operator=(const ObaraSaikaTwoCenterVIDerivRecursion &);

   public:
    ObaraSaikaTwoCenterVIDerivRecursion(int max_am1, int max_am2);
    ~ObaraSaikaTwoCenterVIDerivRecursion() override;

    double ***vx() const override { return vx_; }
    double ***vy() const override { return vy_; }
    double ***vz() const override { return vz_; }

    void compute(double PA[3], double PB[3], double PC[3], double zeta, int am1, int am2) override;
};

/*! \ingroup MINTS
 *  \class ObaraSaikaTwoCenterVIDerivRecursion
 *  \brief Obara and Saika recursion object for computing potential derivatives.
 */
class ObaraSaikaTwoCenterVIDeriv2Recursion : public ObaraSaikaTwoCenterVIDerivRecursion {
   protected:
    double ***vxx_;
    double ***vxy_;
    double ***vxz_;
    double ***vyy_;
    double ***vyz_;
    double ***vzz_;

   private:
    // No default constructor();
    ObaraSaikaTwoCenterVIDeriv2Recursion();
    // No assignment operator
    ObaraSaikaTwoCenterVIDeriv2Recursion &operator=(const ObaraSaikaTwoCenterVIDeriv2Recursion &);

   public:
    ObaraSaikaTwoCenterVIDeriv2Recursion(int max_am1, int max_am2);
    ~ObaraSaikaTwoCenterVIDeriv2Recursion() override;

    double ***vxx() const override { return vxx_; }
    double ***vxy() const override { return vxy_; }
    double ***vxz() const override { return vxz_; }
    double ***vyy() const override { return vyy_; }
    double ***vyz() const override { return vyz_; }
    double ***vzz() const override { return vzz_; }

    void compute(double PA[3], double PB[3], double PC[3], double zeta, int am1, int am2) override;
};

/*! \ingroup MINTS
 *  \class ObaraSaikaTwoCenterElectricField
 *  \brief Obara and Saika recursion object for computing electric field integrals.
 */
class ObaraSaikaTwoCenterElectricField : public ObaraSaikaTwoCenterVIRecursion {
   protected:
    double ***q_;
    double ***x_;
    double ***y_;
    double ***z_;

   private:
    // No default constructor
    ObaraSaikaTwoCenterElectricField();
    // No assignment operator
    ObaraSaikaTwoCenterElectricField &operator=(const ObaraSaikaTwoCenterElectricField &);

   public:
    ObaraSaikaTwoCenterElectricField(int max_am1, int max_am2);
    ~ObaraSaikaTwoCenterElectricField() override;

    // We could also add the getter for the q_ (potential ints here)
    double ***x() const { return x_; }
    double ***y() const { return y_; }
    double ***z() const { return z_; }

    void compute(double PA[3], double PB[3], double PC[3], double zeta, int am1, int am2) override;
};

/*! \ingroup MINTS
 *  \class ObaraSaikaTwoCenterElectricFieldGradient
 *  \brief Obara and Saika recursion object for computing electric field gradient integrals.
 */
class ObaraSaikaTwoCenterElectricFieldGradient : public ObaraSaikaTwoCenterElectricField {
   protected:
    double ***exx_;
    double ***eyy_;
    double ***ezz_;
    double ***exy_;
    double ***exz_;
    double ***eyz_;

   private:
    // No default constructor
    ObaraSaikaTwoCenterElectricFieldGradient();
    // No assignment operator
    ObaraSaikaTwoCenterElectricFieldGradient &operator=(const ObaraSaikaTwoCenterElectricFieldGradient &);

   public:
    ObaraSaikaTwoCenterElectricFieldGradient(int max_am1, int max_am2);
    ~ObaraSaikaTwoCenterElectricFieldGradient() override;

    double ***ex() const { return x_; }
    double ***ey() const { return y_; }
    double ***ez() const { return z_; }

    void compute(double PA[3], double PB[3], double PC[3], double zeta, int am1, int am2) override;
};

/*! \ingroup MINTS
 *  \class ObaraSaikaTwoCenterEFPRecursion
 *  \brief Obara and Saika recursion object for potential integrals, and electric derivatives thereof, required in EFP.
 */
class ObaraSaikaTwoCenterEFPRecursion {
   protected:
    int max_am1_;
    int max_am2_;
    int size_;

    double ***q_;
    double ***x_;
    double ***y_;
    double ***z_;
    double ***xx_;
    double ***xy_;
    double ***xz_;
    double ***yy_;
    double ***yz_;
    double ***zz_;
    double ***xxx_;
    double ***xxy_;
    double ***xxz_;
    double ***xyy_;
    double ***xyz_;
    double ***xzz_;
    double ***yyy_;
    double ***yyz_;
    double ***yzz_;
    double ***zzz_;

    // Forms Fm(U) from A20 (OS 1986)
    void calculate_f(double *F, int n, double t);

   private:
    // No default constructor
    ObaraSaikaTwoCenterEFPRecursion();
    // No assignment operator
    ObaraSaikaTwoCenterEFPRecursion &operator=(const ObaraSaikaTwoCenterEFPRecursion &);

   public:
    /// Constructor, max_am1 and max_am2 are the max angular momentum on center 1 and 2.
    /// Needed to allocate enough memory.
    ObaraSaikaTwoCenterEFPRecursion(int max_am1, int max_am2);
    virtual ~ObaraSaikaTwoCenterEFPRecursion();

    /// Returns the potential integral 3D matrix
    double ***q() const { return q_; }
    double ***x() const { return x_; }
    double ***y() const { return y_; }
    double ***z() const { return z_; }
    double ***xx() const { return xx_; }
    double ***yy() const { return yy_; }
    double ***zz() const { return zz_; }
    double ***xy() const { return xy_; }
    double ***xz() const { return xz_; }
    double ***yz() const { return yz_; }
    double ***xxx() const { return xxx_; }
    double ***yyy() const { return yyy_; }
    double ***zzz() const { return zzz_; }
    double ***xxy() const { return xxy_; }
    double ***xxz() const { return xxz_; }
    double ***xyy() const { return xyy_; }
    double ***yyz() const { return yyz_; }
    double ***xzz() const { return xzz_; }
    double ***yzz() const { return yzz_; }
    double ***xyz() const { return xyz_; }

    /// Computes the potential integral 3D matrix using the data provided.
    virtual void compute(double PA[3], double PB[3], double PC[3], double zeta, int am1, int am2);
};

/*! \ingroup MINTS
    \class ObaraSaikaThreeCenterRecursion
 */
class ObaraSaikaThreeCenterRecursion {
    int max_am1_;
    int max_am2_;
    int max_am3_;

    double ***x_;
    double ***y_;
    double ***z_;

    // No default constructor
    ObaraSaikaThreeCenterRecursion();

    // No assignment operator
    ObaraSaikaThreeCenterRecursion &operator=(const ObaraSaikaThreeCenterRecursion &);

   public:
    ObaraSaikaThreeCenterRecursion(int max_am1, int max_am2, int max_am3);
    ~ObaraSaikaThreeCenterRecursion();

    /// Returns the x recursion matrix.
    double ***x() const { return x_; }
    /// Returns the y recursion matrix.
    double ***y() const { return y_; }
    /// Returns the z recursion matrix.
    double ***z() const { return z_; }

    /// Computes the recursion matrices for the data provided.
    void compute(double GA[3], double GB[3], double GC[3], double gammac, int amA, int amB, int amC);
};

}  // namespace psi

#endif
