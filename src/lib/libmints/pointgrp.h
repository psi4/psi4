//
// pointgrp.h
//
// Additional modifications made by Justin Turney <jturney@ccqc.uga.edu>
// for use in PSI4.
//
// Modifications are
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Edward Seidl <seidl@janed.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

/* pointgrp.h -- definition of the point group classes
 *
 *      THIS SOFTWARE FITS THE DESCRIPTION IN THE U.S. COPYRIGHT ACT OF A
 *      "UNITED STATES GOVERNMENT WORK".  IT WAS WRITTEN AS A PART OF THE
 *      AUTHOR'S OFFICIAL DUTIES AS A GOVERNMENT EMPLOYEE.  THIS MEANS IT
 *      CANNOT BE COPYRIGHTED.  THIS SOFTWARE IS FREELY AVAILABLE TO THE
 *      PUBLIC FOR USE WITHOUT A COPYRIGHT NOTICE, AND THERE ARE NO
 *      RESTRICTIONS ON ITS USE, NOW OR SUBSEQUENTLY.
 *
 *  Author:
 *      E. T. Seidl
 *      Bldg. 12A, Rm. 2033
 *      Computer Systems Laboratory
 *      Division of Computer Research and Technology
 *      National Institutes of Health
 *      Bethesda, Maryland 20892
 *      Internet: seidl@alw.nih.gov
 *      June, 1993
 */

#ifndef _math_symmetry_pointgrp_h
#define _math_symmetry_pointgrp_h

#include <psi4-dec.h>
#include <libmints/vector3.h>
#include <string>
#include <cstdio>
#include <stdint.h>
#include <vector>
#include <map>

namespace psi {

extern FILE *outfile;

#define NUM_TO_OPERATOR_ID(x) ((x) ? 1<<((x)-1) : 0)
#define SKIP_THIS_OPERATOR(num,bit) ((bit) ? !((1<<((bit)-1)) & (num)) : 0 )

namespace SymmOps {
   enum Operations { E = 0, C2_z = 1, C2_y = 2, C2_x = 4, i = 8, Sigma_xy = 16, Sigma_xz = 32, Sigma_yz = 64, ID = 128 };
}

namespace PointGroups {
enum Groups {
    C1    = SymmOps::E,
    Ci    = SymmOps::E | SymmOps::i,
    C2X   = SymmOps::E | SymmOps::C2_x ,
    C2Y   = SymmOps::E | SymmOps::C2_y ,
    C2Z   = SymmOps::E | SymmOps::C2_z ,
    CsZ   = SymmOps::E | SymmOps::Sigma_xy ,
    CsY   = SymmOps::E | SymmOps::Sigma_xz ,
    CsX   = SymmOps::E | SymmOps::Sigma_yz ,
    D2    = SymmOps::E | SymmOps::C2_x | SymmOps::C2_y | SymmOps::C2_z ,
    C2vX  = SymmOps::E | SymmOps::C2_x | SymmOps::Sigma_xy | SymmOps::Sigma_xz ,
    C2vY  = SymmOps::E | SymmOps::C2_y | SymmOps::Sigma_xy | SymmOps::Sigma_yz ,
    C2vZ  = SymmOps::E | SymmOps::C2_z | SymmOps::Sigma_xz | SymmOps::Sigma_yz ,
    C2hX  = SymmOps::E | SymmOps::C2_x | SymmOps::Sigma_yz | SymmOps::i ,
    C2hY  = SymmOps::E | SymmOps::C2_y | SymmOps::Sigma_xz | SymmOps::i ,
    C2hZ  = SymmOps::E | SymmOps::C2_z | SymmOps::Sigma_xy | SymmOps::i ,
    D2h   = SymmOps::E | SymmOps::C2_x | SymmOps::C2_y | SymmOps::C2_z | SymmOps::i |
            SymmOps::Sigma_xy | SymmOps::Sigma_xz | SymmOps::Sigma_yz
};

void similar(unsigned char bits, unsigned char* sim, char& cnt);

}

// //////////////////////////////////////////////////////////////////

/** The SymmetryOperation class provides a 3 by 3 matrix
    representation of a symmetry operation, such as a rotation or reflection.
*/
class SymmetryOperation {
  private:
    double d[3][3];
    unsigned short bits_;

    void analyze_d();

  public:
    SymmetryOperation();
    SymmetryOperation(const SymmetryOperation &);
    ~SymmetryOperation();

    SymmetryOperation & operator = (SymmetryOperation const & a); // Assignment operator

    /// returns the trace of the transformation matrix
    double trace() const { return d[0][0]+d[1][1]+d[2][2]; }

    /// returns the i'th row of the transformation matrix
    double* operator[](int i) { return d[i]; }

    /// const version of the above
    const double* operator[](int i) const { return d[i]; }

    /** returns a reference to the (i,j)th element of the transformation
        matrix */
    double& operator()(int i, int j) { return d[i][j]; }

    /// const version of the above
    double operator()(int i, int j) const { return d[i][j]; }

    /// zero out the symop
    void zero() { memset(d,0,sizeof(double)*9); }

    /// This operates on this with r (i.e. return r * this).
    SymmetryOperation operate(const SymmetryOperation& r) const;

    /// This performs the transform r * this * r~
    SymmetryOperation transform(const SymmetryOperation& r) const;

    /// Get the bit value.
    unsigned char bit() const { return bits_; }

    /// Set equal to a unit matrix
    void unit() { zero(); d[0][0] = d[1][1] = d[2][2] = 1.0;}

    /// Set equal to E
    void E() { unit(); bits_ = SymmOps::E; }

    /// Set equal to an inversion
    void i() { zero(); d[0][0] = d[1][1] = d[2][2] = -1.0; bits_ = SymmOps::i; }

    /// Set equal to reflection in xy plane
    void sigma_xy() { unit(); d[2][2] = -1.0; bits_ = SymmOps::Sigma_xy; }

    /// Set equal to reflection in xz plane
    void sigma_xz() { unit(); d[1][1] = -1.0; bits_ = SymmOps::Sigma_xz; }

    /// Set equal to reflection in yz plane
    void sigma_yz() { unit(); d[0][0] = -1.0; bits_ = SymmOps::Sigma_yz; }

    /// Set equal to a clockwise rotation by 2pi/n
    void rotation(int n);
    void rotation(double theta);

    /// Set equal to C2 about the x axis
    void c2_x() { i(); d[0][0] = 1.0; bits_ = SymmOps::C2_x; }

    /// Set equal to C2 about the y axis
    void c2_y() { i(); d[1][1] = 1.0; bits_ = SymmOps::C2_y; }

    /// Set equal to C2 about the z axis
    void c2_z() { i(); d[2][2] = 1.0; bits_ = SymmOps::C2_z; }

    void transpose();

    /// print the matrix
    void print(FILE *out);
};

// //////////////////////////////////////////////////////////////////

/** The SymRep class provides an n dimensional matrix representation of a
    symmetry operation, such as a rotation or reflection.  The trace of a
    SymRep can be used as the character for that symmetry operation.  d is
    hardwired to 5x5 since the H irrep in Ih is 5 dimensional.
*/
class SymRep {
  private:
    int n;
    double d[5][5];

  public:
    SymRep(int =0);
    SymRep(const SymmetryOperation&);
    ~SymRep();

    /// Cast to a SymmetryOperation.
    operator SymmetryOperation() const;

    /// returns the trace of the transformation matrix
    inline double trace() const;

    /// set the dimension of d
    void set_dim(int i) { n=i; }

    /// returns the i'th row of the transformation matrix
    double* operator[](int i) { return d[i]; }
    /// const version of the above
    const double* operator[](int i) const { return d[i]; }

    /** returns a reference to the (i,j)th element of the transformation
        matrix */
    double& operator()(int i, int j) { return d[i][j]; }
    /// const version of double& operator()(int i, int j)
    double operator()(int i, int j) const { return d[i][j]; }

    /// zero out the symop
    void zero() { memset(d,0,sizeof(double)*25); }

    /// This operates on this with r (i.e. return r * this).
    SymRep operate(const SymRep& r) const;

    /// This performs the transform r * this * r~
    SymRep transform(const SymRep& r) const;

    /// Set equal to a unit matrix
    void unit() {
      zero(); d[0][0] = d[1][1] = d[2][2] = d[3][3] = d[4][4] = 1.0;
    }

    /// Set equal to the identity
    void E() { unit(); }

    /// Set equal to an inversion
    void i() { zero(); d[0][0] = d[1][1] = d[2][2] = d[3][3] = d[4][4] = -1.0;}

    /// Set equal to reflection in xy plane
    void sigma_h();

    /// Set equal to reflection in xz plane
    void sigma_xz();

    /// Set equal to reflection in yz plane
    void sigma_yz();

    /// Set equal to a clockwise rotation by 2pi/n
    void rotation(int n);
    void rotation(double theta);

    /// Set equal to C2 about the x axis
    void c2_x();

    /// Set equal to C2 about the y axis
    void c2_y();

    /// Set equal to C2 about the z axis
    void c2_z();

    /// print the matrix
    // void print(std::ostream& =ExEnv::out0()) const;
};

inline double
SymRep::trace() const
{
  double r=0;
  for (int i=0; i < n; i++)
    r += d[i][i];
  return r;
}

// //////////////////////////////////////////////////////////////////


class CharacterTable;

/** The IrreducibleRepresentation class provides information associated
    with a particular irreducible representation of a point group.  This
    includes the Mulliken symbol for the irrep, the degeneracy of the
    irrep, the characters which represent the irrep, and the number of
    translations and rotations in the irrep.  The order of the point group
    is also provided (this is equal to the number of characters in an
    irrep).  */
class IrreducibleRepresentation {
  friend class CharacterTable;

  private:
    int g;         // the order of the group
    int degen;     // the degeneracy of the irrep
    int nrot_;     // the number of rotations in this irrep
    int ntrans_;   // the number of translations in this irrep
    int complex_;  // true if this irrep has a complex representation
    char *symb;    // mulliken symbol for this irrep
    char *csymb;    // mulliken symbol for this irrep w/o special characters

    SymRep *rep;   // representation matrices for the symops

  public:
    IrreducibleRepresentation();
    IrreducibleRepresentation(const IrreducibleRepresentation&);
    /** This constructor takes as arguments the order of the point group,
     the degeneracy of the irrep, and the Mulliken symbol of the irrep.
     The Mulliken symbol is copied internally. */
    IrreducibleRepresentation(int,int,const char*,const char* =0);

    ~IrreducibleRepresentation();

    IrreducibleRepresentation& operator=(const IrreducibleRepresentation&);

    /// Initialize the order, degeneracy, and Mulliken symbol of the irrep.
    void init(int =0, int =0, const char* =0, const char* =0);

    /// Returns the order of the group.
    int order() const { return g; }

    /// Returns the degeneracy of the irrep.
    int degeneracy() const { return degen; }

    /// Returns the value of complex_.
    int complex() const { return complex_; }

    /// Returns the number of projection operators for the irrep.
    int nproj() const { return degen*degen; }

    /// Returns the number of rotations associated with the irrep.
    int nrot() const { return nrot_; }

    /// Returns the number of translations associated with the irrep.
    int ntrans() const { return ntrans_; }

    /// Returns the Mulliken symbol for the irrep.
    const char * symbol() const { return symb; }

    /** Returns the Mulliken symbol for the irrep without special
        characters.
    */
    const char * symbol_ns() const { return (csymb?csymb:symb); }

    /** Returns the character for the i'th symmetry operation of the point
     group. */
    double character(int i) const {
      return complex_ ? 0.5*rep[i].trace() : rep[i].trace();
    }

    /// Returns the element (x1,x2) of the i'th representation matrix.
    double p(int x1, int x2, int i) const { return rep[i](x1,x2); }

    /** Returns the character for the d'th contribution to the i'th
     representation matrix. */
    double p(int d, int i) const {
      int dc=d/degen; int dr=d%degen;
      return rep[i](dr,dc);
    }

    /** This prints the irrep to the given file, or stdout if none is
     given.  The second argument is an optional string of spaces to offset
     by. */
     void print(FILE *out=outfile) const;
};

// ///////////////////////////////////////////////////////////
/** The CharacterTable class provides a workable character table
 for all of the non-cubic point groups.  While I have tried to match the
 ordering in Cotton's book, I don't guarantee that it is always followed.
 It shouldn't matter anyway.  Also note that I don't lump symmetry
 operations of the same class together.  For example, in C3v there are two
 distinct C3 rotations and 3 distinct reflections, each with a separate
 character.  Thus symop has 6 elements rather than the 3 you'll find in
 most published character tables. */
class CharacterTable {
    int nt;                              //< order of the princ rot axis
    PointGroups::Groups pg;              //< the class of the point group
    int nirrep_;                         //< the number of irreps in this pg
    IrreducibleRepresentation *gamma_;   //< an array of irreps
    SymmetryOperation *symop;            //< the matrices describing sym ops
    int *_inv;                           //< index of the inverse symop
    std::string symb;                    //< the Schoenflies symbol for the pg
    unsigned char bits_;                 //< Bitwise representation of the symmetry operations

    /// this fills in the irrep and symop arrays.
    int make_table();
    void common_init();

  public:
    CharacterTable();
    /** This constructor takes the Schoenflies symbol of a point group as
        input. */
    CharacterTable(const std::string&);
    /** This constructor takes the bitswise representation of a point group as
        input. */
    CharacterTable(unsigned char);

    CharacterTable(const CharacterTable&);
    ~CharacterTable();

    CharacterTable& operator=(const CharacterTable&);

    /// Returns the number of irreps.
    int nirrep() const { return nirrep_; }
    /// Returns the order of the point group
    int order() const { return nirrep_; }
    /// Returns the Schoenflies symbol for the point group
    const std::string& symbol() const { return symb; }
    /// Returns the i'th irrep.
    IrreducibleRepresentation& gamma(int i) { return gamma_[i]; }
    /// Returns the i'th symmetry operation.
    SymmetryOperation& symm_operation(int i) { return symop[i]; }

    /** Cn, Cnh, Sn, T, and Th point groups have complex representations.
        This function returns 1 if the point group has a complex
        representation, 0 otherwise. */
    int complex() const {
      return 0;
    }

    /// Returns the index of the symop which is the inverse of symop[i].
    int inverse(int i) const { return _inv[i]; }

    int ncomp() const {
      int ret=0;
      for (int i=0; i < nirrep_; i++) {
        int nc = (gamma_[i].complex()) ? 1 : gamma_[i].degen;
        ret += nc;
      }
      return ret;
    }

    /// Returns the irrep component i belongs to.
    int which_irrep(int i) {
      for (int ir=0, cn=0; ir < nirrep_; ir++) {
        int nc = (gamma_[ir].complex()) ? 1 : gamma_[ir].degen;
        for (int c=0; c < nc; c++,cn++)
          if (cn==i)
            return ir;
      }
      return -1;
    }

    /// Returns which component i is.
    int which_comp(int i) {
      for (int ir=0, cn=0; ir < nirrep_; ir++) {
        int nc = (gamma_[ir].complex()) ? 1 : gamma_[ir].degen;
        for (int c=0; c < nc; c++,cn++)
          if (cn==i)
            return c;
      }
      return -1;
    }

    unsigned char bits();

    /// This prints the irrep to the given file, or stdout if none is given.
     void print(FILE *out=outfile) const;
};

// ///////////////////////////////////////////////////////////

/** The PointGroup class is really a place holder for a CharacterTable.  It
 contains a string representation of the Schoenflies symbol of a point
 group, a frame of reference for the symmetry operation transformation
 matrices, and a point of origin.  The origin is not respected by the
 symmetry operations, so if you want to use a point group with a nonzero
 origin, first translate all your coordinates to the origin and then set
 the origin to zero.  */
class PointGroup {
  private:
    std::string symb;
    Vector3 origin_;
    unsigned char bits_;

  public:
    PointGroup();

    // These 2 constructors do not work right now.
    /** This constructor takes a string containing the Schoenflies symbol
        of the point group as its only argument. */
    PointGroup(const std::string&);
    /** Like the above, but this constructor also takes a point of origin
        as an argument. */
    PointGroup(const std::string&,const Vector3&);

    /** Using the bitwise representation constructor the point group object. */
    PointGroup(unsigned char bits);
    /** Using the bitwise representation constructor the point group object. */
    PointGroup(unsigned char bits, const Vector3&);

    /** The PointGroup KeyVal constructor looks for three keywords:
       symmetry, symmetry_frame, and origin. symmetry is a string
       containing the Schoenflies symbol of the point group.  origin is an
       array of doubles which gives the x, y, and z coordinates of the
       origin of the symmetry frame.  symmetry_frame is a 3 by 3 array of
       arrays of doubles which specify the principal axes for the
       transformation matrices as a unitary rotation.

       For example, a simple input which will use the default origin and
       symmetry_frame ((0,0,0) and the unit matrix, respectively), might
       look like this:

       <pre>
       pointgrp<PointGroup>: (
         symmetry = "c2v"
       )
       </pre>

       By default, the principal rotation axis is taken to be the z axis.
       If you already have a set of coordinates which assume that the
       rotation axis is the x axis, then you'll have to rotate your frame
       of reference with symmetry_frame:

       <pre>
       pointgrp<PointGroup>: (
         symmetry = "c2v"
         symmetry_frame = [
           [ 0 0 1 ]
           [ 0 1 0 ]
           [ 1 0 0 ]
         ]
       )
       </pre>
     */
    // PointGroup(const Ref<KeyVal>&);

    // PointGroup(StateIn&);
    PointGroup(const PointGroup&);
    PointGroup(const boost::shared_ptr<PointGroup>&);
    ~PointGroup();

    PointGroup& operator=(const PointGroup&);

    /// Returns 1 if the point groups are equivalent, 0 otherwise.
    int equiv(const boost::shared_ptr<PointGroup> &, double tol = 1.0e-6) const;

    /// Returns the CharacterTable for this point group.
    CharacterTable char_table() const;
    /// Returns the Schoenflies symbol for this point group.
    std::string symbol() const { return symb; }
    /// Returns the frame of reference for this point group.
    /// Returns the origin of the symmetry frame.
    Vector3& origin() { return origin_; }
    const Vector3& origin() const { return origin_; }

    /// Sets (or resets) the Schoenflies symbol.
    void set_symbol(const std::string&);

    /// Returns the bitwise representation of the point group
    unsigned char bits() const { return bits_; }

    static const char* bits_to_full_name(unsigned char bits);
    static const char* bits_to_basic_name(unsigned char bits);
    static bool full_name_to_bits(const std::string& pg, unsigned char& bits);

    void print(FILE *out = outfile) const;
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
