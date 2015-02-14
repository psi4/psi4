/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */
#ifndef SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_ORIENTATION_H_
#define SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_ORIENTATION_H_

#include "LibMoleculeBase.h"
#include "Molecule.h"
#include "GeomManipulator.h"
namespace psi{
namespace LibMolecule{
/** \brief The base class for defining a way to orient a molecule
 *
 *   Every orientation can be thought of as involving three main
 *   steps:
 *   - Find an origin
 *   - Define the Z-axis
 *   - Define the +/- directions on each axis
 *   - Rotate into either the XZ or the YZ plane
 *
 *   Of the above the list, the latter two may require comparing double
 *   precision values, and so we require the user to give us a threshold
 *   value to determine if two doubles are equivalent
 */
class Orientation:public LibMoleculeBase{
   protected:
      ///The translation that takes us to the origin
      double Trans2Origin_[3];
      ///The overall rotation
      double Rotation_[9];
      ///The value of the tolerance
      double Tolerance_;
      ///Derived classes should make this function set-up Trans2Origin
      virtual void CalcOrigin(const Molecule& Mol)=0;
      /** \brief Function that defines how to find the z-axis
       *
       *  This function should set the rotation matrix up, such that
       *  when it is transposed and right-applied to the molecule
       *  the molecule's Z-axis is aligned (the new orientation should
       *  be on the rows, and the old-orientation on the columns)
       */
      virtual void SetZAxis(const Molecule& Mol,
                     const GeomManipulator& Manip)=0;
      ///Fxn to set Rotation such that the proper axis parity is obtained
      virtual void SetDirection(const Molecule& Mol,
                                const GeomManipulator& Manip,
                                double* Rotation)=0;
      ///Defines the angle to rotate by to align XZ or YZ plane
      virtual double GetTheta(const Molecule& Mol,
                              const GeomManipulator& Manip)const=0;
   public:
      ///No memory to free up
      virtual ~Orientation(){}
      void TranslateToOrigin(GeomManipulator& Manip)const;
      Orientation(double Tolerance):Tolerance_(Tolerance){
         memset(Trans2Origin_,0.0,3*sizeof(double));
         memset(Rotation_,0.0,9*sizeof(double));
      }
      void Reorient(Molecule& Mol);
};

/** The base class for determining orientations from eigenvectors
 *
 *  I couldn't think of a better description, but basically the idea is
 *  that this class works by diagonalizing a matrix to get a set of
 *  eigenvectors, which in turn are used as the coordinate frame.  The
 *  most common example is the reference frame obtained from diagonalizing
 *  the inertia tensor.
 */
class EigenOrientation:public Orientation{
   protected:
      ///The eigen values of the rotation
      double Omega_[3];
      ///Returns the value for atom i, that we are using to weight
      virtual double Value(const Atom* I)const=0;
      ///Computes Trans2Origin_
      void CalcOrigin(const Molecule& Mol);
      ///Computes Rotation so that it aligns the Z-axis
      void SetZAxis(const Molecule& Mol, const GeomManipulator& Manip);
      ///Returns the angle necessary to align the plane of choice
      double GetTheta(const Molecule& Mol,const GeomManipulator& Manip)const;
      ///Sets + direction such that the larger 4 moment
      void SetDirection(const Molecule &Mol,
            const GeomManipulator& Manip,double* Rotation);

   public:
      EigenOrientation(const double Tol):Orientation(Tol){}
};

class MassOrientation:public EigenOrientation{
   protected:
      double Value(const Atom* I)const;
   public:
      MassOrientation(const double Tol=1e-4):EigenOrientation(Tol){}
};

class ChargeOrientation:public EigenOrientation{
   protected:
      double Value(const Atom* I)const;
   public:
      ChargeOrientation(const double Tol=1e-4):EigenOrientation(Tol){}
};

}}



#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_ORIENTATION_H_ */
