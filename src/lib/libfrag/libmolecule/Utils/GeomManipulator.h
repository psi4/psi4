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
#ifndef SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_GEOMMANIPULATOR_H_
#define SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_GEOMMANIPULATOR_H_
#include <boost/shared_ptr.hpp>
#include "LibMoleculeBase.h"
namespace psi{
namespace LibMolecule{
class Molecule;

enum Plane{XY,XZ,YZ};
/** \brief A class that can manipulate the Cartesian coordinates of a
 *         molecule.  Once an atom is made, this is the only class that
 *         can change its geometry.
 *
 *  For our intents and purposes the Cartesian coordinates form a matrix,
 *  \f$X\f$, with elements such that \f$X_{ij}\f$ is the \f$j\f$-th
 *  Cartesian component of atom \f$i\f$.  Given a translation vector, \f$T\f$,
 *  our new coordinates, \f$X^\prime\f$, are given by:
 *  \f[
 *  X^{\prime}_{ij}=X_{ij}+T_{j}
 *  \f]
 *  and given a rotation matrix \f$U: v^\prime=Uv\f$, where \f$v^\prime\f$ is
 *  our new vector, and \f$v\f$ is our original vector, a rotation of all the
 *  coordinates is given by:
 *  \f[
 *  X^{\prime}=XU^\dagger
 *  \f]
 */
class GeomManipulator: protected LibMoleculeBase{
   private:
      ///Function that plucks the cartesians from the molecule
      boost::shared_ptr<double[]> GetCarts()const;

      ///Function for writing the Cartesians to the molecule
      void SetCarts(boost::shared_ptr<double[]> NewCarts);

      ///The Cartesians we are manipulating
      boost::shared_ptr<double[]> Carts_;

      ///The molecule this class is working on
      Molecule* Mol_;
   public:
      ///Gets ready to manipulate the molecule given by Mol
      GeomManipulator(Molecule* Mol);

      ///Gets cart j of atom i
      double operator()(const int i,const int j)const;

      ///Scales the carts
      void Scale(const double ScaleFac);

      ///Performs: X'=UX , Trans=true means you gave U not U^dagger
      void Rotate(const double* RotMat,const bool Trans=true);
      ///Rotates the molecule by an angle theta in the given plane
      void Rotate(const double Theta,const Plane RotPlan);

      ///Translates, and possibly scales all the carts
      void Translate(const double* TransVec, const double ScaleFac=1.0);

      ///Translates, and possibly scales the carts of atom i
      void Translate(const double* TransVec, const int i, const double ScaleFac=1.0){
         for(int j=0;j<3;j++)Carts_[i*3+j]+=ScaleFac*TransVec[j];
      }

      ///Call this when you are done manipulating the carts to set them
      void Set();
};


}}//End namespaces



#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_GEOMMANIPULATOR_H_ */
