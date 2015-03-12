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
#ifndef SRC_LIB_LIBFRAG_LIBMOLECULE_MOLECULETYPES_H_
#define SRC_LIB_LIBFRAG_LIBMOLECULE_MOLECULETYPES_H_

#include "LibFragMolecule.h"
namespace psi{
namespace LibMolecule{

/** \brief A specialization of Molecule that removes all ghosts,
 *          point charges, etc.
 *
 *   This is the "real" projection of a molecule
 */

class RealProjMolecule:public Molecule{
   public:
      RealProjMolecule(const Molecule& other);
      RealProjMolecule(const RealProjMolecule& other):Molecule(other){}
      const RealProjMolecule& operator=(const RealProjMolecule& other){
         Molecule::operator=(other);
         return *this;
      }
};

/** \brief A specialization that reorients the molecule into standard nuclear
 *         orientation.
 *
 *  Given a set of \f$N\f$ atoms, the \f$i\f$-th one with nuclear charge
 *  \f$Z_i\f$, and \f$j\f$-th Cartesian component, \f$X_{ij}\f$,
 *  Standard nuclear orientation is defined by:
 *  \f{eqnarray}{
 *     X^{SNO}=&X^\prime U^\dagger\\
 *     X^\prime_{ij}=&\left(X_{ij}-T_j\right)\\
 *     T_{j}=&\sum_{k=1}^N\frac{Z_kX_{kj}}{Z_k}\\
 *     U:\ I^\prime=&UIU^\dagger\\
 *     I^\prime_{ij}=&\omega_{ij}\delta_{i,j},
 *  \f}
 *  and
 *  \f[I=\left(
 *   \begin{array}{ccc}
 *  \sum_{i=1}^N\left({X_{i2}}^2+{X_{i3}}^2\right)Z_i &
 *  -\sum_{i=1}^NX_{i1}X_{i2}Z_i & -\sum_{i=1}^NX_{i1}X_{i3}Z_i\\
 *  -\sum_{i=1}^NX_{i1}X_{i2}Z_i &
 *   \sum_{i=1}^N\left({X_{i1}}^2+{X_{i3}}^2\right)Z_i &
 *  -\sum_{i=1}^NX_{i2}X_{i3}Z_i\\
 *  -\sum_{i=1}^NX_{i1}X_{i3}Z_i &
 *  -\sum_{i=1}^NX_{i2}X_{i3}Z_i  &
 *  \sum_{i=1}^N\left({X_{i1}}^2+{X_{i2}}^2\right)Z_i
 *  \end{array}\right),
 *  \f]
 *  That is to say the origin of SNO is set to the center
 *  of nuclear charge, and then rotated such that the x,y, and z axes lie
 *  along the principle moments of charge.  The principle moments of charge
 *  ordering depends on the \f$\omega_{i}\f$.
 *  - If the three \f$\omega_{i}\f$ are all different the largest lies along
 *    z, the second largest on y, and the smallest on x
 *  - If two are the same, the unique axis is z.  The atoms must then lie
 *    symmetrically around the Z axis in "circular sets".  Note circular
 *    sets contain atoms that are equivalent by spatial symmetry.  So
 *    if you have say 4 atoms lying in the XY plane like:
 *    <pre>
 *           A
 *           |
 *           |
 *    B---<Z-axis>----B
 *           |
 *           |
 *           A
 *    </pre>
 *    you have two circular sets, one with atoms A, and one with B. The
 *    first circular set to be distinguished from the others, by the following
 *    criteria, is then the "key circular set" and an arbitrary member is
 *    chosen to lie in the YZ plane.
 *      - Smallest absolute Z component
 *      - Smallest, positive Z component
 *      - Smallest distance to the Z axis
 *      - Lowest atomic number
 *  - For symmetric molecules
 *  - For linear molecules the Z-axis is the molecular axis
 *
 */
class SNOMolecule:public Molecule{
   public:
      ///Orients a molecule to SNO subject to the tolerance given
      SNOMolecule(const Molecule& other,const double Tolerance=1e-4);
      SNOMolecule(const SNOMolecule& other):Molecule(other){}
      const SNOMolecule& operator=(const SNOMolecule& other){
         Molecule::operator=(other);
         return *this;
      }
};

}}//End namespaces



#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_MOLECULETYPES_H_ */
