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
#ifndef SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_FRAGMENTER_H_
#define SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_FRAGMENTER_H_
#include <vector>
#include <boost/shared_ptr.hpp>
#include "Geometry.h"
#include "OrganicGeom.h"
namespace psi{
namespace LibMolecule{
class Molecule;
class Fragment;
class Fragmenter{
   protected:
      boost::shared_ptr<const Molecule> Mol_;
      std::vector<boost::shared_ptr<Fragment> > FoundFrags_;
   public:
      Fragmenter(boost::shared_ptr<const Molecule> Mol):Mol_(Mol){}
      virtual ~Fragmenter(){}
      virtual std::vector<boost::shared_ptr<Fragment> > MakeFrags()=0;
};

class BondFragmenter:public Fragmenter{
   private:
      int NBonds_;
      boost::shared_ptr<OrganicGeom> Geom_;
      void Recurse(
        std::vector<boost::shared_ptr<const FxnalGroup> >& FoundGroups,
          const Connections& Conns,const ConnGroups& FxnGroups,
          long int& value);
      ///Makes fragment number "value" comprised of FoundGroups
      void AddFragment(const std::vector<
                              boost::shared_ptr<const FxnalGroup> >&
                            FoundGroups,const long int value);
   public:
      std::vector<boost::shared_ptr<Fragment> > MakeFrags();
      BondFragmenter(boost::shared_ptr<const Molecule> Mol, const int NBonds=3);
};

}}//End namespaces

#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_FRAGMENTER_H_ */
