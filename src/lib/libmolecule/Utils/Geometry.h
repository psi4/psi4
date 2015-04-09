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
#ifndef SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_GEOMETRY_H_
#define SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_GEOMETRY_H_
#include <vector>
#include "LibMoleculeBase.h"
#include "Parameter.h"

namespace psi{
namespace LibMolecule{
class Molecule;
enum GeomQuant{BONDS,ANGLES,PAIRS,TORSIONS,FXNGROUPS};

class Connections:public std::vector<std::vector<int> >{
   private:
      typedef std::vector<std::vector<int> > Base_t;
   public:
      std::string PrintOut()const;
      Connections(const int NAtoms):Base_t(NAtoms){}
};


class Geometry:public LibMoleculeBase{
   private:
      void FormDistance();
      std::vector<boost::shared_ptr<const Pair> > Pairs_;
      std::vector<boost::shared_ptr<const Bond> > Bonds_;

   protected:
      Connections Connections_;
      const Molecule& Mol_;
      boost::shared_ptr<double[]> Distance_;
   public:
      ///Returns the number of bonds atom i makes
      int NBonds(const int i)const{return Connections_[i].size();}
      typedef const std::vector<boost::shared_ptr<const Bond> > Bond_t;
      typedef const std::vector<boost::shared_ptr<const Pair> > Pair_t;
      Bond_t& GetBonds()const{return Bonds_;}
      Pair_t& GetPairs()const{return Pairs_;}
      const Connections& GetConns()const{return Connections_;}
      virtual ~Geometry(){}
      Geometry(const Molecule& Mol);
};



}}//End namespaces



#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_UTILS_GEOMETRY_H_ */
