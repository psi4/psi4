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
#ifndef SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_RINGNODETYPES_HH_
#define SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_RINGNODETYPES_HH_
#include "RingFinder.h"

namespace psi{
namespace LibMolecule{

class Benzene: public RingFinder<DBCC2,DBCC2,DBCC2>{
   private:
      typedef RingFinder<DBCC2,DBCC2,DBCC2> Base_t;
   public:
      Benzene():Base_t(F_t::Bz0,F_t::CBz0,F_t::HBz0,F_t::CBz0,F_t::HBz0,
                      F_t::CBz0,F_t::HBz0,F_t::CBz0,F_t::HBz0,F_t::CBz0,
                      F_t::HBz0,F_t::CBz0,F_t::HBz0){}
};

class Phenyl: public RingFinder<Alkenyl3,Alkenyl2,Alkenyl2,
                                Alkenyl2,Alkenyl2,Alkenyl2>{
   private:
      typedef RingFinder<Alkenyl3,Alkenyl2,Alkenyl2,
                         Alkenyl2,Alkenyl2,Alkenyl2> Base_t;
   public:
      Phenyl():Base_t(F_t::Bz1,
                F_t::C1Bz1,F_t::C2Bz1,F_t::H1Bz1,
                F_t::C3Bz1,F_t::H2Bz1,F_t::C4Bz1,F_t::H3Bz1,
                F_t::C3Bz1,F_t::H2Bz1,F_t::C2Bz1,F_t::H1Bz1){}
};

class OrthoBenzene: public RingFinder<Alkenyl3,Alkenyl3,Alkenyl2,
                                      Alkenyl2,Alkenyl2,Alkenyl2>{
      private:
      typedef RingFinder<Alkenyl3,Alkenyl3,Alkenyl2,
                         Alkenyl2,Alkenyl2,Alkenyl2> Base_t;
      public:
      OrthoBenzene():Base_t(F_t::Bz2o,
            F_t::C1Bz2o,F_t::C1Bz2o,
            F_t::C2Bz2o,F_t::H1Bz2o,F_t::C3Bz2o,F_t::H2Bz2o,
            F_t::C2Bz2o,F_t::H1Bz2o,F_t::C3Bz2o,F_t::H2Bz2o){}
};

class MetaBenzene: public RingFinder<Alkenyl3,Alkenyl2,Alkenyl3,
                                      Alkenyl2,Alkenyl2,Alkenyl2>{
      private:
      typedef RingFinder<Alkenyl3,Alkenyl2,Alkenyl3,
                         Alkenyl2,Alkenyl2,Alkenyl2> Base_t;
      public:
      MetaBenzene():Base_t(F_t::Bz2m,
            F_t::C1Bz2m,F_t::C2Bz2m,F_t::H1Bz2m,
            F_t::C1Bz2m,F_t::C3Bz2m,F_t::H2Bz2m,
            F_t::C4Bz2m,F_t::H3Bz2m,F_t::C3Bz2m,F_t::H2Bz2m){}
};

class ParaBenzene: public RingFinder<Alkenyl3,Alkenyl2,Alkenyl2,
                                      Alkenyl3,Alkenyl2,Alkenyl2>{
      private:
      typedef RingFinder<Alkenyl3,Alkenyl2,Alkenyl2,
                         Alkenyl3,Alkenyl2,Alkenyl2> Base_t;
      public:
      ParaBenzene():Base_t(F_t::Bz2p,
            F_t::C1Bz2p,F_t::C2Bz2p,F_t::HBz2p,
            F_t::C2Bz2p,F_t::HBz2p,F_t::C1Bz2p,
            F_t::C2Bz2p,F_t::HBz2p,F_t::C2Bz2p,F_t::HBz2p){}
};

class Benzene6: public RingFinder<DBCC4,DBCC4,DBCC4,DBCC4,DBCC4,DBCC4>{
   private:
      typedef RingFinder<DBCC4,DBCC4,DBCC4,DBCC4,DBCC4,DBCC4> Base_t;
   public:
      Benzene6():Base_t(F_t::Bz6,F_t::CBz6,F_t::CBz6,F_t::CBz6,
                        F_t::CBz6,F_t::CBz6,F_t::CBz6){}
};

//This header file is all 256 possible substitutions of an Indole ring
#include "Indole.hh"
}}//End namespaces

#endif /* SRC_LIB_LIBMOLECULE_AUTOFXNALGROUP_RINGNODETYPES_HH_ */
