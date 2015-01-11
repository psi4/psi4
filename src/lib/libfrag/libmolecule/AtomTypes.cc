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

#include "AtomTypes.h"

namespace psi{
namespace LibMolecule{

DummyAtom::DummyAtom(const Carts_t& Carts, const bool IsBohr):
      Atom(Carts,0,IsBohr){}

PointCharge::PointCharge(const Carts_t& Carts,
      const double Charge,const bool IsBohr):
            Atom(Carts,0,IsBohr,"q",Charge){}

GhostAtom::GhostAtom(const Carts_t& Carts,const int Z,
      const bool IsBohr):Atom(Carts,-Z,IsBohr){}


}}//End namespaces

