/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

#ifndef _psi_src_bin_psimrcc_updater_h_
#define _psi_src_bin_psimrcc_updater_h_

#include "psi4/liboptions/liboptions.h"

/**
 *  @file updater.h
 *  @ingroup (PSIMRCC)
*/

namespace psi{ namespace psimrcc{

class Hamiltonian;

/**
 *  @class Updater
 *  @brief Containts the procedure for updating the amplitudes
*/
class Updater{
public:
  Updater(Options &options);
  virtual ~Updater();
  virtual void update(int cycle,Hamiltonian* heff) = 0;
  void zero_internal_amps();
  void zero_t1_internal_amps();
  void zero_internal_delta_amps();
protected:
  Options &options_;
};

class MkUpdater : public Updater{
public:
  MkUpdater(Options &options);
  virtual ~MkUpdater();
  virtual void update(int cycle,Hamiltonian* heff);
};

class BWUpdater : public Updater{
public:
  BWUpdater(Options &options);
  virtual ~BWUpdater();
  virtual void update(int cycle,Hamiltonian* heff);
};

}} /* End Namespaces */

#endif // _psi_src_bin_psimrcc_updater_h_
