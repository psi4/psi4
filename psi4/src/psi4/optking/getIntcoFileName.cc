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

#include <string>

#include "package.h"

#if defined(OPTKING_PACKAGE_PSI)
  #include "psi4/psi4-dec.h"
  #include "psi4/libmints/writer_file_prefix.h"
  #include "psi4/libmints/molecule.h"

#endif

#if defined(OPTKING_PACKAGE_QCHEM)
 #include "qcsys.h"
 extern void getPrefix(char*& pref);
#endif

const char* getIntcoFileName()
{
   static std::string strintco("");
   if (strintco.empty() ) {
#if defined(OPTKING_PACKAGE_QCHEM)
      char* pref=NULL;
      getPrefix(pref);
      strintco = std::string(pref) + "intco.dat";
#elif defined(OPTKING_PACKAGE_PSI)
      strintco = psi::get_writer_file_prefix(psi::Process::environment.legacy_molecule()->name()) + ".intco";
#endif
   }
   return strintco.c_str();
}

const char* getOptdataFileName()
{
   static std::string stroptdata("");
   if (stroptdata.empty() ) {
#if defined(OPTKING_PACKAGE_QCHEM)
      char* pref=NULL;
      getPrefix(pref);
      stroptdata = std::string(pref) + "opt_data.1";
#elif defined(OPTKING_PACKAGE_PSI)
      // In PSI, the opt data file is file 1, with the name set by input and io_start
#endif
   }
   return stroptdata.c_str();
}
