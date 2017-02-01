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

/*! \file    binary_io.cc
    \ingroup optking
    \brief   functions to read and write to optking's binary file
*/

#include "io.h"

#include "print.h"
#define EXTERN
#include "globals.h"

// PSI unit number for opt_data binary file
#if defined (OPTKING_PACKAGE_PSI)
 #define PSI_OPTDATA_FILE_NUM 1
 #include "psi4/libparallel/parallel.h"
 #include "psi4/libpsio/psio.h"
 #include "psi4/libpsio/psio.hpp"
 using namespace psi;
#endif

// Name of opt_data file for QCHEM; PSI uses file 1 in the default name
const char* getOptdataFileName();
#define QCHEM_OPTDATA_FILENAME getOptdataFileName()

// binary QCHEM data file
#if defined (OPTKING_PACKAGE_QCHEM)
#include <fstream>
namespace opt_io {
  std::fstream opt_data_stream;
}
#endif

namespace opt {

// returns true if the binary file exists and is not empty
bool opt_io_is_present(void) {
  bool file_present = false;

#if defined(OPTKING_PACKAGE_PSI)
  psio_open(PSI_OPTDATA_FILE_NUM, PSIO_OPEN_OLD);
  if (psio_rd_toclen(PSI_OPTDATA_FILE_NUM) > 0)
    file_present = true;
  psio_close(PSI_OPTDATA_FILE_NUM, 1);

#elif defined(OPTKING_PACKAGE_QCHEM)
  using opt_io::opt_data_stream;
  using namespace std;

  opt_data_stream.open(QCHEM_OPTDATA_FILENAME, fstream::in | fstream::binary);
  if (opt_data_stream.is_open()) {
    if (opt_data_stream.good())
      file_present = true;
    opt_data_stream.close();
  }
#endif

  return file_present;
}

void opt_io_remove(void) {
#if defined(OPTKING_PACKAGE_PSI)
  // check retention setting in .psi4rc - maybe the user likes file 1 !
  if (! psi::_default_psio_manager_->get_specific_retention(1)) {
    if (!psio_open_check(PSI_OPTDATA_FILE_NUM)) // if not open, open it
      psio_open(PSI_OPTDATA_FILE_NUM, PSIO_OPEN_OLD);
    psio_close(PSI_OPTDATA_FILE_NUM, 0);        // close and delete it
  }

#elif defined(OPTKING_PACKAGE_QCHEM)
  using opt_io::opt_data_stream;

  if (opt_data_stream.is_open())       // if open, close it
    opt_data_stream.close();
  std::remove(QCHEM_OPTDATA_FILENAME); // remove file
#endif
}

void opt_intco_dat_remove(void) {
  std::remove(FILENAME_INTCO_DAT); // rm intco definitions
}

void opt_clean(void) {
  opt_io_remove();        // remove file1
  if (!Opt_params.keep_intcos)
    opt_intco_dat_remove(); // remove intco.dat
  oprintf_out("\tCleaning optimization helper files.\n");
}

// if OPT_IO_OPEN_OLD, open old file or new one
// if OPT_IO_OPEN_NEW, open new file, deleting any existing file
void opt_io_open(OPT_IO_FILE_STATUS status) {
#if defined(OPTKING_PACKAGE_PSI)
  // if file is already open, then close it
  // delete it if NEW is requested
  if (psio_open_check(PSI_OPTDATA_FILE_NUM)) {
    if (status == OPT_IO_OPEN_OLD)
      psio_close(PSI_OPTDATA_FILE_NUM, 1);
    else if (status == OPT_IO_OPEN_NEW)
      psio_close(PSI_OPTDATA_FILE_NUM, 0);
  }

  psio_open(PSI_OPTDATA_FILE_NUM, PSIO_OPEN_OLD);

#elif defined(OPTKING_PACKAGE_QCHEM)
  using opt_io::opt_data_stream;
  using namespace std;

  if ( opt_data_stream.is_open() && (status == OPT_IO_OPEN_OLD))
    return;

  if ( opt_data_stream.is_open() && (status == OPT_IO_OPEN_NEW) )
    opt_data_stream.close();

  if (status == OPT_IO_OPEN_NEW)
    remove(QCHEM_OPTDATA_FILENAME);

  if (opt_io_is_present())
    opt_data_stream.open(QCHEM_OPTDATA_FILENAME, fstream::in | fstream::out | fstream::binary);
  else
    opt_data_stream.open(QCHEM_OPTDATA_FILENAME, fstream::out | fstream::binary);

#endif
  return;
}


void opt_io_close(int keep) {

#if defined(OPTKING_PACKAGE_PSI)
  psio_close(PSI_OPTDATA_FILE_NUM, 1);
#elif defined(OPTKING_PACKAGE_QCHEM)
  opt_io::opt_data_stream.close();
#endif

  if (!keep) opt_io_remove();
  return;
}


// key    = char * ; label for entry ; not used by QChem
// buffer = char * ; stream from which to read
// size   = unsigned long int ; number of bytes to read
void opt_io_read_entry(const char *key, char *buffer, ULI size) {
#if defined(OPTKING_PACKAGE_PSI)
  psio_read_entry(PSI_OPTDATA_FILE_NUM, key, buffer, size);
#elif defined(OPTKING_PACKAGE_QCHEM)
  opt_io::opt_data_stream.read(buffer, size);
#endif
}


// key    = char * ; label for entry ; not used by QChem
// buffer = char * ; stream from which to read
// size   = unsigned long int ; number of bytes to read
void opt_io_write_entry(const char *key, char *buffer, ULI size) {
#if defined(OPTKING_PACKAGE_PSI)
  psio_write_entry(PSI_OPTDATA_FILE_NUM, key, buffer, size);
#elif defined(OPTKING_PACKAGE_QCHEM)
  opt_io::opt_data_stream.write(buffer,size);
#endif
}


} // namespace opt
