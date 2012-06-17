/*! \file    binary_io.cc
    \ingroup optking
    \brief   functions to read and write to optking's binary file
*/

#include "io.h"

#define EXTERN
#include "globals.h"

// PSI unit number for opt_data binary file
#if defined (OPTKING_PACKAGE_PSI)
#include <libparallel/parallel.h>
#define PSI_OPTDATA_FILE_NUM 1
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#endif

#if defined (OPTKING_PACKAGE_PSI)
using namespace psi;
#endif

// Name of opt_data file for QCHEM
#define QCHEM_OPTDATA_FILENAME "opt_data.1"

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
    Communicator::world->sync();
    Communicator::world->sync();
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

  Communicator::world->sync();
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
  printf("\tCleaning optimization helper files.\n");
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

  Communicator::world->sync();
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

    Communicator::world->sync();
#if defined(OPTKING_PACKAGE_PSI)
  psio_close(PSI_OPTDATA_FILE_NUM, 1);
#elif defined(OPTKING_PACKAGE_QCHEM)
  Communicator::world->sync();
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

