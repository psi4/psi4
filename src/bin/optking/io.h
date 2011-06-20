/*! \file    binary_io.h
    \ingroup optking
    \brief   header for binary reading and writing functions
*/

#ifndef _opt_io_h_
#define _opt_io_h_

#include "package.h"

namespace opt {

typedef unsigned long int ULI;

enum OPT_IO_FILE_STATUS {OPT_IO_OPEN_NEW, OPT_IO_OPEN_OLD} ;

bool opt_io_is_present(void);
void opt_io_remove(void);
void opt_io_open(OPT_IO_FILE_STATUS status);
void opt_io_close(int keep);
void opt_io_read_entry(const char *key, char *buffer, ULI size);
void opt_io_write_entry(const char *key, char *buffer, ULI size);
void opt_intco_dat_remove(void);
void opt_clean(void);

}

#endif

