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

/*!
 \file
 \ingroup PSIO
 */

#include <cstdio>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <cstring>
#include <cstdlib>
#include <unistd.h>
#include <string>
#include <map>
#include <sstream>
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libparallel/parallel.h"

namespace psi {

void PSIO::rename_file(unsigned int old_unit, unsigned int new_unit) {
  char*old_name,*new_name;
  /* Get the file name prefix */
  get_filename(old_unit, &old_name);
  get_filename(new_unit, &new_name);

  /* Get the path */
  std::string sold_path = PSIOManager::shared_object()->get_file_path(old_unit).c_str();
  std::string snew_path = PSIOManager::shared_object()->get_file_path(new_unit).c_str();
  const char* old_path = sold_path.c_str();
  const char* new_path = snew_path.c_str();

  /* build the full path */
  char*old_full_path =
      (char*)malloc((strlen(old_path)+strlen(old_name)+80)*sizeof(char));
  char*new_full_path =
      (char*)malloc((strlen(new_path)+strlen(new_name)+80)*sizeof(char));

  sprintf(old_full_path, "%s%s.%u", old_path, old_name, old_unit);
  sprintf(new_full_path, "%s%s.%u", new_path, new_name, new_unit);

  /* move the file.  i don't know how to do this without a system call */
  char*systemcall =
      (char*)malloc((strlen(old_full_path)+strlen(new_full_path)+100)*sizeof(char));
  sprintf(systemcall,"mv %s %s",old_full_path,new_full_path);
  int returnvalue=system(systemcall);

  free(systemcall);
  free(old_name);
  free(new_name);
  free(old_full_path);
  free(new_full_path);
}

}
