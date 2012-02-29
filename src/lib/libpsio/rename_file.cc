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
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#include <libparallel/parallel.h>

namespace psi {

void PSIO::rename_file(unsigned int old_unit, unsigned int new_unit) {
  char*old_name,*new_name;
  /* Get the file name prefix */
  get_filename(old_unit, &old_name);
  get_filename(new_unit, &new_name);

  /* Get the path */
  const char* new_path = 
      PSIOManager::shared_object()->get_file_path(new_unit).c_str(); 
  const char* old_path = 
      PSIOManager::shared_object()->get_file_path(old_unit).c_str(); 
    
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
  system(systemcall);

  free(systemcall);
  free(old_name);
  free(new_name);
  free(old_full_path);
  free(new_full_path);
}

}
