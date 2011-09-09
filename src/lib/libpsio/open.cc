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

void PSIO::open(unsigned int unit, int status) {
  unsigned int i;
  char *name, *path;
  psio_ud *this_unit;
  
  //std::cout << "proc " << Communicator::world->me() << "    status = " << status << std::endl;
  /* check for too large unit */
  if (unit > PSIO_MAXUNIT)
    psio_error(unit, PSIO_ERROR_MAXUNIT);
  
  this_unit = &(psio_unit[unit]);
  
  /* Get number of volumes to stripe across */
  this_unit->numvols = get_numvols(unit);
  if (this_unit->numvols > PSIO_MAXVOL)
    psio_error(unit, PSIO_ERROR_MAXVOL);
  if (!(this_unit->numvols))
    this_unit->numvols = 1;
  
  /* Check to see if this unit is already open */
  for (i=0; i < this_unit->numvols; i++) {
    if (this_unit->vol[i].stream != -1)
      psio_error(unit, PSIO_ERROR_REOPEN);
  }
  
  /* Get the file name prefix */
  get_filename(unit, &name);
  //printf("%s\n",name);
  
  // Check if any files will have the same name
  {
    using std::string;
    typedef std::map<string,int> Names;
    Names names;
    for (i=0; i < this_unit->numvols; i++) {
      std::ostringstream oss;
      get_volpath(unit, i, &path);
      oss << path << name << "." << unit;
      const std::string fullpath = oss.str();
      typedef Names::const_iterator citer;
      citer n = names.find(fullpath);
      if (n != names.end())
        psio_error(unit, PSIO_ERROR_IDENTVOLPATH);
      names[fullpath] = 1;
      free(path);
    }
  }
  
  /* Build the name for each volume and open the file */
  for (i=0; i < this_unit->numvols; i++) {
    char* fullpath;
    get_volpath(unit, i, &path);

    #pragma warn A bit of a hack in psio open at the moment, breaks volumes and some error checking
    const char* path2 = PSIOManager::shared_object()->get_file_path(unit).c_str(); 
    
    fullpath = (char*) malloc( (strlen(path2)+strlen(name)+80)*sizeof(char));
    sprintf(fullpath, "%s%s.%u", path2, name, unit);
    this_unit->vol[i].path = strdup(fullpath);
    free(fullpath);
    
    /* Register the file */
    PSIOManager::shared_object()->open_file(std::string(this_unit->vol[i].path), unit);

    /* Now open the volume */
    if (status == PSIO_OPEN_OLD) {
      if (Communicator::world->me() == 0) {
        this_unit->vol[i].stream = ::open(this_unit->vol[i].path,O_CREAT|O_RDWR,0644);
      }
      Communicator::world->bcast(&(this_unit->vol[i].stream), 1, 0);
      //Communicator::world->raw_bcast(&(this_unit->vol[i].stream), sizeof(int), 0);
      if(this_unit->vol[i].stream == -1)
        psio_error(unit,PSIO_ERROR_OPEN);
    }
    else if(status == PSIO_OPEN_NEW) {
      if (Communicator::world->me() == 0) {
        this_unit->vol[i].stream = ::open(this_unit->vol[i].path,O_CREAT|O_RDWR|O_TRUNC,0644);
      }
      Communicator::world->bcast(&(this_unit->vol[i].stream), 1, 0);
      //Communicator::world->raw_bcast(&(this_unit->vol[i].stream), sizeof(int), 0);
      if(this_unit->vol[i].stream == -1)
        psio_error(unit,PSIO_ERROR_OPEN);
    }
    else psio_error(unit,PSIO_ERROR_OSTAT);

    free(path);
  }

  if (status == PSIO_OPEN_OLD) tocread(unit);
  else if (status == PSIO_OPEN_NEW) {
    /* Init the TOC stats and write them to disk */
    this_unit->toclen = 0;
    this_unit->toc = NULL;
    wt_toclen(unit, 0);
  }
  else psio_error(unit,PSIO_ERROR_OSTAT);

  free(name);
}

void
PSIO::rehash(unsigned int unit)
{
  if (open_check(unit)) {
    close(unit,1);
    open(unit,PSIO_OPEN_OLD);
  }
}

  int psio_open(unsigned int unit, int status) {
    _default_psio_lib_->open(unit, status);
    return 1;
  }   

}

