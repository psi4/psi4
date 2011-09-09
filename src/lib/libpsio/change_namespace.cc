/*!
 \file
 \ingroup PSIO
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#include <libparallel/parallel.h>

namespace psi {

void PSIO::change_file_namespace(unsigned int unit, const std::string & ns1, const std::string & ns2) {
    char *old_name, *new_name, *old_fullpath, *new_fullpath;
    _default_psio_lib_->get_filename(unit, &old_name, true);
    _default_psio_lib_->get_filename(unit, &new_name, true);
    //_default_psio_lib_->get_volpath(unit, 0, &path);  
    const char* path = PSIOManager::shared_object()->get_file_path(unit).c_str();

    old_fullpath = (char*) malloc( (strlen(path)+strlen(old_name)+80)*sizeof(char));
    new_fullpath = (char*) malloc( (strlen(path)+strlen(new_name)+80)*sizeof(char));
    sprintf(old_fullpath, "%s%s.%s.%u", path, old_name, ns1.c_str(), unit);
    sprintf(new_fullpath, "%s%s.%s.%u", path, new_name, ns2.c_str(), unit);

    //printf("%s\n",old_fullpath);
    //printf("%s\n",new_fullpath);

    PSIOManager::shared_object()->move_file(std::string(old_fullpath), std::string(new_fullpath)); 

    if (Communicator::world->me() == 0)
        ::rename(old_fullpath,new_fullpath);
}

}

