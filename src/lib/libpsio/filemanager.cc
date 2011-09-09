#include <boost/thread.hpp>

#include "psio.hpp"
#include "psio.h"
#include <unistd.h>
#include <cstdio>
#include "exception.h"

namespace psi{

PSIOManager::PSIOManager() : default_path_("/tmp/")
{
}
PSIOManager::~PSIOManager()
{
}
boost::shared_ptr<PSIOManager> PSIOManager::shared_object()
{
    return _default_psio_manager_;
}

void PSIOManager::set_specific_path(int fileno, const std::string& path)
{
    specific_paths_[fileno] = path;
}
std::string PSIOManager::get_file_path(int fileno)
{
    if (specific_paths_.count(fileno) != 0)
        return specific_paths_[fileno];
    else
        return default_path_;
}
void PSIOManager::set_specific_retention(int fileno, bool retain)
{
    if (retain)
        specific_retains_.insert(fileno);
    else 
        specific_retains_.erase(fileno);
    mirror_to_disk();
}

void PSIOManager::write_scratch_file(const std::string & full_path, const std::string &text)
{
    files_[full_path] = true;
    FILE* fh = fopen(full_path.c_str(),"w");
    if(!fh)
        throw PSIEXCEPTION("Unable to write to " + full_path);
    fprintf(fh, "%s",text.c_str());
    fclose(fh);
    mirror_to_disk();
}

void PSIOManager::open_file(const std::string& full_path, int fileno)
{
    files_[full_path] = true;
    if (specific_retains_.count(fileno) != 0)
        retained_files_.insert(full_path);
    mirror_to_disk();
}
void PSIOManager::close_file(const std::string& full_path, int fileno, bool keep)
{
    if (keep)
        files_[full_path] = false;
    else
        files_.erase(full_path);
    mirror_to_disk();
}
void PSIOManager::move_file(const std::string& old_full_path, const std::string& new_full_path)
{
    files_[new_full_path] = files_[old_full_path];
    files_.erase(old_full_path);
    mirror_to_disk();
}
void PSIOManager::print(FILE* out)
{
    fprintf(out, "                    --------------------------------\n");
    fprintf(out, "                    ==> PSI4 Current File Status <==\n");
    fprintf(out, "                    --------------------------------\n");
    fprintf(out, "\n");

    fprintf(out, "  Default Path: %s\n\n", default_path_.c_str());
    
    fprintf(out, "  Specific File Paths:\n\n");
    fprintf(out, "  %-6s %-50s\n", "FileNo", "Path");
    fprintf(out, "  ----------------------------------------------------------------------\n");
    for (std::map<int, std::string>::iterator it = specific_paths_.begin(); it != specific_paths_.end(); it++) {
        fprintf(out, "  %-6d %-50s\n", (*it).first, (*it).second.c_str()); 
    }
    fprintf(out, "\n");
    
    fprintf(out, "  Specific File Retentions:\n\n");
    fprintf(out, "  %-6s \n", "FileNo");
    fprintf(out, "  -------\n");
    for (std::set<int>::iterator it = specific_retains_.begin(); it != specific_retains_.end(); it++) {
        fprintf(out, "  %-6d\n", (*it)); 
    }
    fprintf(out, "\n");
    
    fprintf(out, "  Current File Retention Rules:\n\n");
    
    fprintf(out, "  %-6s \n", "Filename");
    fprintf(out, "  --------------------------------------------------\n");
    for (std::set<std::string>::iterator it = retained_files_.begin(); it != retained_files_.end(); it++) {
        fprintf(out, "  %-50s\n", (*it).c_str()); 
    }
    fprintf(out, "\n");

    fprintf(out, "  Current Files:\n\n");

    fprintf(out, "  %-50s%-9s%-13s\n", "Filename", "Status", "Fate");
    fprintf(out, "  ----------------------------------------------------------------------\n");
    for (std::map<std::string, bool>::iterator it = files_.begin(); it != files_.end(); it++) {
        fprintf(out, "  %-50s%-9s%-13s\n", (*it).first.c_str(), ((*it).second ? "OPEN": "CLOSED"), \
            (retained_files_.count((*it).first) == 0 ? "DEREZZ" : "SAVE"));
    }
    fprintf(out, "\n");
    fflush(out);
}
void PSIOManager::mirror_to_disk()
{
 
    //if (Communicator::world->me() == 0) {
      FILE* fh = fopen("psi.clean","w");
      if (fh == NULL) throw PSIEXCEPTION("PSIOManager cannot get a mirror file handle\n");
    
      for (std::map<std::string, bool>::iterator it = files_.begin(); it != files_.end(); it++) {
          if (retained_files_.count((*it).first) == 0) {
              fprintf(fh, "%s\n", (*it).first.c_str()); 
          }
      }
    
      fclose(fh);
    //}
}
void PSIOManager::build_from_disk()
{
    //if (Communicator::world->me() == 0) {

      FILE* fh = fopen("psi.clean","r");
      if (fh == NULL) throw PSIEXCEPTION("PSIOManager cannot get a mirror file handle. Is there a psi.clean file there?\n");
    
      files_.clear();
      retained_files_.clear();   

      char* in = new char[1000];
   
      while (fgets(in, 1000, fh) != NULL) {
          std::string str(in);
          str.resize(str.size()-1); // crush the newline
          files_[str] = false; 
      } 
      delete[] in;

      fclose(fh);
    //}

}
void PSIOManager::crashclean()
{
    build_from_disk();
    psiclean();
}
void PSIOManager::mark_file_for_retention(const std::string& file, bool retain)
{
    if (retain)
        retained_files_.insert(file);
    else
        retained_files_.erase(file);
    mirror_to_disk();
}
void PSIOManager::psiclean()
{
    std::map<std::string, bool> temp;
    for (std::map<std::string, bool>::iterator it = files_.begin(); it != files_.end(); it++) {
        if (retained_files_.count((*it).first) == 0) {
            //Safe to delete
            //if (Communicator::world->me() == 0)
                unlink((*it).first.c_str());
        } else {
            temp[(*it).first] = (*it).second;
        }
    }
    files_.clear();
    files_ = temp;
    //if (Communicator::world->me() == 0)
        unlink("psi.clean");
}

}
