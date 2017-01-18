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

#include <unistd.h>
#include <cstdio>
#include <cstdlib>
#include <string>

#include "psio.hpp"
#include "psio.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/psi4-dec.h"
#include "psi4/libparallel/ParallelPrinter.h"

namespace psi{

PSIOManager::PSIOManager()
{
    pid_ = psio_getpid();

    // set the default to /tmp unless one of the
    // TMP environment variables is set
    if(std::getenv("TMPDIR"))
        set_default_path(std::getenv("TMPDIR"));
    else if(std::getenv("TEMP"))
        set_default_path(std::getenv("TEMP"));
    else if(std::getenv("TMP"))
        set_default_path(std::getenv("TMP"));
    else
        set_default_path("/tmp");
}

PSIOManager::~PSIOManager()
{
}

std::shared_ptr<PSIOManager> PSIOManager::shared_object()
{
    return _default_psio_manager_;
}

void PSIOManager::set_default_path(const std::string& path)
{
    default_path_ = path + "/";
}
void PSIOManager::set_specific_path(int fileno, const std::string& path)
{
    specific_paths_[fileno] = path + "/";
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
    if (retain) {
        specific_retains_.insert(fileno);
    }
    else {
        specific_retains_.erase(fileno);
        std::string filenum = std::to_string((long long) fileno);
        retained_files_.erase((get_file_path(fileno) + "psi." + pid_ + "."+ PSIO::get_default_namespace() + "." + filenum).c_str());
    }
    mirror_to_disk();
}

bool PSIOManager::get_specific_retention(int fileno)
{
  bool retaining = false;

  for (std::set<int>::iterator it = specific_retains_.begin(); it != specific_retains_.end(); it++) {
    if (fileno == (*it))
      retaining = true;
  }
  return retaining;
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
void PSIOManager::print(std::string out)
{
   std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
            std::shared_ptr<OutFile>(new OutFile(out)));
    printer->Printf("                    --------------------------------\n");
    printer->Printf("                    ==> Psi4 Current File Status <==\n");
    printer->Printf( "                    --------------------------------\n");
    printer->Printf( "\n");

    printer->Printf( "  Default Path: %s\n\n", default_path_.c_str());

    printer->Printf( "  Specific File Paths:\n\n");
    printer->Printf( "  %-6s %-50s\n", "FileNo", "Path");
    printer->Printf( "  ----------------------------------------------------------------------\n");
    for (std::map<int, std::string>::iterator it = specific_paths_.begin(); it != specific_paths_.end(); it++) {
        printer->Printf( "  %-6d %-50s\n", (*it).first, (*it).second.c_str());
    }
    printer->Printf( "\n");

    printer->Printf( "  Specific File Retentions:\n\n");
    printer->Printf( "  %-6s \n", "FileNo");
    printer->Printf( "  -------\n");
    for (std::set<int>::iterator it = specific_retains_.begin(); it != specific_retains_.end(); it++) {
        printer->Printf( "  %-6d\n", (*it));
    }
    printer->Printf( "\n");

    printer->Printf( "  Current File Retention Rules:\n\n");

    printer->Printf( "  %-6s \n", "Filename");
    printer->Printf( "  --------------------------------------------------\n");
    for (std::set<std::string>::iterator it = retained_files_.begin(); it != retained_files_.end(); it++) {
        printer->Printf( "  %-50s\n", (*it).c_str());
    }
    printer->Printf( "\n");

    printer->Printf( "  Current Files:\n\n");

    printer->Printf( "  %-50s%-9s%-13s\n", "Filename", "Status", "Fate");
    printer->Printf( "  ----------------------------------------------------------------------\n");
    for (std::map<std::string, bool>::iterator it = files_.begin(); it != files_.end(); it++) {
        printer->Printf( "  %-50s%-9s%-13s\n", (*it).first.c_str(), ((*it).second ? "OPEN": "CLOSED"), \
            (retained_files_.count((*it).first) == 0 ? "DEREZZ" : "SAVE"));
    }
    printer->Printf( "\n");
}
void PSIOManager::mirror_to_disk()
{


//      FILE* fh = fopen("psi.clean","w");
    FILE* fh = fopen(("psi." + pid_ + ".clean").c_str(), "w");
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

                unlink((*it).first.c_str());
        } else {
            temp[(*it).first] = (*it).second;
        }
    }
    files_.clear();
    files_ = temp;

//        unlink("psi.clean");
    unlink(("psi." + pid_ + ".clean").c_str());
}

}
