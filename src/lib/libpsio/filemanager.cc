#include "psio.hpp"
#include "psio.h"
#include <unistd.h>

namespace psi{

PSIOManager::PSIOManager()
{
}
PSIOManager::~PSIOManager()
{
}
void PSIOManager::open_file(const std::string& full_path)
{
    files_[full_path] = true;
}
void PSIOManager::close_file(const std::string& full_path, bool keep)
{
    if (keep)
        files_[full_path] = false;
    else
        files_.erase(full_path);
}
void PSIOManager::move_file(const std::string& old_full_path, const std::string& new_full_path)
{
    files_[new_full_path] = files_[old_full_path];
    files_.erase(old_full_path);
}
void PSIOManager::print(FILE* out)
{
    fprintf(out, "                    --------------------------------\n");
    fprintf(out, "                    ==> PSI4 Current File Status <==\n");
    fprintf(out, "                    --------------------------------\n");
    fprintf(out, "\n");

    fprintf(out, "%-50s%-9s%-13s\n", "Filename", "Status", "Fate");
    fprintf(out, "----------------------------------------------------------------------\n");
    for (std::map<std::string, bool>::iterator it = files_.begin(); it != files_.end(); it++) {
        fprintf(out, "%-50s%-9s%-13s\n", (*it).first.c_str(), ((*it).second ? "OPEN": "CLOSED"), \
            (retained_files_.count((*it).first) == 0 ? "DELETE" : "SAVE"));
    }
    fflush(outfile);
}
void PSIOManager::mark_file_for_retention(const std::string& file)
{
    retained_files_.insert(file);
}
void PSIOManager::mark_file_for_deletion(const std::string& file)
{
    retained_files_.erase(file);
}
void PSIOManager::psiclean()
{
    for (std::map<std::string, bool>::iterator it = files_.begin(); it != files_.end(); it++) {
        if (retained_files_.count((*it).first) == 0) {
            //Safe to delete
            unlink((*it).first.c_str());
        }
    }
    files_.clear();
}

}
