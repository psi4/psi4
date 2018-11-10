/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include <map>
#include <vector>
#include <string>
#include <cstdio>
#include <cstdlib>

#include "psi4/liboptions/liboptions.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/psi4-dec.h"
#include "memory_manager.h"
#include "psi4/libpsi4util/PsiOutStream.h"
namespace psi {

double bytes_to_MiB(size_t n) {
    // 1 byte =  1/ 1048576 MiB
    return (static_cast<double>(n) / static_cast<double>(1048576));
}

MemoryManager::MemoryManager(size_t maxcor) {
    CurrentAllocated = 0;
    MaximumAllocated = 0;
    MaximumAllowed = maxcor;
}

MemoryManager::~MemoryManager() {}

void MemoryManager::RegisterMemory(void *mem, AllocationEntry &entry, size_t size) {
    AllocationTable[mem] = entry;
    CurrentAllocated += size;
    if (CurrentAllocated > MaximumAllocated) MaximumAllocated = CurrentAllocated;
    //  if(options_get_int("DEBUG") > 1){
    //    outfile->Printf( "\n  ==============================================================================");
    //    outfile->Printf( "\n  MemoryManager Allocated   %12ld bytes (%8.1f Mb)",size,double(size)/1048576.0);
    //    outfile->Printf( "\n  %-15s allocated   at %s:%d", entry.variableName.c_str(), entry.fileName.c_str(),
    //    entry.lineNumber);
    //    outfile->Printf( "\n  Currently used            %12ld bytes (%8.1f Mb)",CurrentAllocated,
    //                 double(CurrentAllocated)/1048576.0);
    //    outfile->Printf( "\n  ==============================================================================");
    //
    //  }
}

void MemoryManager::UnregisterMemory(void *mem, size_t size, const char *fileName, size_t lineNumber) {
    CurrentAllocated -= size;
    //  AllocationEntry& entry = AllocationTable[mem];
    //  if(options_get_int("DEBUG") > 1){
    //    outfile->Printf( "\n  ==============================================================================");
    //    outfile->Printf( "\n  MemoryManager Deallocated %12ld bytes (%8.1f Mb)",size,double(size)/1048576.0);
    //    outfile->Printf( "\n  %-15s allocated   at %s:%d", entry.variableName.c_str(), entry.fileName.c_str(),
    //    entry.lineNumber);
    //    outfile->Printf( "\n  %-15s deallocated at %s:%d", entry.variableName.c_str(), fileName, lineNumber);
    //    outfile->Printf( "\n  Currently used            %12ld bytes (%8.1f Mb)",CurrentAllocated,
    //                 double(CurrentAllocated)/1048576.0);
    //    outfile->Printf( "\n  ==============================================================================");
    //
    //  }
    AllocationTable.erase(mem);
}

void MemoryManager::MemCheck(std::string out) {
    std::shared_ptr<psi::PsiOutStream> printer = (out == "outfile" ? outfile : std::make_shared<PsiOutStream>(out));
    static bool alreadyChecked = false;

    printer->Printf("\n\n");
    printer->Printf("  ==============================================================================\n");
    printer->Printf("  Memory Usage Report\n\n");
    printer->Printf("  Maximum memory used: %8.1f Mb \n", double(MaximumAllocated) / 1048576.0);
    printer->Printf("  Number of objects still in memory: %-6lu  Current bytes used: %-14lu",
                    (long unsigned)CurrentAllocated, (long unsigned)AllocationTable.size());

    if (AllocationTable.size() > 0) {
        if (alreadyChecked == false)
            printer->Printf("\n\n  Attempting to free the following objects:\n");
        else
            printer->Printf("\n\n  Unable to delete the following objects:\n");

        std::map<void *, AllocationEntry>::iterator it;

        for (it = AllocationTable.begin(); it != AllocationTable.end(); it++)
            printer->Printf("  %15s allocated at %s:%lu\n", (*it).second.variableName.c_str(),
                            (*it).second.fileName.c_str(), (long unsigned)(*it).second.lineNumber);
        //
        //    it = AllocationTable.begin();
        //    while (it != AllocationTable.end()) {
        //      if ((*it).second.type == "double") {
        //        if ((*it).second.argumentList.size() == 1) {
        //          double *m = (double*)(*it).second.variable;
        //          release_one(m,__FILE__,__LINE__);
        //        }
        //        else if ((*it).second.argumentList.size() == 2) {
        //          double **m = (double**)(*it).second.variable;
        //          release_two(m,__FILE__,__LINE__);
        //        }
        //        else if ((*it).second.argumentList.size() == 3) {
        //          double ***m = (double***)(*it).second.variable;
        //          release_three(m,__FILE__,__LINE__);
        //        }
        //      }
        //      else if ((*it).second.type == "int") {
        //        if ((*it).second.argumentList.size() == 1) {
        //          int *m = (int*)(*it).second.variable;
        //          release_one(m,__FILE__,__LINE__);
        //        }
        //        else if ((*it).second.argumentList.size() == 2) {
        //          int **m = (int**)(*it).second.variable;
        //          release_two(m,__FILE__,__LINE__);
        //        }
        //        else if ((*it).second.argumentList.size() == 3) {
        //          int ***m = (int***)(*it).second.variable;
        //          release_three(m,__FILE__,__LINE__);
        //        }
        //      }
        //      else if ((*it).second.type == "char") {
        //        if ((*it).second.argumentList.size() == 1) {
        //          char *m = (char*)(*it).second.variable;
        //          release_one(m,__FILE__,__LINE__);
        //        }
        //        else if ((*it).second.argumentList.size() == 2) {
        //          char **m = (char**)(*it).second.variable;
        //          release_two(m,__FILE__,__LINE__);
        //        }
        //        else if ((*it).second.argumentList.size() == 3) {
        //          char ***m = (char***)(*it).second.variable;
        //          release_three(m,__FILE__,__LINE__);
        //        }
        //      }
        //      else if ((*it).second.type == "float") {
        //        if ((*it).second.argumentList.size() == 1) {
        //          float *m = (float*)(*it).second.variable;
        //          release_one(m,__FILE__,__LINE__);
        //        }
        //        else if ((*it).second.argumentList.size() == 2) {
        //          float **m = (float**)(*it).second.variable;
        //          release_two(m,__FILE__,__LINE__);
        //        }
        //        else if ((*it).second.argumentList.size() == 3) {
        //          float ***m = (float***)(*it).second.variable;
        //          release_three(m,__FILE__,__LINE__);
        //        }
        //      }
        //      else if ((*it).second.type == "size_t") {
        //        if ((*it).second.argumentList.size() == 1) {
        //          size_t *m = (size_t*)(*it).second.variable;
        //          release_one(m,__FILE__,__LINE__);
        //        }
        //        else if ((*it).second.argumentList.size() == 2) {
        //          size_t **m = (size_t**)(*it).second.variable;
        //          release_two(m,__FILE__,__LINE__);
        //        }
        //        else if ((*it).second.argumentList.size() == 3) {
        //          size_t ***m = (size_t***)(*it).second.variable;
        //          release_three(m,__FILE__,__LINE__);
        //        }
        //      }
        //      else if ((*it).second.type == "unsigned char") {
        //        if ((*it).second.argumentList.size() == 1) {
        //          unsigned char *m = (unsigned char*)(*it).second.variable;
        //          release_one(m,__FILE__,__LINE__);
        //        }
        //        else if ((*it).second.argumentList.size() == 2) {
        //          unsigned char **m = (unsigned char**)(*it).second.variable;
        //          release_two(m,__FILE__,__LINE__);
        //        }
        //        else if ((*it).second.argumentList.size() == 3) {
        //          unsigned char ***m = (unsigned char***)(*it).second.variable;
        //          release_three(m,__FILE__,__LINE__);
        //        }
        //      }
        //      it = AllocationTable.begin();
        //    }

        if (alreadyChecked == false && AllocationTable.size() > 0) {
            alreadyChecked = true;
            printer->Printf("\nRechecking memory.\n");
            MemCheck("output");
        }
    }
    printer->Printf("\n  ==============================================================================\n");
}

} /* End Namespace */
