#include <map>
#include <vector>
#include <string>
#include <cstdio>
#include <cstdlib>

#include <liboptions/liboptions.h>
#include <libciomr/libciomr.h>

#include "memory_manager.h"

namespace psi{

using namespace std;

double bytes_to_MiB(size_t n)
{
  // 1 byte =  1/ 1048576 MiB
  return(static_cast<double>(n) / static_cast<double>(1048576));
}

MemoryManager::MemoryManager(size_t maxcor){
  CurrentAllocated    = 0;
  MaximumAllocated    = 0;
  MaximumAllowed      = maxcor;
}

MemoryManager::~MemoryManager()
{
}

void MemoryManager::RegisterMemory(void *mem, AllocationEntry& entry, size_t size)
{
  AllocationTable[mem] = entry;
  CurrentAllocated += size;
  if (CurrentAllocated > MaximumAllocated)
    MaximumAllocated = CurrentAllocated;
//  if(options_get_int("DEBUG") > 1){
//    fprintf(outfile, "\n  ==============================================================================");
//    fprintf(outfile, "\n  MemoryManager Allocated   %12ld bytes (%8.1f Mb)",size,double(size)/1048576.0);
//    fprintf(outfile, "\n  %-15s allocated   at %s:%d", entry.variableName.c_str(), entry.fileName.c_str(), entry.lineNumber);
//    fprintf(outfile, "\n  Currently used            %12ld bytes (%8.1f Mb)",CurrentAllocated,
//                 double(CurrentAllocated)/1048576.0);
//    fprintf(outfile, "\n  ==============================================================================");
//    fflush(outfile);
//  }
}

void MemoryManager::UnregisterMemory(void *mem, size_t size, const char *fileName, size_t lineNumber)
{
  CurrentAllocated -= size;
//  AllocationEntry& entry = AllocationTable[mem];
//  if(options_get_int("DEBUG") > 1){
//    fprintf(outfile, "\n  ==============================================================================");
//    fprintf(outfile, "\n  MemoryManager Deallocated %12ld bytes (%8.1f Mb)",size,double(size)/1048576.0);
//    fprintf(outfile, "\n  %-15s allocated   at %s:%d", entry.variableName.c_str(), entry.fileName.c_str(), entry.lineNumber);
//    fprintf(outfile, "\n  %-15s deallocated at %s:%d", entry.variableName.c_str(), fileName, lineNumber);
//    fprintf(outfile, "\n  Currently used            %12ld bytes (%8.1f Mb)",CurrentAllocated,
//                 double(CurrentAllocated)/1048576.0);
//    fprintf(outfile, "\n  ==============================================================================");
//    fflush(outfile);
//  }
  AllocationTable.erase(mem);
}

void MemoryManager::MemCheck(FILE *output)
{
  static bool alreadyChecked = false;

  fprintf(output, "\n\n");
  fprintf(output, "  ==============================================================================\n");
  fprintf(output, "  Memory Usage Report\n\n");
  fprintf(output, "  Maximum memory used: %8.1f Mb \n",double(MaximumAllocated)/1048576.0);
  fprintf(output, "  Number of objects still in memory: %-6lu  Current bytes used: %-14lu",(long unsigned)CurrentAllocated,(long unsigned)AllocationTable.size());

  fflush(output);
  if (AllocationTable.size() > 0) {
    if (alreadyChecked == false)
      fprintf(output, "\n\n  Attempting to free the following objects:\n");
    else
      fprintf(output, "\n\n  Unable to delete the following objects:\n");
    fflush(output);

    std::map<void*, AllocationEntry>::iterator it;

    for (it=AllocationTable.begin(); it != AllocationTable.end(); it++)
      fprintf(output, "  %15s allocated at %s:%lu\n", (*it).second.variableName.c_str(), (*it).second.fileName.c_str(), (long unsigned)(*it).second.lineNumber);
      fflush(output);
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
//      else if ((*it).second.type == "unsigned int") {
//        if ((*it).second.argumentList.size() == 1) {
//          unsigned int *m = (unsigned int*)(*it).second.variable;
//          release_one(m,__FILE__,__LINE__);
//        }
//        else if ((*it).second.argumentList.size() == 2) {
//          unsigned int **m = (unsigned int**)(*it).second.variable;
//          release_two(m,__FILE__,__LINE__);
//        }
//        else if ((*it).second.argumentList.size() == 3) {
//          unsigned int ***m = (unsigned int***)(*it).second.variable;
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
            fprintf(output, "\nRechecking memory.\n");
            MemCheck(output);
    }
  }
  fprintf(output, "\n  ==============================================================================\n");
}

} /* End Namespace */
