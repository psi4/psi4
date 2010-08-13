/*! \file execute_bp.cc
\defgroup PSI4 Function to execute basic procedures in psi.dat
*/
#include <iostream>
#include <fstream>              // file I/O support
#include <libipv1/ip_lib.h>
#include <libparallel/parallel.h>
#include <psi4-dec.h>
#include <string.h>
#include <algorithm>
#include <ctype.h>
#include "psi4.h"

#define MAX_ARGS (20)

using namespace std;

namespace psi {

extern void psiclean(void);
extern int read_options(const std::string &jobName, Options &options,
  bool call_ipv1 = true, bool suppress_printing = false);

PsiReturnType execute_bp(std::string & bp, Options & options) {
    std::string calcType = bp;

    char *jobList = const_cast<char*>(calcType.c_str());

    // This version assumes that the array contains only module names, not
    // macros for other job types, like $SCF
    int numTasks = 0;
    int errcod;

    if(ip_exist(const_cast<char*>("EXEC"), 0)) {
      // Override psi.dat with the commands in the exec array in input
      errcod = ip_count(const_cast<char*>("EXEC"), &numTasks, 0);
      jobList = const_cast<char*>("EXEC");
    }
    else {
      errcod = ip_count(jobList, &numTasks, 0);

      if (!ip_exist(jobList, 0)){
	std::string err("Error: jobtype ");
	err += jobList;
	err += " is not defined in psi.dat";
	throw PsiException(err, __FILE__, __LINE__);
      }

      if (errcod != IPE_OK){
	std::string err("Error: trouble reading ");
	err += jobList;
	err += " array from psi.dat";
	throw PsiException(err, __FILE__, __LINE__);
      }
    }
/*    else {
      errcod = ip_count(jobList, &numTasks, 0);
      if (!ip_exist(jobList, 0)){
	std::string err("Error: jobtype ");
	err += jobList;
	err += " is not defined in psi.dat";
	throw PsiException(err, __FILE__, __LINE__);
      }
      if (errcod != IPE_OK){
	std::string err("Error: trouble reading ");
	err += jobList;
	err += " array from psi.dat";
	throw PsiException(err, __FILE__, __LINE__);
      }
    } */

    fprintf(outfile, "    List of tasks to execute:\n");
    for(int n = 0; n < numTasks; ++n) {
      char *thisJob;
      ip_string(jobList, &thisJob, 1, n);
      fprintf(outfile, "    %s\n", thisJob);
      free(thisJob);
    }

    if(check_only)
      return Success;

    // variables to parse string in psi.dat
    string thisJobWithArguments;
    stringstream ss;
    vector<string> tokens;
    string buf;

    for(int n = 0; n < numTasks; ++n){
      char *thisJob;
      errcod = ip_string(jobList, &thisJob, 1, n);

      // tokenize string in psi.dat
      thisJobWithArguments.assign(thisJob);
      ss.clear();
      ss << thisJobWithArguments;
      while (ss >> buf)
	tokens.push_back(buf);

      free(thisJob);
      thisJob = const_cast<char *>(tokens[0].c_str()); // module name for dispatch table

      // Make sure the job name is all upper case
      int length = strlen(thisJob);
      std::transform(thisJob, thisJob + length, thisJob, ::toupper);
      if(Communicator::world->me() == 0) {
	fprintf(outfile, "\n  Job %d is %s\n", n, thisJob); fflush(outfile);
      }

      // Read the options for thisJob.
      read_options(thisJob, options);

      // Handle MODE (command line argument passed in task list
      if (tokens.size() > 1) {
	// Convert token to upper case
	std::transform(tokens[1].begin(), tokens[1].end(), tokens[1].begin(), ::toupper);
	// Set the option overriding anything the user said.
	options.set_str("MODE", tokens[1]);
      }

      // If the function call is LMP2, run in parallel
      if(strcmp(thisJob, "LMP2") == 0 || strcmp(thisJob, "DFMP2") == 0) {
	// Needed a barrier before the functions are called
	Communicator::world->sync();

	if (dispatch_table[thisJob](options) != Success) {
	  // Good chance at this time that an error occurred.
	  // Report it to the user.
	  fprintf(stderr, "%s did not return a Success code.\n", thisJob);
	  throw PsiException("Module failed.", __FILE__, __LINE__);
	}
      }
      // If any other functions are called only process 0 runs the function
      else {
	if (dispatch_table.find(thisJob) != dispatch_table.end()) {
	  // Needed a barrier before the functions are called
	  Communicator::world->sync();
	  if(Communicator::world->me() == 0) {
	    if (dispatch_table[thisJob](options) != Success) {
	      // Good chance at this time that an error occurred.
	      // Report it to the user.
	      fprintf(stderr, "%s did not return a Success code.\n", thisJob);
	      throw PsiException("Module failed.", __FILE__, __LINE__);
	    }
	  }
	}
	else {
	  std::transform(thisJob, thisJob + length, thisJob, ::tolower);

	  // Close the output file, allowing the external program to write to it.
	  if (!outfile_name.empty())
	    fclose(outfile);

	  // Attempt to run the external program
	  int ret = ::system(thisJob);

	  // Re-open the output file, allowing psi4 to take output control again.
	  if (!outfile_name.empty())
	    fopen(outfile_name.c_str(), "a");

	  if (ret == -1 || ret == 127) {
	    std::string err = "Module ";
	    err += thisJob;
	    err += " is not known to PSI4.  Please update the driver\n";
	    throw PsiException(err, __FILE__, __LINE__);
	  }
	}
      }

      tokens.clear();

      fflush(outfile);
    }

    return Success;
  }

} // Namespace

