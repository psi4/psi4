/*!
** \file
** \brief Initialize input, output, file prefix, etc.
** \ingroup
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <boost/regex.hpp>
#include <getopt.h>
#include <psifiles.h>
//#include <libipv1/ip_lib.h>
#include <psiconfig.h>
#include <libparallel/parallel.h>
#include <libplugin/plugin.h>
#include "psi4.h"

using namespace std;

namespace psi {

void print_version(FILE *);
void print_usage();

/*!
** psi_start():
** This function initializes the input, output files, file prefix, etc.,
** by checking command line arguments and environmental variables. It also
** initializes the input parsing library.
**
** \param argc       = number of command-line arguments passed
** \param argv       = command-line arguments
**
** Returns: one of standard PSI error codes
** \ingroup PSI4
*/
int psi_start(int argc, char *argv[])
{
  int i, errcod;
  int next_option;
  // Defaults:
  std::string ifname;
  std::string ofname;
  std::string fprefix;
  bool append         = false;
  check_only          = false;
  clean_only          = false;
  verbose             = false;
  script              = true;

  // A string listing of valid short option letters
  const char* const short_options = "ahvVcwo:p:i:sm";
  const struct option long_options[] = {
    { "append",  0, NULL, 'a' },
    { "help",    0, NULL, 'h' },
    { "verbose", 0, NULL, 'v' },
    { "version", 0, NULL, 'V' },
    { "check",   0, NULL, 'c' },
    { "wipe",    0, NULL, 'w' },
    { "output",  1, NULL, 'o' },
    { "prefix",  1, NULL, 'p' },
    { "input",   1, NULL, 'i' },
//    { "script",  0, NULL, 's' },
    { "messy",   0, NULL, 'm' },
    { NULL,      0, NULL,  0  }
  };

  // Check the command line arguments
  do {
    next_option = getopt_long(argc, argv, short_options, long_options, NULL);

    switch (next_option) {
      case 'a': // -a or --append
      append = true;
      break;

      case 'h': // -h or --help
      print_usage();
      exit(EXIT_SUCCESS);
      break;

      case 'i': // -i or --input
      ifname = optarg;
      break;

      case 'o': // -o or --output
      ofname = optarg;
      break;

      case 'p': // -p or --prefix
      fprefix = optarg;
      break;

      case 'v': // -v or --verbose
      verbose = true;
      break;

      case 'm': // -m or --messy
      messy = true;
      break;

      case 'V': // -V or --version
      print_version(stdout);
      exit(EXIT_SUCCESS);
      break;

      case 'c': // -c or --check
      check_only = true;
      outfile = stdout;
      append = true;
      break;

//      case 's': // -s or --script
//      script = true;
//      break;

      case 'w': // -w or --wipe
      clean_only = true;
      append = true; // to avoid deletion of an existing output file
      break;

      case -1: // done with options
      break;

      default: // Something else, unexpected
      printf("Unexpected argument given: -%c\n", next_option);
      exit(EXIT_FAILURE);
      break;
    }
  } while (next_option != -1);

  char *cfname = NULL;
  char *userhome;
  FILE *psirc;
  char *arg;
  char *tmpstr1;

  if (optind <= argc-2) {
      ifname = argv[optind];
      ofname = argv[optind+1];
  }
  else if (optind <= argc-1) {
      ifname = argv[optind];
  }

  /* if some args were not specified on command-line - check the environment */
  if (ifname.empty() && Process::environment("PSI_INPUT").size())
      ifname = Process::environment("PSI_INPUT");
  if (ofname.empty() && Process::environment("PSI_OUTPUT").size())
      ofname = Process::environment("PSI_OUTPUT");
  if (fprefix.empty() && Process::environment("PSI_PREFIX").size())
      fprefix = Process::environment("PSI_PREFIX");

  /* if some arguments still not defined - assign default values */
  if (ifname.empty()) ifname = "input.dat";
  if (ofname.empty()) ofname = "output.dat";
  /* default prefix is not assigned here yet because need to check input file first */

  /* open input and output files */
  if(ifname == "stdin") infile=stdin;
  else infile = fopen(ifname.c_str(), "r");
  if(infile == NULL) {
    fprintf(stderr, "Error: could not open input file %s\n",ifname.c_str());
    return(PSI_RETURN_FAILURE);
  }

  // Check the input for "psi:(" if found assume IPV1, if not assumes Python
//  boost::regex commentLine("\\s*(%|#).*", boost::regbase::normal | boost::regbase::icase);
//  boost::regex psicommand("^\\s*psi\\s*:\\s*\\(.*", boost::regbase::normal | boost::regbase::icase);
//  char line[256];

//  script = true;
//  smatch  matches;
//  while(fgets(line, sizeof(line), infile)) {
//      string fileline = line;
//      string nocomment = boost::regex_replace(fileline, commentLine, "");
//      if (boost::regex_match(nocomment, matches, psicommand)) {
//          script = false;
//          break;
//      }
//  }

//  // Rewind the file since we've been reading from it.
//  fseek(infile, 0L, SEEK_SET);

  if(append) {
    if(ofname == "stdout") outfile=stdout;
    else outfile = fopen(ofname.c_str(), "a");
  }
  else {
    if(ofname == "stdout") outfile=stdout;
    else outfile = fopen(ofname.c_str(), "w");
  }
  if(outfile == NULL) {
    fprintf(stderr, "Error: could not open output file %s\n",ofname.c_str());
    return(PSI_RETURN_FAILURE);
  }

  // Initialize Yeti's env
  yetiEnv.init(Communicator::world->me(), ofname.c_str());
  // This seems a bit daft, but it's necessary to make sure the compiler doesn't
  // nuke it in the mistaken belief that it's unused; plugins may use it.
  yeti::Env::out0() << "";

  /* if prefix still NULL - check input file */
  if (fprefix.empty()) {
    fprefix = PSI_DEFAULT_FILE_PREFIX;
  }

  /* copy over file prefix, etc. into their appropriate variables */
  psi_file_prefix = strdup(fprefix.c_str());

  // If check_only, force output to stdout because we don't need anything more
  if(check_only) outfile = stdout;

  outfile_name = ofname;

  return(PSI_RETURN_SUCCESS);
}

/*! Print command-line usage information. */
void print_usage(void)
{
  printf("Usage:  psi4 [options] inputfile\n");
  printf(" -h  --help               Display this usage information.\n");
  printf(" -v  --verbose            Print a lot of information.\n");
  printf(" -V  --version            Print version information.\n");
  printf(" -o  --output filename    Redirect output elsewhere. Default: output.dat\n");
  printf(" -i  --input filename     Input file name. Default: input.dat\n");
  printf(" -p  --prefix prefix      Prefix name for psi files. Default: psi\n");
  printf(" -a  --append             Append results to output file. Default: Truncate first\n");
  printf(" -s  --script             Assume input file is a Python script to execute\n");

  exit(EXIT_FAILURE);
}

}

