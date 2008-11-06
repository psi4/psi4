/*! \file psi4.cc
\defgroup PSI4 The new PSI4 driver.
*/

//
//  PSI4 driver
//  Justin Turney
//

#include <getopt.h>
#include <stdio.h>
#include <ruby.h>
#include <psiconfig.h>
#include <libciomr/libciomr.h>
#include <liboptions/liboptions.h>
#include <physconst.h>

// #include <Molecular_system.h>

#include <psi4-def.h>
#define MAIN
#include "psi4.h"

FILE *outfile = NULL;
  
namespace psi { 
  // namespace input    { int input(Options &, char **atom_basis, Molecular_system & ); }
  namespace CINTS    { int cints(Options &, int argc, char *argv[]); }
  namespace cscf     { int cscf(int argc, char *argv[]); }
  namespace psiclean { int psiclean(int argc, char *argv[]); }

  int read_options(std::string name, Options & options);
  void read_atom_basis(char ** & atom_basis, int num_atoms);
}

namespace psi {
  // Functions defined in ruby.c that are only needed here
  extern bool initialize_ruby();
  extern void load_input_file_into_ruby();
  extern void process_input_file();
  extern void finalize_ruby();
  extern int run_interactive_ruby();
  extern bool create_global_task();

    // Functions defined here
  void psi_start_and_parse_command_line(int argc, char *argv[]);
  void print_version();
  void print_usage();
  void redirect_output(const std::string& szFilename, bool append = true);
}

// This is the ONLY main function in PSI
int main(int argc, char *argv[])
{
  using namespace psi;
  // bool run_modules = true;
  // int errcod;

  // int overwrite_output = 1;
  // int num_extra_args = 0;
  // char **extra_args; 
  // extra_args = (char **) malloc(argc*sizeof(char *));
  // for (int i=1; i<argc; i++)
  //   extra_args[num_extra_args++] = argv[i];

  // No need for psi_start to be in a library, as NO module should be calling it 
  // in a single-executable. It's functionality has been moved to
  // psi_start_and_parse_command_line
  // errcod = psi_start(&infile,&outfile,&psi_file_prefix,num_extra_args, 
  //   extra_args,overwrite_output);
  
  // Parse the command-line arguments
  psi_start_and_parse_command_line(argc, argv);
  
  // Sends version information to outfile.
  print_version();

  // Initialize the interpreter
  initialize_ruby();
  // At this point the following has happened:
  //  1. Where output goes has been decided.
  //  2. Input filename has been determined.
  //  3. Ruby has been initialized

  if (create_global_task()) {
    // register_and_enable_modules();
    
    if (g_bIRB == false) { // Are running Interactive Ruby?
      // Process the input file.
      process_input_file();
    }
    else {
      run_interactive_ruby();
    }
  }
    // Test a quick call to input. This is to ensure things can link.
    // If we don't make calls to the code it isn't linked in.
  // if (run_modules) {
  // 
  //   Options options;
  // 
  //   read_options("PSI4", options);
  //     // if read from input.dat, scale.  Elsewhere always use au internally
  //   std::string units = options.get_str_option("UNITS");
  //   double conv_factor;
  //   if (units == "BOHR" || units == "AU")
  //     conv_factor = 1.0;
  //   else if (units == "ANGSTROMS" || units == "ANGSTROM")
  //     conv_factor = 1.0 / _bohr2angstroms;
  // 
  //   Molecular_system molecules(conv_factor);  // by default, reads geometry from input.dat
  //   molecules.print(); fflush(outfile);
  //   options.clear();
  // 
  //   try { read_atom_basis(atom_basis, molecules.get_num_atoms()); }
  //   catch (const char * s) {
  //     fprintf(outfile,"Unable to determine basis set:\n\t %s\n",s);
  //     fprintf(stderr, "Unable to determine basis set:\n\t %s\n",s);
  //     abort();
  //   }
  // 
  //   for (int i=0; i<molecules.get_num_atoms(); ++i)
  //     strcpy(atom_basis[i],"DZ");
  // 
  //   read_options("INPUT", options);
  //   module.set_prgid("INPUT");
  //   input::input(options,atom_basis,molecules);
  //   options.clear();
  // 
  //   read_options("CINTS", options);
  //   module.set_prgid("CINTS");
  //   CINTS::cints(options,argc, argv);
  // 
  //   module.set_prgid("CSCF");
  //   cscf::cscf(argc, argv);
  // 
  //   options.clear();
  //   module.set_prgid("PSICLEAN");
  //   psiclean::psiclean(argc, argv);
  // 
  //   for (int i=0; i<molecules.get_num_atoms(); ++i)
  //     strcpy(atom_basis[i],"DZP");
  // 
  //   read_options("INPUT", options);
  //   module.set_prgid("INPUT");
  //   input::input(options,atom_basis,molecules);
  //   options.clear();
  // 
  //   read_options("CINTS", options);
  //   module.set_prgid("CINTS");
  //   CINTS::cints(options,argc, argv);
  // 
  //   module.set_prgid("CSCF");
  //   cscf::cscf(argc, argv);
  // 
  //   options.clear();
  //   module.set_prgid("PSICLEAN");
  //   psiclean::psiclean(argc, argv);
  // }

  // psi_stop(infile, outfile, psi_file_prefix);

    // Close the interpreter
   finalize_ruby();

    // This needs to be changed a return value from the processed script
  return EXIT_SUCCESS;
}

namespace psi {
  /*! Handles the command line arguments and stores them in 
      global variables. In some cases default globals variables 
      are set. Handles the routine that psi_start would have
      done in individual modules.
        \param argc Number of command-line arguments received.
        \param argv Value of the command-line arguments.
  */
  void psi_start_and_parse_command_line(int argc, char *argv[])
  {
    int next_option;
    // Defaults:
    std::string output_filename = "output.dat";
    std::string input_filename  = "input.dat";
    std::string file_prefix     = "psi";
    bool append                 = false;
    
    // Set the default verbosity value
    g_bVerbose = false;

    // A string listing of valid short option letters
    const char* const short_options = "aIhvVo:p:i:";
    const struct option long_options[] = {
      { "append",  0, NULL, 'a' },
      { "irb",     0, NULL, 'I' },
      { "help",    0, NULL, 'h' },
      { "verbose", 0, NULL, 'v' },
      { "version", 0, NULL, 'V' },
      { "output",  1, NULL, 'o' },
      { "prefix",  1, NULL, 'p' },
      { "input",   1, NULL, 'i' },
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

        case 'I': // -I or --irb
        g_bIRB = true;
        break;

        case 'i': // -i or --input
        input_filename = optarg;
        break;
        
        case 'o': // -o or --output
        output_filename = optarg;
        break;

        case 'v': // -v or --verbose
        g_bVerbose = true;
        break;

        case 'p': // -p or --prefix
        file_prefix = optarg;
        break;
        
        case 'V': // -V or --version
        outfile = stdout;
        print_version();
        exit(EXIT_SUCCESS);
        break;

        case -1: // done with options
        break;

        default: // Something else, unexpected
        exit(EXIT_FAILURE);
        break;
      }
    } while (next_option != -1);

    // Where the output from the modules are to go
    redirect_output(output_filename, append);
    if (!output_filename.empty())
      g_szOutputFilename = output_filename;

    // Check the input file name now
    // Default to input.dat if not given
    if (optind == argc) {
      g_szInputFile = input_filename;
    }
    else {
      g_szInputFile = argv[optind];
    }
  }

  /*! Redirects the screen output from the input file and this module to the file given 
    in the argument otherwise it is sent to the screen. This doesn't quite work right 
    in practice.
    \param szFilename File to redirect to.
    \param append Append to it? default false
  */
  void redirect_output(const std::string& szFilename, bool append)
  {
    if (!szFilename.empty()) {
      // Create the file, truncate if necessary
      outfile = fopen(szFilename.c_str(), append ? "w+" : "w");

      // Error?
      if (outfile == NULL) {
        fprintf(stderr, "Unable to open: %s\n", szFilename.c_str());
        exit(EXIT_FAILURE);
      }
    }
    else {
      // Redirect to the screen
      outfile = stdout;
    }
  }

  /*! Print PSI version information that is was set in configure.ac */
  void print_version()
  {
    const char *PSI_VERSION = "0.1";
    //fprintf(Globals::g_fOutput, "PSI %s Driver\n", PSI_VERSION);
    fprintf(outfile,"\n");
    fprintf(outfile,
      "    -----------------------------------------------------------------------    \n");
    fprintf(outfile,
      "            PSI4: An Open-Source Ab Initio Electronic Structure Package \n");
    fprintf(outfile,
      "                            PSI %s Driver\n", PSI_VERSION);
    fprintf(outfile,
      "    T. D. Crawford, C. D. Sherrill, E. F. Valeev, J. T. Fermann, R. A. King,\n");
    fprintf(outfile,
      "    M. L. Leininger, S. T. Brown, C. L. Janssen, E. T. Seidl, J. P. Kenny,\n");
    fprintf(outfile,
      "    and W. D. Allen, J. Comput. Chem. 28, 1610-1616 (2007)\n");
    fprintf(outfile,
      "    -----------------------------------------------------------------------    \n");
  }

  /*! Print command-line usage information. */
  void print_usage()
  {
    printf("Usage:  psi4 [options] inputfile\n");
    printf(" -h  --help               Display this usage information.\n");
    printf(" -v  --verbose            Print a lot of information.\n");
    printf(" -V  --version            Print version information.\n");
    printf(" -o  --output filename    Redirect output elsewhere. Default: output.dat\n");
    printf(" -i  --input filename     Input file name. Default: input.dat\n");
    printf(" -p  --prefix prefix      Prefix name for psi files. Default: psi\n");
    printf(" -a  --append             Append results to output file. Default: Truncate first\n");
    printf(" -I  --irb                Run psi4 in interactive mode.\n");

    exit(EXIT_FAILURE);
  }

  void psi_abort(void)
  {
    if (outfile)
      fprintf(outfile,"\nPSI4 aborting.\n");
    else
      fprintf(stderr, "\nPSI4 aborting.\n");

    abort();
  }
} // end namespace psi::psi4

