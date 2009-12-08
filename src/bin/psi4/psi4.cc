/*! \file psi4.cc
\defgroup PSI4 The new PSI4 driver.
*/

//  PSI4 driver
//  Justin Turney
//  Rollin King

//#include <mpi.h>
#include <getopt.h>
#include <stdio.h>
//#include <ruby.h>
#include <psiconfig.h>
#include <libciomr/libciomr.h>
#include <liboptions/liboptions.h>
#include <physconst.h>

#include <molecular_system.h>

#include <psi4-def.h>
#define MAIN
#include "psi4.h"

  
namespace psi { 

  void set_memory(FILE *infile, FILE *outfile);

  int psi4_driver(Options & options, int argc, char *argv[]);
  //int psi3_simulator(Options & options, int argc, char *argv[]);

  int read_options(const std::string &name, Options & options);
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
  extern void enable_modules();

    // Functions defined here
  void psi_start_and_parse_command_line(int argc, char *argv[]);
  void print_version();
  void print_usage();
  void redirect_output(const std::string& szFilename, bool append = true);
  
  PSIO *psio = NULL;
  std::map<std::string, PsiReturnType(*)(Options &, int argc, char *argv[])> dispatch_table;
}

// This is the ONLY main function in PSI
int main(int argc, char *argv[])
{
  using namespace psi;

//int err = MPI_Init(&argc, &argv);
//printf("MPI_INIT = %d", err);
  fflush(stdout);

//  int nprocs, myid;
//  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
//  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);


//  const int nproc = MPI::COMM_WORLD.Get_size();
//  const int rank = MPI::COMM_WORLD.Get_rank();


  bool run_modules = true;

  if (run_modules) {

  int num_unparsed, i;
  char *argv_unparsed[100];

  for (i=1, num_unparsed=0; i<argc; ++i)
    argv_unparsed[num_unparsed++] = argv[i];

  psi_start(&infile,&outfile,&psi_file_prefix,num_unparsed, argv_unparsed, 0);

    print_version();

    set_memory(infile, outfile);

    //ip_cwk_add(":OPT09");
    psio_init();
    psio_ipv1_config();
   
     Options options;
   
     //psi3_simulator(options, argc, argv);
    psi4_driver(options, argc, argv);
/*

       // if read from input.dat, scale.  Elsewhere always use au internally
     std::string units = options.get_str("UNITS");
     double conv_factor;
     if (units == "BOHR" || units == "AU")
       conv_factor = 1.0;
     else if (units == "ANGSTROMS" || units == "ANGSTROM")
       conv_factor = 1.0 / _bohr2angstroms;
   
     Molecular_system molecules(conv_factor);  // by default, reads geometry from input.dat
     molecules.print(); fflush(outfile);
     options.clear();
   
     try { read_atom_basis(atom_basis, molecules.get_natoms()); }
     catch (const char * s) {
       fprintf(outfile,"Unable to determine basis set:\n\t %s\n",s);
       fprintf(stderr, "Unable to determine basis set:\n\t %s\n",s);
       abort();
     }
   
     for (int i=0; i<molecules.get_num_atoms(); ++i)
       strcpy(atom_basis[i],"DZ");
     read_options("INPUT", options);
     module.set_prgid("INPUT");
     input::input(options,atom_basis,molecules);
     options.clear();

     read_options("CINTS", options);
     module.set_prgid("CINTS");
     CINTS::cints(options,argc, argv);
     module.set_prgid("CSCF");
     cscf::cscf(argc, argv);
     options.clear();
     module.set_prgid("PSICLEAN");
     psiclean::psiclean(argc, argv);
*/

   }

   psi_stop(infile, outfile, psi_file_prefix);

// MPI_Finalize();

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
    fprintf(outfile,
      "    -----------------------------------------------------------------------    \n");
    fprintf(outfile,
      "            PSI4: An Open-Source Ab Initio Electronic Structure Package \n");
    fprintf(outfile,
      "                            PSI %s Driver\n", PSI_VERSION);
    //fprintf(outfile,
    //  "                      Using Ruby %d.%d.%d interpreter.\n", RUBY_MAJOR, RUBY_MINOR, RUBY_TEENY);
    fprintf(outfile,
      "    T. D. Crawford, C. D. Sherrill, E. F. Valeev, J. T. Fermann, R. A. King,\n");
    fprintf(outfile,
      "    M. L. Leininger, S. T. Brown, C. L. Janssen, E. T. Seidl, J. P. Kenny,\n");
    fprintf(outfile,
      "    and W. D. Allen, J. Comput. Chem. 28, 1610-1616 (2007)\n");
    fprintf(outfile,
      "    -----------------------------------------------------------------------    \n\n");
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

