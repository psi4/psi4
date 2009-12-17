/*! \file psi4.cc
\defgroup PSI4 The new PSI4 driver.
*/

//  PSI4 driver
//  Justin Turney
//  Rollin King

#include "mpi.h"
#include <getopt.h>
#include <stdio.h>
//#include <ruby.h>
#include <psiconfig.h>
#include <libciomr/libciomr.h>
#include <liboptions/liboptions.h>
#include <libparallel/parallel.h>
#include <physconst.h>

#include <molecular_system.h>

#include <psi4-def.h>
#define MAIN
#include "psi4.h"


namespace psi { 

  int psi_start(int argc, char *argv[]);
  int psi_stop(FILE* infile, FILE* outfile, char* psi_file_prefix);
  void print_version(FILE *);
  void set_memory(FILE *infile, FILE *outfile);
  int psi4_driver(Options & options, int argc, char *argv[]);

  int read_options(const std::string &name, Options & options);
  void read_atom_basis(char ** & atom_basis, int num_atoms);

  // Functions defined in ruby.c that are only needed here
  extern bool initialize_ruby();
  extern void load_input_file_into_ruby();
  extern void process_input_file();
  extern void finalize_ruby();
  extern int run_interactive_ruby();
  extern bool create_global_task();
  extern void enable_modules();

  PSIO *psio = NULL;
  std::map<std::string, PsiReturnType(*)(Options &, int argc, char *argv[])> dispatch_table;

  // These are global variable for the number of processes and
  // id for each process
  int nprocs;
  int myid;
}

// This is the ONLY main function in PSI
int main(int argc, char *argv[])
{
    using namespace psi;

    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Create the global world communicator
    Communicator::world = shared_ptr<Communicator>(new MPICommunicator(MPI_COMM_WORLD));
    
    // Needed until codes are converted to using Communicator::world
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    psi_start(argc, argv);

    if(myid==0) print_version(outfile);

    set_memory(infile, outfile);

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

   psi_stop(infile, outfile, psi_file_prefix);

   MPI_Finalize();

  // This needs to be changed a return value from the processed script
  return EXIT_SUCCESS;
}
