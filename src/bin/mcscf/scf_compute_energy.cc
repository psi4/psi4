#include <iostream>
#include <cstdio>

#include <libchkpt/chkpt.hpp>
#include "scf.h"

extern FILE* outfile;

namespace psi{ namespace mcscf{

double SCF::compute_energy()
{
  fprintf(outfile,"\n\n  Running an SCF calculation");

  startup();

  // Read the one electron integrals
  read_so_oei();

  // Read the two electron integrals
  // and construct the PK and K matrices
  read_so_tei();

  // Construct the S^-1/2 Matrix
  construct_S_inverse_sqrt();

  // Guess C
  initial_guess();

  // Iterate the SCF equations
  iterate_scf_equations();

  // Check the orthonormality of the MOs
  check_orthonormality();

  // Canonicalize MOs
  canonicalize_MO();

  // Print eigenvectors and MOs
  print_eigenvectors_and_MO();

  // Canonicalize MOs
  save_info();

  return(0.0);
}

}} /* End Namespaces */
