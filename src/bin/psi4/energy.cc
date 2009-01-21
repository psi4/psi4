/*! \file energy_methods.cc
    \defgroup defines methods for PSI energies
*/

#include <psi4-dec.h>
#include <Molecular_system.h>
#include "psi4.h"

namespace psi {

double energy(energy_type etype)
{
  switch(etype) 
  {
    case SCF : 
      read_options("INPUT", options);
      module.set_prgid("INPUT");
      input::input(options,atom_basis,molecules);
      options.clear();

      read_options("CINTS", options);
      module.set_prgid("CINTS");
      CINTS::cints(options,argc, argv);
      options.clear();

      read_options("CSCF", options);
      module.set_prgid("CSCF");
      cscf::cscf(argc, argv);
      options.clear();

      module.set_prgid("PSICLEAN");
      psiclean::psiclean(argc, argv);
      return 1;
   case CISD :
      fprintf(outfile,"Do CISD energy\n");
      return 1;
   case CISD :
      fprintf(outfile,"Do CISD energy\n");
      return 1;
   case G3 :
      energy(SCF);
      energy(CISD);
      energy(CCSD);
      return 1;
    default :
      throw("invalid energy type");
  }
}

//double energy_cbs(energy_type etype)
} // psi

