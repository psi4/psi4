/*! \file
    \ingroup INTDER
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include "displacements.h"
#include "params.h"
#define EXTERN
#include "globals.h"

#include <libchkpt/chkpt.h>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>
#include <psifiles.h>
#include <physconst.h>

using namespace psi::intder;

namespace psi { namespace intder {
extern Params gParams;
}}

void Displacements::addDisplacement(Molecule& mol)
{
  vectorMolecules.push_back(mol);
}

Molecule* Displacements::displacement(int disp)
{
  return &(vectorMolecules[disp]);
}

void Displacements::printGeometries()
{
  int index;
  
  for (index = 0; index < vectorMolecules.size(); index++)
    {
      if (index == 0)
	fprintf(outfile, "Disp #%d (Reference Geometry)\n", index);
      else
	fprintf(outfile, "Disp #%d\n", index);
      vectorMolecules[index].printGeometry();
    }
}

void Displacements::moveToCenterOfMass()
{
  int index;
  
  for (index = 0; index < vectorMolecules.size(); index++)
    vectorMolecules[index].moveToCenterOfMass();
}

void Displacements::useMasses(double *mass)
{
  int index;

  for (index = 0; index < vectorMolecules.size(); index++)
    vectorMolecules[index].useMasses(mass);
}

