/*! \file
    \ingroup DBOC
    \brief Enter brief description of file here 
*/

#ifndef _psi3_bin_DBOC_molecule_h_
#define _psi3_bin_DBOC_molecule_h_

namespace psi { namespace DBOC {

typedef struct {
  int natom;
  int nuniques;       // Number of unique atoms
  int nirreps;
  double **geom;
  double *zvals;
  double *masses;
  int **ict;          // Transformation table for atoms
  int *ua2a;         // Maps unique atom index -> atom index
  int *ua_degen;      // Degeneracy (length of the orbit) of unique atom
  double **cartrep;   // Group representation in cartesian basis
} Molecule_t;

}} /* namespace psi::DBOC */

#endif

