/*! \file
    \ingroup INPUT
    \brief Enter brief description of file here 
*/
#define EXTERN
#include <cstdio>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <psi4-dec.h>

namespace psi {

  extern int num_atoms;

void punt(const char *mess);

// establish basis label for each atom
void read_atom_basis(char ** & atom_basis) {

  int i,j,errcod, depth = 0;
  char *basis_type, *buf;

  atom_basis = (char **) malloc(num_atoms*sizeof(char *));

  errcod = ip_count("BASIS",&depth,0);
  if (depth == 0) { // the same basis for all atoms, basis = dzp
   errcod = ip_string("BASIS",&basis_type,0);
   if (errcod != IPE_OK)
     throw("There is a problem with the BASIS keyword!");
   atom_basis[0] = basis_type;
   for(i=0;i<num_atoms;i++) {
     atom_basis[i] = (char *) malloc((1+strlen(basis_type))*sizeof(char));
     strcpy(atom_basis[i],basis_type);
   }
  }
  else {
    errcod = ip_count("BASIS",&j,1,0);
    if (errcod == IPE_NOT_AN_ARRAY) { // basis set for each atom, basis = (dzp dz ...)
      if (depth == num_atoms) {
        for(i=0;i<num_atoms;i++) {
          errcod = ip_string("BASIS",&basis_type,1,i);
          if (errcod != IPE_OK)
            throw("There is a problem with the BASIS array!");
          atom_basis[i] = basis_type;
        }
      }
      else
        throw("Number of entries in the BASIS array not the same as num_atoms!");
    }
    /* add this functionality later - it requires atomic numbers
    else {
       // Basis sets for each element type is specified, e.g.
       // basis = ( (h dz) (c dzp) )

     for(i=0;i<depth;i++) {
       errcod = ip_count("BASIS",&j,1,i);
       if (errcod != IPE_OK || j != 2) {
         fprintf(outfile,"  There is a problem with line %d of the BASIS array!\n",i+1);
         throw("Invalid basis set file");
       }
       errcod = ip_string("BASIS",&elem_label,2,i,0);
       atom_num(elem_label,&Z);
       free(elem_label);
       errcod = ip_string("BASIS",&basis_type,2,i,1);
       for(k=0;k<num_atoms;k++)
         if (!strcmp(elem_name[(int)Z],element[k]))
            atom_basis[k] = basis_type;
     }
     for(i=0;i<num_atoms;i++)
       if (atom_basis[i] == NULL) {
         fprintf(outfile,"  Missing basis set for %s\n",element[i]);
         throw("Missing basis set");
       }
    }
  */
  }
  return;
}

}
