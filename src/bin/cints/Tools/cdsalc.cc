/*! \file cdsalc.cc
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <libint/libint.h>
#include "defines.h"
#define  EXTERN
#include "global.h"

#include "cdsalc.h"

namespace psi { namespace CINTS {

void init_cdsalc()
{
  CDSALCs.nsalcs = Molecule.num_atoms*3;
  const int num_cd = CDSALCs.nsalcs;

  CDSALCs.salc2irrep = new int[CDSALCs.nsalcs];
  CDSALC_t::cd2salc_map_t* cd2salc_map = new CDSALC_t::cd2salc_map_t[num_cd];
  for(int cd=0; cd<num_cd; cd++)
    cd2salc_map[cd].nsalcs = 0;

  CDSALCs.atom_irreps = new char[Molecule.num_atoms];
  memset(CDSALCs.atom_irreps,0,sizeof(char)*Molecule.num_atoms);

  for(int irrep=0; irrep<Symmetry.nirreps; irrep++) {
    const int nsalcs_in_irrep = Symmetry.cdsalcpi[irrep];
    const int salc_offset = Symmetry.cdsalc_ioffset[irrep];

    for(int salc=0; salc<nsalcs_in_irrep; salc++) {
      CDSALCs.salc2irrep[salc+salc_offset] = irrep;

      for(int cd=0; cd<CDSALCs.nsalcs; cd++) {
        if (Symmetry.cdsalc2cd[cd][salc+salc_offset] != 0.0) {
          CDSALCs.atom_irreps[cd/3] |= (1 << irrep);
          ++cd2salc_map[cd].nsalcs;
        }
      }
    }
  }
  
  for(int cd=0; cd<num_cd; cd++) {
    cd2salc_map[cd].salcs = new int[cd2salc_map[cd].nsalcs];
    int salc_count = 0;
    for(int salc=0; salc<CDSALCs.nsalcs; salc++) {
      if (Symmetry.cdsalc2cd[cd][salc] != 0.0) {
        cd2salc_map[cd].salcs[salc_count] = salc;
        ++salc_count;
      }
    }
  }
  CDSALCs.cd2salc_map = cd2salc_map;
}


void cleanup_cdsalc()
{
  delete[] CDSALCs.atom_irreps;
  delete[] CDSALCs.salc2irrep;
  for(int cd=0; cd<CDSALCs.nsalcs; cd++) {
    delete[] CDSALCs.cd2salc_map[cd].salcs;
  }
  delete[] CDSALCs.cd2salc_map;
}
}}
