/*!
  \file
  \brief Filter out unneeded frozen core/virt integrals
  \ingroup QT
*/

namespace psi {
	
/*!
** filter(): Filter out undesired (frozen core/virt) integrals
**
** Given a lower-triangle array of integrals in the full
** space of orbitals as well as numbers of frozen core and virtual
** orbitals, this function returns a list of integrals involving only
** active orbitals.
**
** TDC, June 2001
**
** Note: Based on the code written by CDS in the original
** iwl_rd_one_all_act() function in LIBIWL.
**
** \ingroup QT
*/

void filter(double *input, double *output, int *ioff, int norbs, 
            int nfzc, int nfzv)
{
  int i, j, ij, IJ;
  int nact;

  nact = norbs - nfzc - nfzv;

  for(i=0,ij=0; i < nact; i++) {
    for(j=0; j <= i; j++,ij++) {
      IJ = ioff[i+nfzc] + (j + nfzc);
      output[ij] = input[IJ];
    }
  }
}

}

