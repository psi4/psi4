/*! \file
    \ingroup OEPROP
    \brief Enter brief description of file here 
*/
#define EXTERN
#include <stdexcept>
#include "includes.h"
#include "prototypes.h"
#include "globals.h"

namespace psi { namespace oeprop {

void get_nmo()
{
  // throw std::runtime_error("Computing natural orbitals has not been tested and should not be used."); turn this back off for testing, CDS 4/08
  
  int i,k,l,count,dim_i;
  double *evals, *evals_symbl;
  double **nmo_mo, **nmo_symbl;
  double **b, **binv;		/* (C+)*C */
  double **cinv;                /* C_{AO}^-1 = C_{AO}^T S_{AO,AO} */
  double **tmat, **Smat;
  double **psq_ao, **p_symbl;
  
  nmo_ao = init_matrix(nbfao,nmo);
  nmo_so = block_matrix(nbfso,nmo); /* need block matrix for chkpt_wt_scf() */
  evals = init_array(nmo);
  nmo_mo = init_matrix(nmo,nmo);
  psq_ao = init_matrix(nbfao,nbfao);

  tri_to_sq(Ptot,psq_ao,nbfao);

  b = init_matrix(nmo,nmo);
  cinv = init_matrix(nmo, nbfao);
  tmat = init_matrix(nmo,nbfao);
  Smat = init_matrix(nbfao,nbfao);
  tri_to_sq(S,Smat,nbfao);

  mmult(scf_evec_ao,1,Smat,0,cinv,0,nmo,nbfao,nbfao,0);

  mmult(cinv,0,psq_ao,0,tmat,0,nmo,nbfao,nbfao,0);
  mmult(tmat,0,cinv,1,b,0,nmo,nbfao,nmo,0);
  free_matrix(cinv,nmo);
  free_matrix(tmat,nmo);
  free_matrix(Smat,nbfao);

  count = 0;
  for(i=0;i<nirreps;i++) {
    dim_i = orbspi[i];
    if (dim_i != 0) {
      p_symbl = init_matrix(dim_i,dim_i);
      evals_symbl = init_array(dim_i);
      nmo_symbl = init_matrix(dim_i,dim_i);
      for(k=0;k<dim_i;k++) {
        for(l=0;l<=k;l++) {
          p_symbl[k][l] = p_symbl[l][k] = b[count+k][count+l];
        }
      }

      if (print_lvl >= PRINTNMOLEVEL) {
        fprintf(outfile, "  Density Matrix for Symmetry Block %s\n",
	        irr_labs[i]);
        print_mat(p_symbl,dim_i,dim_i,outfile);
        fprintf(outfile,"\n");
      }

      sq_rsp(dim_i,dim_i,p_symbl,evals_symbl,3,nmo_symbl,1.0E-14);

      for(k=0;k<dim_i;k++) {
        evals[count+k] = evals_symbl[k];
        for(l=0;l<dim_i;l++)
          nmo_mo[count+k][count+l] = nmo_symbl[k][l];
      }
      
      if (print_nos) {
        fprintf(outfile, 
          "  Natural Orbital Occupation Numbers for Symmetry Block %s\n\n",
          irr_labs[i]);
        for (k=0; k<dim_i; k++) {
          fprintf(outfile, "%5d  %f\n",k+1,evals_symbl[k]);
        }
        fprintf(outfile, "\n");
      }

      if (print_lvl >= PRINTNMOLEVEL) {
        fprintf(outfile, "  Natural Orbitals for Symmetry Block %s\n",
	        irr_labs[i]);
        eivout(nmo_symbl,evals_symbl,dim_i,dim_i,outfile);
	fprintf(outfile, "\n");
      }
      free_matrix(p_symbl,dim_i);
      free_matrix(nmo_symbl,dim_i);
      free(evals_symbl);
      count += dim_i;
    }
  }
                            
  mmult(scf_evec_ao,0,nmo_mo,0,nmo_ao,0,nbfao,nmo,nmo,0);
  mmult(scf_evec_so,0,nmo_mo,0,nmo_so,0,nbfso,nmo,nmo,0);


  if (print_lvl >= PRINTNMOLEVEL) {
    fprintf(outfile,
            "  Natural orbitals in SO basis computed from density in file%d:\n",
            opdm_file);
    eivout(nmo_so,evals,nbfso,nmo,outfile);
    fprintf(outfile,"\n\n");
  }

  if (wrtnos) {
    chkpt_wt_scf(nmo_so);
    if (print_lvl >= 1)
      fprintf(outfile,
        "  Natural Orbitals have just been written to the checkpoint file\n\n");
  }

  free(evals);
  free_matrix(b,nmo);
  free_matrix(nmo_mo,nmo);
  free_block(nmo_so);
}

}} // namespace psi::oeprop
