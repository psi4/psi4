/*! \file
    \ingroup OEPROP
    \brief Enter brief description of file here 
*/
#include "includes.h"
#include "prototypes.h"
#include "globals.h"

int main(int argc, char* argv[]) {
 using namespace psi::oeprop;
 
 int i,j,k,l,count;
 char buffer[80];	/* buffer string */

	/* Setting defaults */
 read_opdm = 0;
 opdm_file = 76;
 asymm_opdm = 0;
 wrtnos = 0;
 spin_prop = 0;
 print_lvl = 1;
 wrt_dipmom = 1;
 corr = 0;
 mpmax = 1;
 mp_ref = 0;
 nuc_esp = 1;
 grid = 0;
 grid3d = 0;
 nix = 10;
 niy = 10;
 niz = 10;
 grid_zmin = 0.0;
 grid_zmax = 3.0;
 edgrad_logscale = 5;
 zvec_file = 86;
 delete_zvec = 1;
 

	/* Initialization and printing intro */
 start_io(argc,argv);
 print_intro();

/*************************** Main Code *******************************/

	/* Reading in basic information from checkpoint file */

 title = chkpt_rd_label();
 natom = chkpt_rd_natom();
 natom3 = natom * 3;
 nmo = chkpt_rd_nmo();
 nbfso = chkpt_rd_nso();
 nbfao = chkpt_rd_nao();
 natri = nbfao * (nbfao+1)/2;
 nstri = nbfso * (nbfso+1)/2;
 nshell = chkpt_rd_nshell();
 nprim = chkpt_rd_nprim();
 iopen = chkpt_rd_iopen();
 nirreps = chkpt_rd_nirreps(); 
 nsym = chkpt_rd_nsymhf();
 orbspi = chkpt_rd_orbspi();
 sopi = chkpt_rd_sopi();
 clsdpi = chkpt_rd_clsdpi();    
 openpi = chkpt_rd_openpi();
 irr_labs = chkpt_rd_irr_labs();
 geom = chkpt_rd_geom();
 zvals = chkpt_rd_zvals();
 usotao = chkpt_rd_usotao();
    
 /* Parsing */
 parsing();

 if (strcmp(ref,"UHF") != 0) {
   scf_evec_so = chkpt_rd_scf();
   scf_evals = chkpt_rd_evals();
   scf_evec_ao = block_matrix(nbfao,nmo);
   mmult(usotao,1,scf_evec_so,0,scf_evec_ao,0,nbfao,nbfso,nmo,0);
 }
 
 /* parse the grid information */
 if (ip_exist("GRID",0))
   grid_parse();
 
 /* Computing total charge of the system */

 double charge = 0;
 for(i=0;i<nirreps;i++)
   charge -= 2*clsdpi[i] + openpi[i];
 for(i=0;i<natom;i++)
   charge += zvals[i];

 /* Setting up an offset array */

 ioff = init_int_array(nbfao);
 ioff[0] = 0;
 for(i=1;i<nbfao;i++) {
   ioff[i] = ioff[i-1] + i;
 }

	/* Computing double factorials df[i] = (i-1)!! */
 df[0] = 1.0;
 df[1] = 1.0;
 df[2] = 1.0;
 for(i=3; i<MAXFACT*2; i++) {
   df[i] = (i-1)*df[i-2];
 }
                 
 /* Printing tasks and parameters */

 if (print_lvl >= PRINTTASKPARAMLEVEL) {
   print_tasks();
   print_params();
 }

 /* Reading in basis set inforamtion */

 read_basset_info();
 init_xyz(); 

 chkpt_close();

 /* RAK and CDS: Determine number of densities to analyze */
 get_opdm_lbl();

for (irho=0;irho<nrho;irho++) {

 if (transdens && irho == 0) continue;

 fprintf(outfile,"\t** Analyzing ");
 if (transdens)
   fprintf(outfile, "transition ");
 fprintf(outfile, "density number %d **\n", irho+1);

 chkpt_init(PSIO_OPEN_OLD);

 /* Obtain a density matrix */
 if (read_opdm)
   read_density();
 else
   compute_density();

  compute_overlap();
   
  /* Obtain natural orbitals */
  if (read_opdm && wrtnos)
    get_nmo(); 
   
  chkpt_close();
   
  /* Mulliken population analysis */
  print_pop_header();
  populate();
 
  /* Compute coordinates of the MP reference point if needed */
  if (mp_ref != -1)
    compute_mp_ref_xyz();

  /* Moving molecule to the reference center!
   Attention, ever since coordinates of atoms and of the grid box 
   are stored in this new coordinate system. All coordinates are 
   transformed back at the moment of printing out. */

  for(i=0;i<natom;i++) {
    geom[i][0] -= mp_ref_xyz[0];
    geom[i][1] -= mp_ref_xyz[1];
    geom[i][2] -= mp_ref_xyz[2];
    Lm_ref_xyz[0] -= mp_ref_xyz[0];
    Lm_ref_xyz[1] -= mp_ref_xyz[1];
    Lm_ref_xyz[2] -= mp_ref_xyz[2];
  }
  if (grid) {
    grid_origin[0] -= mp_ref_xyz[0];
    grid_origin[1] -= mp_ref_xyz[1];
    grid_origin[2] -= mp_ref_xyz[2];
  }
 
	/* Computing one-electron integrals in 
	   terms of Cartesian Gaussians, 
	   electric first (dipole), second and third 
	   moments W.R.T. origin, electrostatic potential,
	   electric field and field gradients,
	   electron and spin density at atomic centers,
	   and various properties over a grid, if neccessary. */

 compute_oeprops();
 if (grid)
     compute_grid();

 print_mp();
 print_lm();
 if (nuc_esp)
   print_esp();
 if (grid) {
   grid_origin[0] += mp_ref_xyz[0];
   grid_origin[1] += mp_ref_xyz[1];
   grid_origin[2] += mp_ref_xyz[2];
   print_grid();
 }
 print_misc();
}

	/* Cleaning up */

 free_block(usotao);
 if (strcmp(ref,"UHF") != 0) {
   free_block(scf_evec_so);
   free_block(scf_evec_ao);
 }
 free(ioff);
 free(Ptot);
 free(phi);
 free(ex); free(ey); free(ez);
 free(dexx); free(deyy); free(dezz);
 free(dexy); free(dexz); free(deyz);
 if (spin_prop) {
   free(Pspin);
   free(ahfsxx);
   free(ahfsyy);
   free(ahfszz);
   free(ahfsxy);
   free(ahfsxz);
   free(ahfsyz);
 }
 free(S);
/*}  end non-UHF routines */
 stop_io();
 exit(PSI_RETURN_SUCCESS);
 
}


extern "C" const char *gprgid()
{
 const char *prgid = "OEPROP";
   
 return(prgid);
}
