/*! \file
    \ingroup CSCF
    \brief Enter brief description of file here
*/

/**************************************************************************/
/*                                                                        */
/*   CSCF:                                                                */
/*      Written by Edward Seidl (a NON-hog)                               */
/*      September 1990                                                    */
/*      Parts liberally ripped off from HFS group SCF                     */
/*        code written in FORTRAN (Boo!!!)                                */
/*                                                                        */
/*      modified April 17, 1991 to use new input format developed by      */
/*        Curtis Janssen                                                  */
/**************************************************************************/
/*                                                                        */
/*   Description of input                                                 */
/*                                                                        */
/*   LABEL = string                                                       */
/*         This is a character string to be included in the output.       */
/*         There is no default.                                           */
/*                                                                        */
/*   WFN = string                                                         */
/*         This is the type of wavefunction which is ultimately desired.  */
/*         The default is SCF.                                            */
/*                                                                        */
/*   DIRECT_SCF = boolean                                                 */
/*         Flag to request the direct formation of the Fock matrix        */
/*                                                                        */
/*   OPENTYPE = string                                                    */
/*         This specifies the state desired.  It can be one of NONE       */
/*         (for a closed shell singlet), SINGLET (for an open shell       */
/*         singlet), HIGHSPIN (for any high spin open shell system),      */
/*         TWOCON (for a two configuration singlet), or SPECIAL.          */
/*         If SPECIAL is given, then alpha and beta coupling              */
/*         coefficients must be given with the ALPHA and BETA keywords.   */
/*         The default is NONE.                                           */
/*                                                                        */
/*   DOCC = integer_vector                                                */
/*         This gives the number of doubly occupied orbitals in each      */
/*         irreducible representation.  There is no default.              */
/*                                                                        */
/*   SOCC = integer_vector                                                */
/*         This gives the number of singly occupied orbitals in each      */
/*         irreducible representation.  If OPENTYPE = NONE this defaults  */
/*         to the zero vector.  Otherwise, there is no default.           */
/*                                                                        */
/*   DERTYPE = string                                                     */
/*         This specifies the order of derivative that is to even-        */
/*         tually  be  done.   It  is  used  by the scf program to        */
/*         determine if certain files are to be written and it  is        */
/*         also  used  to determine the default convergence of the        */
/*         wavefunction.  The default is FIRST.                           */
/*                                                                        */
/*   MAXITER = integer                                                    */
/*         This gives the maximum number of iterations.  The default      */
/*         is 40.                                                         */
/*                                                                        */
/*   CONVERGENCE = integer                                                */
/*         The convergence criterion is 10**(-integer).  The default is   */
/*         7 if both DERTYPE = NONE and WFN = SCF are given and 10        */
/*         otherwise.                                                     */
/*                                                                        */
/*   LEVELSHIFT = real                                                    */
/*        This specifies the level shift. The default is 1.0.             */
/*                                                                        */
/*                                                                        */
/*   There are also a large number of less  commonly  used  input         */
/*   parameters.   If  you  do  not understand what the following         */
/*   options mean, then make sure that they do not appear in your         */
/*   input.   The defaults will work in the overwhelming majority         */
/*   of cases.  These are specified with the following keywords:          */
/*                                                                        */
/*                                                                        */
/*   REORDER = boolean                                                    */
/*        The molecular orbitals will be  reordered  if  this  is         */
/*        true,  in  which  case,  the  MOORDER parameter must be         */
/*        present.  The default is false.                                 */
/*                                                                        */
/*   MOORDER = integer_vector                                             */
/*        This specifies a molecular orbital  reordering  vector.         */
/*        It  will  only  be  used if REORDER = YES.  This vector         */
/*        contains first the ordering for  the  orbitals  in  the         */
/*        first  irreducible  representation  and then the second         */
/*        and so on.   The  first  orbital  of  each  irreducible         */
/*        representation is numbered 1.  There is no default.             */
/*                                                                        */
/*   ALPHA = real_vector                                                  */
/*        If OPENTYPE = SPECIAL, then this  parameter  gives  the         */
/*        alpha coupling coefficients.  The number of elements in         */
/*        this vector is MM(MM+1)/2, where MM is  the  number  of         */
/*        irreducible  representations containing singly occupied         */
/*        molecular orbitals.  There is no default.                       */
/*                                                                        */
/*   BETA = real_vector                                                   */
/*        If OPENTYPE = SPECIAL, then this  parameter  gives  the         */
/*        beta  coupling coefficients.  The number of elements in         */
/*        this vector is MM(MM+1)/2, where MM is  the  number  of         */
/*        irreducible  representations containing singly occupied         */
/*        molecular orbitals.  There is no default.                       */
/*                                                                        */
/*   RESTART = boolean                                                    */
/*        The calculation will restart from the old  wavefunction         */
/*        if RESTART is true.  If the old wavefunction does not           */
/*        exist, then the cscf program will generate its own  ini-        */
/*        tial  guess  automatically.   Possible  values for this         */
/*        parameter are TRUE, YES, 1,  FALSE,  NO,  and  0.   The         */
/*        default is true.                                                */
/*                                                                        */
/*   IPRINT = integer                                                     */
/*        This is a print option.  The default is 0.                      */
/*                                                                        */
/*   ROTATE = boolean                                                     */
/*        The molecular orbitals will not be rotated if  this  is         */
/*        false.   The rotation only affects the virtual orbitals         */
/*        for open shell systems.  This parameter  must  be  true         */
/*        for  correlated  gradients  and  it  must  be false for         */
/*        second and higher derivatives.  The default is false if         */
/*        WFN = SCF and true otherwise.                                   */
/*                                                                        */
/*   DIIS = boolean                                                       */
/*        This determines whether diis will be used.  The default is      */
/*        false for OPENTYPE = TWOCON and true otherwise.                 */
/*                                                                        */
/*   NDIIS = integer                                                      */
/*        This gives the number of error matrices to use in the diis      */
/*        procedure.  The default is 6 for closed shell, 4 for open       */
/*        shell, and 3 for tcscf.                                         */
/*                                                                        */
/*   DIISSTART = integer                                                  */
/*        This gives the first iteration for which DIIS  will  be         */
/*        used.  The default is 0.                                        */
/*                                                                        */
/*   DIISDAMP = real                                                      */
/*        This gives the damping factor for the diis procedure.  The      */
/*        default is 0.0 for closed shell, 0.02 for open shell, and       */
/*        0.01 for tcscf.                                                 */
/*                                                                        */
/*   INCR = real                                                          */
/*        This is used in tcscf to determine how often the ci             */
/*        coefficients are recalculated.  A small number (~0.25)          */
/*        will cause them to be recalculated nearly every scf             */
/*        iteration.  The default is 0.5.                                 */
/*                                                                        */
/*                                                                        */
/*   FOCK_TYPE = integer                                                  */
/*        Only used for tcscf and open shell calculations.                */
/*        If FOCK_TYPE = 0 use a simple form for fock_eff                 */
/*         "     "     = 1 use form of fock_eff suitable for              */
/*                         high-spin cases                                */
/*         "     "     > 1 experimental fock matrices                     */
/*                                                                        */
/**************************************************************************/



/*-------------------------------------------------------------------------
  READ THIS FIRST: This code is a hack of the original CSCF which employed
  symmetric orthogonalization procedure. This code uses the canonical
  orthogonalization procedure. To decrease the amount of
  rewriting I had to do I left the initialization part (init_scf and
  init_scf2) intact. Thus, although I have to allocate more space than
  I need, I should be fine otherwise.

                                             - Edward Valeev, August'99
 -------------------------------------------------------------------------*/


static char *rcsid = "$Id: cscf.cc 3955 2008-06-07 09:04:04Z rking $";

#include "includes.h"
#include "common.h"
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>
#include <psi4-dec.h>

namespace psi { namespace cscf {
  void print_initial_vec();
  extern void write_scf_matrices(void);

PsiReturnType cscf(Options & options)
{
  int i,j,nn;
  char *prog_name="CSCF3.0: An SCF program written in C";
  char *output="APPEND  ";
  struct symm *s;
  ip_value_t *ipvalue=NULL;
  int errcod, orthog_only, mo_print;
  char *wfn;

  twocut = 1.0;
  readflgc = readflgo = num_bufs_o = num_bufs_c = ibl = iblc = iblo = 0;
  where = wherec = whereo = 0;
  intmx = old_ninta = 25920;

  //errcod = psi_start(&infile,&outfile,&psi_file_prefix,argc-1, argv+1, 0);
  //if (errcod != PSI_RETURN_SUCCESS)
  //  exit(PSI_RETURN_FAILURE);
  ip_cwk_clear();
  ip_cwk_add(":DEFAULT");
  ip_cwk_add(":PSI");
  ip_cwk_add(":SCF");
  tstart();
//  psio_init();
//  psio_ipv1_config();
//
  fprintf(outfile,"%13c------------------------------------------\n",' ');
  fprintf(outfile,"\n%16c%s\n",' ',prog_name);
  fprintf(outfile,"\n%14cWritten by too many people to mention here\n",' ');
  fprintf(outfile,"\n%13c------------------------------------------\n",' ');

  itap30 = 30;
  itap33 = PSIF_SO_TEI;
  itap34 = 34;
  itapS  = PSIF_OEI;
  itapT  = PSIF_OEI;
  itapV  = PSIF_OEI;
  /*
    itapS  = PSIF_SO_S;
    itapT  = PSIF_SO_T;
    itapV  = PSIF_SO_V;
  */
  itapDSCF = PSIF_DSCF;
  itap92 = PSIF_SO_PKSUPER1;
  itap93 = PSIF_SO_PKSUPER2;

  /* JPK 6/1/00 integral accuracy: dynamic(default)=1, static=0 */
  dyn_acc = 1;
  eri_cutoff = 1.0E-14;
  ip_boolean("DYN_ACC",&dyn_acc,0);
  tight_ints=0;
  delta = 1.0;

  /* CDS 3/6/02 add flag to do only orthogonalization */
  errcod = ip_string("WFN",&wfn,0);
  if (strcmp(wfn,"DETCAS")==0 || strcmp(wfn,"CASSCF")==0 ||
      strcmp(wfn,"RASSCF")==0) orthog_only = 1;
  else orthog_only = 0;
  ip_boolean("ORTHOG_ONLY",&orthog_only,0);
  free(wfn);
  wfn = NULL;

  /* STB (6/30/99) - Function added because in order to initialize things
     one must know whether you are doing UHF or restricted */

  chkpt_init(PSIO_OPEN_OLD);

  occ_init();

  /* initialize some constants and arrays */
  if (uhf)
    init_uhf();
  else
    init_scf();

  /* read input.dat, get occupations, flags, etc. */

  scf_input(ipvalue);

  /* we can't just orthogonalize the orbitals if there aren't any */
  if (inflg != 1) orthog_only = 0;

  /* set up other useful arrays */

  init_scf2();

  /* get one electron integrals */

  rdone_iwl();

  if(print & 1) {
    for (i=0; i < num_ir; i++) {
      s = &scf_info[i];
      if (nn=s->num_so) {
    fprintf(outfile,"\nsmat for irrep %s\n",s->irrep_label);
    print_array(s->smat,nn,outfile);
    fprintf(outfile,"\ntmat for irrep %s\n",s->irrep_label);
    print_array(s->tmat,nn,outfile);
    fprintf(outfile,"\nhmat for irrep %s\n",s->irrep_label);
    print_array(s->hmat,nn,outfile);
      }
    }
  }

  /* form S-1/2 matrix sahalf */

  shalf();

  if (print & 1) {
    for (i=0; i < num_ir ; i++) {
      s = &scf_info[i];
      if (nn=s->num_so) {
    fprintf(outfile,"\nsahalf for irrep %s\n",s->irrep_label);
    print_mat(s->sahalf,nn,s->num_mo,outfile);
      }
    }
  }

  /* if no initial vector, form one from core hamiltonian */

  if (inflg == 2) form_vec();
  fflush(outfile);

  /* guess or designate orbital occupations*/
  guess();

  /* orthogonalize old vector and form first density matrix */
  if (inflg == 1) {
    if(uhf)
      schmit_uhf(1);
    else
      schmit(1);
  }

  /* Print out the first vector */
  if (print & 2)
    print_initial_vec();

  /* if we are only orthogonalizing, then quit here */
  if (orthog_only) {
    fprintf(outfile, "Only orbital orthogonalization has been performed\n");
    mo_print = 0;
    errcod = ip_boolean("PRINT_MOS",&mo_print,0);
    if (mo_print) print_mos("Alpha",scf_info);
    write_scf_matrices();
    chkpt_close();
    //psio_done();
    tstop();
    //psi_stop(infile,outfile,psi_file_prefix);
    //exit(PSI_RETURN_SUCCESS);
    return(Success);
  }

  if (!twocon){
    if(!uhf)
      dmat();
    else{
      /*--- Prepare alpha and beta eigenvectors from
    the core Hamiltonian guess ---*/
      if(inflg == 2)
    cmatsplit();
      dmatuhf();

    }
  }

  /* Decide how to form the Fock matrix */

  if (!direct_scf) {

    /* read in tei's and form supermatrix */

    num_ints = 0;
    num_bufs = 0;
    Pmat.unit = itap92;
    Pmat.key = strdup("P-supermatrix");
    Pmat.bufpos = PSIO_ZERO;
    PKmat.unit = itap93;
    PKmat.key = strdup("PK-supermatrix");
    PKmat.bufpos = PSIO_ZERO;
    psio_open(Pmat.unit,PSIO_OPEN_NEW);
    psio_open(PKmat.unit,PSIO_OPEN_NEW);
    rdtwo();

  }
  else {
    /* form the Fock matrix directly */
    /*check for rohf singlet...doesn't work direct*/
    if(!uhf && singlet) {
      fprintf(outfile,"\n  rohf open shell singlet doesn't work direct\n");
      fprintf(outfile,"  remove 'direct_scf = true' from input\n");
      fprintf(stderr,"rohf open shell singlet doesn't work direct\n");
      fprintf(stderr,"remove 'direct_scf = true' from input\n");
      chkpt_close();
      //psio_done();
      //exit(PSI_RETURN_FAILURE);
      return(Success);
    }

    formg_direct();
    if(dyn_acc)  fprintf(outfile,"\n  Using inexpensive integrals");
  }

  /* iterate */

  iter = 0;
  converged = 0;

  // these now don't call cleanup when converged
  if(twocon) scf_iter_2();
  else if(uhf) uhf_iter();
  else scf_iter();

  // free diis stuff - need to add UHF code later
  free(btemp); btemp = NULL;
  free_block(bold); bold = NULL;
  free_block(bmat);  bmat = NULL;
  if (refnum == 0) { //RHF
    for (i=0; i<ndiis; ++i) {
      for(j=0; j<num_ir; ++j) {
        if(scf_info[j].num_so) {
          free_block(diism[i].fock_c[j]);
          free_block(diism[i].error[j]);
          if (iopen) free_block(diism[i].fock_o[j]);
        }
      }
      free(diism[i].fock_c);
      if (iopen) free(diism[i].fock_o);
      free(diism[i].error);
    }
    free(diism); diism = NULL;
  }
  else if (refnum == 1) { //UHF
    for (i=0; i<ndiis; ++i) {
      for(j=0; j<num_ir; ++j) {
        if(scf_info[j].num_so) {
          free_block(udiism[i].fock[0][j]);
          free_block(udiism[i].fock[1][j]);
          free_block(udiism[i].error[0][j]);
          free_block(udiism[i].error[1][j]);
        }
      }
      free(udiism[i].fock[0]);
      free(udiism[i].fock[1]);
      free(udiism[i].fock);
      free(udiism[i].error[0]);
      free(udiism[i].error[1]);
      free(udiism[i].error);
    }
    free(udiism); udiism = NULL;
  }

  cleanup();

  JK = gmat = diis_out = NULL;

  dampsv = repnuc = etot = exc = exch_energy = corr_energy = coulomb_energy = 0.0;
  den_trace = lshift = diiser = save_ci1 = save_ci2 = dampd = dampo = eri_cutoff = 0.0;
  delta = alpha1 = alpha2 = alpha3 = eelec = dconv = 0.0;
  twocut = 1.0;

  stop_lshift = direct_scf = diisflg = scf_conv = iopen = inflg = hcore_guess = print = 0;
  fock_typ = ndiis = it_diis = itmax = use_iwl = delete_ints = delete_1e = delete_2e = 0;
  reset_occ = multp = mflag = charge = natom = nelec = nbfso = nmo = 0;
  phase_check = tight_ints = ok_ints = dyn_acc = acc_switch = 0;
  readflgc = readflgo = 0;
  num_bufs_o = num_bufs_c = 0;
  intmx = old_ninta = 25920;
  last = lasto = lastc = 0;

  if (reference != NULL) { free(reference); reference = NULL; }
  if (functional != NULL) { free(reference); reference = NULL; }

  exitflag = mo_out = n_so_typs = nbasis = nsfmax = n_closed = n_open = a_elec = 0;
  b_elec = num_ir = mxcoef2 = maxbuf = num_bufs = num_ints = iter = 0;
  converged = hsos = singlet = uhf = special = twocon = ksdft = mixing = cscf_nint = 0;
  opshl1 = opshl2 = opblk1 = opblk2 = second_root = icheck_rot = check_mo_orthonormality = 0;
  ibl = iblc = iblo = 0;

  itap30 = itap34 = itapS = itapT = itapV = itap33 = itap92 = itap93 = itapDSCF = 0;

  if (alpha != NULL) {free(alpha); alpha = NULL;}
  if (beta != NULL) {free(beta); beta = NULL;}
  if (zvals != NULL) {free(zvals); zvals = NULL;}
  if (symm_tot != NULL) {free(symm_tot); symm_tot = NULL;}
  if (ener_tot != NULL) {free(ener_tot); ener_tot = NULL;}
  if (i10 != NULL) {free(i10); i10 = NULL;}
  if (so2symblk != NULL) {free(so2symblk); so2symblk = NULL;}

  // Pmat PKmat?
  // need psi_buffer?
  //scf_info

  if (uhf && spin_info != NULL) {
    free(spin_info[0].scf_spin);  free(spin_info[1].scf_spin);
    free(spin_info[0].spinlabel); free(spin_info[1].spinlabel);
    free(spin_info);
  }

  if (c_outbuf != NULL) { free(c_outbuf); c_outbuf = NULL; }
  if (o_outbuf != NULL) { free(o_outbuf); o_outbuf = NULL; }


  // formg stuff
  where = wherec = whereo = 0;

  if (int_nums != NULL) { free(int_nums); int_nums = NULL; }
  if (int_nums_c != NULL) { free(int_nums_c); int_nums_c = NULL; }
  if (int_nums_o != NULL) { free(int_nums_o); int_nums_o = NULL; }
  if (inext != NULL) { free(inext); inext = NULL; }
  if (lbij != NULL) { free(lbij); lbij = NULL; }
  if (lbkl != NULL) { free(lbkl); lbkl = NULL; }

  if (gtmp != NULL) { free(gtmp); gtmp = NULL; }
  if (gtmp2 != NULL) { free(gtmp2); gtmp2 = NULL; }
  if (gtmpo != NULL) { free(gtmpo); gtmpo = NULL; }
  if (gtmpo2 != NULL) { free(gtmpo2); gtmpo2 = NULL; }
  if (ptmp != NULL) { free(ptmp); ptmp = NULL; }
  if (ptmp2 != NULL) { free(ptmp2); ptmp2 = NULL; }
  if (ptmpo != NULL) { free(ptmpo); ptmpo = NULL; }
  if (ptmpo2 != NULL) { free(ptmpo2); ptmpo2 = NULL; }

  if (pa != NULL) { free(pa); pa = NULL; }
  if (pb != NULL) { free(pb); pb = NULL; }
  if (testj != NULL) { free(testj); testj = NULL; }
  if (testk != NULL) { free(testk); testk = NULL; }

  chkpt_close();
  //psio_done();
  return(Success);
}

void print_initial_vec()
{
  int irrep, nso, nmo;
  struct symm *s;

  for(irrep=0;irrep<num_ir;irrep++) {
      s = &scf_info[irrep];
      if (nso=s->num_so)
      nmo = s->num_mo;
      if (uhf) {
          fprintf(outfile,"\n  Initial alpha vector for irrep %s\n",s->irrep_label);
          print_mat(spin_info[0].scf_spin[irrep].cmat,nso,nmo,outfile);
          fprintf(outfile,"\n  Initial beta vector for irrep %s\n",s->irrep_label);
          print_mat(spin_info[1].scf_spin[irrep].cmat,nso,nmo,outfile);
      }
      else {
          fprintf(outfile,"\nInitial vector for irrep %s\n",s->irrep_label);
          print_mat(s->cmat,nso,nmo,outfile);
      }
  }

  return;
}

}} // namespace psi::cscf
