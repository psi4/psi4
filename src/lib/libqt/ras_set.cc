/*!
  \file
  \brief Obtain orbital space and reordering for CI/MCSCF wavefunctions
  \ingroup QT
*/

#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include "qt.h"
#include <psifiles.h>
#include <psi4-dec.h>

namespace psi {

/*!
** ras_set(): Deprecated
**
** This function sets up the number of orbitals per irrep for each of the
** RAS subspaces [frozen core, RAS I, RAS II, RAS III, RAS IV, frozen virts].
** It also obtains the appropriate orbital reordering array.  The
** reordering array takes a basis function in Pitzer ordering (orbitals
** grouped according to irrep) and gives the corresponding index
** in the RAS numbering scheme.  Orbitals are numbered according to
** irrep within each of the subspaces.
**
** Formerly, docc, socc, frdocc, and fruocc were read in this function.
** Now docc and socc will be left as-is if they are not present in input.
**
** C. David Sherrill
** Center for Computational Quantum Chemistry
** University of Georgia, 25 June 1995
**
**  \param nirreps     =  num of irreps in computational point group
**  \param nbfso       =  num of basis functions in symmetry orbitals (num MOs)
**  \param freeze_core =  1 to remove frozen core orbitals from ras_opi
**  \param orbspi      =  array giving num symmetry orbitals (or MOs) per irrep
**  \param docc        =  array of doubly occupied orbitals per irrep
**  \param socc        =  array of singly occupied orbitals per irrep
**  \param frdocc      =  array of frozen core per irrep
**  \param fruocc      =  array of frozen virtuals per irrep
**  \param ras_opi     =  matrix giving num of orbitals per irrep per ras space,
**                        addressed as ras_opi[ras_space][irrep]
**  \param order       =  array nbfso big which maps Pitzer to Correlated order
**  \param ras_type    =  if 1, put docc and socc together in same RAS space
**                        (RAS I), as appropriate for DETCI.  If 0, put socc
**                        in its own RAS space (RAS II), as appropriate for CC.
**
** Returns: 1 for success, 0 otherwise
** \ingroup QT
*/
int ras_set(int nirreps, int nbfso, int freeze_core, int *orbspi,
            int *docc, int *socc, int *frdocc, int *fruocc,
            int **ras_opi, int *order, int ras_type)
{
    fprintf(outfile, "libqt::ras_set: Not converted to liboptions yet...beware!");
    return 1;
#if 0
  int i, irrep, point, tmpi, cnt=0;
  int errcod, errbad, parsed_ras1=0, parsed_ras2=0, do_ras4;
  int *used, *offset, **tras;
  int *tmp_frdocc, *tmp_fruocc;

  used = init_int_array(nirreps);
  offset = init_int_array(nirreps);

  /* if we have trouble reading DOCC and SOCC, we'll take them as
   * provided to this routine.  Zero out everything else
   */
  for (i=0; i<4; i++) {
    zero_int_array(ras_opi[i], nirreps);
  }
  zero_int_array(order, nbfso);


  /* now use the parser to get the arrays we require */
  tmp_frdocc = Process::environment.reference_wavefunction()->get_frzcpi();
  tmp_fruocc = Process::environment.reference_wavefunction()->get_frzvpi();
  for (i=0; i<nirreps; i++) {
    frdocc[i] = tmp_frdocc[i];
    fruocc[i] = tmp_fruocc[i];
  }
  free(tmp_frdocc);
  free(tmp_fruocc);

  /* replace DOCC and SOCC only if they are in input */
  if (ip_exist("DOCC",0))
    errcod = ip_int_array("DOCC",docc,nirreps);
  if (ip_exist("DOCC",0))
    errcod = ip_int_array("SOCC",socc,nirreps);

  errbad=0; do_ras4=1;
  errcod = ip_int_array("RAS1", ras_opi[0], nirreps);
  if (errcod == IPE_OK) parsed_ras1 = 1;
  else if (errcod == IPE_KEY_NOT_FOUND) {
    errcod = ip_int_array("ACTIVE_CORE", ras_opi[0], nirreps);
    if (errcod == IPE_OK) parsed_ras1 = 1;
    else if (errcod != IPE_KEY_NOT_FOUND) errbad = 1;
  }
  else errbad = 1;
  errcod = ip_int_array("RAS2", ras_opi[1], nirreps);
  if (errcod == IPE_OK) parsed_ras2 = 1;
  else if (errcod == IPE_KEY_NOT_FOUND) {
    errcod = ip_int_array("MODEL_SPACE", ras_opi[1], nirreps);
    if (errcod == IPE_OK) parsed_ras2 = 1;
    else if (errcod != IPE_KEY_NOT_FOUND) errbad = 1;
  }
  else errbad = 1;
  errcod = ip_int_array("RAS3", ras_opi[2], nirreps);
  if (errcod != IPE_OK && errcod != IPE_KEY_NOT_FOUND) errbad=1;
  if (errcod == IPE_KEY_NOT_FOUND) do_ras4=0;

  if (errbad == 1) {
    fprintf(stderr, "(ras_set): trouble parsing RAS keyword\n");
    return(0);
  }

  /* if the user has not specified RAS I, we must deduce it.
   * RAS I does not include any FZC orbs but does include COR orbs
   */

  if (!parsed_ras1) {
    for (irrep=0; irrep<nirreps; irrep++) {
      if (ras_type==1) ras_opi[0][irrep] = docc[irrep] + socc[irrep];
      if (ras_type==0) ras_opi[0][irrep] = docc[irrep];
      ras_opi[0][irrep] -= frdocc[irrep]; /* add back later for COR */
    }
  }


  /* if the user hasn't specified RAS II, look for val_orb             */
  /* val_orb should typically be RAS I + RAS II, so subtract out RAS I */
  if (!parsed_ras2) {
    errcod = ip_int_array("VAL_ORB",ras_opi[1],nirreps);
    if (errcod != IPE_OK) {
      for (irrep=0; irrep<nirreps; irrep++) {
    if (ras_type==1) ras_opi[1][irrep] = 0;
    if (ras_type==0) ras_opi[1][irrep] = socc[irrep];
      }
    }
    else {
      for (irrep = 0; irrep<nirreps; irrep++) {
         ras_opi[1][irrep] -= ras_opi[0][irrep];
         if (ras_opi[1][irrep] < 0) {
           fprintf(stderr, "(ras_set): val_orb must be larger than RAS I\n");
           return(0);
         }
      }
    }
  }

  /* set up the RAS III or IV array: if RAS IV is used, RAS III must
   * be specified and then RAS IV is deduced.
   */

  for (irrep=0; irrep<nirreps; irrep++) {
    tmpi = frdocc[irrep] + fruocc[irrep] + ras_opi[0][irrep] +
      ras_opi[1][irrep];
    if (do_ras4) tmpi += ras_opi[2][irrep];
    if (tmpi > orbspi[irrep]) {
      fprintf(stderr, "(ras_set): orbitals don't add up for irrep %d\n",
          irrep);
      return(0);
    }
    if (do_ras4) ras_opi[3][irrep] = orbspi[irrep] - tmpi;
    else ras_opi[2][irrep] = orbspi[irrep] - tmpi;
  }

  /* copy everything to the temporary RAS arrays: */
  /* add subspaces for frozen orbitals            */
  tras = init_int_matrix(6, nirreps);
  for (irrep=0; irrep<nirreps; irrep++) {
    tras[0][irrep] = frdocc[irrep];
    tras[5][irrep] = fruocc[irrep];
  }
  for (i=0; i<4; i++) {
    for (irrep=0; irrep<nirreps; irrep++) {
      tras[i+1][irrep] = ras_opi[i][irrep];
    }
  }

  /* construct the offset array */
  offset[0] = 0;
  for (irrep=1; irrep<nirreps; irrep++) {
    offset[irrep] = offset[irrep-1] + orbspi[irrep-1];
  }

  for (i=0; i<6; i++) {
    for (irrep=0; irrep<nirreps; irrep++) {
      while (tras[i][irrep]) {
    point = used[irrep] + offset[irrep];
    if (point < 0 || point >= nbfso) {
      fprintf(stderr, "(ras_set): Invalid point value\n");
      exit(PSI_RETURN_FAILURE);
    }
    order[point] = cnt++;
    used[irrep]++;
    tras[i][irrep]--;
      }
    }
  }


  /* do a final check */

  for (irrep=0; irrep<nirreps; irrep++) {
    if (used[irrep] > orbspi[irrep]) {
      fprintf(stderr, "(ras_set): on final check, used more orbitals");
      fprintf(stderr, "   than were available (%d vs %d) for irrep %d\n",
          used[irrep], orbspi[irrep], irrep);
      errbad = 1;
    }
    if (used[irrep] < orbspi[irrep]) {
      fprintf(stderr, "(ras_set): on final check, used fewer orbitals");
      fprintf(stderr, "   than were available (%d vs %d) for irrep %d\n",
          used[irrep], orbspi[irrep], irrep);
      errbad = 1;
    }
  }

  /* for restricted COR orbitals */
  if (!freeze_core) {
    for (irrep=0; irrep<nirreps; irrep++) {
      ras_opi[0][irrep] += frdocc[irrep];
    }
  }

  free(used);  free(offset);
  free_int_matrix(tras);

  return(!errbad);
#endif
}


/*!
** ras_set2()
**
** This function sets up the number of orbitals per irrep for each of the
** RAS subspaces [frozen core (FZC), restricted core (COR), RAS I, RAS II,
** RAS III, RAS IV, restricted virts (VIR), frozen virts (FZV)].
** It also obtains the appropriate orbital reordering array.  The
** reordering array takes a basis function in Pitzer ordering (orbitals
** grouped according to irrep) and gives the corresponding index
** in the RAS numbering scheme.  Orbitals are numbered according to
** irrep within each of the subspaces.
**
** Formerly, docc, socc, frdocc, and fruocc were read in this function.
** Now docc and socc will be left as-is if they are not present in input.
**
** Assume we always want integrals (at least some of them...) involving
** restricted orbitals, but we may not need them for frozen orbitals unless
** perhaps it's a gradient.  The frozen orbitals will never enter explicitly
** in DETCI, but restricted orbitals may or may not, depending on the
** type of computation (CI vs CAS vs CASPT2, etc).  For conventional
** CI, FCI, MRCI, RASCI, CASSCF, restricted orbitals will not participate
** explicitly in DETCI.  For CASPT2, perhaps they will.  CLAG and DETCAS
** will still need to have some indices in the restricted (and possibly
** frozen) orbital subspaces?
**
** C. David Sherrill
**
** Updated June 2002 to distinguish between frozen and restricted spaces
**
**  \param nirreps     =  num of irreps in computational point group
**  \param nmo         =  number of MO's
**  \param delete_fzdocc    = 1 to remove frozen core orbitals from ras_opi[0]
**  \param delete_restrdocc = 1 to remove restricted core orbs from ras_opi[0]
**  \param orbspi      =  array giving num symmetry orbitals (or MOs) per irrep
**  \param docc        =  array of doubly occupied orbitals per irrep
**                        (guess should be provided)
**  \param socc        =  array of singly occupied orbitals per irrep
**                        (guess should be provided)
**  \param frdocc      =  array of frozen core per irrep
**                        (returned by function, but allocate before call)
**  \param fruocc      =  array of frozen virtuals per irrep
**                        (returned by function, but allocate before call)
**  \param rstrdocc    =  array of restricted core per irrep
**                        (returned by function, but allocate before call)
**  \param rstruocc    =  array of restricted core per irrep
**                        (returned by function, but allocate before call)
**  \param ras_opi     =  matrix giving num of orbitals per irrep per ras space,
**                        addressed as ras_opi[ras_space][irrep]
**                        (returned by function, but allocate before call)
**  \param order       =  array nmo big which maps Pitzer to Correlated order
**                        (returned by function, but allocate before call)
**  \param ras_type    =  if 1, put docc and socc together in same RAS space
**                        (RAS I), as appropriate for DETCI.  If 0, put socc
**                        in its own RAS space (RAS II), as appropriate for CC.
**  \param hoffmann    =  if > 0, order orbitals according to Mark Hoffmann.
**                        hoffmann==1:
**                        ras1, ras2, ..., rasn, COR, FZC, VIR, FZV.
**                        hoffmann==2:
**                        VIR, ras1, ras2, ..., rasn, COR, FZC, FZV.
**                        Note odd placement of FZC in middle!
**
** Returns: 1 for success, 0 otherwise
** \ingroup QT
*/
int ras_set2(int nirreps, int nmo, int delete_fzdocc,
             int delete_restrdocc, int *orbspi,
             int *docc, int *socc, int *frdocc, int *fruocc,
             int *restrdocc, int *restruocc, int **ras_opi, int *order,
             int ras_type, int hoffmann)
{
    fprintf(outfile, "libqt::ras_set2: has not be converted to liboptions...beware!");
    return 1;
#if 0
  int i, irrep, point, tmpi, cnt=0;
  int errcod, errbad, parsed_ras1=0, parsed_ras2=0, do_ras4;
  int parsed_restr_uocc=0;
  int *used, *offset, **tras;
  int *tmp_frdocc, *tmp_fruocc;

  /* Hoffmann's code never wants FZC or COR lumped in with RAS I,
     even if those orbitals need to be transformed also */
  if (hoffmann > 0) {
    delete_fzdocc = 1;
    delete_restrdocc = 1;
  }

  used = init_int_array(nirreps);
  offset = init_int_array(nirreps);

  /* if we have trouble reading DOCC and SOCC, we'll take them as
   * provided to this routine.  Zero out everything else
   */
  zero_int_array(restrdocc, nirreps);
  zero_int_array(restruocc, nirreps);

  for (i=0; i<MAX_RAS_SPACES; i++) {
    zero_int_array(ras_opi[i], nirreps);
  }
  zero_int_array(order, nmo);

  /* now use the parser to get the arrays we require */
  tmp_frdocc = get_frzcpi();
  tmp_fruocc = get_frzvpi();
  for (i=0; i<nirreps; i++) {
    frdocc[i] = tmp_frdocc[i];
    fruocc[i] = tmp_fruocc[i];
  }
  free(tmp_frdocc);
  free(tmp_fruocc);


  /* replace DOCC and SOCC only if they are in input */
  /* this fills existing DOCC and SOCC arrays        */
  if (ip_exist("DOCC",0))
    errcod = ip_int_array("DOCC",docc,nirreps);
  if (ip_exist("SOCC",0))
    errcod = ip_int_array("SOCC",socc,nirreps);

  /* now use the parser to get the arrays we require */
  errcod = ip_int_array("RESTRICTED_DOCC",restrdocc,nirreps);

  errbad=0; do_ras4=1;
  errcod = ip_int_array("RAS1", ras_opi[0], nirreps);
  if (errcod == IPE_OK) parsed_ras1 = 1;
  else if (errcod != IPE_KEY_NOT_FOUND) errbad = 1;
  errcod = ip_int_array("RAS2", ras_opi[1], nirreps);
  if (errcod == IPE_OK) parsed_ras2 = 1;
  else if (errcod != IPE_KEY_NOT_FOUND) errbad = 1;
  errcod = ip_int_array("RAS3", ras_opi[2], nirreps);
  if (errcod != IPE_OK && errcod != IPE_KEY_NOT_FOUND) errbad=1;
  if (errcod == IPE_KEY_NOT_FOUND) do_ras4=0;

  if (errbad == 1) {
    fprintf(stderr, "(ras_set): trouble parsing RAS keyword\n");
    return(0);
  }

  /*
   * If the user has not specified RAS I, we must deduce it.
   * As of 2004, RAS I will not include any FZC ("frozen core")
   * orbitals *OR* any COR ("restricted core") orbitals.  However,
   * we'll add these back into RAS I later if the user requests it
   * by setting delete_fzdocc or delete_rstrdocc = 0.
   */

  if (!parsed_ras1) {
    for (irrep=0; irrep<nirreps; irrep++) {
      if (ras_type==1) ras_opi[0][irrep] = docc[irrep] + socc[irrep];
      if (ras_type==0) ras_opi[0][irrep] = docc[irrep];
      ras_opi[0][irrep] -= (frdocc[irrep] + restrdocc[irrep]);
    }
  }


  /*
   * if the user hasn't specified RAS II, look for ACTIVE (replaces
   * VAL_ORB).  Assume this always goes after FZC + COR.
   */
  if (!parsed_ras2) {
    errcod = ip_int_array("ACTIVE",ras_opi[1],nirreps);
    if (errcod != IPE_OK) { /* we didn't find ACTIVE keyword */
      for (irrep=0; irrep<nirreps; irrep++) {
    if (ras_type==1) ras_opi[1][irrep] = 0;
    if (ras_type==0) ras_opi[1][irrep] = socc[irrep];
      }
    }
    else { /* we found a ACTIVE keyword */
      if (parsed_ras1)
        fprintf(outfile, "ras_set(): Warning: ACTIVE overrides RAS 1\n");

      for (irrep=0; irrep<nirreps; irrep++)
        ras_opi[0][irrep] = 0; /* ACTIVE overrides RAS 1 */

      if (!ip_exist("RESTRICTED_UOCC",0)) { /* default restrict other virs */
        for (irrep=0; irrep<nirreps; irrep++) {
          ras_opi[0][irrep] = 0; /* ACTIVE overrides RAS 1 */
          restruocc[irrep] = orbspi[irrep] - frdocc[irrep] -
                             restrdocc[irrep] - ras_opi[0][irrep] -
                             ras_opi[1][irrep] - fruocc[irrep];
        }
        parsed_restr_uocc = 1;
      }

      /* do a quick check */
      for (irrep=0; irrep<nirreps; irrep++) {
        if (frdocc[irrep]+restrdocc[irrep]+ras_opi[0][irrep]+ras_opi[1][irrep]
            < docc[irrep] + socc[irrep])
          fprintf(outfile,
                  "ras_set():Warning:Occupied electrons beyond ACTIVE orbs!\n");
      }
    }
  } /* end looking for "ACTIVE" */

  if (!parsed_restr_uocc)
    errcod = ip_int_array("RESTRICTED_UOCC",restruocc,nirreps);

  /* set up the RAS III or IV array: if RAS IV is used, RAS III must
   * be specified and then RAS IV is deduced.
   */

  for (irrep=0; irrep<nirreps; irrep++) {
    tmpi = frdocc[irrep] + restrdocc[irrep] + fruocc[irrep] + restruocc[irrep]
           + ras_opi[0][irrep] + ras_opi[1][irrep];
    if (do_ras4) tmpi += ras_opi[2][irrep];
    if (tmpi > orbspi[irrep]) {
      fprintf(stderr, "(ras_set): orbitals don't add up for irrep %d\n",
          irrep);
      return(0);
    }
    if (do_ras4) ras_opi[3][irrep] = orbspi[irrep] - tmpi;
    else ras_opi[2][irrep] = orbspi[irrep] - tmpi;
  }

  /* copy everything to the temporary RAS arrays: */
  /* add subspaces for frozen orbitals            */
  tras = init_int_matrix(MAX_RAS_SPACES+4, nirreps);

  /* for usual DETCI, DETCAS, etc */
  if (hoffmann==0) {
    for (irrep=0; irrep<nirreps; irrep++) {
      tras[0][irrep] = frdocc[irrep];
      tras[1][irrep] = restrdocc[irrep];
      tras[MAX_RAS_SPACES+2][irrep] = restruocc[irrep];
      tras[MAX_RAS_SPACES+3][irrep] = fruocc[irrep];
    }
    for (i=0; i<MAX_RAS_SPACES; i++) {
      for (irrep=0; irrep<nirreps; irrep++) {
        tras[i+2][irrep] = ras_opi[i][irrep];
      }
    }
  }
  /* for Mark Hoffmann's UND order for GVVPT2, etc. */
  else if (hoffmann == 1) {
    for (i=0; i<MAX_RAS_SPACES; i++) {
      for (irrep=0; irrep<nirreps; irrep++) {
        tras[i][irrep] = ras_opi[i][irrep];
      }
    }
    for (irrep=0; irrep<nirreps; irrep++) {
      tras[MAX_RAS_SPACES+0][irrep] = restrdocc[irrep];
      tras[MAX_RAS_SPACES+1][irrep] = frdocc[irrep];
      tras[MAX_RAS_SPACES+2][irrep] = restruocc[irrep];
      tras[MAX_RAS_SPACES+3][irrep] = fruocc[irrep];
    }
  }
  /* for Mark Hoffmann's UND GUGA code */
  else if (hoffmann == 2) {
    for (i=0; i<MAX_RAS_SPACES; i++) {
      for (irrep=0; irrep<nirreps; irrep++) {
        tras[i+1][irrep] = ras_opi[i][irrep];
      }
    }
    for (irrep=0; irrep<nirreps; irrep++) {
      tras[0][irrep] = restruocc[irrep];
      tras[1+MAX_RAS_SPACES][irrep] = restrdocc[irrep];
      tras[2+MAX_RAS_SPACES][irrep] = frdocc[irrep];
      tras[MAX_RAS_SPACES+3][irrep] = fruocc[irrep];
    }
  }

  /* construct the offset array */
  offset[0] = 0;
  for (irrep=1; irrep<nirreps; irrep++) {
    offset[irrep] = offset[irrep-1] + orbspi[irrep-1];
  }

  for (i=0; i<MAX_RAS_SPACES+4; i++) {
    for (irrep=0; irrep<nirreps; irrep++) {
      while (tras[i][irrep]) {
    point = used[irrep] + offset[irrep];
    if (point < 0 || point >= nmo) {
      fprintf(stderr, "(ras_set): Invalid point value\n");
      abort();
    }
    order[point] = cnt++;
    used[irrep]++;
    tras[i][irrep]--;
      }
    }
  }


  /* do a final check */
  for (irrep=0; irrep<nirreps; irrep++) {
    if (used[irrep] > orbspi[irrep]) {
      fprintf(stderr, "(ras_set): on final check, used more orbitals");
      fprintf(stderr, "   than were available (%d vs %d) for irrep %d\n",
          used[irrep], orbspi[irrep], irrep);
      errbad = 1;
    }
    if (used[irrep] < orbspi[irrep]) {
      fprintf(stderr, "(ras_set): on final check, used fewer orbitals");
      fprintf(stderr, "   than were available (%d vs %d) for irrep %d\n",
          used[irrep], orbspi[irrep], irrep);
      errbad = 1;
    }
  }

  /* add back FZC and COR orbitals to RAS I if they are to be included */
  for (irrep=0; irrep<nirreps; irrep++) {
    if (!delete_restrdocc) ras_opi[0][irrep] += restrdocc[irrep];
    if (!delete_fzdocc)    ras_opi[0][irrep] += frdocc[irrep];
  }

  free(used);  free(offset);
  free_int_matrix(tras);

  return(!errbad);
#endif
}

}

