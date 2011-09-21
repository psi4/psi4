/*!
  \file
  \brief Obtain the QT orbital reordering array between Pitzer and correlated
    order
  \ingroup QT
*/

#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>

namespace psi {

/*!
** reorder_qt()
**
** This function constructs a reordering array according to the
** "Quantum Trio" standard ordering, in which the orbitals are divided
** into the following sets: frozen core, then doubly occupied, then singly
** occupied, then virtuals, then deleted (frozen) virtuals.
** The reordering array takes a basis function in
** Pitzer ordering (orbitals grouped according to irrep) and gives the
** corresponding index in the Quantum Trio numbering scheme.
**
** Should give the same reordering array as in the old libread30 routines.
**
** C. David Sherrill
** Center for Computational Quantum Chemistry
** University of Georgia, 1995
**
** \param docc_in        = doubly occupied orbitals per irrep
** \param socc_in        = singly occupied orbitals per irrep
** \param frozen_docc_in = frozen occupied orbitals per irrep
** \param frozen_uocc_in = frozen unoccupied orbitals per irrep
** \param order          = reordering array (Pitzer->QT order)
** \param nirreps        = number of irreducible representations
**
** \ingroup QT
*/
void reorder_qt(int *docc_in, int *socc_in, int *frozen_docc_in,
      int *frozen_uocc_in, int *order, int *orbs_per_irrep, int nirreps)
{

   int cnt=0, irrep, point, tmpi;
   int *used, *offset;
   int *docc, *socc, *frozen_docc, *frozen_uocc;
   int *uocc;

   used = init_int_array(nirreps);
   offset = init_int_array(nirreps);

   docc = init_int_array(nirreps);
   socc = init_int_array(nirreps);
   frozen_docc = init_int_array(nirreps);
   frozen_uocc = init_int_array(nirreps);
   uocc = init_int_array(nirreps);

   for (irrep=0; irrep<nirreps; irrep++) {
      docc[irrep] = docc_in[irrep];
      socc[irrep] = socc_in[irrep];
      frozen_docc[irrep] = frozen_docc_in[irrep];
      frozen_uocc[irrep] = frozen_uocc_in[irrep];
      }

   /* construct the offset array */
   offset[0] = 0;
   for (irrep=1; irrep<nirreps; irrep++) {
      offset[irrep] = offset[irrep-1] + orbs_per_irrep[irrep-1];
      }

   /* construct the uocc array */
   for (irrep=0; irrep<nirreps; irrep++) {
      tmpi = frozen_uocc[irrep] + docc[irrep] + socc[irrep];
      if (tmpi > orbs_per_irrep[irrep]) {
         fprintf(stderr, "(reorder_qt): orbitals don't add up for irrep %d\n",
            irrep);
         return;
         }
      else
         uocc[irrep] = orbs_per_irrep[irrep] - tmpi;
      }

   /* do the frozen core */
   for (irrep=0; irrep<nirreps; irrep++) {
      while (frozen_docc[irrep]) {
         point = used[irrep] + offset[irrep];
         order[point] = cnt++;
         used[irrep]++;
         frozen_docc[irrep]--;
         docc[irrep]--;
         }
      }

   /* do doubly occupied orbitals */
   for (irrep=0; irrep<nirreps; irrep++) {
      while (docc[irrep]) {
         point = used[irrep] + offset[irrep];
         order[point] = cnt++;
         used[irrep]++;
         docc[irrep]--;
         }
      }

   /* do singly-occupied orbitals */
   for (irrep=0; irrep<nirreps; irrep++) {
      while (socc[irrep]) {
         point = used[irrep] + offset[irrep];
         order[point] = cnt++;
         used[irrep]++;
         socc[irrep]--;
         }
      }

   /* do virtual orbitals */
   for (irrep=0; irrep<nirreps; irrep++) {
      while (uocc[irrep]) {
         point = used[irrep] + offset[irrep];
         order[point] = cnt++;
         used[irrep]++;
         uocc[irrep]--;
         }
      }

   /* do frozen uocc */
   for (irrep=0; irrep<nirreps; irrep++) {
      while (frozen_uocc[irrep]) {
         point = used[irrep] + offset[irrep];
         order[point] = cnt++;
         used[irrep]++;
         frozen_uocc[irrep]--;
         }
      }


   /* do a final check */
   for (irrep=0; irrep<nirreps; irrep++) {
      if (used[irrep] > orbs_per_irrep[irrep]) {
         fprintf(stderr, "(reorder_qt): on final check, used more orbitals");
         fprintf(stderr, "   than were available (%d vs %d) for irrep %d\n",
            used[irrep], orbs_per_irrep[irrep], irrep);
         }
      }

   free(used);  free(offset);
   free(docc);  free(socc);  free(frozen_docc);  free(frozen_uocc);
   free(uocc);
}

/*!
** reorder_qt_uhf()
**
** Generalization of reorder_qt() for UHF case
**
** \param docc        = doubly occupied orbitals per irrep
** \param socc        = singly occupied orbitals per irrep
** \param frozen_docc = frozen occupied orbitals per irrep
** \param frozen_uocc = frozen unoccupied orbitals per irrep
** \param order_alpha = reordering array for alpha (Pitzer->QT order)
** \param order_beta  = reordering array for beta  (Pitzer->QT order)
** \param nirreps     = number of irreducible representations
**
** \ingroup QT
*/
void reorder_qt_uhf(int *docc, int *socc, int *frozen_docc,
                    int *frozen_uocc, int *order_alpha, int *order_beta,
                    int *orbspi, int nirreps)
{
  int p, nmo;
  int cnt_alpha, cnt_beta, irrep, tmpi;
  int *offset, this_offset;
  int *uocc;

  offset = init_int_array(nirreps);

  uocc = init_int_array(nirreps);

  /* construct the offset array */
  offset[0] = 0;
  for (irrep=1; irrep<nirreps; irrep++) {
    offset[irrep] = offset[irrep-1] + orbspi[irrep-1];
  }

  /* construct the uocc array */
  nmo = 0;
  for (irrep=0; irrep<nirreps; irrep++) {
    nmo += orbspi[irrep];
    tmpi = frozen_uocc[irrep] + docc[irrep] + socc[irrep];
    if (tmpi > orbspi[irrep]) {
      fprintf(stderr, "(reorder_qt_uhf): orbitals don't add up for irrep %d\n",
              irrep);
      return;
    }
    else
      uocc[irrep] = orbspi[irrep] - tmpi;
  }

  cnt_alpha = cnt_beta = 0;

  /* do the frozen core */
  for(irrep=0; irrep<nirreps; irrep++) {
    this_offset = offset[irrep];
    for(p=0; p < frozen_docc[irrep]; p++) {
      order_alpha[this_offset+p] = cnt_alpha++;
      order_beta[this_offset+p] = cnt_beta++;
    }
  }

  /* alpha occupied orbitals */
  for(irrep=0; irrep<nirreps; irrep++) {
    this_offset = offset[irrep] + frozen_docc[irrep];
    for(p=0; p < docc[irrep] + socc[irrep] - frozen_docc[irrep]; p++) {
      order_alpha[this_offset + p] = cnt_alpha++;
    }
  }

  /* beta occupied orbitals */
  for(irrep=0; irrep<nirreps; irrep++) {
    this_offset = offset[irrep] + frozen_docc[irrep];
    for(p=0; p < docc[irrep] - frozen_docc[irrep]; p++) {
      order_beta[this_offset + p] = cnt_beta++;
    }
  }

  /* alpha unoccupied orbitals */
  for(irrep=0; irrep<nirreps; irrep++) {
    this_offset = offset[irrep] + docc[irrep] + socc[irrep];
    for(p=0; p < uocc[irrep]; p++) {
      order_alpha[this_offset + p] = cnt_alpha++;
    }
  }

  /* beta unoccupied orbitals */
  for(irrep=0; irrep<nirreps; irrep++) {
    this_offset = offset[irrep] + docc[irrep];
    for(p=0; p < uocc[irrep] + socc[irrep]; p++) {
      order_beta[this_offset + p] = cnt_beta++;
    }
  }

  /* do the frozen uocc */
  for (irrep=0; irrep<nirreps; irrep++) {
    this_offset = offset[irrep] + docc[irrep] + socc[irrep] + uocc[irrep];
    for(p=0; p < frozen_uocc[irrep]; p++) {
      order_alpha[this_offset + p] = cnt_alpha++;
      order_beta[this_offset + p] = cnt_beta++;
    }
  }

  /* do a final check */
  for (irrep=0; irrep<nirreps; irrep++) {
    if (cnt_alpha > nmo) {
      fprintf(stderr, "(reorder_qt_uhf): on final check, used more orbitals");
      fprintf(stderr, "   than were available (%d vs %d) for irrep %d\n",
              cnt_alpha, nmo, irrep);
    }
    if (cnt_beta > nmo) {
      fprintf(stderr, "(reorder_qt_uhf): on final check, used more orbitals");
      fprintf(stderr, "   than were available (%d vs %d) for irrep %d\n",
              cnt_beta, nmo, irrep);
    }
  }

  free(offset);
  free(uocc);
}

}

