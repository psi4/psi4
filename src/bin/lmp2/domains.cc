/*! \file
    \ingroup LMP2
    \brief localized the SCF MO's
 */
//#include <cstdio>
//#include <cstdlib>
//#include <cstring>
#include <cmath>
//#include <iostream>
//#include <fstream>              // file I/O support
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libmints/mints.h>
#include <libparallel/parallel.h>
//#include <libchkpt/chkpt.h>
//#include <libpsio/psio.h>
#include <libqt/qt.h>
//#include <psifiles.h>
//#include <libiwl/iwl.h>
//#include <libint/libint.h>
//#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

namespace psi {

namespace lmp2 {

/*!
 ** domain(): Set up parameters of localization domains.
 **
 ** The orbital domains constructed here are based on those described
 ** in Boughton and Pulay, J. Comp. Chem. 14, 736-740 (1993).  The
 ** localization of the occupied orbitals is done elsewhere (see the
 ** program "localized Pair domains are defined as the union of pairs
 ** of single occupied orbital domains.  "Weak pairs", which are
 ** defined as pair domains whose individual occupied orbital domains
 ** have no atoms in common, are identified (cf. int *weak_pairs).
 **
 ** TDC, Jan-June 2002
 */

void LMP2::domains() {

    int i, j, k, l, ij, m, n, a;
    double xcoord, ycoord, zcoord;
    int next_atom, row, col, errcod;
    int max, nvir, nfzc;
    int atom, am, offset, shell_length;
    int *rank, *boolean, *ipiv;
    int **domain_bp, *domain_len_bp;
    int *distant_pairs, weak;
    int num_entries, entry_len, orbital;
    int print_test;
    double cutoff, dist;
    double **X, **Y, *fR;
    double *eps_all; // All MO energies
    double *charge, *SR, *Z, tmp;
    double **geometry;

    geometry = block_matrix(natom,3);
    geometry = get_geom();
    //for(i=0; i < natom; i++) {
    //    std::cout << "Cartesian coordinates for atom " << i << std::endl;
    //    std::cout << "x = " << geometry[i][0] << "\ty = " << geometry[i][1] << "\tZ = " << geometry[i][2] << std::endl;
    //}

    nfzc = get_frdocc();

    cutoff = 0.02;
    ip_data("LOCAL_CUTOFF", "%lf", &(cutoff), 0);

    nvir = nso - nocc;

    /* A couple of scratch arrays */
    X = block_matrix(nso, nso);
    Y = block_matrix(nso, nso);

    //  fprintf(outfile, "\n\tAO Overlap (S)\n");
    //  print_mat(aoovlp, nso, nso, outfile);

    //  fprintf(outfile, "\n\tAO Density Matrix (D)\n");
    //  print_mat(D, nso, nso, outfile);

    /************* Build the orbital domains ************/

    domain = init_int_matrix(nocc, natom);
    domain_len = init_int_array(nocc);
    domain_bp = init_int_matrix(nocc, natom);
    domain_len_bp = init_int_array(nocc);
    charge = init_array(natom);
    rank = init_int_array(natom);
    boolean = init_int_array(natom);
    SR = init_array(nso);
    Z = init_array(nso);
    ipiv = init_int_array(nso);
    fR = init_array(nocc);

    for (i = 0; i < nocc; i++) {

        // Compute the contribution of each atom to this orbital's charge/population
        for (j = 0; j < natom; j++) {
            charge[j] = 0.0;
            for (k = aostart[j]; k <= aostop[j]; k++) {
                tmp = 0.0;
                for (l = 0; l < nso; l++) tmp += aoovlp[k][l] * C[l][i];
                tmp *= C[k][i];
                charge[j] += tmp;
            }
        }

        //    for(j=0; j < natom; j++) {
        //      fprintf(outfile, "charge[%d] = %20.12f\n", j, charge[j]);
        //    }


        // Rank the atomic contributions to the orbital's charge
        for (j = 0; j < natom; j++) {
            rank[j] = 0;
            boolean[j] = 0;
        }
        for (j = 0, max = 0; j < natom; j++) // find the overall maximum
            if (fabs(charge[j]) >= fabs(charge[max])) max = j;
        rank[0] = max;
        boolean[max] = 1;
        for (j = 1; j < natom; j++) {
            max = 0;
            while (boolean[max]) max++; // find an unused max
            for (k = 0; k < natom; k++)
                if ((fabs(charge[k]) >= fabs(charge[max])) && !boolean[k]) max = k;
            rank[j] = max;
            boolean[max] = 1;
        }

        //    for(j=0; j < natom; j++) {
        //      fprintf(outfile, "rank[%d] = %d\n", j, rank[j]);
        //    }


        // Build the orbital's domain starting in order of decreasing charge contribution
        for (j = 0; j < nso; j++) {
            SR[j] = 0.0;
            for (k = 0; k < nso; k++)
                SR[j] += aoovlp[j][k] * C[k][i];
        }

        domain[i - nfzc][rank[0]] = 1; // at least one atom must be in the domain
        domain_len[i - nfzc] = 1;

        fR[i - nfzc] = 1.0;
        next_atom = 1;
        while (fabs(fR[i - nfzc]) > cutoff) {

            // Completeness check
            for (j = 0, row = 0; j < natom; j++) {
                if (domain[i - nfzc][j]) {
                    for (k = aostart[j]; k <= aostop[j]; k++, row++) {

                        Z[row] = SR[k];

                        for (l = 0, col = 0; l < natom; l++) {
                            if (domain[i - nfzc][l]) {

                                for (m = aostart[l]; m <= aostop[l]; m++, col++)
                                    X[row][col] = aoovlp[k][m];

                            }
                        } // l

                    } // k
                }
            } // j

            errcod = C_DGESV(row, 1, &(X[0][0]), nso, &(ipiv[0]), &(Z[0]), nso);
            if (errcod) {
                //        if(Communicator::world->me() == 0)
                //          fprintf(outfile, "\nError in DGESV return in orbital domain construction.\n");
                //        exit(PSI_RETURN_FAILURE);
                throw PsiException("Error in DGESV in LMP2 orbital domain construction", __FILE__, __LINE__);
            }
            fR[i - nfzc] = 1.0;
            for (j = 0, row = 0; j < natom; j++) {
                if (domain[i - nfzc][j]) {
                    for (k = aostart[j]; k <= aostop[j]; k++, row++) {
                        for (l = 0; l < nso; l++) fR[i - nfzc] -= Z[row] * aoovlp[k][l] * C[l][i];
                    }
                }
            }

            // Augment the domain if necessary
            if (fabs(fR[i - nfzc]) > cutoff) {
                domain[i - nfzc][rank[next_atom++]] = 1;
                domain_len[i - nfzc]++;
            }
        } // cutoff check
    } // i

    for (i = 0; i < nocc; i++) {
        domain_len_bp[i] = domain_len[i];
        for (k = 0; k < natom; k++)
            domain_bp[i][k] = domain[i][k];
    }
    /* Print the orbital domains */
    if (Communicator::world->me() == 0) {
        fprintf(outfile, "\n   ****** Boughton-Pulay Occupied Orbital Domains ******\n");
        //    domain_print(nocc, natom, domain_len_bp, domain_bp, fR);
        print_domains(fR);
    }

    /* Allow user input of selected domains */
    if (ip_exist("DOMAINS", 0)) {
        ip_count("DOMAINS", &num_entries, 0);
        for (i = 0; i < num_entries; i++) {
            ip_count("DOMAINS", &entry_len, 1, i);
            ip_data("DOMAINS", "%d", &orbital, 2, i, 0);

            /* Clear out the current domain for this orbital */
            for (j = 0; j < natom; j++) domain[orbital][j] = 0;
            domain_len[orbital] = 0;

            for (j = 1; j < entry_len; j++) {
                errcod = ip_data("DOMAINS", "%d", &atom, 2, i, j);
                domain[orbital][atom] = 1;
                domain_len[orbital]++;
            }
        }
    }

    /* Recheck Completeness */
    for (i = 0; i < nocc; i++) {

        /* Build the orbital's domain starting in order of decreasing charge contribution */
        for (j = 0; j < nso; j++) {
            SR[j] = 0.0;
            for (k = 0; k < nso; k++)
                SR[j] += aoovlp[j][k] * C[k][i];
        }

        for (j = 0, row = 0; j < natom; j++) {
            if (domain[i - nfzc][j]) {
                for (k = aostart[j]; k <= aostop[j]; k++, row++) {

                    Z[row] = SR[k];

                    for (l = 0, col = 0; l < natom; l++) {
                        if (domain[i - nfzc][l]) {

                            for (m = aostart[l]; m <= aostop[l]; m++, col++)
                                X[row][col] = aoovlp[k][m];

                        }
                    } /* l */

                } /* k */
            }
        } /* j */

        errcod = C_DGESV(row, 1, &(X[0][0]), nso, &(ipiv[0]), &(Z[0]), nso);
        if (errcod) {
            //        if(Communicator::world->me() == 0)
            //          fprintf(outfile, "\nError in DGESV return in orbital domain construction.\n");
            //        exit(PSI_RETURN_FAILURE);
            throw PsiException("Error in DGESV in LMP2 orbital domain construction", __FILE__, __LINE__);
        }

        /*    fR[i-nfzc] = 1.0;
            for(j=0,row=0; j < natom; j++) {
              if(domain[i-nfzc][j]) {
                for(k=aostart[j]; k <= aostop[j]; k++,row++) {
                  for(l=0; l < nso; l++) fR[i-nfzc] -= Z[row] * aoovlp[k][l] * C[l][i];
                }
              }
            }
         */
    } /* i */

    /* Print the orbital domains */
    if (Communicator::world->me() == 0) {
        fprintf(outfile, "\n   ****** Final Occupied Orbital Domains ******\n");
        //    if(!domain_sep)
        //    domain_print(nocc, natom, domain_len, domain, fR);
        print_domains(fR);
    }

    //  for(i=0; i < natom; i++) {
    //    fprintf(outfile, "aostart[%d] = %d\n", i, aostart[i]);
    //    fprintf(outfile, "aostop[%d] = %d\n", i, aostop[i]);
    //  }

    int pairs = (nocc * (nocc + 1)) / 2;


    /* Determine if pairs are distant or not */
    pairdom_exist = init_int_array(pairs);

    ij_pairs = 0;
    for (i = 0; i < nocc; i++) {
        for (m = 0; m < natom; m++) {
            for (j = 0; j <= i; j++) {
                ij = (i * (i + 1)) / 2 + j;
                for (n = 0; n < natom; n++) {
                    if (pairdom_exist[ij] == 0) {
                        if(neglectdp) {
                            if (domain[i][m] && domain[j][n]) {
                                if(m == n || i == j) {
                                    pairdom_exist[ij] = 1;
                                    ij_pairs++;
				    fprintf(outfile,"  (%d,%d) = %d, Significant at atom %d, atmo %d\n", i,j,ij,m,n);
                                }
                                else {
                                    xcoord = (geometry[n][0] - geometry[m][0]) * (geometry[n][0] - geometry[m][0]);
                                    ycoord = (geometry[n][1] - geometry[m][1]) * (geometry[n][1] - geometry[m][1]);
                                    zcoord = (geometry[n][2] - geometry[m][2]) * (geometry[n][2] - geometry[m][2]);
                                    dist = sqrt((xcoord + ycoord + zcoord)) * 0.529177249;
                                    if (dist < dpcutoff) {
                                        pairdom_exist[ij] = 1;
                                        ij_pairs++;
                                    }
                                }
                            }
                        }
                        else {
                            pairdom_exist[ij] = 1;
                            ij_pairs++;
                        }
                    }
                }
            }
        }
    }

    int **ij_map = get_ij_map();
    if(Communicator::world->me() == 0) {
      // if(neglectdp && print > 3) {
      if(neglectdp) {
        fprintf(outfile, "\n   The following ij pairs are significant and will not be negleted\n");
        for(ij=0; ij < ij_pairs; ij++) {
          //if (pairdom_exist[ij]==0) {
            fprintf(outfile, " i = %d \t j = %d\n",  ij_map[ij][0], ij_map[ij][1]);
          //}
        }
      }
    }

    free_int_matrix(ij_map);
    free_int_matrix(domain_bp);
    free(domain_len_bp);

    free(ao2atom);
    free(l_length);
    free(charge);
    free(rank);
    free(boolean);
    free(SR);
    free(Z);
    free(ipiv);
    free(fR);

    print_test = 0;
    ip_boolean("DOMAIN_PRINT", &(print_test), 0);
    if (print_test) {
        //    if(Communicator::world->me() == 0);
        //    fprintf(outfile, "Printing of orbital domains requested...exiting.\n\n");
        //    exit(PSI_RETURN_FAILURE);
        throw PsiException("Printing of orbital domains requested...exiting", __FILE__, __LINE__);
    }

    /************* Orbital Domains Complete ***************/

}

}
} // namespace psi::lmp2
