/*! \file
    \ingroup LMP2
    \brief localized the SCF MO's
 */
#include "mpi.h"
#include <iostream>
#include <fstream>              // file I/O support
//#include <cstdio>
//#include <cstdlib>
//#include <cstring>
//#include <cmath>
//#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
//#include <libchkpt/chkpt.h>
//#include <libpsio/psio.h>
//#include <libqt/qt.h>
//#include <libint/libint.h>
#include <libmints/basisset.h>
//#include <libmints/onebody.h>
#include <libmints/twobody.h>
#include <libmints/integral.h>
#include <libmints/factory.h>
#include <libmints/symmetry.h>
#include <libmints/wavefunction.h>
#include <libparallel/parallel.h>

//#include <psifiles.h>
#define EXTERN
#include "globals.h"

namespace psi {

extern int myid;
extern int nprocs;

namespace lmp2 {

extern int myid_lmp2;
extern int nprocs_lmp2;

void LMP2::direct_transformation() {

    using namespace std;

    int i, j, a, b, ab, t, ij;
    int k, v;
    int count;
    int *ij_owner, *ij_local, *mr_owner;
    int **MR_shell;
    int M, R, N, S;
    double ****eri_1, ***eri_2, ***eri_2_mn, ***eri_3, ***eri_4;
    // MPI_Status stat;

    //  ****  These are required to utilize libmints  ****

    // Required for libmints, allocates and computes the following:
    // ioff, fac, df, bc
    Wavefunction::initialize_singletons();

    // Create a basis set object and initialize it using the checkpoint file.
    shared_ptr<BasisSet> basis(new BasisSet(chkpt));

    // Initialize an integral factory
    shared_ptr<IntegralFactory> integral(new IntegralFactory(basis, basis, basis, basis));

    // Create an integral object for ERIs
    shared_ptr<TwoBodyInt> eri(integral->eri());

    // Get the storage buffer from the eri object
    const double *buffer = eri->buffer();

    // **** Integral direct transformation ****

    //  **** Allocate memory for the first quarter integral transformation ****
    //  double *B1 = (double *) calloc(nso*nocc*nso*nocc , sizeof(double));

    ij_owner = get_ij_owner();
    ij_local = get_ij_local();

    if (ij_pairs % nprocs == 0) {
        eri_2 = (double ***) malloc((ij_pairs / nprocs) * sizeof (double **));
        for (ij = 0; ij < ij_pairs; ij++) {
            if (myid == ij_owner[ij])
                eri_2[ij_local[ij]] = block_matrix(nso, nso);
        }
    } else {
        if (myid < ij_pairs % nprocs) {
            eri_2 = (double ***) malloc(((ij_pairs / nprocs) + 1) * sizeof (double **));
            for (ij = 0; ij < ij_pairs; ij++) {
                if (myid == ij_owner[ij])
                    eri_2[ij_local[ij]] = block_matrix(nso, nso);
            }
        } else {
            eri_2 = (double ***) malloc(ij_pairs / nprocs * sizeof (double **));
            for (ij = 0; ij < ij_pairs; ij++) {
                if (myid == ij_owner[ij])
                    eri_2[ij_local[ij]] = block_matrix(nso, nso);
            }
        }
    }


    num_unique_shells = 0;
    for (M = 0; M < nshell; M++) {
        for (R = 0; R <= M; R++) {
            num_unique_shells++;
        }
    }

  MR_shell = init_int_matrix(4, num_unique_shells);

    count = 0;
    for (M = 0; M < nshell; M++) {
        int numm = basis->shell(M)->nfunction();
        for (R = 0; R <= M; R++, count++) {
            int numr = basis->shell(R)->nfunction();
            MR_shell[0][count] = M;
            MR_shell[1][count] = numm;
            MR_shell[2][count] = R;
            MR_shell[3][count] = numr;
        }
    }

    sort_shell(MR_shell, num_unique_shells);

    int maxshell = MR_shell[1][0];

    eri_1 = (double ****) calloc(maxshell, sizeof (double ***));
    for (i = 0; i < maxshell; i++) {
        eri_1[i] = (double ***) calloc(maxshell, sizeof (double **));
        for (j = 0; j < maxshell; j++) {
            eri_1[i][j] = block_matrix(nso, nocc);
        }
    }

    eri_2_mn = (double ***) calloc(ij_pairs, sizeof (double **));
    for (ij = 0; ij < ij_pairs; ij++) {
        eri_2_mn[ij] = block_matrix(maxshell*2, nso);
    }

    double **rec;
    rec = block_matrix(maxshell*2, nso);

    mr_owner = get_mr_owner(num_unique_shells);

    int *map_m = (int *) calloc(maxshell * 2, sizeof (int));
    int *map_buf = (int *)  calloc(maxshell*2 , sizeof(int));

    //  **** Integral direct half AO-MO integral transformation
    v = 0;
    for (count = 0; count < num_unique_shells; count++) {
        M = MR_shell[0][count];
        int numm = MR_shell[1][count];
        R = MR_shell[2][count];
        int numr = MR_shell[3][count];

        if (myid == mr_owner[count]) {

            for (i = 0; i < maxshell; i++) {
                for (j = 0; j < maxshell; j++) {
                    zero_mat(eri_1[i][j], nso, nocc);
                }
            }

            for (N = 0; N < nshell; N++) {
                int numn = basis->shell(N)->nfunction();
                for (S = 0; S <= N; S++) {
                    int nums = basis->shell(S)->nfunction();

                    eri->compute_shell(M, R, N, S);

                    for (j = 0; j < nocc; j++) {
                        //  **** First quarter integral transformation ****
                        int index = 0;
                        for (int m = 0; m < numm; m++) {
                            int om = basis->shell(M)->function_index() + m;
                            for (int r = 0; r < numr; r++) {
                                int oor = basis->shell(R)->function_index() + r;
                                for (int n = 0; n < numn; n++) {
                                    int oon = basis->shell(N)->function_index() + n;
                                    for (int s = 0; s < nums; s++, index++) {
                                        int os = basis->shell(S)->function_index() + s;

                                        if (fabs(buffer[index]) > 1.0e-14) {
                                            if (N > S) {
                                                eri_1[m][r][oon][j] += C[os][j] * buffer[index];
                                                eri_1[m][r][os][j] += C[oon][j] * buffer[index];
                                            } else {
                                                eri_1[m][r][oon][j] += C[os][j] * buffer[index];
                                            }
                                        }
                                    } // End of s loop
                                } // End of n loop
                            } // End of r loop
                        } // End of m loop
                    } // End of j loop
                } // End of S loop
            } // Enf of N loop

            for (ij = 0; ij < ij_pairs; ij++) {
                zero_mat(eri_2_mn[ij], maxshell * 2, nso);
            }
            for (i = 0; i < maxshell*2; i++) {
                map_m[i] = 0;
            }

            //  **** Second quarter integral transformation ****
            for (i = 0, ij = 0; i < nocc; i++) {
                for (j = 0; j <= i; j++, ij++) {
                    for (int m = 0; m < numm; m++) {
                        int om = basis->shell(M)->function_index() + m;
                        for (int r = 0; r < numr; r++) {
                            int oor = basis->shell(R)->function_index() + r;
                            for (N = 0; N < nshell; N++) {
                                int numn = basis->shell(N)->nfunction();
                                for (int n = 0; n < numn; n++) {
                                    int oon = basis->shell(N)->function_index() + n;

                                    // This is a map for mapping m and r to om,
                                    // which is needed after the communication
                                    map_m[m] = om;
                                    map_m[r + numm] = oor;

                                    if (fabs(eri_1[m][r][oon][j]) > 1.0e-14) {

                                        if ((myid == ij_owner[ij]) && (myid == mr_owner[count])) {
                                            if (om > oor) {
                                                eri_2[ij_local[ij]][om][oon] += C[oor][i] * eri_1[m][r][oon][j];
                                                eri_2[ij_local[ij]][oor][oon] += C[om][i] * eri_1[m][r][oon][j];
                                            } else if (om == oor) {
                                                eri_2[ij_local[ij]][om][oon] += C[oor][i] * eri_1[m][r][oon][j];
                                            }
                                        } else {
                                            if (om > oor) {
                                                eri_2_mn[ij][m][oon] += C[oor][i] * eri_1[m][r][oon][j];
                                                eri_2_mn[ij][r + numm][oon] += C[om][i] * eri_1[m][r][oon][j];
                                            } else if (om == oor) {
                                                eri_2_mn[ij][m][oon] += C[oor][i] * eri_1[m][r][oon][j];
                                            }
                                        }
                                    }
                                } // End of n loop
                            } // End of N loop
                        } // End of r loop
                    } // End of m loop
                } // End of j loop
            } // End of i loop
        } // End if myid loop

        for (ij = 0; ij < ij_pairs; ij++) {
            if ((myid == ij_owner[ij]) && (myid != mr_owner[count])) {
                zero_mat(rec, maxshell * 2, nso);
                for (i = 0; i < maxshell * 2; i++) {
                    map_buf[i] = 0;
                }
                Communicator::world->recv(mr_owner[count], rec[0],  (numm+numr) * nso);
                Communicator::world->recv(mr_owner[count], map_buf,  (numm+numr));
                for (int m = 0; m < (numm + numr); m++) {
                    for (int N = 0; N < nshell; N++) {
                        int numn = basis->shell(N)->nfunction();
                        for (int n = 0; n < numn; n++) {
                            int oon = basis->shell(N)->function_index() + n;

                            if (fabs(rec[m][oon]) > 1.0e-14)
                                eri_2[ij_local[ij]][map_buf[m]][oon] += rec[m][oon];
                        }
                    }
                }
                //}
            } else if ((myid != ij_owner[ij]) && (myid == mr_owner[count])) {
                Communicator::world->send(ij_owner[ij], eri_2_mn[ij][0], (numm+numr) * nso);
                Communicator::world->send(ij_owner[ij], map_m, (numm+numr));
            }
        }

        v++;
    } // End of M loop

    //  **** Free the memory used for the first quarter integral transformation ****
    for (i = 0; i < maxshell; i++) {
        for (j = 0; j < maxshell; j++) {
            free_block(eri_1[i][j]);
        }
        free(eri_1[i]);
    }
    free(eri_1);

    // Freeing the memory used by eri_2_mn, the map_m, and the buffers used in the communication
    for(ij=0; ij < ij_pairs; ij++)
        free_block(eri_2_mn[ij]);
    free(eri_2_mn);
    free(map_m);
    free(map_buf);
    free_block(rec);

    // Allocating the memory need by eri_3
    if (ij_pairs % nprocs == 0) {
        eri_3 = (double ***) malloc((ij_pairs / nprocs) * sizeof (double **));
        for (ij = 0; ij < ij_pairs; ij++) {
            if (myid == ij_owner[ij])
                eri_3[ij_local[ij]] = block_matrix(nso, pairdom_len[ij]);
        }
    }
    else {
        if (myid < ij_pairs % nprocs) {
            eri_3 = (double ***) malloc(((ij_pairs / nprocs) + 1) * sizeof (double **));
            for (ij = 0; ij < ij_pairs; ij++) {
                if (myid == ij_owner[ij])
                    eri_3[ij_local[ij]] = block_matrix(nso, pairdom_len[ij]);
            }
        }
        else {
            eri_3 = (double ***) malloc(ij_pairs / nprocs * sizeof (double **));
            for (ij = 0; ij < ij_pairs; ij++) {
                if (myid == ij_owner[ij])
                    eri_3[ij_local[ij]] = block_matrix(nso, pairdom_len[ij]);
            }
        }
    }

    //  **** Third quarter integral transformation ****
    v = 0;
    for (i = 0, ij = 0; i < nocc; i++) {
        for (j = 0; j <= i; j++, ij++) {
            if (v % nprocs == myid) {
                for (k = 0, b = 0; k < natom; k++) {
                    if (pairdomain[ij][k]) {
                        for (t = aostart[k]; t <= aostop[k]; t++, b++) {
                            for (M = 0; M < nshell; M++) {
                                int numm = basis->shell(M)->nfunction();
                                for (int m = 0; m < numm; m++) {
                                    int om = basis->shell(M)->function_index() + m;
                                    for (N = 0; N < nshell; N++) {
                                        int numn = basis->shell(N)->nfunction();
                                        for (int n = 0; n < numn; n++) {
                                            int oon = basis->shell(N)->function_index() + n;

                                            if (fabs(eri_2[ij_local[ij]][om][oon]) > 1.0e-14) {

                                                eri_3[ij_local[ij]][om][b] += Rt_full[oon][t] * eri_2[ij_local[ij]][om][oon];

                                            }
                                        } // End of n loop
                                    } // End of N loop
                                } // End of m loop
                            } // End of M loop
                        } // End of b loop
                    } // End of if pairdomain
                } // End of k loop
            } // End of if myid loop
            v++;
        } // End of j loop
    } // End of i loop

    // Freeing the memory used by eri_2
    if (ij_pairs % nprocs == 0) {
        for (ij = 0; ij < ij_pairs/nprocs; ij++)
            free_block(eri_2[ij]);
    }
    else {
        if (myid < ij_pairs % nprocs) {
            for (ij = 0; ij < (ij_pairs/nprocs) + 1; ij++)
                free_block(eri_2[ij]);
        }
        else {
            for (ij = 0; ij < ij_pairs/nprocs; ij++)
                free_block(eri_2[ij]);
        }
    }
    free(eri_2);

    // Allocating the memory needed by Ktilde
    if (ij_pairs % nprocs == 0) {
        Ktilde = (double ***) malloc((ij_pairs / nprocs) * sizeof (double **));
        for (ij = 0; ij < ij_pairs; ij++) {
            if (myid == ij_owner[ij])
                Ktilde[ij_local[ij]] = block_matrix(pairdom_len[ij], pairdom_len[ij]);
        }
    }
    else {
        if (myid < ij_pairs % nprocs) {
            Ktilde = (double ***) malloc(((ij_pairs / nprocs) + 1) * sizeof (double **));
            for (ij = 0; ij < ij_pairs; ij++) {
                if (myid == ij_owner[ij])
                    Ktilde[ij_local[ij]] = block_matrix(pairdom_len[ij], pairdom_len[ij]);
            }
        }
        else {
            Ktilde = (double ***) malloc(ij_pairs / nprocs * sizeof (double **));
            for (ij = 0; ij < ij_pairs; ij++) {
                if (myid == ij_owner[ij])
                    Ktilde[ij_local[ij]] = block_matrix(pairdom_len[ij], pairdom_len[ij]);
            }
        }
    }

    //  **** Fourth quarter integral transformation ****
    v = 0;
    for (i = 0, ij = 0; i < nocc; i++) {
        for (j = 0; j <= i; j++, ij++) {
            if (v % nprocs == myid) {
                for (k = 0, a = 0; k < natom; k++) {
                    if (pairdomain[ij][k]) {
                        for (t = aostart[k]; t <= aostop[k]; t++, a++) {
                            for (b = 0; b < pairdom_len[ij]; b++) {
                                //                Ktilde[ij_local[ij]][a][b] = 0;
                                for (M = 0; M < nshell; M++) {
                                    int numm = basis->shell(M)->nfunction();
                                    for (int m = 0; m < numm; m++) {
                                        int om = basis->shell(M)->function_index() + m;

                                        if (fabs(eri_3[ij_local[ij]][om][b]) > 1.0e-14) {
                                            //                      Ktilde[ij][a][b] += Rt_full[om][t] * eri_3[ij][om][b];
                                            Ktilde[ij_local[ij]][a][b] += Rt_full[om][t] * eri_3[ij_local[ij]][om][b];
                                        }

                                    } // End of m loop
                                } // End of M loop
                            } // End of b loop
                        } // End of a loop
                    } // End of if pairdomain
                } // End of k loop
            } // End of if myid loop
            v++;
        } // End of j loop
    } // End of i loop

    if (print > 2) {
        v = 0;
        for (i = 0, ij = 0; i < nocc; i++) {
            for (j = 0; j <= i; j++, ij++) {
                if (v % nprocs == myid) {
                    fprintf(outfile, "Ktilde[%d] matrix", ij);
                    print_mat(Ktilde[ij_local[ij]], pairdom_len[ij], pairdom_len[ij], outfile);
                }
                v++;
            }
        }
    }

    if (ij_pairs % nprocs == 0) {
        for (ij = 0; ij < ij_pairs/nprocs; ij++)
            free_block(eri_3[ij]);
    }
    else {
        if (myid < ij_pairs % nprocs) {
            for (ij = 0; ij < (ij_pairs/nprocs) + 1; ij++)
                free_block(eri_3[ij]);
        }
        else {
            for (ij = 0; ij < ij_pairs/nprocs; ij++)
                free_block(eri_3[ij]);
        }
    }
    free(eri_3);
    free(ij_owner);
    free(mr_owner);
    free(ij_local);
    free(MR_shell);

}

}} // namespace psi::lmp2
