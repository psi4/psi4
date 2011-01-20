/*! \file
    \ingroup LMP2
    \brief localized the SCF MO's
 */
//#include "mpi.h"
#include <iostream>
#include <fstream>              // file I/O support
#include <libciomr/libciomr.h>
#include <libmints/mints.h>
#include <libparallel/parallel.h>

#define EXTERN
#include "globals.h"

namespace psi {

namespace lmp2 {

void LMP2::direct_transformation() {

    using namespace std;

    int i, j, a, b, ab, t, ij;
    int k;
    int count;
    int *ij_owner, *ij_local, *mn_owner;
    int **MN_shell;
    int M, R, N, S;
    int *R_shells, *S_shells;
    int **ij_map;
    int **pairdomain, *pairdom_len;
    double temp_val;
    double **T_NS, *T_N;
    double **P_Ns_max, *P_N_max, **PP_max;
    double **L_Nj_max, *L_N_max, Lmax;
    double **LL_MN_max, *LL_M_max, LLmax;
    double ****eri_1, ***eri_2, ***eri_2_mn, ***eri_3, ***eri_4;
    // MPI_Status stat;

    //  ****  These are required to utilize libmints  ****

    // Create a basis set object and initialize it.
    shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    shared_ptr<BasisSet> basis = BasisSet::construct(parser, Process::environment.molecule(), orbital_basis);

    // Initialize an integral factory
    shared_ptr<IntegralFactory> integral(new IntegralFactory(basis, basis, basis, basis));

    // Create an integral object for ERIs
    shared_ptr<TwoBodyAOInt> eri(integral->eri());

    // Get the storage buffer from the eri object
    const double *buffer = eri->buffer();

    // **** Integral direct transformation ****

    //  **** Allocate memory for the first quarter integral transformation ****
    //  double *B1 = (double *) calloc(nso*nocc*nso*nocc , sizeof(double));

    ij_owner = get_ij_owner();
    ij_local = get_ij_local();
    ij_map = get_ij_map();
    pairdomain = compute_pairdomain(ij_map);
    pairdom_len = compute_pairdomlen(ij_map);



    eri_2 = (double ***) malloc(pairs_per_proc * sizeof (double **));
    for(ij=0; ij < ij_pairs; ij++) {
        if(Communicator::world->me() == ij_owner[ij])
            eri_2[ij_local[ij]] = block_matrix(nso, nso);
    }

/*    if (ij_pairs % nprocs == 0) {
        eri_2 = (double ***) malloc((ij_pairs / nprocs) * sizeof (double **));
        for (ij = 0; ij < ij_pairs; ij++) {
            if (Communicator::world->me() == ij_owner[ij])
                eri_2[ij_local[ij]] = block_matrix(nso, nso);
        }
    } else {
        if (Communicator::world->me() < ij_pairs % nprocs) {
            eri_2 = (double ***) malloc(((ij_pairs / nprocs) + 1) * sizeof (double **));
            for (ij = 0; ij < ij_pairs; ij++) {
                if (Communicator::world->me() == ij_owner[ij])
                    eri_2[ij_local[ij]] = block_matrix(nso, nso);
            }
        } else {
            eri_2 = (double ***) malloc(ij_pairs / nprocs * sizeof (double **));
            for (ij = 0; ij < ij_pairs; ij++) {
                if (Communicator::world->me() == ij_owner[ij])
                    eri_2[ij_local[ij]] = block_matrix(nso, nso);
            }
        }
    }
*/

    num_unique_shells = 0;
    for (M = 0; M < nshell; M++) {
        for (N = 0; N < nshell; N++) {
            num_unique_shells++;
        }
    }

    MN_shell = init_int_matrix(4, num_unique_shells);

    count = 0;
    for (M = 0; M < nshell; M++) {
        int numm = basis->shell(M)->nfunction();
        for (N = 0; N < nshell; N++, count++) {
            int numn = basis->shell(N)->nfunction();
            MN_shell[0][count] = M;
            MN_shell[1][count] = numm;
            MN_shell[2][count] = N;
            MN_shell[3][count] = numn;
        }
    }

    /* Computing the required integral prescreening quantities */
    if(screen_int) {
        T_NS = block_matrix(num_unique_shells, num_unique_shells);
        T_N = init_array(nshell);
        for (M = 0; M < nshell; M++) {
            double maxval = 0.0;
            int numm = basis->shell(M)->nfunction();
            for (N = 0; N < nshell; N++) {
                int numn = basis->shell(N)->nfunction();
                R = M;
                S = N;
                int numr = numm;
                int nums = numn;

                eri->compute_shell(M, N, R, S);

                int index = 0;
                for (int m = 0; m < numm; m++) {
                    int om = basis->shell(M)->function_index() + m;
                    for (int r = 0; r < numr; r++) {
                        int oor = basis->shell(R)->function_index() + r;
                        for (int n = 0; n < numn; n++) {
                            int oon = basis->shell(N)->function_index() + n;
                            for (int s = 0; s < nums; s++, index++) {
                                int os = basis->shell(S)->function_index() + s;

                                if(om == oor && oon == os) {

                                    temp_val = sqrt(fabs(buffer[index]));

                                    if(temp_val > T_NS[M][N])
                                        T_NS[M][N] = temp_val;
                                    if(temp_val > maxval)
                                        maxval = temp_val;
                                }
                            }
                        }
                    }
                }
            }
            T_N[M] = maxval;
        }

        //for (M = 0; M < nshell; M++) {
        //    for (N = 0; N < nshell; N++) {
        //        fprintf(outfile, "T_NS[%d][%d] = %20.12f\n", M, N, T_NS[M][N]);
        //        fflush(outfile);
        //    }
        //}

        P_Ns_max = block_matrix(nshell, nso);
        for(int M = 0; M < nshell; M++) {
            double maxval = 0.0;
            int numm = basis->shell(M)->nfunction();
            for (int m = 0; m < numm; m++) {
                double maxval2 = 0.0;
                int om = basis->shell(M)->function_index() + m;
                for(int s=0; s < nso; s++) {

                        temp_val = fabs(Rt_full[om][s]);

                        if(temp_val > P_Ns_max[M][s])
                            P_Ns_max[M][s] = temp_val;
                }
            }
        }

        P_N_max = init_array(nshell);
        for(int N = 0; N < nshell; N++) {
            for(int s=0; s < nso; s++) {
                if(P_Ns_max[N][s] > P_N_max[N])
                    P_N_max[N] = P_Ns_max[N][s];
            }
        }

        PP_max = block_matrix(nshell, nshell);
        for (int M = 0; M < nshell; M++) {
            for (int N = 0; N < nshell; N++) {
                for (int ij = 0; ij < ij_pairs; ij++) {
                    for (int k = 0, a = 0; k < natom; k++) {
                        if (pairdomain[ij][k]) {
                            for (int t = aostart[k]; t <= aostop[k]; t++, a++) {
                                for (int l = 0, b = 0; l < natom; l++) {
                                    if (pairdomain[ij][l]) {
                                        for (int u = aostart[l]; u <= aostop[l]; u++, b++) {
                                            temp_val = P_Ns_max[M][a] * P_Ns_max[N][b];
                                            if(temp_val > PP_max[M][N])
                                                PP_max[M][N] = temp_val;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        L_Nj_max = block_matrix(nshell, nocc);
        for(int N = 0; N < nshell; N++) {
            int numn = basis->shell(N)->nfunction();
            for (int n = 0; n < numn; n++) {
                int oon = basis->shell(N)->function_index() + n;
                for(int j=0; j < nocc; j++) {

                        temp_val = fabs(C[oon][j]);

                        if(temp_val > L_Nj_max[N][j])
                            L_Nj_max[N][j] = temp_val;

                }
            }
        }

        L_N_max = init_array(nshell);
        for(int N = 0; N < nshell; N++) {
            for(int j=0; j < nocc; j++) {
                if(L_Nj_max[N][j] > L_N_max[N])
                    L_N_max[N] = L_Nj_max[N][j];
            }
        }

        for(int N = 0; N < nshell; N++) {
            if(L_N_max[N] > Lmax)
                Lmax = L_N_max[N];
        }

        LL_MN_max = block_matrix(nshell, nshell);
        for (int M = 0; M < nshell; M++) {
            for (int N = 0; N < nshell; N++) {
                for (int ij = 0; ij < ij_pairs; ij++) {
                    i = ij_map[ij][0];
                    j = ij_map[ij][1];

                    temp_val = L_Nj_max[M][i] * L_Nj_max[N][j];
                    if(temp_val > LL_MN_max[M][N])
                        LL_MN_max[M][N] = temp_val;

                }
            }
        }

        LL_M_max = init_array(nshell);
        for (int M = 0; M < nshell; M++) {
            for (int N = 0; N < nshell; N++) {
                if(LL_MN_max[M][N] > LL_M_max[M])
                    LL_M_max[M] = LL_MN_max[M][N];

            }
        }

        for (int M = 0; M < nshell; M++) {
            if(LL_M_max[M] > LLmax)
                LLmax = LL_M_max[M];
        }
    }
    /* End of computing prescreening quantities */

    sort_shell(MN_shell, num_unique_shells);

    int maxshell = MN_shell[1][0];

    eri_1 = (double ****) malloc(maxshell * sizeof (double ***));
    for (i = 0; i < maxshell; i++) {
        eri_1[i] = (double ***) malloc(nso * sizeof (double **));
        for (j = 0; j < nso; j++) {
            eri_1[i][j] = block_matrix(maxshell, nocc);
        }
    }

    eri_2_mn = (double ***) malloc(ij_pairs * sizeof (double **));
    for (ij = 0; ij < ij_pairs; ij++) {
        eri_2_mn[ij] = block_matrix(nso, nso);
    }

    double **rec;
    rec = block_matrix(nso, nso);

    mn_owner = get_mn_owner(num_unique_shells);

    //  **** Integral direct half AO-MO integral transformation
    for (count = 0; count < num_unique_shells; count++) {
        M = MN_shell[0][count];
        int numm = MN_shell[1][count];
        N = MN_shell[2][count];
        int numn = MN_shell[3][count];

        if (Communicator::world->me() == mn_owner[count]) {

            if(screen_int)
                if(T_N[M] * T_N[N] * PP_max[M][N] * LLmax < escreen) continue;

            for (i = 0; i < maxshell; i++) {
                for (j = 0; j < nso; j++) {
                    zero_mat(eri_1[i][j], maxshell, nocc);
                }
            }

            if(screen_int) R_shells = sort_schwarz_MN(T_NS, M, nshell);
            for (int R_loop = 0; R_loop < nshell; R_loop++) {
                if(screen_int) {
                    R = R_shells[R_loop];

                    if(T_N[N] * T_NS[M][R] * PP_max[M][N] * LLmax < escreen) {
                        continue;
                    }
                    if(T_N[N] * T_NS[M][R] * PP_max[M][N] * LL_M_max[R] < escreen) {
                        continue;
                    }
                }
                else R = R_loop;

                int numr = basis->shell(R)->nfunction();

                if(screen_int) S_shells = sort_schwarz_MN(T_NS, N, nshell);
                for (int S_loop = 0; S_loop < nshell; S_loop++) {
                    if(screen_int) {
                        S = S_shells[S_loop];

                        if(T_N[N] * T_NS[M][R] * PP_max[M][N] * LL_M_max[R] < escreen) continue;
                        if(T_N[N] * T_NS[M][R] * PP_max[M][N] * LL_MN_max[R][S] < escreen) continue;
                    }
                    else S = S_loop;

                    int nums = basis->shell(S)->nfunction();

                    eri->compute_shell(M, R, N, S);

                    for (j = 0; j < nocc; j++) {
                        //  **** First quarter integral transformation ****
                        int index = 0;
                        for (int m = 0; m < numm; m++) {
                            //int om = basis->shell(M)->function_index() + m;
                            for (int r = 0; r < numr; r++) {
                                int oor = basis->shell(R)->function_index() + r;
                                for (int n = 0; n < numn; n++) {
                                    //int oon = basis->shell(N)->function_index() + n;
                                    for (int s = 0; s < nums; s++, index++) {
                                        int os = basis->shell(S)->function_index() + s;

                                        if (fabs(buffer[index]) > 1.0e-14) {

                                            eri_1[m][oor][n][j] += C[os][j] * buffer[index];

                                        }
                                    } // End of s loop
                                } // End of n loop
                            } // End of r loop
                        } // End of m loop
                    } // End of j loop
                } // End of S loop
            } // Enf of N loop


            //  **** Second quarter integral transformation ****
            for (ij = 0; ij < ij_pairs; ij++) {
                i = ij_map[ij][0];
                j = ij_map[ij][1];
                for (int m = 0; m < numm; m++) {
                    int om = basis->shell(M)->function_index() + m;
                    for (int n = 0; n < numn; n++) {
                        int oon = basis->shell(N)->function_index() + n;
                        for (R = 0; R < nshell; R++) {
                            int numr = basis->shell(R)->nfunction();
                            for (int r = 0; r < numr; r++) {
                                int oor = basis->shell(R)->function_index() + r;

                                if (fabs(eri_1[m][oor][n][j]) > 1.0e-14) {

                                    if ((Communicator::world->me() == ij_owner[ij]) && (Communicator::world->me() == mn_owner[count])) {
                                        eri_2[ij_local[ij]][om][oon] += C[oor][i] * eri_1[m][oor][n][j];
                                    }
                                    else {
                                        eri_2_mn[ij][om][oon] += C[oor][i] * eri_1[m][oor][n][j];
                                    }
                                }
                            } // End of n loop
                        } // End of N loop
                    } // End of r loop
                } // End of m loop
            } // End of ij loop
        } // End if Communicator::world->me() loop
    } // End of MN loop

    //  **** Free the memory used for the first quarter integral transformation ****
    for (i = 0; i < maxshell; i++) {
        for (j = 0; j < nso; j++) {
            free_block(eri_1[i][j]);
        }
        free(eri_1[i]);
    }
    free(eri_1);

    /* Redistribute the MN integrals over ij pairs */
    for (count = 0; count < num_unique_shells; count++) {
        M = MN_shell[0][count];
        int numm = MN_shell[1][count];
        N = MN_shell[2][count];
        int numn = MN_shell[3][count];

        for (ij = 0; ij < ij_pairs; ij++) {
            if ((Communicator::world->me() == ij_owner[ij]) && (Communicator::world->me() != mn_owner[count])) {
                zero_mat(rec, nso, nso);
                Communicator::world->recv(mn_owner[count], rec[0], nso * nso);
                for (int m = 0; m < numm; m++) {
                    int om = basis->shell(M)->function_index() + m;
                    for (int n = 0; n < numn; n++) {
                        int oon = basis->shell(N)->function_index() + n;

                        if (fabs(rec[om][oon]) > 1.0e-14)
                            eri_2[ij_local[ij]][om][oon] = rec[om][oon];

                    }
                }
            }
            else if ((Communicator::world->me() != ij_owner[ij]) && (Communicator::world->me() == mn_owner[count])) {
                Communicator::world->send(ij_owner[ij], eri_2_mn[ij][0], nso * nso);
            }
        }
    }

    // Freeing the memory used by eri_2_mn, the map_m, and the buffers used in the communication
    for (ij = 0; ij < ij_pairs; ij++)
        free_block(eri_2_mn[ij]);
    free(eri_2_mn);
    free_block(rec);

    // Allocating the memory need by eri_3
    eri_3 = (double ***) malloc(pairs_per_proc * sizeof (double **));
    for (ij = 0; ij < ij_pairs; ij++) {
        if (Communicator::world->me() == ij_owner[ij])
            eri_3[ij_local[ij]] = block_matrix(nso, pairdom_len[ij]);
    }

/*    if (ij_pairs % nprocs == 0) {
        eri_3 = (double ***) malloc((ij_pairs / nprocs) * sizeof (double **));
        for (ij = 0; ij < ij_pairs; ij++) {
            if (Communicator::world->me() == ij_owner[ij])
                eri_3[ij_local[ij]] = block_matrix(nso, pairdom_len[ij]);
        }
    } else {
        if (Communicator::world->me() < ij_pairs % nprocs) {
            eri_3 = (double ***) malloc(((ij_pairs / nprocs) + 1) * sizeof (double **));
            for (ij = 0; ij < ij_pairs; ij++) {
                if (Communicator::world->me() == ij_owner[ij])
                    eri_3[ij_local[ij]] = block_matrix(nso, pairdom_len[ij]);
            }
        } else {
            eri_3 = (double ***) malloc(ij_pairs / nprocs * sizeof (double **));
            for (ij = 0; ij < ij_pairs; ij++) {
                if (Communicator::world->me() == ij_owner[ij])
                    eri_3[ij_local[ij]] = block_matrix(nso, pairdom_len[ij]);
            }
        }
    }
*/
    
    //  **** Third quarter integral transformation ****
    for (ij = 0; ij < ij_pairs; ij++) {
        i = ij_map[ij][0];
        j = ij_map[ij][1];
        if (Communicator::world->me() == ij_owner[ij]) {
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
        } // End of if Communicator::world->me() loop
    } // End of i loop

    // Freeing the memory used by eri_2
    for(ij=0; ij < pairs_per_proc; ij++) {
        free_block(eri_2[ij]);
    }
/*    if (ij_pairs % nprocs == 0) {
        for (ij = 0; ij < ij_pairs / nprocs; ij++)
            free_block(eri_2[ij]);
    } else {
        if (Communicator::world->me() < ij_pairs % nprocs) {
            for (ij = 0; ij < (ij_pairs / nprocs) + 1; ij++)
                free_block(eri_2[ij]);
        } else {
            for (ij = 0; ij < ij_pairs / nprocs; ij++)
                free_block(eri_2[ij]);
        }
    }
 */
    free(eri_2);

    // Allocating the memory needed by Ktilde
    Ktilde = (double ***) malloc(pairs_per_proc * sizeof (double **));
    for (ij = 0; ij < ij_pairs; ij++) {
        if (Communicator::world->me() == ij_owner[ij])
            Ktilde[ij_local[ij]] = block_matrix(pairdom_len[ij], pairdom_len[ij]);
    }

/*    if (ij_pairs % nprocs == 0) {
        Ktilde = (double ***) malloc((ij_pairs / nprocs) * sizeof (double **));
        for (ij = 0; ij < ij_pairs; ij++) {
            if (Communicator::world->me() == ij_owner[ij])
                Ktilde[ij_local[ij]] = block_matrix(pairdom_len[ij], pairdom_len[ij]);
        }
    } else {
        if (Communicator::world->me() < ij_pairs % nprocs) {
            Ktilde = (double ***) malloc(((ij_pairs / nprocs) + 1) * sizeof (double **));
            for (ij = 0; ij < ij_pairs; ij++) {
                if (Communicator::world->me() == ij_owner[ij])
                    Ktilde[ij_local[ij]] = block_matrix(pairdom_len[ij], pairdom_len[ij]);
            }
        } else {
            Ktilde = (double ***) malloc(ij_pairs / nprocs * sizeof (double **));
            for (ij = 0; ij < ij_pairs; ij++) {
                if (Communicator::world->me() == ij_owner[ij])
                    Ktilde[ij_local[ij]] = block_matrix(pairdom_len[ij], pairdom_len[ij]);
            }
        }
    }
 */

    //  **** Fourth quarter integral transformation ****
    for (ij = 0; ij < ij_pairs; ij++) {
        i = ij_map[ij][0];
        j = ij_map[ij][1];
        if (Communicator::world->me() == ij_owner[ij]) {
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
        } // End of if Communicator::world->me() loop
    } // End of i loop

    if (print > 2) {
        for (ij = 0; ij < ij_pairs; ij++) {
            i = ij_map[ij][0];
            j = ij_map[ij][1];
            if (Communicator::world->me() == ij_owner[ij]) {
                fprintf(outfile, "Ktilde[%d] matrix", ij);
                print_mat(Ktilde[ij_local[ij]], pairdom_len[ij], pairdom_len[ij], outfile);
            }
        }
    }

    for(ij=0; ij < pairs_per_proc; ij++) {
        free_block(eri_3[ij]);
    }

/*    if (ij_pairs % nprocs == 0) {
        for (ij = 0; ij < ij_pairs / nprocs; ij++)
            free_block(eri_3[ij]);
    } else {
        if (Communicator::world->me() < ij_pairs % nprocs) {
            for (ij = 0; ij < (ij_pairs / nprocs) + 1; ij++)
                free_block(eri_3[ij]);
        } else {
            for (ij = 0; ij < ij_pairs / nprocs; ij++)
                free_block(eri_3[ij]);
        }
    }
*/
    free(eri_3);
    free(ij_owner);
    free(mn_owner);
    free(ij_local);
    free(MN_shell);

}

}
} // namespace psi::lmp2
