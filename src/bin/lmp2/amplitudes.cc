/*! \file
    \ingroup LMP2
    \brief compute the mp2 amplitudes
*/

#include "mpi.h"
#include <iostream>
//#include <cstdio>
//#include <cstdlib>
//#include <cstring>
//#include <cmath>
//#include <pthread.h>
//#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
//#include <libchkpt/chkpt.h>
//#include <libpsio/psio.h>
#include <libqt/qt.h>
#include <libparallel/parallel.h>
//#include <psifiles.h>
#define EXTERN
#include "globals.h"
extern pthread_mutex_t compute_mutex;

namespace psi {

extern int myid;
extern int nprocs;

namespace lmp2 {

extern int myid_lmp2;
extern int nprocs_lmp2;

void LMP2::amplitudes() {

    int i, j, kj, ik, kj_, ik_, l, k, m, n, K;
    int a, b, c, d, r, s, ij, L, M;
    int t, N, u, o, O;
    int nvir;
    int *ij_owner, *ij_local;
    int **ij_map, *abs_ij_map;
    int **pairdomain, *pairdom_len;
    double **temp, **F_sum, **temp1, **X, **Sij;
    double **Rtilde;
    double **Rbar, **Tbar;
    double **St, **Ft;
    double ****Tempkj, ****Tempik;
    MPI_Win *T2win;
    MPI_Status stat;

    nvir = nso - nocc;

    ij_owner = get_ij_owner();
    ij_local = get_ij_local();

    ij_map = get_ij_map();
    pairdomain = compute_pairdomain(ij_map);
    pairdom_len = compute_pairdomlen(ij_map);

    /* This maps the absolute ij value to the value corresponding ij value
     * after the removal of distant pairs */
    abs_ij_map = original_ij_map();

    
    // Replicate the amplitudes that each process will need to
    // compute the new amplitudes if iter > 0
    if(iter > 0) {

        Tempkj = (double ****) malloc(pairs_per_proc * sizeof (double ***));
        Tempik = (double ****) malloc(pairs_per_proc * sizeof (double ***));
        for(ij=0; ij < ij_pairs; ij++) {
            Tempkj[ij_local[ij]] = (double ***) malloc(nocc * sizeof (double **));
            Tempik[ij_local[ij]] = (double ***) malloc(nocc * sizeof (double **));
        }

        for (ij = 0; ij < ij_pairs; ij++) {
            i = ij_map[ij][0];
            j = ij_map[ij][1];

            // Get the required amplitudes from their respective proc and put them in a temporary array
            // If ij == ik or kj then do simple DCOPY
            for (k = 0; k < nocc; k++) {
                if (k > j) kj_ = (k * (k + 1)) / 2 + j;
                else kj_ = (j * (j + 1)) / 2 + k;

                if (i > k) ik_ = (i * (i + 1)) / 2 + k;
                else ik_ = (k * (k + 1)) / 2 + i;

                // If we are neglecting distant pairs then kj and ik will need to be mapped to the
                // new values that correspond to the removal of the distant pairs
                kj = abs_ij_map[kj_];
                ik = abs_ij_map[ik_];

                if (pairdom_exist[kj_]) {
                    if ((myid == ij_owner[ij]) && (myid == ij_owner[kj])) {
                        Tempkj[ij_local[ij]][k] = block_matrix(pairdom_len[kj], pairdom_len[kj]);
                        C_DCOPY(pairdom_len[kj] * pairdom_len[kj], &(T[dmat1][ij_local[kj]][0][0]), 1, &(Tempkj[ij_local[ij]][k][0][0]), 1);
                    }
                    else if ((myid == ij_owner[ij]) && (myid != ij_owner[kj])) {
                        Tempkj[ij_local[ij]][k] = block_matrix(pairdom_len[kj], pairdom_len[kj]);
                        Communicator::world->recv(ij_owner[kj], Tempkj[ij_local[ij]][k][0], pairdom_len[kj] * pairdom_len[kj]);
                        //MPI::COMM_WORLD.Recv(&(Tempkj[k][0][0]), pairdom_len[kj]*pairdom_len[kj], MPI::DOUBLE, ij_owner[kj], kj);
                    }
                    else if ((myid != ij_owner[ij]) && (myid == ij_owner[kj])) {
                        Communicator::world->send(ij_owner[ij], T[dmat1][ij_local[kj]][0], pairdom_len[kj] * pairdom_len[kj]);
                        //MPI::COMM_WORLD.Send(&(T[dmat1][ij_local[kj]][0][0]), pairdom_len[kj]*pairdom_len[kj], MPI::DOUBLE, ij_owner[ij], kj);
                    }
                }

                if (pairdom_exist[ik_]) {
                    if ((myid == ij_owner[ij]) && (myid == ij_owner[ik])) {
                        Tempik[ij_local[ij]][k] = block_matrix(pairdom_len[ik], pairdom_len[ik]);
                        C_DCOPY(pairdom_len[ik] * pairdom_len[ik], &(T[dmat1][ij_local[ik]][0][0]), 1, &(Tempik[ij_local[ij]][k][0][0]), 1);
                    }
                    else if ((myid == ij_owner[ij]) && (myid != ij_owner[ik])) {
                        Tempik[ij_local[ij]][k] = block_matrix(pairdom_len[ik], pairdom_len[ik]);
                        Communicator::world->recv(ij_owner[ik], Tempik[ij_local[ij]][k][0], pairdom_len[ik] * pairdom_len[ik]);
                        //MPI::COMM_WORLD.Recv(&(Tempik[k][0][0]), pairdom_len[ik]*pairdom_len[ik], MPI::DOUBLE, ij_owner[ik], ik);
                    }
                    else if ((myid != ij_owner[ij]) && (myid == ij_owner[ik])) {
                        Communicator::world->send(ij_owner[ij], T[dmat1][ij_local[ik]][0], pairdom_len[ik] * pairdom_len[ik]);
                        //MPI::COMM_WORLD.Send(&(T[dmat1][ij_local[ik]][0][0]), pairdom_len[ik]*pairdom_len[ik], MPI::DOUBLE, ij_owner[ij], ik);
                    }
                }
            }
        }
    }

    for (ij = 0; ij < ij_pairs; ij++) {
        i = ij_map[ij][0];
        j = ij_map[ij][1];

        if (myid != ij_owner[ij])
            continue;

        zero_mat(T[div][ij_local[ij]], pairdom_len[ij], pairdom_len[ij]);

        Rtilde = block_matrix(pairdom_len[ij], pairdom_len[ij]);
        temp = block_matrix(pairdom_len[ij], pairdom_len[ij]);

        if (iter > 0) {

            /* Build the virtual space overlap and Fock matrices for this pair */
            St = block_matrix(pairdom_len[ij], pairdom_len[ij]);
            Ft = block_matrix(pairdom_len[ij], pairdom_len[ij]);

            for (r = 0, K = 0; r < natom; r++) {
                if (pairdomain[ij][r]) {
                    for (k = aostart[r]; k <= aostop[r]; k++, K++) {
                        for (s = 0, L = 0; s < natom; s++) {
                            if (pairdomain[ij][s]) {
                                for (l = aostart[s]; l <= aostop[s]; l++, L++) {

                                    St[K][L] = paoovlp[k][l];
                                    Ft[K][L] = paoF[k][l];

                                }
                            }
                        }
                    }
                }
            }

            C_DCOPY(pairdom_len[ij] * pairdom_len[ij], &(Ktilde[ij_local[ij]][0][0]), 1, &(Rtilde[0][0]), 1);


            C_DGEMM('n', 'n', pairdom_len[ij], pairdom_len[ij], pairdom_len[ij], 1, &(Ft[0][0]),
                    pairdom_len[ij], &(T[dmat1][ij_local[ij]][0][0]), pairdom_len[ij], 0, &(temp[0][0]), pairdom_len[ij]);
            C_DGEMM('n', 'n', pairdom_len[ij], pairdom_len[ij], pairdom_len[ij], 1, &(temp[0][0]),
                    pairdom_len[ij], &(St[0][0]), pairdom_len[ij], 1, &(Rtilde[0][0]), pairdom_len[ij]);

            //std::cout << "procid = " << myid << "   made it here in iter = " << iter << "   for ij = " << ij <<
            //"       ij_local[" << ij << "] = " << ij_local[ij] << " pairdom_len[" << ij << "] = " << pairdom_len[ij] << std::endl;


            C_DGEMM('n', 'n', pairdom_len[ij], pairdom_len[ij], pairdom_len[ij], 1, &(St[0][0]),
                    pairdom_len[ij], &(T[dmat1][ij_local[ij]][0][0]), pairdom_len[ij], 0, &(temp[0][0]), pairdom_len[ij]);
            C_DGEMM('n', 'n', pairdom_len[ij], pairdom_len[ij], pairdom_len[ij], 1, &(temp[0][0]),
                    pairdom_len[ij], &(Ft[0][0]), pairdom_len[ij], 1, &(Rtilde[0][0]), pairdom_len[ij]);

            F_sum = block_matrix(nso, nso);


            for (k = 0; k < nocc; k++) {
                if (k > j) kj_ = (k * (k + 1)) / 2 + j;
                else kj_ = (j * (j + 1)) / 2 + k;

                if (i > k) ik_ = (i * (i + 1)) / 2 + k;
                else ik_ = (k * (k + 1)) / 2 + i;

                kj = abs_ij_map[kj_];
                ik = abs_ij_map[ik_];

                if (pairdom_exist[kj_]) {
                    for (r = 0, L = 0; r < natom; r++) {
                        if (pairdomain[kj][r]) {
                            for (l = aostart[r]; l <= aostop[r]; l++, L++) {
                                for (s = 0, M = 0; s < natom; s++) {
                                    if (pairdomain[kj][s]) {
                                        for (m = aostart[s]; m <= aostop[s]; m++, M++) {
                                            if (k <= j)
                                                F_sum[l][m] -= loF[i][k] * Tempkj[ij_local[ij]][k][M][L];
                                            else
                                                F_sum[l][m] -= loF[i][k] * Tempkj[ij_local[ij]][k][L][M];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                if (pairdom_exist[ik_]) {
                    for (r = 0, L = 0; r < natom; r++) {
                        if (pairdomain[ik][r]) {
                            for (l = aostart[r]; l <= aostop[r]; l++, L++) {
                                for (s = 0, M = 0; s < natom; s++) {
                                    if (pairdomain[ik][s]) {
                                        for (m = aostart[s]; m <= aostop[s]; m++, M++) {
                                            if (k >= i)
                                                F_sum[l][m] -= loF[k][j] * Tempik[ij_local[ij]][k][M][L];
                                            else
                                                F_sum[l][m] -= loF[k][j] * Tempik[ij_local[ij]][k][L][M];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

            }

            temp1 = block_matrix(pairdom_len[ij], nso);
            Sij = block_matrix(pairdom_len[ij], nso);
            X = block_matrix(pairdom_len[ij], pairdom_len[ij]);

            for (r = 0, L = 0; r < natom; r++) {
                if (pairdomain[ij][r]) {
                    for (l = aostart[r]; l <= aostop[r]; l++, L++) {
                        for (m = 0; m < nso; m++) {
                            Sij[L][m] = paoovlp[l][m];
                        }
                    }
                }
            }

            C_DGEMM('n', 'n', pairdom_len[ij], nso, nso, 1, &(Sij[0][0]),
                    nso, &(F_sum[0][0]), nso, 0, &(temp1[0][0]), nso);
            C_DGEMM('n', 't', pairdom_len[ij], pairdom_len[ij], nso, 1, &(temp1[0][0]),
                    nso, &(Sij[0][0]), nso, 0, &(X[0][0]), pairdom_len[ij]);
            free_block(Sij);
            free_block(F_sum);
            free_block(temp1);


            //        fprintf(outfile, "\nX[%d]:\n", ij);
            //        fprintf(outfile,   "======\n");
            //        print_mat(X, pairdom_len[ij], pairdom_len[ij], outfile);

            add_mat(X, Rtilde, Rtilde, pairdom_len[ij], pairdom_len[ij]);

            free_block(X);
            free_block(St);
            free_block(Ft);
        }
        else {
            C_DCOPY(pairdom_len[ij] * pairdom_len[ij], &(Ktilde[ij_local[ij]][0][0]), 1, &(Rtilde[0][0]), 1);
        }



        if (print > 4) {
            fprintf(outfile, "\nRtilde[%d]:\n", ij);
            fprintf(outfile, "======\n");
            print_mat(Rtilde, pairdom_len[ij], pairdom_len[ij], outfile);
        }


        Rbar = block_matrix(pairdom_nrlen[ij], pairdom_nrlen[ij]);
        temp = block_matrix(pairdom_nrlen[ij], pairdom_len[ij]);

        // *** Transform the Rij matrices to an orthogonal basis
        C_DGEMM('t', 'n', pairdom_nrlen[ij], pairdom_len[ij], pairdom_len[ij], 1, &(W[ij_local[ij]][0][0]),
                pairdom_nrlen[ij], &(Rtilde[0][0]), pairdom_len[ij], 0, &(temp[0][0]), pairdom_len[ij]);
        C_DGEMM('n', 'n', pairdom_nrlen[ij], pairdom_nrlen[ij], pairdom_len[ij], 1, &(temp[0][0]),
                pairdom_len[ij], &(W[ij_local[ij]][0][0]), pairdom_nrlen[ij], 0, &(Rbar[0][0]), pairdom_nrlen[ij]);

        //        fprintf(outfile, "\nRbar[%d]:\n", ij);
        //        fprintf(outfile,   "======\n");
        //        print_mat(Rbar, pairdom_nrlen[ij], pairdom_nrlen[ij], outfile);

        free_block(Rtilde);

        Tbar = block_matrix(pairdom_nrlen[ij], pairdom_nrlen[ij]);
        // *** Apply the denominator
        for (a = 0; a < pairdom_nrlen[ij]; a++) {
            for (b = 0; b < pairdom_nrlen[ij]; b++) {
                //            if(fabs(Rbar[a][b]) > 1e-14)
                Tbar[a][b] = -Rbar[a][b] / (loevals[ij][a] + loevals[ij][b] - loF[i][i] - loF[j][j]);
            }
        }
        //        fprintf(outfile, "Tbar[%d] Matrix", ij);
        //        print_mat(Tbar, pairdom_nrlen[ij], pairdom_nrlen[ij], outfile);
        //        fprintf(outfile, "F[%d][%d] = %20.12f\n", loF[i][i]);

        free_block(Rbar);

        if (iter > 0) {
            C_DCOPY(pairdom_len[ij] * pairdom_len[ij], &(T[dmat1][ij_local[ij]][0][0]), 1, &(T[div][ij_local[ij]][0][0]), 1);
            C_DGEMM('n', 'n', pairdom_len[ij], pairdom_nrlen[ij], pairdom_nrlen[ij], 1, &(W[ij_local[ij]][0][0]),
                    pairdom_nrlen[ij], &(Tbar[0][0]), pairdom_nrlen[ij], 0, &(temp[0][0]), pairdom_nrlen[ij]);
            C_DGEMM('n', 't', pairdom_len[ij], pairdom_len[ij], pairdom_nrlen[ij], 1, &(temp[0][0]),
                    pairdom_nrlen[ij], &(W[ij_local[ij]][0][0]), pairdom_nrlen[ij], 1, &(T[div][ij_local[ij]][0][0]),
                    pairdom_len[ij]);
        } else {
            // *** Back Transform the amplitudes into the non-orthogonal basis ***
            C_DGEMM('n', 'n', pairdom_len[ij], pairdom_nrlen[ij], pairdom_nrlen[ij], 1, &(W[ij_local[ij]][0][0]),
                    pairdom_nrlen[ij], &(Tbar[0][0]), pairdom_nrlen[ij], 0, &(temp[0][0]), pairdom_nrlen[ij]);
            C_DGEMM('n', 't', pairdom_len[ij], pairdom_len[ij], pairdom_nrlen[ij], 1, &(temp[0][0]),
                    pairdom_nrlen[ij], &(W[ij_local[ij]][0][0]), pairdom_nrlen[ij], 0, &(T[div][ij_local[ij]][0][0]),
                    pairdom_len[ij]);
        }

        if (print > 3) {
            fprintf(outfile, "\nT[%d]:\n", ij);
            fprintf(outfile, "======\n");
            print_mat(T[div][ij_local[ij]], pairdom_len[ij], pairdom_len[ij], outfile);
        }

        free_block(Tbar);
        free_block(temp);

    }
    free(abs_ij_map);

    // Free the memory used by Temp amplitudes
    if(iter > 0) {

        for (ij = 0; ij < ij_pairs; ij++) {
            if (myid == ij_owner[ij]) {
                for (k = 0; k < nocc; k++) {
                    free_block(Tempkj[ij_local[ij]][k]);
                    free_block(Tempik[ij_local[ij]][k]);
                }
                free(Tempkj[ij_local[ij]]);
                free(Tempik[ij_local[ij]]);
            }
        }
        free(Tempkj);
        free(Tempik);
    }

}


}} // namespace psi::lmp2
