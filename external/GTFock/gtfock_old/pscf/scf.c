#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <assert.h>
#include <mpi.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <unistd.h>
#include <sys/time.h>
#include <libgen.h>

#include "pfock.h"
#include "CInt.h"
#include "purif.h"


#define MAX_NUM_D    5
#define NUM_D        3
#define USE_D_ID     2
#define IS_SYMM      0

static void usage(char *call)
{
    printf("Usage: %s <basis> <xyz> "
           "<#nodes per row> "
           "<#nodes per column> "
           "<np for purification> "
           "<ntasks> " "<#iters>\n", call);
}


/// compute initial guesses for the density matrix
static void initial_guess(PFock_t pfock, BasisSet_t basis, int ispurif,
                          int rowstart, int rowend,
                          int colstart, int colend,
                          double *D_block, int ldD)
{
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    PFock_fillDenMat(0.0, USE_D_ID, pfock);
    
    // load initial guess, only process 0
    if (myrank == 0) {
        int num_atoms = CInt_getNumAtoms(basis);
        for (int i = 0; i < num_atoms; i++) {
            double *guess;
            int spos;
            int epos;
            CInt_getInitialGuess(basis, i, &guess, &spos, &epos);
            int ld = epos - spos + 1;
            PFock_putDenMat(spos, epos, spos, epos, ld, guess, USE_D_ID, pfock);
        }
    }
    PFock_sync(pfock);

    if (1 == ispurif) {
        PFock_getMat(pfock, PFOCK_MAT_TYPE_D, USE_D_ID,
                     rowstart, rowend, colstart, colend,
                     ldD, D_block);
    }
}


/// compute Hartree-Fock energy
static double compute_energy(purif_t * purif, double *F_block, double *D_block)
{
    double etmp = 0.0;
    double energy = 0.0;
    double *H_block = purif->H_block;
    int nrows = purif->nrows_purif;
    int ncols = purif->ncols_purif;
    int ldx = purif->ldx;

    if (1 == purif->runpurif) {
        #pragma omp parallel for reduction(+: etmp)
        for (int i = 0; i < nrows; i++) {
            for (int j = 0; j < ncols; j++) {
                F_block[i * ldx + j] += H_block[i * ldx + j];
                etmp += D_block[i * ldx + j] *
                    (H_block[i * ldx + j] + F_block[i * ldx + j]);
            }
        }
    }
    MPI_Allreduce(&etmp, &energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return energy;
}


/// build a Fock matrix
static void fock_build(PFock_t pfock, BasisSet_t basis,
                       int ispurif, int rowstart, int rowend,
                       int colstart, int colend, int stride,
                       double *D_block, double *F_block)
{
    // put density matrix
    if (1 == ispurif) {
        PFock_putDenMat(rowstart, rowend, colstart, colend,
                        stride, D_block, USE_D_ID, pfock);
    }
    PFock_commitDenMats(pfock);

    // compute Fock matrix
    PFock_computeFock(basis, pfock);
    
    // get Fock matrix
    if (1 == ispurif) {
        PFock_getMat(pfock, PFOCK_MAT_TYPE_F, USE_D_ID,
                     rowstart, rowend, colstart, colend,
                     stride, F_block);
    }
}



static void init_oedmat(BasisSet_t basis, PFock_t pfock,
                        purif_t * purif, int nprow, int npcol)
{
    int myrank;
    double t1;
    double t2;
    int ldx = purif->ldx;

    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
    if (myrank == 0) {
        printf ("Preprocessing one electron matrices ...\n");
    }
    int srow_purif = purif->srow_purif;
    int scol_purif = purif->scol_purif;
    int nrows_purif = purif->nrows_purif;
    int ncols_purif = purif->ncols_purif;
    int erow_purif = srow_purif + nrows_purif - 1;
    int ecol_purif = scol_purif + ncols_purif - 1;

    // compute S and X
    if (myrank == 0) {
        printf("  computing H\n");
    }
    t1 = MPI_Wtime();
    PFock_createOvlMat(pfock, basis);
    if (purif->runpurif == 1) {
        PFock_getOvlMat(pfock, srow_purif, erow_purif, scol_purif, ecol_purif,
                        ldx, purif->S_block);
        PFock_getOvlMat2(pfock, srow_purif, erow_purif, scol_purif, ecol_purif,
                         ldx, purif->X_block);
    }
    PFock_destroyOvlMat(pfock);
    t2 = MPI_Wtime();
    if (myrank == 0) {
        printf("  takes %.3lf secs\n", t2 - t1);
        printf("  Done\n");
    }
    
    // Compute H
    if (myrank == 0) {
        printf("  computing H\n");
    }
    t1 = MPI_Wtime();
    PFock_createCoreHMat(pfock, basis);
    if (purif->runpurif == 1) {
        PFock_getCoreHMat(pfock, srow_purif, erow_purif,
                          scol_purif, ecol_purif, ldx, purif->H_block);
    }
    PFock_destroyCoreHMat(pfock);
    t2 = MPI_Wtime();
    if (myrank == 0) {
        printf("  takes %.3lf secs\n", t2 - t1);
        printf("  Done\n");
    }
}


/// main for SCF
int main (int argc, char **argv)
{
    // init MPI
    int myrank;
    int nprocs;
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if (myrank == 0)  {
        printf("MPI level: %d\n", provided);
    }
#if 0
    char hostname[1024];
    gethostname (hostname, 1024);
    printf ("Rank %d of %d running on node %s\n", myrank, nprocs, hostname);
#endif

    // create basis set    
    BasisSet_t basis;
    CInt_createBasisSet(&basis);

    // input parameters and load basis set
    int nprow_fock;
    int npcol_fock;
    int nblks_fock;
    int nprow_purif;
    int nshells;
    int natoms;
    int nfunctions;
    int niters;
    if (myrank == 0) {
        if (argc != 8) {
            usage(argv[0]);
            MPI_Finalize();
            exit(0);
        }
        // init parameters
        nprow_fock = atoi(argv[3]);
        npcol_fock = atoi(argv[4]);
        nprow_purif = atoi(argv[5]);
        nblks_fock = atoi(argv[6]);
        niters = atoi(argv[7]);
        assert(nprow_fock * npcol_fock == nprocs);
        assert(nprow_purif * nprow_purif * nprow_purif  <= nprocs);
        assert(niters > 0);       
        CInt_loadBasisSet(basis, argv[1], argv[2]);
        nshells = CInt_getNumShells(basis);
        natoms = CInt_getNumAtoms(basis);
        nfunctions = CInt_getNumFuncs(basis);
        assert(nprow_fock <= nshells && npcol_fock <= nshells);
        assert(nprow_purif <= nfunctions && nprow_purif <= nfunctions);
        printf("Job information:\n");
        char *fname;
        fname = basename(argv[2]);
        printf("  molecule:  %s\n", fname);
        fname = basename(argv[1]);
        printf("  basisset:  %s\n", fname);
        printf("  #atoms     = %d\n", natoms);
        printf("  #shells    = %d\n", nshells);
        printf("  #functions = %d\n", nfunctions);
        printf("  fock build uses   %d (%dx%d) nodes\n",
               nprow_fock * npcol_fock, nprow_fock, npcol_fock);
        printf("  purification uses %d (%dx%dx%d) nodes\n",
               nprow_purif * nprow_purif * nprow_purif,
               nprow_purif, nprow_purif, nprow_purif);
        printf("  #tasks = %d (%dx%d)\n",
               nblks_fock * nblks_fock * nprow_fock * nprow_fock,
               nblks_fock * nprow_fock, nblks_fock * nprow_fock);
        int nthreads = omp_get_max_threads();
        printf("  #nthreads_cpu = %d\n", nthreads);   
    }
    int btmp[8];
    btmp[0] = nprow_fock;
    btmp[1] = npcol_fock;
    btmp[2] = nprow_purif;
    btmp[3] = nblks_fock;
    btmp[4] = niters;
    btmp[5] = natoms;
    btmp[6] = nshells;
    btmp[7] = nfunctions;
    MPI_Bcast(btmp, 8, MPI_INT, 0, MPI_COMM_WORLD);
    nprow_fock = btmp[0];
    npcol_fock = btmp[1];
    nprow_purif = btmp[2];
    nblks_fock = btmp[3];
    niters = btmp[4];
    natoms = btmp[5];
    nshells = btmp[6];
    nfunctions = btmp[7];

    // broadcast basis set
    void *bsbuf;
    int bsbufsize;
    if (myrank == 0) {
        CInt_packBasisSet(basis, &bsbuf, &bsbufsize);
        MPI_Bcast(&bsbufsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(bsbuf, bsbufsize, MPI_CHAR, 0, MPI_COMM_WORLD);
    }
    else {
        MPI_Bcast(&bsbufsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
        bsbuf = (void *)malloc(bsbufsize);
        assert(bsbuf != NULL);
        MPI_Bcast(bsbuf, bsbufsize, MPI_CHAR, 0, MPI_COMM_WORLD);
        CInt_unpackBasisSet(basis, bsbuf);  
        free(bsbuf);
    }

    // init PFock
    if (myrank == 0) {
        printf("Initializing pfock ...\n");
    }
    PFock_t pfock;
    PFock_create(basis, nprow_fock, npcol_fock, nblks_fock, 1e-10,
                 MAX_NUM_D, IS_SYMM, &pfock);
    if (myrank == 0) {
        double mem_cpu;
        PFock_getMemorySize(pfock, &mem_cpu);
        printf("  CPU uses %.3lf MB\n", mem_cpu / 1024.0 / 1024.0);
        printf("  Done\n");
    }

    // init purif
    purif_t *purif = create_purif(basis, nprow_purif, nprow_purif, nprow_purif);
    init_oedmat(basis, pfock, purif, nprow_fock, npcol_fock);

    // compute SCF
    if (myrank == 0) {
        printf("Computing SCF ...\n");
    }
    int rowstart = purif->srow_purif;
    int rowend = purif->nrows_purif + rowstart - 1;
    int colstart = purif->scol_purif;
    int colend = purif->ncols_purif + colstart - 1;
    double energy0 = -1.0;
    double totaltime = 0.0;
    double purif_flops = 2.0 * nfunctions * nfunctions * nfunctions;
    double diis_flops;

    // set initial guess
    if (myrank == 0) {
        printf("  initialing D ...\n");
    }
    PFock_setNumDenMat(NUM_D, pfock);
    initial_guess(pfock, basis, purif->runpurif,
                  rowstart, rowend, colstart, colend,
                  purif->D_block, purif->ldx);

    // compute nuc energy
    double ene_nuc = CInt_getNucEnergy(basis);
    if (myrank == 0) {
        printf("  nuc energy = %.8lf\n", ene_nuc);

    }

    MPI_Barrier(MPI_COMM_WORLD);
    // main loop
    double t1, t2, t3, t4;
    for (int iter = 0; iter < niters; iter++) {
        if (myrank == 0) {
            printf("  iter %d\n", iter);
        }
        t3 = MPI_Wtime();

        // fock matrix construction
        t1 = MPI_Wtime();
        fock_build(pfock, basis, purif->runpurif,
                   rowstart, rowend, colstart, colend,
                   purif->ldx, purif->D_block, purif->F_block);
        // compute energy
        double energy = compute_energy(purif, purif->F_block, purif->D_block);
        t2 = MPI_Wtime();
        if (myrank == 0) {
            printf("    fock build takes %.3lf secs\n", t2 - t1);
            if (iter > 0) {
                printf("    energy %.8lf (%.8lf), %le\n",
                       energy + ene_nuc, energy, fabs (energy - energy0));
            }
            else {
                printf("    energy %.8lf (%.8lf)\n", energy + ene_nuc,
                       energy);
            }
        }
        if (iter > 0 && fabs (energy - energy0) < 1e-8) {
            niters = iter + 1;
            break;
        }
        energy0 = energy;

        // compute DIIS
        t1 = MPI_Wtime();
        compute_diis(pfock, purif, purif->D_block, purif->F_block, iter);
        t2 = MPI_Wtime();

        if (myrank == 0) {
            if (iter > 1) {
                diis_flops = purif_flops * 6.0;
            } else {
                diis_flops = purif_flops * 2.0;
            }
            printf("    diis takes %.3lf secs, %.3lf Gflops\n",
                   t2 - t1, diis_flops / (t2 - t1) / 1e9);
        }
        
    #ifdef __SCF_OUT__
        if (myrank == 0) {
            double outbuf[nfunctions];
            char fname[1024];
            sprintf(fname, "XFX_%d_%d.dat", nfunctions, iter);
            FILE *fp = fopen(fname, "w+");
            assert(fp != NULL);
            for (int i = 0; i < nfunctions; i++) {
                PFock_getMat(pfock, PFOCK_MAT_TYPE_F, USE_D_ID,
                             i, i, USE_D_ID, nfunctions - 1,
                             outbuf, nfunctions);
                for (int j = 0; j < nfunctions; j++) {
                    fprintf(fp, "%le\n", outbuf[j]);
                }
            }
            fclose(fp);
        }
    #endif
    
        // purification
        MPI_Barrier(MPI_COMM_WORLD);
        t1 = MPI_Wtime();
        int it = compute_purification(purif, purif->F_block, purif->D_block);
        t2 = MPI_Wtime();
        MPI_Barrier(MPI_COMM_WORLD);
        if (myrank == 0) {
            printf("    purification takes %.3lf secs,"
                   " %d iterations, %.3lf Gflops\n",
                   t2 - t1, it,
                   (it * 2.0 + 4.0) * purif_flops / (t2 - t1) / 1e9);
        }

        t4 = MPI_Wtime ();
        totaltime += t4 - t3;

#ifdef __SCF_TIMING__
        PFock_getStatistics(pfock);
        double purif_timedgemm;
        double purif_timepass;
        double purif_timetr;
        MPI_Reduce(&purif->timedgemm, &purif_timedgemm,
                   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&purif->timepass, &purif_timepass,
                   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&purif->timetr, &purif_timetr,
                   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (myrank == 0) {
            printf("    Purification Statistics:\n");
            printf("      average totaltime = %.3lf\n"
                   "      average timetr    = %.3lf\n"
                   "      average timedgemm = %.3lf, %.3lf Gflops\n",
                   purif_timepass / purif->np_purif,
                   purif_timetr / purif->np_purif,
                   purif_timedgemm / purif->np_purif,
                   (it * 2.0 + 4.0) *
                   purif_flops / (purif_timedgemm / purif->np_purif) / 1e9);
        }
#endif
    } /* for (iter = 0; iter < NITERATIONS; iter++) */

    if (myrank == 0) {
        printf("  totally takes %.3lf secs: %.3lf secs/iters\n",
               totaltime, totaltime / niters);
        printf("  Done\n");
    }

    destroy_purif(purif);
    PFock_destroy(pfock);
    CInt_destroyBasisSet(basis);

    MPI_Finalize();

    return 0;
}
