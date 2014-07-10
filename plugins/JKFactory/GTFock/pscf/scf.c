#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <assert.h>
#include <mpi.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <unistd.h>
#include <mkl.h>
#include <mkl_trans.h>
#include <ga.h>
#include <macdecls.h>
#include <sys/time.h>

#include "PFock.h"
#include "CInt.h"
#include "purf.h"
#include "putils.h"


#define NUM_D    1
#define USE_D    0
#define ISSYMM   1

#define TESTSYMM 0

static void usage (char *call)
{
    printf ("Usage: %s <basis> <xyz> "
            "<#nodes per row> "
            "<#nodes per column> "
            "<nprow for purification> "
            "<npcol for purification> "
            "<ntasks> "
            "<#iters>\n", call);
}


#ifdef _MYGA_
static void gather_stat (mystat_t *mystat, int nprocs)
{
    long tmpl;
    double tmpd;
  
    MPI_Reduce (&(mystat->numcre), &tmpl, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    mystat->numcre = (long)((double)tmpl/nprocs);
    
    MPI_Reduce (&(mystat->numdes), &tmpl, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    mystat->numdes = (long)((double)tmpl/nprocs);
    
    MPI_Reduce (&(mystat->numget), &tmpl, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    mystat->numget = (long)((double)tmpl/nprocs);
    
    MPI_Reduce (&(mystat->numput), &tmpl, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    mystat->numput = (long)((double)tmpl/nprocs);
    
    MPI_Reduce (&(mystat->numacc), &tmpl, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    mystat->numacc = (long)((double)tmpl/nprocs);
    
    MPI_Reduce (&(mystat->numsca), &tmpl, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    mystat->numsca = (long)((double)tmpl/nprocs);
    
    MPI_Reduce (&(mystat->numgat), &tmpl, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    mystat->numgat = (long)((double)tmpl/nprocs);
    
    MPI_Reduce (&(mystat->numrdi), &tmpl, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    mystat->numrdi = (long)((double)tmpl/nprocs);
    
    MPI_Reduce (&(mystat->numser), &tmpl, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    mystat->numser = (long)((double)tmpl/nprocs);
    
    MPI_Reduce (&(mystat->curmem), &tmpl, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    mystat->curmem = (long)((double)tmpl/nprocs);
    
    MPI_Reduce (&(mystat->maxmem), &tmpl, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    mystat->maxmem = (long)((double)tmpl/nprocs);
    
    MPI_Reduce (&(mystat->numget_procs), &tmpl, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    mystat->numget_procs = (long)((double)tmpl/nprocs);
    
    MPI_Reduce (&(mystat->numput_procs), &tmpl, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    mystat->numput_procs = (long)((double)tmpl/nprocs);
    
    MPI_Reduce (&(mystat->numacc_procs), &tmpl, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    mystat->numacc_procs = (long)((double)tmpl/nprocs);
    
    MPI_Reduce (&(mystat->numsca_procs), &tmpl, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    mystat->numsca_procs = (long)((double)tmpl/nprocs);
    
    MPI_Reduce (&(mystat->numgat_procs), &tmpl, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    mystat->numgat_procs = (long)((double)tmpl/nprocs);
    
    MPI_Reduce (&(mystat->acctot), &tmpd, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    mystat->acctot = tmpd/nprocs;
    
    MPI_Reduce (&(mystat->accloc), &tmpd, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    mystat->accloc = tmpd/nprocs;
    
    MPI_Reduce (&(mystat->gettot), &tmpd, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    mystat->gettot = tmpd/nprocs;
    
    MPI_Reduce (&(mystat->getloc), &tmpd, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    mystat->getloc = tmpd/nprocs;
    
    MPI_Reduce (&(mystat->puttot), &tmpd, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    mystat->puttot = tmpd/nprocs;
    
    MPI_Reduce (&(mystat->putloc), &tmpd, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    mystat->putloc = tmpd/nprocs;
    
    MPI_Reduce (&(mystat->rditot), &tmpd, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    mystat->rditot = tmpd/nprocs;
    
    MPI_Reduce (&(mystat->rdiloc), &tmpd, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    mystat->rdiloc = tmpd/nprocs;
    
    MPI_Reduce (&(mystat->gattot), &tmpd, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    mystat->gattot = tmpd/nprocs;
    
    MPI_Reduce (&(mystat->gatloc), &tmpd, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    mystat->gatloc = tmpd/nprocs;
    
    MPI_Reduce (&(mystat->scatot), &tmpd, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    mystat->scatot = tmpd/nprocs;
    
    MPI_Reduce (&(mystat->scaloc), &tmpd, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    mystat->scaloc = tmpd/nprocs;
}
#endif /* #ifdef _MYGA_ */


static double compute_energy (purf_t *purf, double *F_block, double *D_block)
{
    double etmp;
    double energy;
    double *H_block;
    int i;
    int j;
    int nrows;
    int ncols;

    etmp = 0.0;
    energy = 0.0;
    H_block = purf->H_block;
    nrows = purf->nrows_purf;
    ncols = purf->ncols_purf;

    if (1 == purf->runpurf)
    {
        for (i = 0; i < nrows; i++)
        {
            for (j = 0; j < ncols; j++)
            {
                F_block[i * ncols + j] += H_block[i * ncols + j];
                etmp += D_block[i * ncols + j] *
                    (H_block[i * ncols + j] + F_block[i * ncols + j]);
            }
        }
    }
    
    MPI_Reduce (&etmp, &energy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    return energy;
}


static void fock_build (PFock_t pfock, BasisSet_t basis,
                        int ispurf, int rowstart, int rowend,
                        int colstart, int colend, int stride,
                        double *D_block, double *F_block)
{
    PFock_setNumDenMat (pfock, NUM_D);       

    // put density matrix
    if (1 == ispurf)
    {
        PFock_putDenMat (pfock, USE_D,
                         rowstart, rowend,
                         colstart, colend,
                         D_block, stride);
    }
    PFock_commitDenMats (pfock);

#ifdef _MYGA_   
    GA_Reset_stats ();
#endif

    PFock_computeFock (pfock, basis);      
        
    // get fock matrix
    if (1 == ispurf)
    {
        PFock_getMat (pfock, USE_D, PFOCK_MAT_TYPE_F,
                      rowstart, rowend,
                      colstart, colend,
                      F_block, stride);
    }
}


int main (int argc, char **argv)
{
    // MPI comm
    int nprocs,myrank;
    //Number of processors per row adn column of fock matrix (nblks_fock is idk, but told to set to 5)
    int nprow_fock,npcol_fock,nblks_fock;
    //Same as for Fock, except now for purification
    int nprow_purf,npcol_purf;
    //Pointer to a BasisSet
    BasisSet_t basis;
    //Pointer to a Fock matrix
    PFock_t pfock;
    //The number of shells, atoms, and basis functions in our molecule
    int nshells, natoms, nfunctions;
    //Pointer to a denisty matrix object
    purf_t *purf;    
    //Iteration number, and max number of iterations
    int iter,niters;
    //Energies before and after the iteration
    double energy,energy0;
    //Where blocks start and end
    int rowstart,rowend;
    int colstart,colend;
    int stride;

    void *bsbuf;
    int bsbufsize;
    double t1;
    double t2;
    int it;
    
    // init MPI
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
    MPI_Comm_size (MPI_COMM_WORLD, &nprocs);
#if 0
    char hostname[1024]; 
    gethostname (hostname, 1024);  
    printf ("Rank %d of %d running on node %s\n", myrank, nprocs, hostname);
#endif

    // create basis set
    CInt_createBasisSet (&basis);
    
    // input parameters and load basis set
    if (myrank == 0)
    {
        if (argc < 8)
        {
            usage (argv[0]);
            MPI_Finalize ();
            exit (0);
        }
        
        // init parameters
        nprow_fock = atoi (argv[3]);
        npcol_fock = atoi (argv[4]);
        nprow_purf = atoi (argv[5]);
        npcol_purf = atoi (argv[6]);
        nblks_fock = atoi (argv[7]);
        niters = atoi (argv[8]);
        printf("NProcs= %d ,nprow= %d, npcol= %d\n",nprocs,nprow_fock,npcol_fock);
        assert (nprow_fock * npcol_fock == nprocs);
        assert (nprow_purf * npcol_purf <= nprocs);
        assert (niters > 0);

        // load basis set
        CInt_loadBasisSet (basis, argv[1], argv[2]);
        nshells = CInt_getNumShells (basis);
        natoms = CInt_getNumAtoms (basis);
        nfunctions = CInt_getNumFuncs (basis);
        assert (nprow_fock <= nshells && npcol_fock <= nshells);
        assert (nprow_purf <= nfunctions && npcol_purf <= nfunctions);
        
        printf ("Job information:\n");
        printf ("  #atoms = %d\n", natoms);
        printf ("  #shells = %d\n", nshells);
        printf ("  #functions = %d\n", nfunctions);  
        printf ("  fock using %d (%dx%d) nodes\n", nprow_fock*npcol_fock, nprow_fock, npcol_fock);
        printf ("  purification using %d (%dx%d) nodes\n", nprow_purf*npcol_purf, nprow_purf, npcol_purf);
        printf ("  #tasks = %d (%dx%d)\n",
            nblks_fock * nblks_fock * nprow_fock * nprow_fock,
            nblks_fock * nprow_fock, nblks_fock * nprow_fock);       
    }

    MPI_Bcast (&nprow_fock, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&npcol_fock, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&nprow_purf, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&npcol_purf, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&nblks_fock, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&niters, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&natoms, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&nshells, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&nfunctions, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // broadcast basis set
    if (myrank == 0)
    {
        CInt_packBasisSet (basis, &bsbuf, &bsbufsize);
        MPI_Bcast (&bsbufsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast (bsbuf, bsbufsize, MPI_CHAR, 0, MPI_COMM_WORLD);
    }
    else    
    {
        MPI_Bcast (&bsbufsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
        bsbuf = (void *)malloc (bsbufsize);
        assert (bsbuf != NULL);    
        MPI_Bcast (bsbuf, bsbufsize, MPI_CHAR, 0, MPI_COMM_WORLD);
        CInt_unpackBasisSet (basis, bsbuf);
    }

    // init PFock
    PFock_GAInit (nfunctions, nprow_fock, npcol_fock, 1, 0, 0);           
    PFock_create (basis, nprow_fock, npcol_fock, nblks_fock, NUM_D, ISSYMM, 1e-10, &pfock);

    // init purf
    purf = create_purf (basis, nprow_purf, npcol_purf);    
    init_oedmat (basis, pfock, purf, nprow_fock, npcol_fock);  
    if (purf->runpurf == 1)  
    {
        memcpy (purf->F_block, purf->H_block,
                purf->nrows_purf * purf->ncols_purf * sizeof(double));
        memset (purf->DD_block, 0,
                sizeof(double) * purf->nrows_purf * purf->ncols_purf);
    }
    energy0 = 0.0;
    
    if (myrank == 0)
    {
        printf ("Computing SCF ...\n");
    }
    t1 = MPI_Wtime ();
    
    rowstart = purf->srow_purf;
    rowend = purf->nrows_purf + rowstart - 1;
    colstart = purf->scol_purf;
    colend = purf->ncols_purf + colstart - 1;
    stride = purf->ncols_purf;
    
    // main loop
    for (iter = 0; iter < niters; iter++)
    {
        if (myrank == 0)
        {
            printf ("  iter %d\n", iter);
        }
        
        // purification
        t1 = MPI_Wtime ();
        it = compute_purification (purf, purf->F_block, purf->D_block);
        t2 = MPI_Wtime ();        
        if (myrank == 0)
        {
            printf ("    purification takes %.3lf secs, %d\n", t2 - t1, it);
        }

        // fock matrix construction
        t1 = MPI_Wtime ();
        fock_build (pfock, basis, purf->runpurf,
                    rowstart, rowend, colstart, colend,
                    stride, purf->D_block, purf->FF_block);
        t2 = MPI_Wtime ();
        if (myrank == 0)
        {
            printf ("    fock takes %.3lf secs\n", t2 - t1);
        }

        energy = compute_energy (purf, purf->FF_block, purf->D_block);
        if (myrank == 0)
        {
            printf ("    energy %.8lf, %le\n", energy, fabs(fabs(energy) - fabs(energy0)));
        }
        energy0 = energy;
                
        // correct F and D
        t1 = MPI_Wtime ();
        correction (purf);
        t2 = MPI_Wtime ();
        if (myrank == 0)
        {
            printf ("    correction takes %.3lf secs\n", t2 - t1);
        }      
    } /* for (iter = 0; iter < NITERATIONS; iter++) */
    
    t2 = MPI_Wtime ();
    if (myrank == 0)
    {
        printf ("  in total takes %.3lf secs: %.3lf secs/iters\n",
            t2 - t1, (t2 - t1) / niters);
        printf ("  Done\n");
    }

#ifndef _MYGA_    
    if (myrank == 0)
    {
        GA_Print_stats ();
    }
#else
    mystat_t *mystat;
    mystat = (mystat_t *)malloc (sizeof(mystat_t));
    GA_Get_stats (mystat);
    gather_stat (mystat, nprow_fock*npcol_fock);   
    if (myrank == 0)
    {
        GA_Print_stats ();        
        GA_Put_stats (mystat);
    }  
    free (mystat);
#endif /* #ifndef _MYGA_ */
  
    destroy_purf (purf);
    PFock_destroy (pfock);
    CInt_destroyBasisSet (basis);
   
    MPI_Finalize ();
    
    return 0;
}
