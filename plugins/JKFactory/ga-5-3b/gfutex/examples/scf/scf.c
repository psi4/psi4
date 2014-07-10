/* C version of the GA SCF code */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <ga.h>
#include <macdecls.h>

#include "cscc.h"
#include "integ.h"
#include "input.h"
#include "output.h"
#include "oneel.h"
#include "twoel.h"
#include "scf.h"
#include "timer.h"

// #include "global.h"
// #include "tcgmsg.h"

#include "sum.h"

//function declaration
void makesz(double* schwmax);
void ininrm(void);
double s(int i,int j);
void makden(void);
int nxtask(int nproc);
int verify_zero_chunk(double chunk[][ichunk]);
void damp(double* fac);
void prnout(int iter, double energy, double deltad, double tester);
double dendif(void);
double testfock(void);
void shiftfock(double shift);
void prnfin(double energy);
void diagon(double* tester, int iter);
void makeob(void);
void denges(void);
void setarrays(void);
void closearrays(void);
void makoverlap(void);
void setoverlap(double *a,int* lo,int* hi, int ld1,int ld2);
void print_GA_block(int g_a);
void print_GA_block_ij(int g_a,int tlo);
void dump_chunk(double *a,int ld1,int ld2);
// double abs_double(double x); 

#ifdef CLOG
FILE *clog = NULL;
#endif

double enrep, q[maxatom], ax[maxatom], ay[maxatom], az[maxatom],
  x[maxnbfn], y[maxnbfn], z[maxnbfn], expnt[maxnbfn], rnorm[maxnbfn];
long long int iky[maxnbfn], nocc, nbfn, nnbfn;
long long int icut1,icut2,icut3,icut4; 
int natom; //long long int --> long 

double eigv[maxnbfn];
int g_counter, g_dens, g_fock, g_tfock, g_schwarz, g_work, g_ident, g_orbs;
int g_tmp, g_proc; //temporay global array for storage major transformation

int main(int argc, char **argv) 
{
  long long int nints, maxint; 

  //     CAUTION: int precision requirements;
  //     nints, maxint, etc. are proportional to the number of basis functions;
  //     to the fourth power! 216^4 is greater than the largest number;
  //     that can be represented as a 32-bit signed interger, so 64-bit;
  //     arithmetic is needed to count integrals when calculating more than;
  //     14 Be atoms with 15 basis functions each. Since integrals are counted;
  //     over all iterations, 18 iterations with 7 atoms can result in precision;
  //     problems. Note that the wave function could be calculated correctly;
  //     for much larger basis sets without 64-bit ints because the required;
  //     indexing is usually proportional to nbfn^2, which is good to 46,340;
  //     basis functions, except that the task counter runs as (nbfn/ichunk)^4,;
  //     so with ichunk = 10, 32-bit ints yield correct wavefunctions out to;
  //     2145 basis functions (maxatom=143), or 4290 (maxatom=286) with ichunk =;
  //     20, ...;
  //;
  //     This warning applies to the Global Arrays implementation as well!;
  //     functions of special concern are GA_igop and NGA_Read_inc.;
  //;
#define USE_TRANSFORM 1
#define MPI 1

  int heap, stack;
  double tinit=0.0, tonel=0.0, ttwoel=0.0, tdiag=0.0, tdens=0.0, tprint=0.0;
  double eone=0.0, etwo=0.0, energy=0.0, deltad=0.0;
      
  //implicit variables in fortran
  int me;
  double frac;
  int iter;
  double scale;
  int nproc = 0;
  double rjunk;
  double totsec, tester;
  double schwmax;
#ifdef CLOG
  char fname[50];
#endif

  //     initalize the parallel message passing environment;
#ifdef MPI
  MPI_Init(&argc,&argv);
#else
  //pbeginf();
  tcg_pbeginf();
#endif
      
  GA_Initialize();
      
  //   Allocate memory;
  heap =  32000000;
  stack = 64000000;
  if (!MA_init(MT_DBL, stack, heap))
    GA_Error("ma_init failed",-1);

#if 1
  GFInitialize();
#endif

  me = GA_Nodeid();
  nproc = GA_Nnodes();
      
  //     initialize a bunch of stuff and initial density matrix;
  rjunk = timer();

#ifdef CLOG
  sprintf(fname, "clog.dat.%d", me);
  clog = fopen(fname, "w");
#endif

  //     get input from file be.inpt;
  input();

  //     create and allocate global arrays;
  setarrays();

  if (me == 0) {
    printf(" bytes of memory used by node 0: %lld\n\n", GA_Inquire_memory());
  }

  ininrm();

  //     create initial guess for density matrix by using single atom;
  //     densities;
  denges();

#if 1
  GA_Print_distribution(g_schwarz);
  GA_Print_distribution(g_dens);
#endif

#if USE_TRANSFORM
  //     make initial orthogonal orbital set for solution method using;
  //     similarity transform;
  makeob();
#endif
  //     make info for sparsity test;
  makesz(&schwmax);

  tinit = timer();

  //     print preliminary data before any long compute segments start;
  if (me == 0)
    printf("\n");

      //     ^* iterate ^*;
  for (iter = 0; iter < mxiter; iter++) {
    double s;

    //     make the one particle contribution to the fock matrix;
    //     and get the partial contribution to the energy;
    oneel(schwmax, &eone, nproc);
    tonel = tonel + timer();

#if 0
    s = sum(g_fock);
    GA_Dgop(&s, 1, "+");

    if (me == 0)
      printf("g_fock sum after oneel: %f\n", s);

    s = sum(g_orbs);
    GA_Dgop(&s, 1, "+");

    if (me == 0)
      printf("g_orbs sum after oneel: %f\n", s);
#endif
        
    //     compute the two particle contributions to the fock matrix and;
    //     get the total energy.;
    twoel(schwmax, &etwo, nproc);
    ttwoel = ttwoel + timer();

#if 0
    s = sum(g_fock);
    GA_Dgop(&s, 1, "+");

    if (me == 0)
      printf("g_fock sum after twoel: %f\n", s);
#endif

    //     Diagonalize the fock matrix. The diagonalizers used in this;
    //     are actually sequential, not parallel.;
    diagon(&tester, iter);
    tdiag = tdiag + timer();

#if 0
    s = sum(g_fock);
    GA_Dgop(&s, 1, "+");

    if (me == 0)
      printf("g_fock sum after diagon: %f\n", s);

    s = sum(g_orbs);
    GA_Dgop(&s, 1, "+");

    if (me == 0)
      printf("g_orbs sum after diagon: %f\n", s);
#endif

    //     make the new density matrix in g_work from orbitals in g_orbs,;
    //     compute the norm of the change in the density matrix and;
    //     then update the density matrix in g_dens with damping.;
    makden();
         
    deltad = dendif();

    if (iter == 0)
      scale = 0.00;
    else if (iter < 5) {
      if (nbfn > 60)
	scale = 0.50;
      else
	scale = 0.00;
    } 
    else
      scale = 0.00;

    damp(&scale);
    tdens = tdens + timer();
          
    //     add up energy and print out convergence information;
    if (me == 0) {
      energy = enrep + eone + etwo;
      prnout(iter, energy, deltad, tester);
      tprint = tprint + timer();
    }

    //     if converged exit iteration loop;
    if (deltad < tol)
      goto L20;

#if GA_VERSION_MAJOR >= 5
    GA_Llgop(&icut4, 1, "+");
#else
    GA_Igop(&icut4, 1, "+");
#endif

    if(icut4 == 0) {
      //     something has gone wrong--print what you know and quit.;
      printf("no two-electron integrals computed!\n");
      goto L20;
    }

  }

  iter = iter - 1;
  if(me == 0)
    printf("SCF failed to converge in %d iters\n", iter);
  //...v....1....v....2....v....3....v....4....v....5....v....6....v....7..;
  //;
  //     finished ... print out eigenvalues and occupied orbitals;
  //;

 L20:

#if GA_VERSION_MAJOR >= 5
  GA_Llgop(&icut1, 1, (char*) "+");
  GA_Llgop(&icut2, 1, (char*) "+");
  GA_Llgop(&icut3, 1, (char*) "+");
#else
  GA_Igop(&icut1, 1, (char*) "+");
  GA_Igop(&icut2, 1, (char*) "+");
  GA_Igop(&icut3, 1, (char*) "+");
#endif

  if (me == 0) {
    //     print out timing information;
    prnfin(energy);
    printf("      init       onel      twoel       diag       dens       print       ncpu    \n");
    printf("      ----       ----      -----       ----       ----       ----        ----    \n");
    totsec = tinit + tonel + ttwoel + tdiag + tdens + tprint;
    printf("%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f        %d",tinit, tonel, ttwoel, tdiag,
	   tdens,tprint,nproc);
    printf("\n elapsed time in seconds %10.2f\n\n",totsec);
    //     print out information on # integrals evaluated each iteration;
    nints = icut1 + icut2 + icut3;
    frac = (double)icut3 / (double)nints;
    printf("No. of integrals screened or computed (all iters) \n\n");
    printf("-------------------------------------\n");
    printf("     failed #ij    test failed #kl    test #compute    #total    fraction\n");
    printf("     ----------    ----------------   -------------    ------    --------\n");
    printf("%15lld %15lld %15lld %15lld %9.6f\n", icut1, icut2, icut3, nints, frac);
         
    maxint = nbfn;
    maxint = pow(maxint, 4) * (iter + 1); 
    if(nints != maxint) {
      printf("Inconsistent number of integrals, should be %ld\n", maxint);
      printf("Note: largest 32-bit int is 2,147,483,647\n");
      printf("Probably due to insufficient int precision in GA.\n");
    }
#ifndef MPI
    tcg_stats();
#endif
  }

  closearrays();

#ifdef CLOG
  fprintf(stderr, "Before GFFinalize()\n");
  fflush(stderr);
#endif

#ifdef CLOG
  fclose(clog);
#endif

#if 1
  GFFinalize();
#endif

#ifdef CLOG
  fprintf(stderr, "After GFFinalize()\n");
  fflush(stderr);
#endif
      
  GA_Terminate();

#ifdef CLOG
  fprintf(stderr, "After GA_Terminate()\n");
  fflush(stderr);
#endif

#ifdef MPI
  MPI_Finalize();
#else
  tcg_pend();
#endif

}

void makesz(double *schwmax)
{
  double work[ichunk][ichunk];
  int lo[2],hi[2],i,j,iloc,jloc,ld;
  int dotask;

  int lo_c[2],hi_c[2]; //column-major inter-patch switch

  //implicit declared variable in fortran
  double gg;
    
  //schwarz(ij) = (ij|ij) for sparsity test;
  icut1 = 0;
  icut2 = 0;
  icut3 = 0;

  GA_Zero(g_schwarz); 
  GA_Zero(g_counter);
  *schwmax = 0.00;
  dotask = next_chunk(lo, hi);
  ld = ichunk;

  while (dotask) {
    for (i = lo[0]; i <= hi[0]; i++) {
      iloc = i - lo[0];
      for (j = lo[1]; j <= hi[1]; j++) {
	jloc = j - lo[1]; 
	g(&gg, i, j, i, j);
	work[iloc][jloc] = sqrt(gg);
	//work[jloc][iloc] = sqrt(gg);
	*schwmax = MAX((*schwmax), work[iloc][jloc]);
	// *schwmax = MAX((*schwmax), work[jloc][iloc]);
      }
    }

    // lo_c[0]=lo[1];
    // lo_c[1]=lo[0];
    // hi_c[0]=hi[1];
    // hi_c[1]=hi[0];

    NGA_Put(g_schwarz, lo, hi, work, &ld);
    // NGA_Put(g_schwarz,lo_c,hi_c,work,&ld); //column-major inter-patch switch
    dotask = next_chunk(lo, hi);
  }
  GA_Dgop(schwmax, 1, "max");

  return;
}

void ininrm(void)  
{
  long long int maxint;
  //implicit declared variables in fortran
  int i;
      
  //     write a little welcome message;
  maxint = nbfn;
  maxint = pow(maxint, 4);

  if (GA_Nodeid()==0) {
    printf(" Example Direct Self Consistent Field Program \n");
    printf(" -------------------------------------------- \n\n");
    printf(" no. of atoms .............. %5d\n", natom);
    printf(" no. of occupied orbitals .. %5ld\n", nocc);
    printf(" no. of basis functions .... %5ld\n", nbfn);
    printf(" basis functions^4 ......... %5ld\n", maxint);
    printf(" convergence threshold ..... %9.4f\n", tol);
    printf(" chunk size .................%5d\n", ichunk);
  }

  //     generate normalisation coefficients for the basis functions;
  //     and the index array iky;
  for (i = 0; i < nbfn; i++) {
    //iky[i] = i*(i-1)/2;
    iky[i] = (i + 1) * i / 2;
  }

  for (i = 0;i < nbfn;i++) {
    rnorm[i] = pow((expnt[i] * 2.00 / pi), 0.750);
  }

  //     initialize common for computing f0;
  setfm();
}

double h(int i,int j) 
{
//vd$r novector;
//vd$r noconcur;
//     generate the one particle hamiltonian matrix element;
//     over the normalized primitive 1s functions i and j;
      double ret;
      
      //implicit declared variables in fortran
      double f0val = 0.00;
      double sum = 0.00;
      double facij,expij,repij;
      int iat;
      //long iat;
      double xp,yp,zp,rpc2;
      double rab2;
     
      rab2 = (x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]) + (z[i]-z[j])*(z[i]-z[j]);
      facij = expnt[i]*expnt[j]/(expnt[i]+expnt[j]);
      expij = exprjh(-facij*rab2);
      repij = (2.00*pi/(expnt[i]+expnt[j])) * expij;

//     first do the nuclear attraction integrals;

      for (iat = 0;iat < natom;iat++) 
      {
         xp = (x[i]*expnt[i] + x[j]*expnt[j])/(expnt[i]+expnt[j]);
         yp = (y[i]*expnt[i] + y[j]*expnt[j])/(expnt[i]+expnt[j]);
         zp = (z[i]*expnt[i] + z[j]*expnt[j])/(expnt[i]+expnt[j]);
         rpc2 = (xp-ax[iat])*(xp-ax[iat]) + (yp-ay[iat])*(yp-ay[iat]) + (zp-az[iat])*(zp-az[iat]);
//;
         f0(&f0val, (expnt[i]+expnt[j])*rpc2);
         sum = sum - repij * q[iat] * f0val;
        }
//     add on the kinetic energy term;

      sum = sum + facij*(3.00-2.00*facij*rab2) *
            pow((pi/(expnt[i]+expnt[j])),1.50) * expij;

//     finally multiply by the normalization constants;
      ret = sum * rnorm[i] * rnorm[j];
      return ret;
}

double s(int i, int j) 
{
  //     generate the overlap matrix element between the normalized;
  //     primitve gaussian 1s functions i and j;
  double ret;

  //implicit declared variables in fortran
  double rab2, facij;

  rab2 = (x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]) +
    (z[i] - z[j]) * (z[i] - z[j]);
  facij = expnt[i] * expnt[j] / (expnt[i] + expnt[j]);
  ret = pow((pi / (expnt[i] + expnt[j])), 1.50) * exprjh(-facij * rab2) * rnorm[i] * rnorm[j];
  return ret;
}

void makden(void) 
{
  double work[ichunk][ichunk], orbsi[maxnbfn][ichunk]; //maxnbfn
  double orbsj[maxnbfn][ichunk]; //maxnbfn
  int lo[2], hi[2], tlo[2], thi[2], i, j, iloc, jloc, ld;
  int dotask;

  //implicit declared variables in fortran
  double p;
  int k;

  //transform between column-major row-major
  // double a[ichunk][ichunk], b[ichunk][ichunk];
  // int lo_c[2],hi_c[2]; 

  //     generate density matrix from orbitals in g_orbs. the first;
  //     nocc orbitals are doubly occupied.;
  // GA_Zero(g_counter);
  // dotask = next_chunk(lo, hi);
  ld = ichunk;

  //transform g_orbs from column-major to row-major 
  /*   while (dotask) { */
  /*     NGA_Get(g_orbs, lo, hi, a, &ld); */
  /*     for (j = 0; j <= hi[1] - lo[1]; j++) { */
  /*       for (i = 0; i <= hi[0] - lo[0]; i++) { */
  /*   	b[j][i] = a[i][j]; //intra-patch transfrom */
  /*       } */
  /*     } */
  /*     lo_c[0] = lo[1]; //inter-patch transform */
  /*     lo_c[1] = lo[0]; */
  /*     hi_c[0] = hi[1]; */
  /*     hi_c[1] = hi[0]; */
  /*     NGA_Put(g_tmp, lo_c, hi_c, b, &ld); */
  /*     dotask = next_chunk(lo, hi); */
  /*   } */
  /*   GA_Copy(g_tmp, g_orbs); */

  // GA_Transpose(g_orbs, g_tmp);
  // GA_Copy(g_tmp, g_orbs);
      
  GA_Zero(g_counter);
  dotask = next_chunk(lo, hi);
  while (dotask) {
    ld = ichunk;

    tlo[0] = 0;
    thi[0] = nocc - 1;

    tlo[1] = lo[0];
    thi[1] = hi[0];
    NGA_Get(g_orbs, tlo, thi, orbsi, &ld);
        
    tlo[1] = lo[1];
    thi[1] = hi[1];
    NGA_Get(g_orbs, tlo, thi, orbsj, &ld);

    for (i = lo[0]; i <= hi[0]; i++) {
      iloc = i - lo[0];// + 1;
      for (j = lo[1]; j <= hi[1]; j++) {
	jloc = j - lo[1];// + 1;
	p = 0.0;
	for (k = 0; k < nocc; k++) {
	  p = p + orbsi[k][iloc] * orbsj[k][jloc];
	}
	work[iloc][jloc] = 2.0 * p;
      }
    }
    ld = ichunk;
    NGA_Put(g_work, lo, hi, work, &ld);
    dotask = next_chunk(lo, hi);
  }

  // GA_Transpose(g_orbs, g_tmp);
  // GA_Copy(g_tmp, g_orbs);

  // GA_Transpose(g_work, g_tmp);
  // GA_Copy(g_tmp, g_work);

  //transform g_orbs from row-major to column major
  /*   ld = ichunk; */
  /*   GA_Zero(g_counter); */
  /*   dotask = next_chunk(lo, hi); */
  /*   while (dotask) { */
  /*     NGA_Get(g_orbs, lo, hi, a, &ld); */
  /*     for (j = 0; j <= hi[1] - lo[1]; j++) { */
  /*       for (i = 0; i <= hi[0] - lo[0]; i++) { */
  /*   	b[j][i] = a[i][j]; //intra-patch transfrom */
  /*       } */
  /*     } */
  /*     lo_c[0] = lo[1]; //inter-patch transform */
  /*     lo_c[1] = lo[0]; */
  /*     hi_c[0] = hi[1]; */
  /*     hi_c[1] = hi[0]; */
  /*     NGA_Put(g_tmp, lo_c, hi_c, b, &ld); */
  /*     dotask = next_chunk(lo, hi); */
  /*   } */
  /*   GA_Copy(g_tmp, g_orbs); */
  // transform g_work from row-major to column major
  /*   ld = ichunk; */
  /*   GA_Zero(g_counter); */
  /*   dotask = next_chunk(lo, hi); */
  /*   while (dotask) { */
  /*     NGA_Get(g_work, lo, hi, a, &ld); */
  /*     for (j = 0; j <= hi[1] - lo[1]; j++) { */
  /*       for (i = 0; i <= hi[0] - lo[0]; i++) { */
  /*   	b[j][i] = a[i][j]; //intra-patch transfrom */
  /*       } */
  /*     } */
  /*     lo_c[0] = lo[1]; //inter-patch transform */
  /*     lo_c[1] = lo[0]; */
  /*     hi_c[0] = hi[1]; */
  /*     hi_c[1] = hi[0]; */
  /*     NGA_Put(g_tmp, lo_c, hi_c, b, &ld); */
  /*     dotask = next_chunk(lo, hi); */
  /*   } */
  /*   GA_Copy(g_tmp, g_work); */
  return;
}

int nxtask(int nproc) //not used
{
      const int ichunk_local = 10; //in the header file, ichunk is defined as 20, so I rename it as ichunk_local
      static double nleft=0, icount=0;
      
      int ret = 0;

      //implicit declared variables in fortran
      int junk;

//     wrapper round nxtval() to increase granularity;
//     and thus reduce no. of requests to shared counter;
#ifndef MPI
      if(nproc>0) 
      {
         if(nleft==0) 
         {
            icount = tcg_nxtval(nproc) * ichunk_local;
            nleft = ichunk_local;
         }
         ret = icount;
         icount = icount + 1;
         nleft = nleft -1;
      } 
      else 
      {
          nleft = 0;
          ret = 0;
          junk = tcg_nxtval(nproc);
      }
#endif
//     following does dumb static load balancing;
//;
//$$$      if(nproc>0) {
//$$$         if (nleft == 0) {
//$$$            icount = GA_Nodeid();
//$$$            nleft = 1;
//$$$         }
//$$$         nxtask = icount;
//$$$         icount = icount + GA_Nnodes();
//$$$      } else {
//$$$          nleft = 0;
//$$$          nxtask = 0;
//$$$      }

     return ret;
}

int next_chunk(int* lo,int* hi) 
{
  //const int one=1;
  int one = 0;
  int imax, ilo, jlo;

  //implicit declared variables in fortran
  int itask;   
      
  int ret;   

  itask = NGA_Read_inc(g_counter, &one, 1);
  imax = nbfn / ichunk;
  if (nbfn - ichunk * imax > 0)
    imax = imax + 1;
  if (itask < imax * imax) {
    ilo = itask % imax;
    jlo = (itask - ilo) / imax;
    lo[0] = ilo * ichunk; 
    lo[1] = jlo * ichunk; 
    hi[0] = MIN((ilo + 1) * ichunk - 1, nbfn); 
    hi[1] = MIN((jlo + 1) * ichunk - 1, nbfn); 
    ret = 1;
  } 
  else 
    ret = 0;

  return ret;
}

long int acquire_tasks(int numTasks)
{
  int one = 0;
  long int itask;

  itask = NGA_Read_inc(g_counter, &one, numTasks);

  return itask;
} // acquire_tasks

int translate_task(long int itask, int *lo, int *hi, int *ilo, int *jlo, int *klo, int *llo)
{
  long int imax;
  int itmp;
  int ret;

  imax = nbfn / ichunk;
  if (nbfn - ichunk * imax > 0)
    imax = imax + 1;
  if (itask < 0) {
    printf("translate_task: itask negative: %ld imax: %ld nbfn: %ld ichunk: %d\n", itask, imax,
	   nbfn, ichunk);
    printf("probable GA int precision problem if imax^4 > 2^31\n");
    printf("\n"); 
    printf("translate_task\n");
    exit(0);
  }
  if (itask < pow(imax, 4)) {
    *ilo = itask % imax;
    itmp = (itask - (*ilo)) / imax;
    *jlo = itmp % imax;
    itmp = (itmp - (*jlo)) / imax;
    *klo = itmp % imax;
    *llo = (itmp - (*klo)) / imax;
    lo[0] = (*ilo) * ichunk;
    lo[1] = (*jlo) * ichunk;
    lo[2] = (*klo) * ichunk;
    lo[3] = (*llo) * ichunk;
    hi[0] = MIN(((*ilo) + 1) * ichunk - 1, nbfn);
    hi[1] = MIN(((*jlo) + 1) * ichunk - 1, nbfn);
    hi[2] = MIN(((*klo) + 1) * ichunk  -1, nbfn);
    hi[3] = MIN(((*llo) + 1) * ichunk - 1, nbfn);
    ret = 1;
  }
  else
    ret = 0;

  return ret;
} // translate_task

int next_4chunk(int *lo,int *hi,int* ilo,int* jlo,int* klo,int* llo, long int *pitask)
{
  int one = 0;
  long int imax;
  long int itask;
  int itmp;

  int ret;

  itask = NGA_Read_inc(g_counter, &one, 1);
  *pitask = itask;

  imax = nbfn / ichunk;
  if (nbfn - ichunk * imax > 0)
    imax = imax + 1;
  if (itask < 0) {
    printf("next_4chunk: itask negative: %ld imax: %ld nbfn: %ld ichunk: %d\n", itask, imax,
	   nbfn, ichunk);
    printf("probable GA int precision problem if imax^4 > 2^31\n");
    printf("\n"); 
    printf("next_4chunk\n");
    exit(0);
  }
  if (itask < pow(imax, 4)) {
    *ilo = itask % imax;
    itmp = (itask - (*ilo)) / imax;
    *jlo = itmp % imax;
    itmp = (itmp - (*jlo)) / imax;
    *klo = itmp % imax;
    *llo = (itmp - (*klo)) / imax;
    lo[0] = (*ilo) * ichunk;
    lo[1] = (*jlo) * ichunk;
    lo[2] = (*klo) * ichunk;
    lo[3] = (*llo) * ichunk;
    hi[0] = MIN(((*ilo) + 1) * ichunk - 1, nbfn);
    hi[1] = MIN(((*jlo) + 1) * ichunk - 1, nbfn);
    hi[2] = MIN(((*klo) + 1) * ichunk  -1, nbfn);
    hi[3] = MIN(((*llo) + 1) * ichunk - 1, nbfn);
    ret = 1;
  }
  else
    ret = 0;

  return ret;
}

void clean_chunk(double chunk[][ichunk])
{
  int i, j;

  for (i = 0; i < ichunk; i++)
    for (j = 0; j < ichunk; j++)
      chunk[i][j] = 0.0;
} // clean_chunk

int verify_zero_chunk(double chunk[][ichunk])
{
  int i, j, flag;

  flag = 1;

  for (i = 0; (i < ichunk) && flag; i++)
    for (j = 0; (j < ichunk) && flag; j++)
      if (chunk[i][j] != 0.0) {
	flag = 0;
	fprintf(stderr, "chunk[%d][%d] != 0.0, %.16f\n", i, j, chunk[i][j]);
	break;
      }

  return flag;
}

void damp(double *fac)
{
  //    create damped density matrix as a linear combination of;
  //    old density matrix and density matrix formed from new orbitals;
      
  //implicit declared variables in fortran
  double ofac;

  ofac = 1.00 - (*fac);
  GA_Add(fac, g_dens, &ofac, g_work, g_dens);
  return;
}

void prnout(int iter, double energy, double deltad, double tester) 
{
  //     printout results of each iteration;
  if (GA_Nodeid() != 0)
    return;

  printf(" iter= %3d, energy=%15.8f, deltad= %9.7f, deltaf=%9.7f\n",
	 iter, energy, deltad, tester);
  return;
}

double dendif(void)
{
  double xdiff;
  double dens_c[ichunk][ichunk],work_c[ichunk][ichunk];
  int lo[2], hi[2], i, j, ld;
  int dotask;

  // int lo_c[2],hi_c[2]; //column-major inter-patch switch

  double ret;

  //implicit declared variables in fortran
  double denmax;
  //     compute largest change in density matrix elements;
  denmax = 0.00;
  GA_Zero(g_counter);
  dotask = next_chunk(lo,hi);
  ld = ichunk;
  while(dotask) {
    //column-major inter-path switch
    // lo_c[0] = lo[1];
    // lo_c[1] = lo[0];
    // hi_c[0] = hi[1];
    // hi_c[1] = hi[0];
        
    NGA_Get(g_dens, lo, hi, dens_c, &ld);
    NGA_Get(g_work, lo, hi, work_c, &ld);

    //column-major inter-patch switch        
    // NGA_Get(g_dens,lo_c,hi_c,dens_c,&ld);
    // NGA_Get(g_work,lo_c,hi_c,work_c,&ld);
        
    for (j = 0; j <= hi[1] - lo[1]; j++) {
      for (i = 0; i <= hi[0] - lo[0]; i++) {
	//xdiff = abs(dens_c[i][j]-work_c[i][j]); 
	//xdiff = abs_double(dens_c[j][i]-work_c[j][i]); //column-major intra-patch switch 
	xdiff = fabs(dens_c[i][j] - work_c[i][j]); //column-major intra-patch switch 
	if (xdiff > denmax)
	  denmax = xdiff;
      }
    }
    dotask = next_chunk(lo, hi);
  }
  GA_Dgop(&denmax, 1, "max");
  ret = denmax;
  return ret;
}

double testfock(void) 
{
  double xmax, xtmp;
  double work[ichunk][ichunk];
  int lo[2], hi[2], i, j, iloc, jloc, ld;
  int dotask;

  // int lo_c[2],hi_c[2]; //column-major inter-patch switch

  double ret;
  //     compute largest change in density matrix elements;
  xmax = 0.00;
  GA_Zero(g_counter);
  dotask = next_chunk(lo, hi);
  ld = ichunk;
  while(dotask) {
    //column-major inter-patch switch
    // lo_c[0]=lo[1];
    // lo_c[1]=lo[0];
    // hi_c[0]=hi[1];
    // hi_c[1]=hi[0];
    
    NGA_Get(g_fock, lo, hi, work, &ld); 
    // NGA_Get(g_fock, lo_c, hi_c, work, &ld);//column-major inter-patch switch

    for (j = lo[1]; j <= hi[1]; j++) {
      jloc = j - lo[1];
      for (i = lo[0]; i <= hi[0]; i++) {
	iloc = i - lo[0];
	if (i != j) {
	  //xtmp = abs(work[iloc][jloc]);
	  // xtmp = abs_double(work[jloc][iloc]); //column-major intra-patch switch
	  xtmp = fabs(work[iloc][jloc]); //column-major intra-patch switch
	  if (xtmp > xmax)
	    xmax = xtmp;
	}
      }
    }
    dotask = next_chunk(lo, hi);
  }
  GA_Dgop(&xmax, 1, "max");
  ret = xmax;

  return ret;
}

void shiftfock(double shift)
{
  double work[ichunk][ichunk];
  int lo[2], hi[2], i, j, iloc, jloc, ld, icnt;
  int dotask;

  // int lo_c[2],hi_c[2]; //column-major inter-patch switch

  //     compute largest change in density matrix elements;
  GA_Zero(g_counter);
  dotask = next_chunk(lo, hi);
  ld = ichunk;
  while (dotask) {
    //column-major inter-path switch
    // lo_c[0]=lo[1];
    // lo_c[1]=lo[0];
    // hi_c[0]=hi[1];
    // hi_c[1]=hi[0];

    // NGA_Get(g_fock, lo_c, hi_c, work, &ld); //column-major inter-patch switch
    NGA_Get(g_fock, lo, hi, work, &ld);
    icnt = 0;
    for (j = lo[1]; j <= hi[1]; j++) {
      jloc = j - lo[1];
      for (i = lo[0]; i <= hi[0]; i++) {
	iloc = i - lo[0];
	if ((i == j) && (i > (nocc - 1))) {
	  work[iloc][jloc] = work[iloc][jloc] + shift;
	  // work[jloc][iloc] = work[jloc][iloc] + shift; //column-major inter-patch switch
	  icnt = icnt + 1;
	}
      }
    }
    if (icnt > 0)
      //NGA_Put(g_fock, lo_c, hi_c, work, &ld);//column-major inter-patch switch
      NGA_Put(g_fock, lo, hi, work, &ld);
    dotask = next_chunk(lo, hi);
  }
  return;
}

void prnfin(double energy) 
{
  double orbs[maxnbfn][maxnbfn];
  int lo[2],hi[2],ld;

  //     printout final results;
  if (GA_Nodeid() != 0)
    return;

  printf("\n\nfinal energy = %18.11f\n",energy);
  printf("\neigenvalues\n\n");
  output(eigv, 0, MIN(nbfn,nocc+5), 0, 1, nbfn, 1, 1);
  //output(eigv, 1, MIN(nbfn,nocc+5), 1, 1, nbfn, 1, 1);

  return;
}

void g(double *value, int i, int j, int k, int l)
{
  //implicit declared variables in fortran
  double f0val, rab2, rcd2, facij,fackl,exijkl,denom, fac, xp,yp,zp,xq,yq,zq,rpq2;

  //     compute the two electon integral (ij|kl) over normalized;
  //     primitive 1s gaussians;
  f0val = 0.00;
  rab2 = (x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]) +
    (z[i] - z[j]) * (z[i] - z[j]);
  rcd2 = (x[k] - x[l]) * (x[k] - x[l]) + (y[k] - y[l]) * (y[k] - y[l]) +
    (z[k] - z[l]) * (z[k] - z[l]);
  facij = expnt[i] * expnt[j] / (expnt[i] + expnt[j]);
  fackl = expnt[k] * expnt[l] / (expnt[k] + expnt[l]);
  exijkl = exprjh(-facij * rab2 - fackl * rcd2);
  denom = (expnt[i] + expnt[j]) * (expnt[k] + expnt[l]) *
    sqrt(expnt[i] + expnt[j] + expnt[k] + expnt[l]);
  fac = (expnt[i] + expnt[j]) * (expnt[k] + expnt[l]) /
    (expnt[i] + expnt[j] + expnt[k] + expnt[l]);

  xp = (x[i] * expnt[i] + x[j] * expnt[j]) / (expnt[i] + expnt[j]);
  yp = (y[i] * expnt[i] + y[j] * expnt[j]) / (expnt[i] + expnt[j]);
  zp = (z[i] * expnt[i] + z[j] * expnt[j]) / (expnt[i] + expnt[j]);
  xq = (x[k] * expnt[k] + x[l] * expnt[l]) / (expnt[k] + expnt[l]);
  yq = (y[k] * expnt[k] + y[l] * expnt[l]) / (expnt[k] + expnt[l]);
  zq = (z[k] * expnt[k] + z[l] * expnt[l]) / (expnt[k] + expnt[l]);
  rpq2 = (xp - xq) * (xp - xq) + (yp - yq) * (yp - yq) + (zp - zq) * (zp - zq);

  f0(&f0val, fac * rpq2);
  *value = (2.00 * pow(pi, 2.50) / denom) * exijkl * f0val *
    rnorm[i] * rnorm[j] * rnorm[k] * rnorm[l];
  return;
}

void diagon(double *tester, int iter)
{
  //      diagon(fock, orbs, evals, work, tester, iter);
  double r_zero, r_one, shift;

  //implicit declared variables in fortran
  int i, me;

  double test = -1.0, s;

  me = GA_Nodeid();

#if USE_TRANSFORM
  //     use similarity transform to solve standard eigenvalue problem;
  //     (overlap matrix has been transformed out of the problem);
  r_one = 1.00;
  r_zero = 0.00;
  GA_Dgemm('n', 'n', nbfn, nbfn, nbfn, r_one, g_fock, g_orbs, r_zero, g_tfock);

#if 0
  s = sum(g_tfock);
  GA_Dgop(&s, 1, "+");

  if (me == 0)
    printf("g_tfock sum after first dgemm(): %f\n", s);
#endif

  GA_Dgemm('t', 'n', nbfn, nbfn, nbfn, r_one, g_orbs, g_tfock, r_zero, g_fock);

#if 0
  s = sum(g_fock);
  GA_Dgop(&s, 1, "+");

  if (me == 0)
    printf("g_fock sum after second dgemm(): %f\n", s);
#endif

  *tester = testfock();
  shift = 0.00;
  if ((*tester) > 0.30) {
    shift = 0.30;
  }
#if 1
  else {
    if (nbfn > 60)
      shift = 0.10;
    else
      shift = 0.00;
  }
#endif
      
  //if (iter>=2 && shift!=0.00) 
  if (iter >= 1 && shift != 0.00) { //iter 2 in Fotran is iter 1 in C
    shiftfock(shift);
  }
  GA_Copy(g_orbs, g_tfock);
  GA_Diag_std_seq(g_fock, g_work, eigv);

  //     Back transform eigenvectors;
  GA_Dgemm('n', 'n', nbfn, nbfn, nbfn, r_one, g_tfock, g_work, r_zero, g_orbs);

  if (iter>= 1 && shift != 0.00) { //>=2 --> >=1
    for (i = nocc; i < nbfn; i++) {
      eigv[i] = eigv[i] - shift;
    }
  }
#else
  //;
  //     Keep remaking overlap matrix since GA_Diag_seq does not;
  //     guarantee that g_ident is preserved.;
  //;
  makoverlap();
  GA_Diag_seq(g_fock, g_ident, g_orbs, eigv); 
  *tester = 0.00;
#endif
  return;
}

void makeob(void) 
{
  double work[ichunk][ichunk],orbs[ichunk][ichunk];
  double eval[maxnbfn];
  int lo[2],hi[2],ld,me,i,j,iloc,jloc;
  int dotask;

  int lo_c[2],hi_c[2]; //for column-major switch 
  //     generate set of orthonormal vectors by creating a random;
  //     symmetric matrix and solving associated generalized eigenvalue;
  //     problem using the correct overlap matrix.;

  me = GA_Nodeid();
  GA_Zero(g_counter);
  dotask = next_chunk(lo, hi);
  ld = ichunk;

  while (dotask) {
    for (i = lo[0]; i <= hi[0]; i++) {
      iloc = i - lo[0];
      for (j = lo[1]; j <= hi[1]; j++) {
	jloc = j - lo[1];
	work[iloc][jloc] = s(i, j);
	//work[jloc][iloc] = s(i,j); //column-major intra-patch switch
	orbs[iloc][jloc] = 0.5;//drand48();
	//orbs[jloc][iloc] = 0.5; //column-major intra-patch switch
      }
    }

    NGA_Put(g_ident, lo, hi, work, &ld);
    NGA_Put(g_fock, lo, hi, orbs, &ld);
    // lo_c[0]=lo[1];
    //lo_c[1]=lo[0];
    // hi_c[0]=hi[1];
    // hi_c[1]=hi[0];
    // NGA_Put(g_ident,lo_c,hi_c,work,&ld);
    // NGA_Put(g_fock,lo_c,hi_c,orbs,&ld);

    dotask = next_chunk(lo, hi);
  }
  GA_Symmetrize(g_fock);
  GA_Diag_seq(g_fock, g_ident, g_orbs, eval);
  return;
}

void denges(void)
{
  //     Form guess density from superposition of atomic densities in the AO;
  //     basis set ... instead of doing the atomic SCF hardwire for this;
  //     small basis set for the Be atom.;
  int one, itask, lo[2], hi[2], ld;

  double atdens[15][15] = {
    {0.000002,0.000027,0.000129,0.000428,0.000950,0.001180,
     0.000457,-0.000270,-0.000271,0.000004,0.000004,0.000004,
     0.000004,0.000004,0.000004},
    {0.000027,0.000102,0.000987,
     0.003269,0.007254,0.009007,0.003492,-0.002099,-0.002108,
     0.000035,0.000035,0.000035,0.000035,0.000035,0.000035},
    {0.000129,0.000987,0.002381,0.015766,0.034988,0.043433,
     0.016835,-0.010038,-0.010082,0.000166,0.000166,0.000166,
     0.000166,0.000166,0.000166},
    {0.000428,0.003269,0.015766,
     0.026100,0.115858,0.144064,0.055967,-0.035878,-0.035990,
     0.000584,0.000584,0.000584,0.000584,0.000584,0.000584},
    {0.000950,0.007254,0.034988,0.115858,0.128586,0.320120,
     0.124539,-0.083334,-0.083536,0.001346,0.001346,0.001346,
     0.001346,0.001346,0.001346},
    {0.001180,0.009007,0.043433,
     0.144064,0.320120,0.201952,0.159935,-0.162762,-0.162267,
     0.002471,0.002471,0.002471,0.002471,0.002471,0.002471},
    {0.000457,0.003492,0.016835,0.055967,0.124539,0.159935,
     0.032378,-0.093780,-0.093202,0.001372,0.001372,0.001372,
     0.001372,0.001372,0.001372},
    {-0.000270,-0.002099,-0.010038,
     -0.035878,-0.083334,-0.162762,-0.093780,0.334488,0.660918,
     -0.009090,-0.009090,-0.009090,-0.009090,-0.009090,-0.009090},
    {-0.000271,-0.002108,-0.010082,-0.035990,-0.083536,-0.162267,
     -0.093202,0.660918,0.326482,-0.008982,-0.008982,-0.008981,
     -0.008981,-0.008981,-0.008982},
    {0.000004,0.000035,0.000166,
     0.000584,0.001346,0.002471,0.001372,-0.009090,-0.008982,
     0.000062,0.000124,0.000124,0.000124,0.000124,0.000124},
    {0.000004,0.000035,0.000166,0.000584,0.001346,0.002471,
     0.001372,-0.009090,-0.008982,0.000124,0.000062,0.000124,
     0.000124,0.000124,0.000124},
    {0.000004,0.000035,0.000166,
     0.000584,0.001346,0.002471,0.001372,-0.009090,-0.008981,
     0.000124,0.000124,0.000062,0.000124,0.000124,0.000124},
    {0.000004,0.000035,0.000166,0.000584,0.001346,0.002471,
     0.001372,-0.009090,-0.008981,0.000124,0.000124,0.000124,
     0.000062,0.000124,0.000124},
    {0.000004,0.000035,0.000166,
     0.000584,0.001346,0.002471,0.001372,-0.009090,-0.008981,
     0.000124,0.000124,0.000124,0.000124,0.000062,0.000124},
    {0.000004,0.000035,0.000166,0.000584,0.001346,0.002471,
     0.001372,-0.009090,-0.008982,0.000124,0.000124,0.000124,
     0.000124,0.000124,0.000062}};
      
  //implicit declared variables in fortran
  int ioff,i,j;
      
  double atdens_f[15][15]; 

  //   Create initial guess for density matrix in global array;
  GA_Zero(g_dens);
  GA_Zero(g_counter);
  //one = 1;
  one = 0;
  ld = 15;

  //   Correct for a factor of two along the diagonal;
  for (i = 0; i < ld; i++) {
    atdens[i][i] = 2.00 * atdens[i][i];
  }

  //to get a fortran based array by switching the indices
  for (i = 0; i < 15; i++) {
    for (j = 0; j < 15; j++) {
      atdens_f[j][i] = atdens[i][j];
    }
  }

  itask = NGA_Read_inc(g_counter, &one, 1);
      
  while(itask < natom) {
    ioff = itask * 15;
    lo[0] = ioff + 0; //1--->0
    lo[1] = ioff + 0; //1--->0
    hi[0] = ioff + 14;//15--->14
    hi[1] = ioff + 14;//15--->14 
    //NGA_Put(g_dens,lo,hi,atdens,&ld);
    NGA_Put(g_dens, lo, hi, atdens_f, &ld);
    itask = NGA_Read_inc(g_counter, &one, 1);
  }
  GA_Sync();
  return;
}

void setarrays(void)
{
  int one, two, dims[2];
  int status;
  int nproc;

  one = 1;
  two = 2;

  nproc = GA_Nnodes();

  g_counter = GA_Create_handle();
  GA_Set_data(g_counter,one,&one,MT_INT);
  status = GA_Allocate(g_counter);
  GA_Zero(g_counter);
  
  g_proc = GA_Create_handle();
  GA_Set_data(g_proc, one, &nproc, MT_INT);
  status = GA_Allocate(g_proc);
  GA_Zero(g_proc);
      
  dims[0] = nbfn;
  dims[1] = nbfn;
  g_dens = GA_Create_handle();
  GA_Set_data(g_dens, two, dims, MT_DBL);
  status = GA_Allocate(g_dens);
  GA_Zero(g_dens);

  g_schwarz = GA_Create_handle();
  GA_Set_data(g_schwarz, two, dims, MT_DBL);
  status = GA_Allocate(g_schwarz);
  GA_Zero(g_schwarz);

  g_fock = GA_Create_handle();
  GA_Set_data(g_fock, two, dims, MT_DBL);
  status = GA_Allocate(g_fock);
  GA_Zero(g_fock);

  g_tfock = GA_Create_handle();
  GA_Set_data(g_tfock, two, dims, MT_DBL);
  status = GA_Allocate(g_tfock);
  GA_Zero(g_tfock);

  g_work = GA_Create_handle();
  GA_Set_data(g_work, two, dims, MT_DBL);
  status = GA_Allocate(g_work);
  GA_Zero(g_work);

  g_ident = GA_Create_handle();
  GA_Set_data(g_ident, two, dims, MT_DBL);
  status = GA_Allocate(g_ident);
  GA_Zero(g_ident);

  g_orbs = GA_Create_handle();
  GA_Set_data(g_orbs, two, dims, MT_DBL);
  status = GA_Allocate(g_orbs);
  GA_Zero(g_orbs);

  //temporay global array for storage major transformation
  g_tmp = GA_Create_handle();
  GA_Set_data(g_tmp, two, dims, MT_DBL);
  status = GA_Allocate(g_tmp);
  GA_Zero(g_tmp);

  return;
}

void closearrays(void)
{
  GA_Destroy(g_counter);
  GA_Destroy(g_proc);
  GA_Destroy(g_dens);
  GA_Destroy(g_schwarz);
  GA_Destroy(g_fock);
  GA_Destroy(g_tfock);
  GA_Destroy(g_work);
  GA_Destroy(g_ident);
  GA_Destroy(g_orbs);
  //temporay global array for storage major transformation
  GA_Destroy(g_tmp);

  return;
}

void makoverlap(void) //not used due to USE_TRANSFORM 
{
      int me, lo[2], hi[2], ld[2];
      int ld1, ld2;
      double *ptr; //int--->double
 
      me = GA_Nodeid();
      NGA_Distribution(g_ident, me, lo, hi);
      NGA_Access(g_ident, lo, hi, &ptr, ld);
      ld1 = hi[0] - lo[0] + 1;
      ld2 = hi[1] - lo[1] + 1;
      
      setoverlap(ptr,lo,hi,ld1,ld2); 
      NGA_Release(g_ident,lo,hi);
      return;
}

void setoverlap(double *a,int* lo,int* hi, int ld1,int ld2) //not used due to USE_TRANSFORM 
{
      int i,j,ii, jj;
      
      for (i = 0;i < ld1;i++) 
      {
        ii = i + lo[0]; 
        for (j = 0;j < ld2;j++) 
        {
          jj = j + lo[1];
#if USE_TRANSFORM
          if (ii==jj)
            //a[i][j] = 1.0;
            a[ld2*i+j]=1.0;
          else
            //a[i][j] = 0.0;
            a[ld2*i+j] = 0.0;
#else
          //a[i][j] = s[ii][jj];
          a[ld1*i+j] = s[ii][jj];
#endif
        }
      }
      return;
}

void print_GA_block(int g_a) //not used 
{
      int lo[2], hi[2], ld1, ld2, ld;
      double *ptr; //int--->double

      //implicit declared variables in fortran
      int me;

      me = GA_Nodeid();
      NGA_Distribution(g_a, me, lo, hi);
      ld1 = hi[0] - lo[0] + 1;
      ld2 = hi[1] - lo[1] + 1;
      NGA_Access(g_a, lo, hi, &ptr, &ld);
      dump_chunk(ptr,ld1,ld2);
      NGA_Release(g_a,lo,hi);
      return;
}

void print_GA_block_ij(int g_a,int tlo) //not used
{
      int lo[2], hi[2], ld1, ld2, ld;
      double *ptr; //int--->double

      //implicit declared variables in fortran
      int me;

      me = GA_Nodeid();
      NGA_Distribution(g_a, me, lo, hi);
      ld1 = hi[0] - lo[0] + 1;
      ld2 = hi[1] - lo[1] + 1;
      NGA_Access(g_a, &tlo, hi, &ptr, &ld);
      dump_chunk(ptr,ld1,ld2);
      NGA_Release(g_a,lo,hi);
      return;
}

void dump_chunk(double *a,int ld1,int ld2) //not used 
{
      //implicit declared variables in fortran
      int i,j;
      double trace;

      for (i = 0; i < MIN(10,ld1);i++) 
      {
         for (j = 0; j < MIN(10,ld2); j++)
         {
            //printf("%10.4f\n",a[i][j]);
            printf("%10.4f\n",a[MIN(10,ld2)*i+j]);
         }
      }
      trace = 0.00;
      for (i=0;i<ld2;i++) 
      {
         //trace = trace +a[i][i];
         trace = trace + a[MIN(10,ld2)*i+j];
      }
      printf("trace= %10.4f\n",trace);
      return;
}

double contract_matrices(int g_a, int g_b)
{
  int lo[2], hi[2], ptr_a, ptr_b, ld, ld1, ld2;
  double a[ichunk][ichunk],b[ichunk][ichunk];
  double value;
  int dotask;

  //implicit declared variables in fortran
  int i,j;

  double ret;

  //   evalute sum_ij a_ij*b_ij;
  value = 0.00;
  GA_Zero(g_counter);
  dotask = next_chunk(lo, hi);
  ld = ichunk;
  while (dotask) {
    NGA_Get(g_a, lo, hi, a, &ld);
    NGA_Get(g_b, lo, hi, b, &ld);

    for (i = 0; i <= hi[0] - lo[0]; i++) {
      for (j = 0; j <= hi[1] - lo[1]; j++) {
	value += a[i][j] * b[i][j];   
      }
    }

    /*     for (j = 0; j <= hi[0] - lo[0]; j++) //column-major switch */
    /*       { */
    /* 	for (i = 0;i <= hi[1]-lo[1];i++)  */
    /*           { */
    /*             value = value + a[j][i]*b[j][i];    */
    /*           } */
    /*       } */

    dotask = next_chunk(lo, hi);
  }
  GA_Dgop(&value, 1, "+");
  ret = value;
  return ret;
}

#if 0
double abs_double(double x)
{
   double ret;
   if (x<0.0)
      ret=-x;
   else
      ret=x;

   return ret;   
  
}
#endif
