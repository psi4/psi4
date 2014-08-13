#if HAVE_CONFIG_H
#   include "config.h"
#endif

/*$Id: grid.c,v 1.2 1995-02-02 23:24:11 d3g681 Exp $*/
#include <stdio.h>
#include <math.h>

#ifdef PLOT
#include <sys/types.h>
#include <netinet/in.h>
FILE *plot_file;

#define PLOT_VALUE    1
#define PLOT_ERROR    2
#define PLOT_RESIDUAL 3

static long plot_type = 0;      /* 0 means no plot */
#endif

#include "sndrcv.h"

#if defined(DELTA) || defined(IPSC)
#define htonl(a) (a)
#endif

#if !defined(AIX) && !defined(DECOSF)
extern char *malloc();
#endif
extern void exit();

#define TCG_MAX(a,b) (((a)>(b)) ? (a) : (b))
#define TCG_MIN(a,b) (((a)<(b)) ? (a) : (b))
#define TCG_ABS(a)   (((a)>=0 ) ? (a) : -(a))

#define LO -3.1415926535       /* Phsyical dimensions   */
#define HI  3.1415926535

double scale;                  /* Physical mesh grain   */

long P, ncol_P, nrow_P;         /* P = No. of processes = ncol_P * nrow_P */
long my_col_P, my_row_P;        /* Co-ords of this proc. in grid of processes*/

double *buffer;                /* Buffer used for exchanging rows */

long north, south, east, west;  /* No. of process in that direction on 
                                  process grid or -1 if no-one there */

long col_low = 0;                /* Grid co-ords of top left corner */
long row_low = 0;

/* Mapping from GLOBAL grid index to physical coordinates */
#define MAPCOL(I) (scale*((double) (I+col_low)) + LO)
#define MAPROW(J) (scale*((double) (J+row_low)) + LO)

long cs_exchange = 0;            /* Timing information */
long cs_global = 0;
long cs_interpolate = 0;
long cs_plot = 0;
long cs_total = 0;

/*
void ZeroGrid(grid, ncols, nrows)
     double **grid;
     long ncols, nrows;
{
  long i,j;

  for (i=0; i<=ncols; i++)
    for (j=0; j<=nrows; j++)
      grid[i][j] = 0.0;
}
*/

double Solution(x,y)
     double x,y;
/*
  Model potential that defines solution and boundary conditions.
  This particular choice makes the discretization error quite apparent.
*/
{
/*
  return cos(x)*exp(-y) + 0.1*sin(2.*x)*exp(-2.*y) + 0.1*cos(3.*x)*exp(-3.*y);
*/

  return 0.5 * (cos(x)*exp(-y) + cos(y)*exp(-x) +
        0.1*(sin(2.*x)*exp(-2.*y) + sin(2.*y)*exp(-2.*x)) +
        0.1*(cos(3.*x)*exp(-3.*y) + cos(3.*y)*exp(-3.*x)));
}

double GridError(grid, ncols, nrows, ngrid)
     long ncols, nrows, ngrid;
     double **grid;
/*
  Compute mean absolute error at grid points relative to
  the analytic solution ... error is due to either lack of
  convergence or discretization error.
*/
{
  long i,j, type=9, length=1;
  long start;
  double error = 0.0;

  for (i=1; i<ncols; i++)
    for (j=1; j<nrows; j++)
      error += fabs(grid[i][j] - Solution(MAPCOL(i), MAPROW(j)));

  start = MTIME_();
  DGOP_(&type, &error, &length, "+");
  cs_global += MTIME_() - start;

  return error / ((ngrid-1)*(ngrid-1));
}
  
void BoundaryConditions(grid, ncols, nrows)
     double **grid;
     long ncols, nrows;
/* 
  Fill in b.c.s on a grid
*/
{
  long i,j;

  for (i=0; i<=ncols; i++) {
    if (my_row_P == 0)
      grid[i][0]  = Solution(MAPCOL(i), LO);
    if (my_row_P == (nrow_P-1))
      grid[i][nrows] = Solution(MAPCOL(i), HI);
  }

  for (j=0; j<=nrows; j++) {
    if (my_col_P == 0)
      grid[0][j]  = Solution(LO, MAPROW(j));
    if (my_col_P == (ncol_P-1))
      grid[ncols][j] = Solution(HI, MAPROW(j));
  }
}

void Initialize(grid, ncols, nrows)
     long ncols, nrows;
     double **grid;
/*
  Fill in boundary values and zero initial guess for interior.
*/
{
  long i,j;

  BoundaryConditions(grid, ncols, nrows);

  for (i=1; i<ncols; i++)
    for (j=1; j<nrows; j++)
      grid[i][j] = 0.0;
}

void Exchange(grid, ncols, nrows)
     double ** grid;
     long ncols, nrows;
/*
  Exchange data with neighboring ndoes. In Operate only need
  to exchange red and black elements separately but do NOT do
  this at the moment ... thus are sending twice as much data
  as necessary. However, are still dominated by the latency
  so 'no big whup'.  Interpolate needs to exchange the full
  boundary information (I think).
*/
{
  long type1=1, type2=2, type3=3, type4=4, type5=5, type6=6, type7=7, type8=8;
  long bncols = (ncols+1)*sizeof(double);
  long bnrows = (nrows+1)*sizeof(double);
  long synch = 1, lenmes, nodefrom, i;
  long start = MTIME_();

#define GATHER(k)  for (i=0; i<=ncols; i++) buffer[i] = grid[i][k]
#define SCATTER(k) for (i=0; i<=ncols; i++) grid[i][k] = buffer[i]

  if (my_col_P%2) {
    if (west >= 0) {
      SND_(&type1, (char *) grid[1], &bnrows, &west, &synch);
      RCV_(&type2, (char *) grid[0], &bnrows, &lenmes, &west, 
       &nodefrom, &synch);
    }
    if (east >= 0) {
      SND_(&type3, (char *) grid[ncols-1], &bnrows, &east, &synch);
      RCV_(&type4, (char *) grid[ncols], &bnrows, &lenmes, &east, 
       &nodefrom, &synch);
    }
  }
  else {
    if (east >= 0) {
      RCV_(&type1, (char *) grid[ncols], &bnrows, &lenmes, &east, 
       &nodefrom, &synch);
      SND_(&type2, (char *) grid[ncols-1], &bnrows, &east, &synch);
    }
    if (west >= 0) {
      RCV_(&type3, (char *) grid[0], &bnrows, &lenmes, &west, 
       &nodefrom, &synch);
      SND_(&type4, (char *) grid[1], &bnrows, &west, &synch);
    }
  }

  if (my_row_P%2) {
    if (north >= 0) {
      GATHER(1);
      SND_(&type5, (char *) buffer, &bncols, &north, &synch);
      RCV_(&type6, (char *) buffer, &bncols, &lenmes, &north, 
       &nodefrom, &synch);
      SCATTER(0);
    }
    if (south >= 0) {
      GATHER(nrows-1);
      SND_(&type7, (char *) buffer, &bncols, &south, &synch);
      RCV_(&type8, (char *) buffer, &bncols, &lenmes, &south, 
       &nodefrom, &synch);
      SCATTER(nrows);
    }
  }
  else {
    if (south >= 0) {
      RCV_(&type5, (char *) buffer, &bncols, &lenmes, &south, 
       &nodefrom, &synch);
      SCATTER(nrows);
      GATHER(nrows-1);
      SND_(&type6, (char *) buffer, &bncols, &south, &synch);
    }
    if (north >= 0) {
      RCV_(&type7, (char *) buffer, &bncols, &lenmes, &north, 
       &nodefrom, &synch);
      SCATTER(0);
      GATHER(1);
      SND_(&type8, (char *) buffer, &bncols, &north, &synch);
    }
  }
  cs_exchange += MTIME_() - start;
}

#ifdef PLOT

void ClosePlotFile()
{
  if (!plot_type)
    return;

  if (NODEID_() == 0)
    (void) fclose(plot_file);
}

void OpenPlotFile(maxgrid)
     long maxgrid;
{
  if (!plot_type)
    return;

  if (NODEID_())
    return;              /* Only node 0 needs to do this */

  if (!(plot_file = fopen("plot","w+")))
    Error("OpenPlotFile: failed to open plot file", (long) -1);

  maxgrid = htonl((long) maxgrid); /* For portability */
  if (fwrite((char *) &maxgrid, sizeof(int), 1, plot_file) != 1)
    Error("OpenPlotFile: failed to write maxgrid", (long) -1);

  (void) fflush(plot_file);
}

void InsertPlot(plot_full, ngrid, plot_grid, nrows, ncols, col_low, row_low)
     unsigned char *plot_full, *plot_grid;
     long ngrid, nrows, ncols, col_low, row_low;
{
  long i, j;
  long i_lo = (col_low == 0) ? 0 : 1;
  long j_lo = (row_low == 0) ? 0 : 1;
  unsigned char *temp;

  /* Dink around with i_lo, j_lo so that the interior edges in the
     parallel case are not overwritten with incorrect values */

  if (i_lo)
    plot_grid += nrows;

  for (i=i_lo; i<ncols; i++) {
    temp = plot_full + (i+col_low)*ngrid + row_low + j_lo;
    if (j_lo)
      plot_grid++;
    for (j=j_lo; j<nrows; j++)
      *temp++ = *plot_grid++;
  }
}

void PlotGrid(grid, ngrid, ncols, nrows)
     double **grid;
     long ngrid, ncols, nrows;
{
  unsigned char *plot_grid, *plot_full;

  double d_63_19 = 63.0 / 19.0;
  double rlog2 = 1.0 / log((double) 2.0);
#define LOG2(x)  ((x != 0.0) ? (rlog2 * log(x)) : -1024.0)
#define VALUE(x) (d_63_19 * (LOG2(fabs(x)) + 9.0))
  register long i, j, value;
  register unsigned char *temp;
  long n;
  double residual, factor;
  long net_ngrid = htonl((long) ngrid);
  long msg[4], synch_type=13, anyone = (-1);
  long msg_len=sizeof(msg), msg_type=120, length, node, zero=0, synch=1;
  long start = MTIME_();

  if (!plot_type)
    return;

  n = nrows*ncols;
  
  plot_grid = (unsigned char *) malloc((unsigned) n);
  if (!plot_grid)
    Error("PlotGrid: failed to allocate plot grid", (long) -1);

  temp = plot_grid;

  switch (plot_type) {
  case PLOT_VALUE:
    for (i=0; i<ncols; i++)
      for (j=0; j<nrows; j++) {
    value = VALUE(grid[i][j]);
    value = TCG_MIN(value, 63);
    value = TCG_MAX(value, 0);
    *temp++ = (unsigned char) value;
      }
    break;

  case PLOT_ERROR:
    for (i=0; i<ncols; i++)
      for (j=0; j<nrows; j++) {
    value = VALUE(grid[i][j] - Solution(MAPCOL(i), MAPROW(j)));
    value = TCG_MIN(value, 63);
    value = TCG_MAX(value, 0);
    *temp++ = (unsigned char) value;
      }
    break;

  case PLOT_RESIDUAL:

    Exchange(grid, ncols, nrows);  /* This to get the edges up date */

    factor = 1.0/(scale*scale);
    for (j=0; j<nrows; j++)        /* Edge has residual zero */
      *temp++ = (unsigned char) VALUE(0.0);
    for (i=1; i<ncols; i++) {
      *temp++ = (unsigned char) VALUE(0.0);
      for (j=1; j<nrows; j++) {
    residual = grid[i+1][j] + grid[i-1][j] + grid[i][j+1] +
      grid[i][j-1] - 4.0*grid[i][j];
    value = VALUE(residual*factor);
    value = TCG_MIN(value, 63);
    value = TCG_MAX(value, 0);
    *temp++ = (unsigned char) value;
      }
    }
    break;
    
  default:
    Error("PlotGrid: unknown plot type requested", (long) plot_type);
  }
  
  /* Now have grid with my portion of the plot. This gets sent
     to node 0 which assembles the full grid before dumping to disc */
  
  if (NODEID_()) {
    msg_type = 120;
    msg[0] = nrows; msg[1] = ncols; msg[2] = col_low; msg[3] = row_low;
    SND_(&msg_type, (char *) msg, &msg_len, &zero, &synch);
    msg_type++;
    SND_(&msg_type, (char *) plot_grid, &n, &zero, &synch);
    msg_type++;
  }
  else {
    plot_full = (unsigned char *) malloc((unsigned) (ngrid*ngrid));
    if (!plot_full)
      Error("PlotGrid: failed to allocate plot full", (long) -1);
    InsertPlot(plot_full, ngrid, plot_grid, nrows, ncols, col_low, row_low);

    for (i=1; i<NNODES_(); i++) {
      msg_type = 120;
      RCV_(&msg_type, (char *) msg, &msg_len, &length, &anyone, &node, &synch);
      msg_type++;
      RCV_(&msg_type, (char *) plot_grid, &n, &length, &node, &node, &synch);
      msg_type++;
      InsertPlot(plot_full, ngrid, plot_grid, msg[0], msg[1], msg[2], msg[3]);
    }

    if (fwrite((char *) &net_ngrid, sizeof(int), 1, plot_file) != 1)
      Error("PlotGrid: failed to write ngrid", (long) -1);
    n = ngrid*ngrid;
    if (fwrite(plot_full, 1, n, plot_file) != n)
      Error("PlotGrid: failed writing to file", (long) n);
    (void) free((char *) plot_full);
    (void) fflush(plot_file);
    (void) fflush(plot_file);
  }
    
  BRDCST_(&synch_type, &zero, &zero, &zero);    /* Synchronization */
  (void) free((char *) plot_grid);

  cs_plot += MTIME_() - start;
}
#endif

void PrintGrid(grid, ncols, nrows)
     double **grid;
     long ncols, nrows;
{
  long i,j,newline;

  for (i=0; i<=ncols; i++) {
    (void) printf("Column %d \n",i+col_low);
    newline = 0;
    for(j=0; j<=nrows; j++) {
      if (grid[i][j] != 0.0) {
    (void) printf("(%3d,%10.4f) ",j+row_low,grid[i][j]);
    if (++newline == 4) {
      (void) printf("\n"); newline = 0;
    }
      }
    }
    if (newline) 
      (void) printf("\n");
    (void) printf("\n");
  }
}

double Operate(grid, ncols, nrows, ngrid, do_sums)
     long ncols, nrows, ngrid, do_sums;
     double **grid;
/*
  Update grid in place according to the simple rule

     new[i][j] = 0.25 * 
       (old[i-1][j] + old[i+1][j] + old[i][j-1] + old[i][j+1]);

  for 0<i<ncols, 0<j<nrows.

  To improve convergence rate for large grids use red/black
  checkerboard -> update red then update black using new reds
  (gauss-seidel red-black w-relaxation).

  Return the mean abs. error (value of the laplacian) on the old grid.
*/
{
  double residual = 0.0;
  double omega = 0.94;          /* Over relaxation parameter */
  long i,j, type=10, length=1, jlo;
  double *gg, *ggm, *ggp;

  Exchange(grid, ncols, nrows);
  
  /* Update red */
  for (i=1; i<ncols; i++) {
    jlo = 1+((i-1+col_low+row_low)%2);
    gg = grid[i]; ggm = grid[i-1]; ggp = grid[i+1];
#ifdef ALLIANT
#pragma safe gg
#pragma safe ggp
#pragma safe ggm
#endif
#ifdef ARDENT
#pragma IVDEP
#endif
    for (j=jlo; j<nrows; j+=2) {
      double new, diff;
      new = 0.25 * (ggp[j] + ggm[j] + gg[j-1] + gg[j+1]);
      new = new + omega*(new - gg[j]);
      diff = new - gg[j];
      residual += TCG_ABS(diff);
      gg[j] = new;
    }
  }

  Exchange(grid, ncols, nrows);

  /* Update black */
  for (i=1; i<ncols; i++) {
    jlo = 1+((i+col_low+row_low)%2);
    gg = grid[i]; ggm = grid[i-1]; ggp = grid[i+1];
#ifdef ALLIANT
#pragma safe gg
#pragma safe ggp
#pragma safe ggm
#endif
#ifdef ARDENT
#pragma IVDEP
#endif
    for (j=jlo; j<nrows; j+=2) {
      double new, diff;
      new = 0.25 * (ggp[j] + ggm[j] + gg[j-1] + gg[j+1]);
      new = new + omega*(new - gg[j]);
      diff = new - gg[j];
      residual += TCG_ABS(diff);
      gg[j] = new;
    }
  }

  if (do_sums) {
    long start = MTIME_();
    DGOP_(&type, &residual, &length, "+");
    cs_global += MTIME_() - start;

    return 4.0 * residual / (scale*scale*(ngrid-1)*(ngrid-1));
  }
  else
    return 999999.0;
}

double **Make2d(ncols, nrows)
     long ncols, nrows;
/*
  C does not allow variably dimensioned multi-dim. arrays.
  Allocate memory to simulate 2-D array with array of pointers.
*/
{
  double **temp, **array;

  array = temp = (double **) malloc((unsigned) (ncols*sizeof(double *)));

  if (!temp)
    Error("Make2d: malloc 1 failed",(long) (ncols*sizeof(double *)));

  while (ncols--) {
    if (!(*temp++ = (double *) malloc((unsigned) (nrows*sizeof(double)))))
      Error("Make2d: malloc 2 failed", (long) (nrows*sizeof(double)));
  }
  return array;
}

void Interpolate(old_grid, old_ncols, old_nrows, old_col_low, old_row_low,
         new_grid, new_ncols, new_nrows, new_col_low, new_row_low)
     double **old_grid, **new_grid;
     long old_ncols, old_nrows, old_col_low, old_row_low;
     long new_ncols, new_nrows, new_col_low, new_row_low;
/*
  Interpolate down from a solution on grid ncols/2 to an initial
  guess on a grid dimension ncols.
  The parallel version is complicated by the fact that rows may
  have shifted between processes on scaling up. Ugh.
*/
{
  long col_shift = 2*old_col_low - new_col_low;
  long row_shift = 2*old_row_low - new_row_low;
  long start = MTIME_();
  long i, i1, i2, j, j1, j2;

  for (i=0; i<=old_ncols; i++)
    for (j=0; j<=old_nrows; j++) {
      i1 = 2*i+col_shift; j1 = 2*j+row_shift;
      if ( (i1>=0) && (i1<=new_ncols) && (j1>=0) && (j1<=new_nrows) )
    new_grid[i1][j1] = old_grid[i][j];
    }

  for (i=1; i<=old_ncols; i++)
    for (j=1; j<=old_nrows; j++) {
      i2 = 2*i-1+col_shift; j2 = 2*j-1+row_shift;
      if ( (i2>=0) && (i2<=new_ncols) && (j2>=0) && (j2<=new_nrows) )
    new_grid[i2][j2] = 0.25 * (old_grid[i  ][j  ] + old_grid[i-1][j-1] +
                   old_grid[i-1][j  ] + old_grid[i  ][j-1]);
    }
  
  BoundaryConditions(new_grid, new_ncols, new_nrows);
  Exchange(new_grid, new_ncols, new_nrows);

  for (i=1; i<new_ncols; i++)
    for (j=1+((i+new_col_low+new_row_low)%2); j<new_nrows; j+=2)
      new_grid[i][j] = 0.25 * (new_grid[i+1][j] + 
                   new_grid[i-1][j] + 
                   new_grid[i][j-1] + 
                   new_grid[i][j+1]);

  Exchange(new_grid, new_ncols, new_nrows);

  cs_interpolate += MTIME_() - start;
  
}

void Solve(grid, ncols, nrows, ngrid, niter, nprint, thresh)
     double **grid, thresh;
     long ncols, nrows, ngrid, niter, nprint;
/*
  Apply iterative procedure to solve current grid.

  Parallel version only does global sums every nsums iterations.
*/
{
  long iter, do_sums, nsums;
  double residual, error;
#ifdef PLOT
  long nplots = TCG_MAX(5, ngrid/10);  /* Only plot when have changed a lot */
#endif

  nsums = TCG_MIN(10, nprint);     /* Need sums whenever we print */
  if (nprint%nsums)
    nprint= nprint + nsums - (nprint%nsums); /* Make nprint a multiple 
                            of nsums */
  for (iter=0; iter<niter; iter++) {

    /* For efficiency only do global sums every 10 iters */
    do_sums = !(iter%nsums);

    /* Actually do the work */
    residual = Operate(grid, ncols, nrows, ngrid, do_sums);

    /* Print the results every now and again or if converged */
    if ((NODEID_()==0) && ((iter%nprint == 0) || (residual < thresh))) {
      (void) printf("ngrid=%d iter=%d residual=%f\n",
            ngrid, iter+1, residual);
      (void) fflush(stdout);
    }
    
    /* Are we converged ? */
    if (do_sums && (residual < thresh)) {
      if (NODEID_() == 0) {
    (void) printf("Converged!\n");
    (void) fflush(stdout);
      }
      break;
    }

#ifdef PLOT
    if (!(iter%nplots))
      PlotGrid(grid, ngrid, ncols, nrows);
#endif

  }
  
  /* Have to do a final exchange to get edges correct */
  Exchange(grid, ncols, nrows);

  error = GridError(grid, ncols, nrows, ngrid);
  if (NODEID_() == 0) {
    (void) printf("Mean abs. error to exact soln. = %f, ngrid=%d\n\n",
          error, ngrid);
    (void) fflush(stdout);
  }

#ifdef PLOT
  PlotGrid(grid, ngrid, ncols, nrows);
#endif
}

void ParseArguments(argc, argv, pngrid, pniter, pnprint, pnlevel, pthresh)
     long argc, *pngrid, *pniter, *pnprint, *pnlevel;
     double *pthresh;
     char **argv;
{
  argc--;
  argv++;
  while (argc--) {
    if (strcmp(*argv, "-niter") == 0)
      *pniter = atoi(*++argv);
    else if (strcmp(*argv, "-ngrid") == 0)
      *pngrid = atoi(*++argv);
    else if (strcmp(*argv, "-nprint") == 0)
      *pnprint = atoi(*++argv);
    else if (strcmp(*argv, "-nlevel") == 0)
      *pnlevel = atoi(*++argv);
    else if (strcmp(*argv, "-thresh") == 0)
      (void) sscanf(*++argv, "%lf", pthresh);
#ifdef PLOT
    else if (strcmp(*argv, "-plot") == 0) {
      argv++;
      if (strcmp(*argv,"value") == 0)
    plot_type = PLOT_VALUE;
      else if (strcmp(*argv, "error") == 0)
    plot_type = PLOT_ERROR;
      else if (strcmp(*argv, "residual") == 0)
    plot_type = PLOT_RESIDUAL;
      else
    Error("Unknown plot type - use error|value|residual",(long) -1);
    }    
#endif
    else if (strcmp(*argv, "-help") == 0) {
      (void) fprintf(stderr,"gridtest [-ngrid #] [-nprint #] [-niter #]\n");
      (void) fprintf(stderr,"         [-thresh #] [-nlevel #] [-help]\n");
#ifdef PLOT
      (void) fprintf(stderr,"         [-plot value|error|residual]\n");
#endif
      PEND_();
      exit(1);
    }
    argv++; argc--;
  }
}

void Factor(N, pn, pm)
  long N, *pn, *pm;
/*
  Factor N into two integers that are as close together as possible
*/
{
  long n = ceil(sqrt((double) N));
  long m = floor(sqrt((double) N));
 
  while (n*m != N) {
    if (n*m < N)
      n++;
    else
      m--;
  }
  *pn = n; *pm = m;
}

void Partition(ngrid, pncols, pnrows)
    long ngrid, *pncols, *pnrows;
/*
  The square mesh actually has (ngrid-1)*(ngrid-1) interior points.
  Divide these up equitably between all the processors.
*/
{
  long col_chunk = (ngrid-1) / ncol_P;
  long col_extra = (ngrid-1) - col_chunk*ncol_P;
  long row_chunk = (ngrid-1) / nrow_P;
  long row_extra = (ngrid-1) - row_chunk*nrow_P;

  col_low = my_col_P*col_chunk + TCG_MIN(my_col_P,col_extra); /* Col of top LHS */

  row_low = my_row_P*row_chunk + TCG_MIN(my_row_P,row_extra); /* Row of top LHS */

  *pncols = col_chunk + ( (my_col_P<col_extra) ? 1 : 0 ) + 1;
  *pnrows = row_chunk + ( (my_row_P<row_extra) ? 1 : 0 ) + 1;
}

int main(argc, argv)
     int argc;
     char **argv;
/*
  Grid test code. Solve Laplaces eqn. on a square grid subject
  to b.c.s on the boundary.  Use 5 point discretization of the
  operator and a heirarchy of grids with red/black gauss seidel
  w-relaxation (omega=0.94 as 1.0 diverges).

  command arguments:

  -ngrid  # = dimension of largest grid (default = 128)
  -niter  # = max. no. of interations   (default = 1000)
  -nprint # = print residual every nprint iterations (default = 100)
  -nlevel # = no. of grid levels        (default = 4)
  -thresh # = convergence criterion on the residual (default = 0.1)
  -help     = print usage and exit with error
  -plot value|error|residual = if compiled with -DPLOT generate plots
                               displaying 
  
*/
{
  double **grid, thresh=0.1,**new_grid;
  long niter=1000, nprint=100, ngrid, ncols=0, nrows=0, maxgrid=128;
  long nlevel=4, level, ngrid1, on=0;
  long old_ncols, old_nrows, old_col_low, old_row_low;
  long start = MTIME_();

  tcg_pbegin(argc, argv);            /* Initialize parallel environment */
  SETDBG_(&on);

  P = NNODES_();
  Factor(P, &ncol_P, &nrow_P);    /* Arrange processes into a grid   */
  my_col_P = NODEID_()/nrow_P;
  my_row_P = NODEID_() - my_col_P*nrow_P; 

  /* Who are my neighbors on the process grid */

  north = (my_row_P > 0)          ? NODEID_()-1      : -1;
  south = (my_row_P < (nrow_P-1)) ? NODEID_()+1      : -1;
  east  = (my_col_P < (ncol_P-1)) ? NODEID_()+nrow_P : -1;
  west  = (my_col_P > 0)          ? NODEID_()-nrow_P : -1;

  ParseArguments(argc, argv, 
         &maxgrid, &niter, &nprint, &nlevel, &thresh);

  ngrid1 = TCG_MAX(ncol_P,maxgrid>>(nlevel-1));
  ngrid1 = TCG_MAX(nrow_P,ngrid1);                 /* Size of first grid */
  maxgrid = ngrid1<<(nlevel-1);                /* Actual size of final grid */

  if (!(buffer = (double *) malloc((unsigned) (sizeof(double)*(maxgrid+1)))))
    Error("failed to allocate comms buffer", (long) (maxgrid+1));

#ifdef PLOT
  OpenPlotFile(maxgrid);
#endif

  /* Loop from coarse to fine grids */

  for (level=0; level<nlevel; level++) {
    
    ngrid = ngrid1<<level;        /* Grid dimension */
    scale = (HI - LO) / ngrid;

    /* Partition grid between processors */

    old_ncols = ncols; old_nrows = nrows;        /* Save info for interp. */
    old_col_low = col_low; old_row_low = row_low;
    Partition(ngrid, &ncols, &nrows);

    if (level == 0) {
      grid = Make2d(ncols+1, nrows+1);    /* Allocate first grid */
      Initialize(grid, ncols, nrows);     /* Intitial guess */
    }
    else {
      new_grid = Make2d(ncols+1, nrows+1);
      Interpolate(grid, old_ncols, old_nrows, old_col_low, old_row_low,
          new_grid, ncols, nrows, col_low, row_low);
      grid = new_grid;
    }

    Solve(grid, ncols, nrows, ngrid, niter, nprint, thresh);
  }

#ifdef PLOT
  ClosePlotFile();
#endif

  cs_total = MTIME_() - start;
  if (NODEID_() == 0) {
    (void) printf("\n  Plot    Exchange   Global-sum   Interpolate    Total ");
    (void) printf("\n ------   --------   ----------   -----------   -------");
    (void) printf("\n %6d    %6d     %6d       %6d      %7d\n",
          cs_plot, cs_exchange, cs_global,
          cs_interpolate, cs_total);
  }
          
  PEND_();                               /* Terminate parallel env. */
  return 0;
}

