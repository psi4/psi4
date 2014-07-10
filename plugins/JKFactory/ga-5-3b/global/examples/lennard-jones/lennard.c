#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_UNISTD_H
#   include <unistd.h>
#endif
#if HAVE_MATH_H
#   include <math.h>
#endif

#include "ga.h"
#include "macdecls.h"
#include "mp3.h"

#define NATOMS            256
#define BLOCK_SIZE        64
#define DENSITY           1.05
#define TEMPERATURE       2.0
#define TIMESTEP          0.004
#define RESCALE_STEPS     10
#define EQUILIBRIUM_TIME  10.0
#define SIMULATION_TIME   20.0

#define STOP_ITERATING  0
#define DEBUG           0
#define PRINT_LEVEL_1   1 /**< dumps result: level 1 */
#define PRINT_LEVEL_2   1 /**< dumps result: level 2..more results */
#define WRITE_TO_FILE   0 /**< dumps the coordinates in molden viz format */
#define NDIM            3 /**< always 3-d */
#define MAX_PROC        256
#define MAX_BLOCKS      256
#define SAFELIMIT       3

#ifdef MPI
#   define CLOCK_ MPI_Wtime
#else
#   define CLOCK_ tcg_time
#endif

/* #define GA_ABS(a) ( ((a) > 0) ? (a) : -(a) ) */

typedef struct {
  int x;
  int y;
} topo_t;

/* box specifications */
typedef struct {
  double length;
  double width;
  double height;
} boxSpec_t;

int gMemHandle;

static int gMe, gNproc;
static int g_X; /* Coordinate array */
static int g_G; /* Gradient Array */
static int g_V; /* Velocity Array */
static int g_T; /* Atomic Task: 1-d array whose size = 1 */
static int gBlockSize, gNblocks;
static double gDensity, gDesiredTemperature, gTimeStep, gComputeTime=0.0;
static topo_t btopo[MAX_BLOCKS];
static boxSpec_t box;

static int measurementStep;
static double temperature, temperatureSum, temperatureSqdSum;
static double pressure, pressureSum, pressureSqdSum;
static double potentialEnergySum, potentialEnergySqdSum;
static double kineticEnergySum, totalEnergySum, totalEnergySqdSum;


/**
 * To get the next task id. It is an atomic operation.
 */
static int nxtask(void) {

  int subscript = 0;
  return NGA_Read_inc(g_T, &subscript, 1);
}

/**
 * Block Topology (of Force Matrix): 
 * Say for example: If there are 4 block and 100 atoms, the size of 
 * the force matrix is 100x100 and each block size is 50x50. 
 * 
 *  -----------
 * |     |     |
 * | 0,0 | 0,1 |
 *  -----------
 * | 1,0 | 1,1 |
 * |     |     |
 *  -----------
 */
void LJ_Setup(int natoms, double **x_i, double **x_j, double **grad) {
  
  int i, j, k=0, n;
  MA_AccessIndex maindex;

  n = natoms/gBlockSize;
   
  /* block topology */
  for(i=0; i<n; i++) 
    for(j=0; j<n; j++, k++) {
      btopo[k].x = i; 
      btopo[k].y = j;
    }
  
  /* create task array */
  n = 1;
  g_T = NGA_Create(C_INT, 1, &n,"Atomic Task", NULL);
  
  
  /* Allocate memory to the buffers used */
  i = gBlockSize * NDIM * 2; /* for gX_i and gX_j used in computeFG() */
  j = natoms * NDIM;         /* for gGrad in computeFG() */
  n = i + j + SAFELIMIT;                 /* total memory required */
  if(MA_push_get(C_DBL, n, "GA LJ bufs", (void *)&gMemHandle, &maindex))
    MA_get_pointer(gMemHandle, x_i);
  else GA_Error("ma_alloc_get failed",n);
  
  *x_j  = *x_i + i/2 + 1;
  *grad = *x_j + i/2 + 1;
}

/**
 * Initial assumption is verified here 
 */
void check(int natoms) {

  int temp = natoms/gBlockSize;
  gNblocks = temp * temp;

  if(gBlockSize <= 0) GA_Error("Invalid Block Size", gBlockSize);
  
  if(natoms%gBlockSize) {
    GA_Error("CHECK: # of atoms should be a multiple of block size. Chose a different block size.", 0L);
  } 
  
  if(gNblocks > MAX_BLOCKS)
    GA_Error("Number of blocks is greater than MAX_BLOCKS(512): Solution is either increase the defined MAX_BLOCKS or increase your block size", 0L); 
  
  if(gNblocks < gNproc)
    GA_Error("Number of blocks should be greater than number of processors",0L);
}




/**
 * The entire Force matrix is sub-divided into many blocks.
 * This function gets the next block to be computed by a 
 * process. Once a process finishes a block it gets a new 
 * block to be computed.
 */
void getBlock(int taskId, int size, double *x_i, double *x_j) { 
  int lo, hi;
  
#if DEBUG
  printf("%d: new task = %d: topo: %d,%d\n", gMe, taskId, 
     btopo[taskId].x,  btopo[taskId].y);
#endif
  
  /** get the coordinates of the atoms in the corresponding rows in the 
      block */
  lo = btopo[taskId].x * size;
  hi = lo + size -1;
  NGA_Get(g_X, &lo, &hi, x_i, &hi);

  /** get the coordinates of the atoms in the corresponding columns in 
      the block */
  lo = btopo[taskId].y * size;
  hi = lo + size -1;
  NGA_Get(g_X, &lo, &hi, x_j, &hi);
}


/**
 * LJ Function Gradient Computation.
 */
void LJ_FG(int taskId, double *x_i, double *x_j, double *f, 
       double *grad) {
  int b_x, b_y; /* block topology */
  int i, j, start_x, start_y, tempA, tempB;
  int start_i=0, end_i=0, start_j=0, *end_j=NULL;
  int sign_x, sign_y, sign_z;
  double xx, yy, zz, rij, temp,r2,r6,r12, xtmp, ytmp, ztmp;
  
  b_x = btopo[taskId].x;
  b_y = btopo[taskId].y;
  start_x = gBlockSize * NDIM * btopo[taskId].x;
  start_y = gBlockSize * NDIM * btopo[taskId].y;
  
  if(b_x == b_y) { /* computer lower triagular matrix only */
    start_i = 1; end_i = gBlockSize;
    start_j = 0; end_j = &i;
  }
  else if(b_x > b_y) { /* compute right half of the block */
    start_i = 0; end_i = gBlockSize;
    start_j = gBlockSize/2; end_j = &gBlockSize;    
  }
  else if(b_x < b_y) { /* compute upper half of the block */
    start_i = 0; end_i = gBlockSize/2;
    start_j = 0; end_j = &gBlockSize;    
  }
  
  /* Calculating Force exerted on 'i' by 'j' */
  for(i=start_i; i< end_i; i++) {
    xtmp = x_i[NDIM*i];
    ytmp = x_i[NDIM*i+1];
    ztmp = x_i[NDIM*i+2];
    for(j=start_j; j< *end_j; j++) {
      
      xx = xtmp - x_j[NDIM*j];
      yy = ytmp - x_j[NDIM*j+1];
      zz = ztmp - x_j[NDIM*j+2];

      sign_x = xx > 0.0 ? 1 : -1;
      sign_y = yy > 0.0 ? 1 : -1;
      sign_z = zz > 0.0 ? 1 : -1;
      /** 
       * Using Nearest Image Approximation in computing the distance between 
       * any 2 atoms: If any component of (Ri-Rj) is greater than L/2, then 
       * there is an image particle located closer which exerts a larger force 
       * (because we use periodic boundary conditions).
       */
      if(xx*sign_x > box.length/2) xx -= sign_x*box.length; 
      if(yy*sign_y > box.width/2)  yy -= sign_y*box.width ; 
      if(zz*sign_z > box.height/2) zz -= sign_z*box.height; 
      rij = xx*xx + yy*yy + zz*zz;
#if DEBUG
      if(rij <= 0.0) GA_Error("Divide by Zero Error\n", 0L);
#endif
      r2 = 1.0/rij;
      r6 = r2*r2*r2;
      r12 = r6*r6;
      *f = *f + 4.0*(r12-r6); /* function */
      temp = r2 * (48.0 * r12 - 24.0 * r6);
      tempA = start_x + NDIM*i;   tempB = start_y + NDIM*j;
      grad[tempA]   += xx*temp; /* gradient */
      grad[tempA+1] += yy*temp;
      grad[tempA+2] += zz*temp;
      grad[tempB]   -= xx*temp;
      grad[tempB+1] -= yy*temp;
      grad[tempB+2] -= zz*temp;
    }
  }
}


/**
 * Compute Function and Gradient
 */
void computeFG(double *force, int natoms,  double *x_i, double *x_j, 
           double *grad) {
  int taskId, size, lo, hi, i;
  double tt;
  
  taskId = gMe;
  size = gBlockSize * NDIM;

  for(i=0; i<natoms*NDIM; i++)
    grad[i] = 0.0;
  
  while(taskId < gNblocks) {
    getBlock(taskId, size, x_i, x_j);
    tt = CLOCK_();
    LJ_FG(taskId, x_i, x_j, force, grad);
    gComputeTime += CLOCK_() - tt;
    taskId = nxtask(); /* atomic */
  }

  /* function */
  GA_Dgop(force, 1, "+");
  
  /* gradient */
  GA_Dgop(grad, natoms*NDIM, "+");
  GA_Sync();
  NGA_Distribution(g_G, gMe, &lo, &hi);
  NGA_Put(g_G, &lo, &hi, &grad[lo], &hi);
  GA_Sync();
}


/**
 * Process the command line arguments
 */
void
commandLine(int argc, char **argv) {
 
#if 0
  int n;
#endif
 
  /* default options */
  gBlockSize   = BLOCK_SIZE; /* size of the sub-block (force matrix) */
  gDensity     = DENSITY;
  gDesiredTemperature = TEMPERATURE;
  gTimeStep    = TIMESTEP;
  
#if 0
  if(argc > 5) {
    printf("Argc = %d\n", argc);
    while((n = getopt(argc, argv, "b:")) != EOF)
      switch (n) {
      case 'b':
    gBlockSize = atoi(optarg); 
      }
  }
#endif
}

void rescaleVelocities (int totalAtoms) {
  int i, p, lo, hi, ld, natms;
  double vSqdSum = 0.0, scale;
  double *v;
  
  GA_Sync();
  
  NGA_Distribution(g_V, gMe, &lo, &hi);
  NGA_Access(g_V, &lo, &hi, &v, &ld);
  
  natms = (hi-lo+1)/NDIM;
  for (p = 0; p < natms; p++)
    for (i = 0; i < NDIM; i++)
      vSqdSum += v[p*NDIM + i] * v[p*NDIM + i];
  
  GA_Dgop(&vSqdSum, 1, "+");
  
  scale = NDIM * (totalAtoms) * gDesiredTemperature / vSqdSum;
  scale = sqrt(scale);
  for (p = 0; p < natms; p++)
    for (i = 0; i < NDIM; i++)
      v[p*NDIM + i] *= scale;
  
  GA_Sync();
}

double gaussianDistribution() {

  static int available = 0;
  static double savedDeviate;
  double r[2], rSqd, factor;
  int i;
  
  if (available) {
    available = 0;
    return savedDeviate;
  }
  do {
    rSqd = 0.0;
    for (i = 0; i < 2; i++) {
      r[i] = (2.0 * (double) rand()) / RAND_MAX - 1.0;
      rSqd += r[i] * r[i];
    }
  } while (rSqd >= 1.0 || rSqd == 0.0);
  factor = sqrt(-2.0 * log(rSqd) / rSqd);
  savedDeviate = r[0] * factor;
  available = 1;
  return r[1] * factor;
  
}

void LJ_Initialize(int natoms) {

  double L;  
  
  /* compute the length, width and height of the box */
  L = pow(natoms/gDensity, 1.0/3.0);
  box.length = box.width = box.height = L;
  
#if PRINT_LEVEL_1
  if(gMe == 0) {
    printf("\n");
    printf(" =========================================================\n\n");
    printf("              Molecular Dynamics Simulation               \n\n");
    printf("                           of                             \n\n");
    printf("                   Lennard Jones System                   \n\n");
    printf(" =========================================================\n\n\n");
    printf(" Number of Atoms/Particles   =   %d\n",    natoms);
    printf(" Block Size                  =   %d\n",    gBlockSize);
    printf(" Number of Blocks            =   %d\n",    gNblocks);
    printf(" Density                     =   %f\n",    gDensity);
    printf(" Temperature                 =   %f\n",    gDesiredTemperature);
    printf(" Time Step                   =   %f\n\n",  gTimeStep);
    printf(" Box Specifications:\n");
    printf("    Size of the Cube (Box)   =   %f\n",     box.length);
    printf("    System Volume            =   %f\n\n\n", pow(L, 3.0));
  }
#endif
  
  if(gMe == 0) {
    
    int c, i, j, k, m, n, p;
    double b, xSum[3] = {0.0, 0.0, 0.0};
    double rFCC[4][3] = {{0.0, 0.0, 0.0}, {0.0, 0.5, 0.5},
             {0.5, 0.0, 0.5}, {0.5, 0.5, 0.0}};
    double rCell[3];
    double *x;
    int lo, hi, handle;  
    MA_AccessIndex maindex;
    
    n = NDIM * natoms + 1;
    if(MA_push_get(C_DBL, n, "GA LJ_Init bufs", (void *)&handle, &maindex))
      MA_get_pointer(handle, &x);
    else GA_Error("ma_alloc_get failed",n);
    
    /* Use face centered cubic (FCC) lattice for initial positions.
       Find number of unit cells (c) needed to place all atoms */
    for (c = 1; ; c++)
      if (4*c*c*c >= natoms)
    break;
    b = L / c;            /* side of unit cell */
    p = 0;            /* atoms placed so far */
    
    for (i = 0; i < c; i++) {
      rCell[0] = i * b;
      for (j = 0; j < c; j++) {
    rCell[1] = j * b;
    for (k = 0; k < c; k++) {
      rCell[2] = k * b;
      for (m = 0; m < 4; m++)    /* 4 particles in cell */
        if (p < natoms) {
          for (n = 0; n < NDIM; n++)  /* 3-dimensions - x, y, z */
        x[p*NDIM + n] = rCell[n] + b * rFCC[m][n];
          ++p;
        }
    }
      }
    }
    
    lo = 0;
    hi = natoms*NDIM-1; 
    NGA_Put (g_X, &lo, &hi, x, &hi);

    /* Random Gaussian distribution of initial velocities */ 
    for(i=0; i<natoms; i++) 
      for(j=0; j<NDIM; j++) 
    xSum[j] += x[i*NDIM + j] = gaussianDistribution();
    
    /* with zero total momentum */
    for(i=0; i<natoms; i++) 
      for(j=0; j<NDIM; j++) 
    x[i*NDIM + j] -= xSum[j]/natoms;

    NGA_Put (g_V, &lo, &hi, x, &hi); /* velocity array */
    
    if(!MA_pop_stack(handle)) GA_Error("LJ_Init:MA_pop_stack failed",0);  
  }
    
  /* rescale to desired temperature */
  rescaleVelocities(natoms);

#if  DEBUG
  GA_Print(g_X);
#endif
    
  GA_Sync(); /* wait till process 0 is done */
}

void initializeProperties () {
  measurementStep = 0;
  temperatureSum = temperatureSqdSum = 0;
  pressureSum = pressureSqdSum = 0;
  potentialEnergySum = potentialEnergySqdSum = 0;
  kineticEnergySum = totalEnergySqdSum = totalEnergySum = 0;
}

void computeProperties(int natoms, double potentialEnergy, 
               double *totalEnergy) {
  
  int i, p, lo, hi, ld, natms;
  double kineticEnergy = 0.0, vir = 0.0;
  double *v, *x, *a;
  
  GA_Sync();
  
  NGA_Distribution(g_V, gMe, &lo, &hi);
  NGA_Access(g_V, &lo, &hi, &v, &ld);
  NGA_Access(g_X, &lo, &hi, &x, &ld);
  NGA_Access(g_G, &lo, &hi, &a, &ld);
  
  natms = (hi-lo+1)/NDIM;
  for (p = 0; p < natms; p++)
    for (i = 0; i < NDIM; i++) {
      kineticEnergy += 0.5 * v[p*NDIM + i] * v[p*NDIM + i];
      vir += x[p*NDIM + i] * a[p*NDIM + i];
    }

  GA_Dgop(&kineticEnergy, 1, "+");
  GA_Dgop(&vir, 1, "+");
  (*totalEnergy) = kineticEnergy + potentialEnergy;

  temperature = 2.0 * kineticEnergy / 3.0 / natoms;
  pressure = temperature + vir / 3.0 / natoms;
  pressure *= gDensity;

  ++measurementStep;
  temperatureSum += temperature;
  temperatureSqdSum += temperature * temperature;
  pressureSum += pressure;
  pressureSqdSum += pressure * pressure;
  potentialEnergySum += potentialEnergy;
  potentialEnergySqdSum += potentialEnergy * potentialEnergy;

  /* for checking the fluctuation */
  kineticEnergySum += kineticEnergy;
  totalEnergySum   += (*totalEnergy);
  totalEnergySqdSum += (*totalEnergy)*(*totalEnergy);
}

void printProperties (int natoms) {
  double average, stdDev;

  average = temperatureSum / measurementStep;
  stdDev = temperatureSqdSum / measurementStep;
  stdDev = sqrt(stdDev - average * average);
  printf("\nTemperature = %f  +-  %f\n", average, stdDev);
    
  average = pressureSum / measurementStep;
  stdDev = pressureSqdSum / measurementStep;
  stdDev = sqrt(stdDev - average * average);
  printf("Pressure      = %f  +-  %f\n", average, stdDev);
    
  average = potentialEnergySum / measurementStep;
  stdDev = potentialEnergySqdSum / measurementStep;
  stdDev = sqrt(stdDev - average * average);
  printf("Potential Energy per particle = %f  +-  %f\n",
     average / natoms, stdDev / natoms);
  
  stdDev  = (totalEnergySqdSum - totalEnergySum);
  printf("Energy Fluctuation            = %f\n\n", 
     sqrt(stdDev)/kineticEnergySum);
}


void LJ_Update() {
  int i, p, lo, hi, ld, natms;
  double *x, *a, *v, *b;
  
  GA_Sync();
  
  NGA_Distribution(g_X, gMe, &lo, &hi);
  NGA_Access(g_X, &lo, &hi, &x, &ld);
  NGA_Access(g_G, &lo, &hi, &a, &ld); /* X and G has the same distribution */
  NGA_Access(g_V, &lo, &hi, &v, &ld); 
  
  natms = (hi-lo+1)/NDIM;
  b = (double *)&box;
  for (p = 0; p < natms; p++)
    for (i = 0; i < NDIM; i++) {
      x[p*NDIM + i] += v[p*NDIM + i] * gTimeStep + 
    0.5 * a[p*NDIM + i] * gTimeStep * gTimeStep;
      /* impose periodic boundary conditions */
      if (x[p*NDIM+i] < 0.0)   x[p*NDIM+i] += b[i];
      if (x[p*NDIM+i] >= b[i]) x[p*NDIM+i] -= b[i];
      v[p*NDIM+i] += 0.5 * a[p*NDIM+i] * gTimeStep;
    }
  GA_Sync();
}

/**
 * Using Velocity-Verlet Algorithm.
 */
void LJ_Update_Velocity() {
  int i, p, lo, hi, ld, natms;
  double *a, *v;
  
  GA_Sync();
  
  NGA_Distribution(g_V, gMe, &lo, &hi);
  NGA_Access(g_V, &lo, &hi, &v, &ld); 
  NGA_Access(g_G, &lo, &hi, &a, &ld); 
  
  natms = (hi-lo+1)/NDIM;
  for (p = 0; p < natms; p++)
    for (i = 0; i < NDIM; i++) 
      v[p*NDIM+i] += 0.5 * a[p*NDIM+i] * gTimeStep;

  GA_Sync();
}

void initializeTaskArray() {
  int lo, hi, n;
  if(gMe ==0) {
    /* Initialize the task array */
    lo = 0; hi = 0; /* only one element */
    n = gNproc; /* initial value of the task counter/id */
    NGA_Put (g_T, &lo, &hi, &n, &hi);
  }
  GA_Sync();
}

void solveOneTimeStep(double *potentialEnergy, int natoms, 
              double *x_i, double *x_j, double *grad) {
  *potentialEnergy = 0.0;

  LJ_Update();    /* Update the coordinates */
  initializeTaskArray();
  computeFG(potentialEnergy, natoms, x_i, x_j, grad);
  LJ_Update_Velocity();
}

void writeToFile(int natoms) {
  static int first_time = 1;
  static FILE *gOutfile;
  int n=0, i=0, lo, hi;
  double *p_data;
  
  if(gMe == 0) {
    
    if(first_time) {
      gOutfile = fopen("output.dat", "w"); 
      first_time = 0;
    }
    else if(natoms == -1) { fclose(gOutfile); return; }

    p_data = (double *) malloc(sizeof(double) * natoms * NDIM);
    lo = 0; hi = natoms * NDIM -1;
    NGA_Get(g_X, &lo, &hi, p_data, &hi);
    
    /* in molden format */
    fprintf(gOutfile, "%d\n\n", natoms) ; /* 2 new lines needed */
    do {
      fprintf(gOutfile, "%s %f %f %f\n", "XX",
          p_data[i], p_data[i+1], p_data[i+2]);
      i+=NDIM;
    }while(++n < natoms);
    
    free(p_data);
  }
}

/**
 * Lennard Jones Molecular Dynamics: starts here...
 */
void LJ_Solve(int natoms) {
  double execTime, potentialEnergy, totalEnergy = 0.0, time=0.0;
  double *x_i, *x_j, *grad, ip;
  int s=0;

  /* Initial Setup */
  LJ_Setup(natoms, &x_i, &x_j, &grad);
  initializeTaskArray();
#if WRITE_TO_FILE
  writeToFile(natoms);
#endif
 
  execTime = CLOCK_();
  potentialEnergy = 0.0;
  computeFG(&potentialEnergy, natoms, x_i, x_j, grad);
  initializeProperties();
  computeProperties(natoms, potentialEnergy, &totalEnergy);
#if WRITE_TO_FILE
  writeToFile(natoms);
#endif
  
  if(gMe == 0) printf("Equilibrium Steps:\n");
  if(gMe == 0) printf("Time = %.3f \tEnergy = %f\n", time, totalEnergy);
  
  /* Equilibrium Steps */
  while(time < EQUILIBRIUM_TIME) {
    solveOneTimeStep(&potentialEnergy, natoms, x_i, x_j, grad);
    computeProperties(natoms, potentialEnergy, &totalEnergy);
#if WRITE_TO_FILE
    writeToFile(natoms);
#endif
    time += gTimeStep;
    if(gMe == 0 && modf(time,&ip) < gTimeStep) 
      printf("Time = %.2f \tEnergy = %f\n", ip, totalEnergy);
    if(++s >= RESCALE_STEPS)  { rescaleVelocities(natoms); s = 0; }    
  }
  
  /* Actual Production Steps */
  if(gMe == 0) printf("Production Steps:\n");
  initializeProperties();
  while(time < EQUILIBRIUM_TIME + SIMULATION_TIME) {
    solveOneTimeStep(&potentialEnergy, natoms, x_i, x_j, grad);
    computeProperties(natoms, potentialEnergy, &totalEnergy);
#if WRITE_TO_FILE
    writeToFile(natoms);
#endif
    time += gTimeStep;
    if(gMe == 0 && modf(time,&ip) < gTimeStep) 
      printf("Time = %.2f \tEnergy = %f\n", ip, totalEnergy);
  }
  
#if WRITE_TO_FILE
  writeToFile(-1);
#endif
  if(gMe == 0) printProperties(natoms);
#if PRINT_LEVEL_2
  if(gMe == 0) {
    execTime = CLOCK_()-execTime;
    printf("%d: Total Elapsed Time  = %f\n", gMe, execTime);
    printf("%d: Computation Time    = %f\n", gMe, gComputeTime);
    printf("%d: Percentage Overhead = %f\n\n", gMe, 
       100*(execTime-gComputeTime)/execTime);
  }
#endif
  if(!MA_pop_stack(gMemHandle)) GA_Error("LJ_Init:MA_pop_stack failed",0); 
}


/**
 * main(int argc, char **argv)
 */
int main(int argc, char **argv) {

  int heap=4000000, stack=4000000;
  int natoms=NATOMS;
  int dims[NDIM];

  MP_INIT(argc, argv);

  /**
   * Initialize Global Arrays.
   */
  GA_Initialize_args(&argc, &argv);
  gMe    = GA_Nodeid();
  gNproc = GA_Nnodes();
  heap  /= gNproc;      
  stack /= gNproc;
  if(! MA_init(C_DBL, stack, heap)) 
    GA_Error("MA_init failed",stack+heap);  /* initialize memory allocator*/

  /** 
   * Parse the command line and check the initial conditions (and assumptions)
   */
  commandLine(argc, argv);
  check(natoms);

  /**
   * Create coordinate(x), gradient(g) and velocity(v) vectors
   */
  dims[0] = natoms*NDIM;
  g_X = NGA_Create(C_DBL, 1, dims, "Coordinate Array - X", NULL);
  g_G = GA_Duplicate (g_X, "Gradient");     
  g_V = GA_Duplicate (g_X, "Velocity");     
    
  /** 
   * Molecular Dynamics of the Lennard Jones(lj) clusters.
   * (molecular conformation problem).
   */
  LJ_Initialize(natoms);
  LJ_Solve(natoms);
  
  /**
   * Deaalocate the arrays and free the resources 
   */
  GA_Destroy(g_X);  GA_Destroy(g_G);  GA_Destroy(g_V); GA_Destroy(g_T);

  /**
   * Termination signal to release the resources, etc.
   */
  GA_Terminate();
  
#ifdef MPI
  MPI_Finalize();
#else
  tcg_pend(); 
#endif
  
  return 0;
}
