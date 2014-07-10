/*                                                                                                       
 * Test Program for GA                                                                                   
 * This is to test GA_Dgop (is a collective operation)                                                 
 * 
 * 
 */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include"mpi.h"
#include"ga.h"
#include"macdecls.h"

addition_operator(int rank, int nprocs, int n)
{
  int i;
  double x[n], temp[n];

  for(i=0; i<n; i++) x[i]=rand()%5;
  for(i=0; i<n; i++) temp[i]=x[i];
  
  if(rank==1)
    for(i=0; i<n; i++) printf("%f \n", x[i]);
  
  
  GA_Dgop(x, n, "+");
  printf("\n");
  
  if(rank==0)
    {
      for(i=0; i<n; i++) printf("%f :: %f\n", x[i], temp[i]);
      
      for(i=0; i<n; i++)
	if(x[i]!=temp[i]*sizeof(int)) printf(" GA ERROR: \n");
    }
}

multiply_operator(int rank, int nprocs, int n)
{
  int i;
  double x[n], temp[n];

  for(i=0; i<n; i++) x[i]=rand()%5;
  for(i=0; i<n; i++) temp[i]=x[i];
  
  if(rank==1)
    for(i=0; i<n; i++) printf("%f \n", x[i]);
  
  
  GA_Dgop(x, n,"*");
  printf("\n");
  
  if(rank==0)
    {
      //      for(i=0; i<n; i++) printf("%f :: %f\n", x[i], temp[i]);
      //      for(i=0; i<n; i++) printf("%f :: %f \n", x[i], (int)pow((int)temp[i],(int)4));
      
      for(i=0; i<n; i++) printf("%f :: %f \n", x[i], pow(temp[i],4));

      //      for(i=0; i<n; i++)
      //if(x[i]!=(int)pow((int)temp[i],(int)4) printf(" GA ERROR: \n");

      for(i=0; i<n; i++)
	if(x[i]!=pow(temp[i],4)) printf(" GA ERROR: \n");
    }
}

max_operator(int rank, int nprocs, int n)
{
  int i;
  double x[n], temp[n];

  for(i=0; i<n; i++) x[i]=rand()%5;
  for(i=0; i<n; i++) temp[i]=x[i];
  
  if(rank==0)
    for(i=0; i<n; i++) printf("%f \n", x[i]);
  
  
  GA_Dgop(x, n, "max");
  printf("\n");
  
  if(rank==1)
    {
      for(i=0; i<n; i++) printf("%f :: %f\n", x[i], temp[i]);
      
      for(i=0; i<n; i++)
	if(x[i]!=temp[i]) printf(" GA ERROR: \n");
      
    }
}

min_operator(int rank, int nprocs, int n)
{
  int i;
  double x[n], temp[n];

  for(i=0; i<n; i++) x[i]=rand()%5;
  for(i=0; i<n; i++) temp[i]=x[i];
  
  if(rank==1)
    for(i=0; i<n; i++) printf("%f \n", x[i]);
  
  
  GA_Dgop(x, n, "min");
  printf("\n");
  
  if(rank==0)
    {
      for(i=0; i<n; i++) printf("%f :: %f\n", x[i], temp[i]);
      
      /*for(i=0; i<n; i++)
	if(x[i]!=temp[i]*sizeof(int)) printf(" GA ERROR: \n");
      */
    }
}

absmax_operator(int rank, int nprocs, int n)
{
  int i;
  double x[n], temp[n];

  for(i=0; i<n; i++) x[i]=rand()%5;
  for(i=0; i<n; i++) temp[i]=x[i];
  
  if(rank==1)
    for(i=0; i<n; i++) printf("%f \n", x[i]);
  
  
  GA_Dgop(x, n, "absmax");
  printf("\n");
  
  if(rank==0)
    {
      for(i=0; i<n; i++) printf("%f :: %f\n", x[i], temp[i]);
      
      /*for(i=0; i<n; i++)
	if(x[i]!=temp[i]*sizeof(int)) printf(" GA ERROR: \n");
      */
    }
}

absmin_operator(int rank, int nprocs, int n)
{
  int i;
  double x[n], temp[n];

  for(i=0; i<n; i++) x[i]=rand()%5;
  for(i=0; i<n; i++) temp[i]=x[i];
  
  if(rank==1)
    for(i=0; i<n; i++) printf("%f \n", x[i]);
  
  
  GA_Dgop(x, n, "absmin");
  printf("\n");
  
  if(rank==0)
    {
      for(i=0; i<n; i++) printf("%f :: %f\n", x[i], temp[i]);
      
      /*for(i=0; i<n; i++)
	if(x[i]!=temp[i]*sizeof(int)) printf(" GA ERROR: \n");
      */
    }
}


main(int argc, char **argv)
{

  int rank, nprocs, n=10;

  MPI_Init(&argc, &argv);
  GA_Initialize();
  MA_init(C_DBL, 1000, 1000);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  //  addition_operator(rank, nprocs, n);
  GA_Sync();
  //multiply_operator(rank, nprocs, n);
  
  max_operator(rank, nprocs, n);
   printf("\n");
   // min_operator(rank, nprocs, n);
  
   //   absmax_operator(rank, nprocs, n);
   //   absmin_operator(rank, nprocs, n);
  
  if(rank==1)
    printf("GA Test Completed \n");
  GA_Terminate();

  MPI_Finalize();

}
