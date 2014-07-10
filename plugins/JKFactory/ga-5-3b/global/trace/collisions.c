#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif

int **Patches,  *Col, *lastCol;
int p;
   
#define hash(ilo,ihi,p) (((ilo)+(ihi))%(p))
#define hash2(ilo,ihi,jlo,jhi,p) ( ((ilo)+(ihi) + (jlo)+(jhi))%(p) )

int LocFound(int *patch);


int **idim2(row,col)
int row,col;
{
   register int **prow, *pdata, i;
   
   pdata = (int*) calloc(row*col, sizeof(int));
   if(pdata == (int*) NULL){
      printf("Memory allocation failed - data: %dx%d\n",row,col);
      exit(1);
   }
   prow = (int **)calloc(row,sizeof(int *));
   if(prow == (int**) NULL){
      printf("memory allocation failed for prow");
      exit(1);
   }
   for(i=0;i<row;i++){
      prow[i]= pdata;
      pdata += col;
   }
   return(prow);
}

float **fdim2(row,col)
int row,col;
{
   int i;
   register float **prow, *pdata;
   
   pdata = (float*) calloc(row*col, sizeof(float));
   if(pdata == (float*) NULL){
      printf("Memory allocation failed - data: %dx%d\n",row,col);
      exit(1);
   }
   prow = (float **)calloc(row,sizeof(float*));
   if(prow == (float**) NULL){
      printf("memory allocation failed for prow");
      exit(1);
   }
   for(i=0;i<row;i++){
      prow[i]= pdata;
      pdata += col;
   }
   return(prow);
}

int main(argc,argv)
int argc;
char **argv;
{
FILE *fin;
char ln[100];
int i,j,k,patch[5],from,to;
int flag, cur_col, loc, max_col=0, last_col;
unsigned long time;

   if(argc<5){
     printf("Usage: collisions <number of processors> <filename> <from> <to>\n");
     printf("(from,to) specifies the range of collision numbers to be printed\n");
     exit(1);
   }
   sscanf(argv[1],"%d",&p);
   if(p>1000||p<1){
      printf("%d processors ? Please be reasonable ...\n",p);
      exit(1);
   }
   fin = fopen(argv[2],"r");
   if(!fin){
          printf("%s: File Not Found, Exiting ...\n",argv[2]);
          exit(1);
   }

   sscanf(argv[3],"%d",&from);
   sscanf(argv[4],"%d",&to);
   if(from<0||from>p||to>p||to<from){
      printf("Wrong values of <from> and <to> \n");
      exit(1);
   }
    if(!(Col = (int*)calloc(p,sizeof(int)))){
                     printf("couldn't allocate memory 2\n");
                     exit(2);
   }

   Patches = idim2(p,5);
   for(i=0;i<p;i++){
      Patches[0][i]=0;
   }

   for(j=0;;j++){
     fgets(ln,100,fin);
     if(feof(fin))break;
     if(9>sscanf(ln,"%d%d%d%d%d%d%d%d%lu",
            &i,&k,patch+1,patch+2,patch+3,patch+4,&i,&flag,&time))continue;

     loc = LocFound(patch);  
     if(!Patches[loc][0])for(i=1;i<5;i++)Patches[loc][i]=patch[i];
     last_col = Patches[loc][0];

     if(flag == 1){ 
       cur_col = ++Patches[loc][0];
       if(cur_col>p){
         printf("%d -- Error in collision number. Record: %d\n",cur_col,j);
         exit(3);
       }
       Col[cur_col-1]++;
       if(max_col<cur_col) max_col=cur_col;
     }else{
       cur_col = Patches[loc][0]--;
       if(cur_col<0){
         printf("%d -- Error in collision number. Record: %d\n",cur_col,j);
         exit(3);
       }
       Col[cur_col-1] --;
       if(!Col[max_col-1])max_col--;
     }

     /* selective output */
     if((last_col>=from && last_col<=to) || (cur_col>=from && cur_col<=to)){
       printf("%f ",1e-6*time);
       for(i=from-1;i<to;i++)printf("%d ",Col[i]);
       printf("\n");
     }
   }

   return 0;
}


int diff(a,b,n)
int *a,*b,n;
{
  int i;
  for (i=0;i<n;i++) if(a[i]!=b[i])return(1);
  return(0);
}


int LocFound(patch)
int *patch;
/* Uses 2 hash functions and exhaustive search to find existing patch *
 * if not found then returns the first found empty slot               */
{
int loc, empty, i;

  empty = -1;
  loc = hash(patch[1],patch[2],p);
  if(!Patches[loc][0]) empty = loc;
       /* check if this is our patch */
  else if(!diff(Patches[loc]+1,patch+1,4)) return(loc);

  /* not yet found so check another location */
  loc = hash2(patch[1],patch[2],patch[3],patch[4],p);

  if(!Patches[loc][0]) empty = (empty != -1)? empty: loc;
  else if(!diff(Patches[loc]+1,patch+1,4)) return(loc);

  /* all is left is to search */ 

  for(loc=0;loc<p;loc++){
    if(!Patches[loc][0]) empty = (empty != -1)? empty: loc;
    else if(!diff(Patches[loc]+1,patch+1,4)) return(loc); 
  }
  if(empty != -1) return(empty);

  printf("Error in  LocFound ");
  for(i=0;i<5; i++)printf("%d ",patch[i]); printf("\n");
  exit(5);
  return(0); /* never gets here */
}
