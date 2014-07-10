#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_MATH_H
#   include <math.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_SYS_TYPES_H
#   include <sys/types.h>
#endif
#if HAVE_SYS_STAT_H
#   include <sys/stat.h>
#endif
#if HAVE_FCNTL_H
#   include <fcntl.h>
#endif
#if HAVE_STRING_H
#   include <string.h>
#endif

#include "armci.h"
#include "message.h"

extern int na;
extern int nz;
extern double *dmvec,*svec,*bvec,*dvec,*amat,*xvec,*axvec,*rvec,*qvec;
extern int *ridx,*cidx;
extern int me, nproc;
extern int myfirstrow,mylastrow;
static int *columnmap,*allfirstrow,*alllastrow;
static FILE *fd;

void generate_random_file(int naa,int nnz){
    fd = fopen("randominput.dat", "w");
}

void read_and_create(int argc, char **argv)
{
int ri,i;
int tmp1,idealelementsperproc;
void **amatptrs,**xvecptrs;

    na = atoi(argv[1]);
    nz = atoi(argv[2]);

    if(strncmp("random",argv[3],6)){
       if(me==0){
         fd = fopen(argv[3], "r");
         if(fd==NULL)ARMCI_Error("unable to open given file",0);
       }
    }
    else{
       if(na==0 || nz==0){
         printf("\nERROR:exiting-no input file given and na or nz is 0");
         fflush(stdout);
         ARMCI_Finalize();
         armci_msg_finalize();
         return;
       }
       if(me==0){
         generate_random_file(na,nz);
         fd = fopen("randominput.dat", "r");
       }
    }
    if(me==0){
       if(na==0)
         fread(&na, sizeof(na), 1, fd);
       if(nz==0)
         fread(&nz, sizeof(nz), 1, fd);
       printf("\nReading CG input\n");
       printf("Number of rows: %d\n", na);
       printf("Number of non-zeros: %d\n", nz);
    }

    armci_msg_bcast(&nz,sizeof(int),0);
    armci_msg_bcast(&na,sizeof(int),0);
    armci_msg_barrier();

    amatptrs = (void **)malloc(sizeof(void *)*nproc); 
    xvecptrs = (void **)malloc(sizeof(void *)*nproc);
    if(xvecptrs==NULL || amatptrs==NULL)
      ARMCI_Error("xvecptrs amatptrs malloc failed",sizeof(void *)*nproc);

    if(ARMCI_Malloc(amatptrs,((me==0)?(sizeof(double)*nz):0)))
      ARMCI_Error("amat malloc failed",sizeof(double)*nz);
    amat = (double *)amatptrs[0];
    
    if(ARMCI_Malloc(amatptrs,((me==0)?(sizeof(int)*(nz+1)):0)))
      ARMCI_Error("icol malloc failed",sizeof(int)*(nz+1));
    cidx = (int *)amatptrs[0];
    
    ARMCI_Malloc(xvecptrs,((me==0)?(sizeof(int)*(na+1)):0)); /*+1 for end of last row*/
    ridx = (int *)xvecptrs[0];

    ARMCI_Malloc(xvecptrs,((me==0)?(sizeof(double)*(na+1)):0));
    xvec = (double *)xvecptrs[0];

    ARMCI_Malloc(xvecptrs,((me==0)?(sizeof(double)*(na+1)):0));
    bvec = (double *)xvecptrs[0];

    if(me==0){

      for (i = 0; i < na + 1; i++)
        xvec[i] = 0.0;

      fread(amat, sizeof(double), nz, fd);
      fread(ridx, sizeof(int), (na+1), fd);
      ridx[na]=nz;
      fread(cidx, sizeof(int), (nz+1), fd);
      fread(bvec, sizeof(double), (na+1), fd);

      /* the c adjustment */
      for (i = 0; i < na; i++)
        ridx[i] -= 1;
         
      for (i = 0; i < nz; i++)
        cidx[i] -= 1;
    }
   
    armci_msg_barrier();
    /*acg_matvecmul(amat,xvec,bvec,ridx,cidx);*/
    if(0){
    for(i=0;i<nz+1;i++)
      printf("\n%d:amat[%d]=%f icol[%d]=%d",me,i,amat[i],i,cidx[i]);
    for(i=0;i<na+1;i++)
      printf("\n%d:irow[%d]=%d bvec[%d]=%f",me,i,ridx[i],i,bvec[i]);
    }
    allfirstrow = (int *)malloc(sizeof(int)*nproc);
    alllastrow = (int *)malloc(sizeof(int)*nproc);
    columnmap = (int *)malloc(sizeof(int)*nproc);
    if(!allfirstrow || !alllastrow || !columnmap)
      ARMCI_Error("malloc failed allfirstrow ",0);
    armci_msg_barrier();
    /* 
     * next decide who works on which rows, this will decide the
     * distribution of a,d,r,q,x,and ax
     */
    /*create the mapping for all vectors, row matrix and column matrix*/
    if(me==0){
       idealelementsperproc = nz/nproc;
       tmp1=0;
       for(i=0;i<nproc;i++){
         int elementsperproc=0;
         allfirstrow[i]=tmp1;
         for(ri=tmp1;ri<na;ri++,tmp1++){
           elementsperproc+=(ridx[ri+1]-ridx[ri]);
       if(elementsperproc>=idealelementsperproc){
             if((elementsperproc-idealelementsperproc) > 
                idealelementsperproc-(elementsperproc-(ridx[ri+1]-ridx[ri]))){
               alllastrow[i] = ri-1;  
           if((ri-1)<0)ARMCI_Error("run on a smaller processor count",0);
               /*tmp1--;*/
             }
             else{
               alllastrow[i] = ri;  
               if(ri<0)ARMCI_Error("run on a smaller processor count",0);
               tmp1++;
             }
             elementsperproc=0;
             break;
       }
         }
       }
       alllastrow[nproc-1]=na-1;
       for(i=0;i<nproc;i++)columnmap[i]=ridx[allfirstrow[i]];
    }
    armci_msg_bcast(columnmap,nproc*sizeof(int),0);
    armci_msg_bcast(allfirstrow,nproc*sizeof(int),0);
    armci_msg_bcast(alllastrow,nproc*sizeof(int),0);
    myfirstrow = allfirstrow[me];
    mylastrow = alllastrow[me];
    if(me==0)for(i=0;i<nproc;i++){
      printf("\nDISTRIBUTION:first row of process\t%d is %d last row of process\t%d is %d",i,allfirstrow[i],i,alllastrow[i]);
    }
    /*
    for(i=myfirstrow;i<mylastrow;i++){
            xvec[i]=0.0;
    }
    */
    ARMCI_Malloc(xvecptrs,((me==0)?(sizeof(double)*na):0));
    rvec = (double *)xvecptrs[0];

    ARMCI_Malloc(xvecptrs,((me==0)?(sizeof(double)*na):0));
    dvec = (double *)xvecptrs[0];

    ARMCI_Malloc(xvecptrs,((me==0)?(sizeof(double)*na):0));
    svec = (double *)xvecptrs[0];

    ARMCI_Malloc(xvecptrs,((me==0)?(sizeof(double)*na):0));
    dmvec = (double *)xvecptrs[0];

    ARMCI_Malloc(xvecptrs,((me==0)?(sizeof(double)*na):0));
    qvec = (double *)xvecptrs[0];

    ARMCI_Malloc(xvecptrs,((me==0)?(sizeof(double)*na):0));
    axvec = (double *)xvecptrs[0];

    if(me==0)fclose(fd);
    /*dont forget to free mallocs*/
    free(allfirstrow);
    free(alllastrow);
    free(columnmap);
}
