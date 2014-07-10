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

#include "ga.h"
#include "macdecls.h"
#include "mp3.h"

extern int na;
extern int nz;
extern int dmvec,svec,bvec,dvec,m_dvec,amat,xvec,axvec,rvec,qvec,ridx,cidx;
extern int me, nproc;
extern int myfirstrow,mylastrow;
static int *columnmap,*allfirstrow,*alllastrow;
extern double *ga_vecptr;
extern int isvectormirrored;
static FILE *fd;

void generate_random_file(int naa,int nnz){
    fd = fopen("randominput.dat", "w");
}

void read_and_create(int argc, char **argv)
{
int ri,i;
int ph;
double d_one=1.0;
double *a,*x;
int *icol, *irow;
#if 0
int dims[2];
#endif
int tmp1,idealelementsperproc;
int lo,hi;

    na = atoi(argv[1]);
    nz = atoi(argv[2]);

    if(strncmp("random",argv[3],6)){
       if(me==0){
         fd = fopen(argv[3], "r");
         if(fd==NULL)GA_Error("unable to open given file",0);
       }
    }
    else{
       if(na==0 || nz==0){
         printf("\nERROR:exiting-no input file given and na or nz is 0");
         fflush(stdout);
         GA_Terminate();
         MP_FINALIZE();
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

       a = (double *)malloc(sizeof(double)*nz);
       icol = (int *)malloc(sizeof(int)*(nz+1));
       x = (double *)malloc(sizeof(double)*(na+1));
       irow = (int *)malloc(sizeof(int)*(na+1));

       for (i = 0; i < na + 1; i++)
         x[i] = 1.0;

       fread(a, sizeof(double), nz, fd);
       fread(irow, sizeof(int), na + 1, fd);
       fread(icol, sizeof(int), nz + 1, fd);
       fread(x, sizeof(double), na + 1, fd);

       /* the c adjustment */
       for (i = 0; i < na + 1; i++)
         irow[i] -= 1;
       for (i = 0; i < nz + 1; i++)
         icol[i] -= 1;
       /* MPI_Bcast(&nz,1,MPI_INT,0,MPI_COMM_WORLD); */
       GA_Brdcst(&nz, sizeof(int), 0);
       /* MPI_Bcast(&na,1,MPI_INT,0,MPI_COMM_WORLD); */
       GA_Brdcst(&na, sizeof(int), 0);
    }
    else {
       /* MPI_Bcast(&nz,1,MPI_INT,0,MPI_COMM_WORLD); */
       GA_Brdcst(&nz, sizeof(int), 0);
       /* MPI_Bcast(&na,1,MPI_INT,0,MPI_COMM_WORLD); */
       GA_Brdcst(&na, sizeof(int), 0);
       /*for now, others dont need to malloc really*/
       a = (double *)malloc(sizeof(double)*nz);
       icol = (int *)malloc(sizeof(int)*(nz+1));
       x = (double *)malloc(sizeof(double)*(na+1));
       irow = (int *)malloc(sizeof(int)*(na+1));
       if(!a || !icol || !x || !irow)GA_Error("malloc failed in ga_cg",0);
    }
    allfirstrow = (int *)malloc(sizeof(int)*nproc);
    alllastrow = (int *)malloc(sizeof(int)*nproc);
    columnmap = (int *)malloc(sizeof(int)*nproc);
    if(!allfirstrow || !alllastrow || !columnmap)GA_Error("malloc failed in ga_cg",0);

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
           elementsperproc+=(irow[ri+1]-irow[ri]);
       if(elementsperproc>=idealelementsperproc){
             if((elementsperproc-idealelementsperproc) > 
                idealelementsperproc-(elementsperproc-(irow[ri+1]-irow[ri]))){
               alllastrow[i] = ri-1;  
           if((ri-1)<0)GA_Error("run on a smaller processor count",0);
               tmp1--;
             }
             else{
               alllastrow[i] = ri;  
               if(ri<0)GA_Error("run on a smaller processor count",0);
             }
             elementsperproc=0;
             break;
       }
         }
       }
       alllastrow[nproc-1]=na-1;
       for(i=0;i<nproc;i++)columnmap[i]=irow[allfirstrow[i]];
    }
    /* MPI_Bcast(columnmap,nproc,MPI_INT,0,MPI_COMM_WORLD); */
    GA_Brdcst(columnmap, sizeof(int)*nproc, 0);
    /* MPI_Bcast(allfirstrow,nproc,MPI_INT,0,MPI_COMM_WORLD); */
    GA_Brdcst(allfirstrow, sizeof(int)*nproc, 0);
    /* MPI_Bcast(alllastrow,nproc,MPI_INT,0,MPI_COMM_WORLD); */
    GA_Brdcst(alllastrow, sizeof(int)*nproc, 0);
    myfirstrow = allfirstrow[me];
    mylastrow = alllastrow[me];

    
    amat = NGA_Create_irreg(MT_C_DBL, 1, &nz , "A", &nproc,columnmap);
    if(!amat) GA_Error("create failed: A",nz); 
    if(me==0){
      lo=0;
      hi=nz-1;
      NGA_Put(amat,&lo,&hi,a,&hi);
    }

    /*now create column matrix */
    cidx = NGA_Create_irreg(MT_C_INT, 1, &nz , "COLUMN",&nproc,columnmap);
    if(!cidx) GA_Error("create cidx failed",nz); 
    if(me==0){
      lo=0;
      hi=nz-1;
      NGA_Put(cidx,&lo,&hi,icol,&hi);
    }
    GA_Print_distribution(cidx);
    GA_Sync();
    /* ************** */

    /*now create row matrix and all vectors*/

    ridx = NGA_Create_irreg(MT_C_INT, 1, &na , "ROW", &nproc,allfirstrow);
    if(!ridx) GA_Error("create ridx failed",na); 
    lo=0;
    hi=na-1;
    if(me==0){
      NGA_Put(ridx,&lo,&hi,irow,&hi);
    }
    GA_Sync();
    GA_Print_distribution(ridx);
    MP_BARRIER();

    xvec = NGA_Create_irreg(MT_C_DBL, 1, &na , "X",&nproc,allfirstrow);
    if(!xvec) GA_Error("create x failed",na); 
    GA_Fill(xvec, &d_one);
    /*GA_Print(xvec);*/

    bvec = NGA_Create_irreg(MT_C_DBL, 1, &na , "B",&nproc,allfirstrow);
    if(!bvec) GA_Error("create b failed",na); 
    if(me==0){
      lo=0;
      hi=na-1;
      NGA_Put(bvec,&lo,&hi,x,&hi);
    }
    /*GA_Print(xvec);*/

    axvec = NGA_Create_irreg(MT_C_DBL, 1, &na , "Ax", &nproc,allfirstrow);
    if(!axvec) GA_Error("create Ax failed",na); 

    rvec = NGA_Create_irreg(MT_C_DBL, 1, &na , "R", &nproc, allfirstrow);
    if(!rvec) GA_Error("create r failed",na); 

    dvec = GA_Duplicate(rvec,"D");                /* d = r */

    svec = NGA_Create_irreg(MT_C_DBL, 1, &na , "S", &nproc, allfirstrow);
    if(!svec) GA_Error("create s failed",na); 

    dmvec = NGA_Create_irreg(MT_C_DBL, 1, &na , "Dm", &nproc, allfirstrow);
    if(!dmvec) GA_Error("create dm failed",na); 

#if 0
    dims[0]=1;dims[1]=na;
    t_rvec = NGA_Create(MT_C_DBL, 2, dims, "T_R", NULL);
    if(!t_rvec) GA_Error("create q failed",na); 
#endif

    qvec = NGA_Create_irreg(MT_C_DBL, 1, &na , "Q", &nproc, allfirstrow);
    if(!qvec) GA_Error("create q failed",na); 
    /*GA_Print_distribution(qvec);*/

    if(isvectormirrored){
       ph= GA_Pgroup_get_mirror();
       m_dvec = 
             NGA_Create_irreg_config(MT_C_DBL,1,&na,"mD",&nproc,allfirstrow,ph);
       if(!m_dvec) GA_Error("create mirrored dvec failed",na);
    }

    /* JAD What is GA_Malloc_local?  Where is it? */
    /* ga_vecptr = (double *)GA_Malloc_local(na*sizeof(double)); */

    
    if(me==0)fclose(fd);
    /*dont forget to free mallocs*/
    free(allfirstrow);
    free(alllastrow);
    free(columnmap);
    free(a);
    free(x);
    free(irow);
    free(icol);
}
