#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_MATH_H
#   include <math.h>
#endif

#include "armci.h"
#include "message.h"

extern int na,nz;
extern int me, nproc;
extern int myfirstrow,mylastrow;

void computeminverse(double *minvptr,double *aptr,int *rowptr,int *colptr)
{
int i,j;
    for(i=myfirstrow;i<=mylastrow;i++){
      for(j=rowptr[i];j<rowptr[i+1];j++){
        if(colptr[j]>=i){
          if(colptr[j]==i){
            /*printf("\n%d:i=%d j=%d aptr=%f",me,i,j,aptr[j]);*/
            minvptr[i]=10.0/aptr[j];
          }
          if(colptr[j]>i){
            minvptr[i]=0.0;
            /*printf("\n%d:l=%d i=%d mycolptr[j]=%d",me,j,i,colptr[j]);*/
          }
          break;
        }
      }
    }
    /*armci_msg_barrier();*/
}

void computeminverser(double *minvptr,double *rvecptr,double *minvrptr)
{
int i;
    for(i=myfirstrow;i<=mylastrow;i++)
       minvrptr[i]=minvptr[i]*rvecptr[i];
    /*armci_msg_barrier();*/
}

void acg_printvec2(char *v, double *vec, char *v1, double *vec1)
{
int i;
    for(i=myfirstrow;i<=mylastrow;i++)
      printf("\n%d:%s[%d]=%f %s[%d]=%f",me,v,i,vec[i],v1,i,vec1[i]);
    fflush(stdout);
    armci_msg_barrier();
}

void acg_printvec(char *v, double *vec)
{
int i;
    for(i=myfirstrow;i<=mylastrow;i++)
      printf("\n%d:%s[%d]=%f",me,v,i,vec[i]);
    fflush(stdout);
    armci_msg_barrier();
}

double acg_ddot(double *vec1,double *vec2)
{
int i;
double dt=0.0;
    for(i=myfirstrow;i<=mylastrow;i++)
      dt+=(vec1[i]*vec2[i]);
    armci_msg_dgop(&dt,1,"+");
    /*armci_msg_barrier();*/
    return(dt);
}


void acg_zero(double *vec1)
{
int i;
    for(i=myfirstrow;i<=mylastrow;i++)
      vec1[i]=0.0;
    armci_msg_barrier();
}

void acg_addvec(double *pscale1,double *vec1,double *pscale2,double *vec2, double *result)
{
int i;
double scale1=*pscale1,scale2=*pscale2;
    for(i=myfirstrow;i<=mylastrow;i++)
      result[i]=(scale1*vec1[i]+scale2*vec2[i]);
    /*armci_msg_barrier();*/
}

void acg_2addvec(double *pscale1a,double *vec1a, double *pscale2a,double *vec2a,
        double *resulta, 
                double *pscale1b, double *vec1b,double *pscale2b, double *vec2b,
        double *resultb, 
        int *rowptr, int *colptr)
{
int i;
double scale1a=*pscale1a,scale2a=*pscale2a, scale1b=*pscale1b,scale2b=*pscale2b;
    for(i=myfirstrow;i<=mylastrow;i++){
      resulta[i]=vec1a[i]*scale1a+vec2a[i]*scale2a;
      resultb[i]=vec1b[i]*scale1b+vec2b[i]*scale2b;
    }
    /*armci_msg_barrier();*/
}

void acg_matvecmul(double *aptr,double *vec, double *result,int *rowptr, int *colptr)
{
int i,j;
double tmprowsum=0.0;
    ARMCI_Barrier();
    for(i=myfirstrow;i<=mylastrow;i++){
       for(j=rowptr[i];j<rowptr[i+1];j++){
         tmprowsum=tmprowsum+aptr[j]*vec[colptr[j]];
         /*printf("\n%d:%d %d %f %f %f",
            me,j,colptr[j],aptr[j],vec[colptr[j]],tmprowsum);*/
       }
       result[i]=tmprowsum;
       tmprowsum=0.0;
    }
    /*ARMCI_Barrier();*/
}
