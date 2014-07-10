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

#include "ga.h"
#include "macdecls.h"
#include "mp3.h"

extern int na,nz;
extern int me, nproc;
extern int myfirstrow,mylastrow;
double *ga_vecptr;
extern double time_get;


void computeminverse(double *minvptr,double *aptr,int *myrowptr,int *mycolptr)
{
    int i,j,k,l=0;
    for(k=0,i=myfirstrow;i<mylastrow;i++,k++){
        for(j=myrowptr[k];j<myrowptr[k+1];j++,l++){
            if(mycolptr[l]>=i){
                if(mycolptr[l]==i){
                    /*printf("\n%d:i=%d j=%d aptr=%d",me,i,j,aptr[l]);*/
                    minvptr[k]=10.0/aptr[l];
                    if(minvptr[k]<0)minvptr[k]=1.0;
                }
                if(mycolptr[l]>i)
                    minvptr[k]=1.0;
                /*printf("\n%d:l=%d i=%d mycolptr[l]=%d",me,l,i,mycolptr[l]);*/
                l+=(myrowptr[k+1]-j);
                break;
            }
        }
    }
}

void computeminverser(double *minvptr,double *rvecptr,double *minvrptr)
{
    int i,k;
    for(k=0,i=myfirstrow;i<mylastrow;i++,k++){
        minvrptr[k]=minvptr[k]*rvecptr[k];
    }
}

void matvecmul(double *aptr,int ga_vec, double *myresultptr,int isvectormirrored, int *myrowptr, int *mycolptr)
{
#if 1
#else
    int lo,hi;
#endif
    int i,j,k,l=0,zero=0;
    double t0;
    double tmprowsum=0.0;
    double *vecptr=ga_vecptr;

    if(isvectormirrored){
        NGA_Access(ga_vec,&zero,&zero,&vecptr,&zero);
    }
    else { /*preliminary, later this will be done inside loop*/
#if 1
        t0=MP_TIMER();
        na--;
        NGA_Get(ga_vec,&zero,&na,vecptr,&na);
        na++;
        time_get+=(MP_TIMER()-t0);
#else /* 1 */
        NGA_Access(ga_vec,&lo,&hi,&vecptr,&zero);
#endif /* 1 */
    }

    for(k=0,i=myfirstrow;i<mylastrow;i++,k++){
        for(j=myrowptr[k];j<myrowptr[k+1];j++,l++){
            tmprowsum=tmprowsum+aptr[l]*vecptr[mycolptr[l]];
        }
        myresultptr[k]=tmprowsum;
        tmprowsum=0.0;
    }
}
