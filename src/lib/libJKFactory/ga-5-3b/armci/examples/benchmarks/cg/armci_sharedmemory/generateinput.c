#if HAVE_CONFIG_H
#   include "config.h"
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
#if HAVE_STDIO_H
#   include <stdio.h>
#endif

void printsparse_(int *irow, int *icol, double *a,  double *b, int n,int nz)
{
int fd;
int i,max,j;

    fd = open("matrix.bin",O_CREAT|O_WRONLY,0777);
    if(fd<=0) { printf("open failed"); return;}

    write (fd, &n, sizeof(int));
    write (fd, &nz, sizeof(int));
    write (fd, a, sizeof(double)*(nz));
    write (fd, irow, sizeof(int)*(n));
    write (fd, icol, sizeof(int)*(nz));
    write (fd, b, sizeof(double)*(n));

    close(fd);

    printf("dumped sparse matrix: dim = %d %d nonzeros\n", n, nz);

}

int main(int argc, char **argv)
{
int NZROW=8,NZSQ=64;
int i,k,irow[NZROW+1],icol[NZSQ+4],n=NZROW,nz=NZSQ;
double    a[NZSQ],  b[NZROW], xvec[NZROW];
    for(i=0;i<NZROW;i++){
      for(k=0;k<NZROW;k++)
        a[i*NZROW+k] =  k+1;
      xvec[i] = i+3.5;
    }
    for(i=0;i<NZROW;i++){
      irow[i]=i*NZROW+1;
      b[i]=0;
      for(k=0;k<NZROW;k++){
        b[i] = b[i]+xvec[k]*a[i*NZROW+k];
        icol[i*NZROW+k]=i+1;
      }
      printf("\nb[%d]=%f",i,b[i]);
    }
    printsparse_(irow,icol,a,b,n,nz);
    return(1);
}
