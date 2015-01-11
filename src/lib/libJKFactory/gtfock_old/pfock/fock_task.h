#ifndef __FOCK_TASK_H__
#define __FOCK_TASK_H__


#include "CInt.h"

void fock_task(BasisSet_t basis, int ncpu_f, int num_dmat,
               int *shellptr, double *shellvalue,
               int *shellid, int *shellrid, int *f_startind,
               int *rowpos, int *colpos, int *rowptr, int *colptr,
               double tolscr2, int startrow, int startcol,
               int startM, int endM, int startP, int endP,
               double **D1, double **D2, double **D3,
               double *F1, double *F2, double *F3,
               double *F4, double *F5, double *F6,
               int ldX1, int ldX2, int ldX3,
               int ldX4, int ldX5, int ldX6,
               int sizeX1, int sizeX2, int sizeX3,
               int sizeX4, int sizeX5, int sizeX6,
               double *nitl, double *nsq);

/*void fock_task(BasisSet_t basis, ERD_t erd, int ncpu_f, int num_dmat,
               int *shellptr, double *shellvalue,
               int *shellid, int *shellrid, int *f_startind,
               int *rowpos, int *colpos, int *rowptr, int *colptr,
               double tolscr2, int startrow, int startcol,
               int startM, int endM, int startP, int endP,
               double **D1, double **D2, double **D3,
               double *F1, double *F2, double *F3,
               double *F4, double *F5, double *F6, 
               int ldX1, int ldX2, int ldX3,
               int ldX4, int ldX5, int ldX6,
               int sizeX1, int sizeX2, int sizeX3,
               int sizeX4, int sizeX5, int sizeX6,
               double *nitl, double *nsq);*/

void reset_F(int numF, int num_dmat, double *F1, double *F2, double *F3,
             double *F4, double *F5, double *F6,
             int sizeX1, int sizeX2, int sizeX3,
             int sizeX4, int sizeX5, int sizeX6);

void reduce_F(int numF, int num_dmat,
              double *F1, double *F2, double *F3,
              double *F4, double *F5, double *F6,
              int sizeX1, int sizeX2, int sizeX3,
              int sizeX4, int sizeX5, int sizeX6,
              int maxrowsize, int maxcolsize,
              int maxrowfuncs, int maxcolfuncs,
              int iX3M, int iX3P,
              int ldX3, int ldX4, int ldX5, int ldX6);


#endif /* #define __FOCK_TASK_H__ */
