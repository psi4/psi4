#define f_matvecmul        F77_FUNC_(f_matvecmul,F_MATVECMUL)
#define f_computeminverse  F77_FUNC_(f_computeminverse,F_COMPUTEMINVERSE)
#define f_computeminverser F77_FUNC_(f_computeminverser,F_COMPUTEMINVERSER)
#define f_addvec           F77_FUNC_(f_addvec,F_ADDVEC)
#define f_2addvec          F77_FUNC_(f_2addvec,F_2ADDVEC)

extern void f_matvecmul(double*,double*,double*,int*,int*,int*,int*,int*);
extern void f_computeminverse(double*,double*,int*,int*,int*,int*); 
extern void f_computeminverser(double*,double*,double*,int*,int*); 
extern void f_addvec(double*,double*,double*,double*,double*,int*,int*);
extern void f_2addvec(double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,int*,int*);
