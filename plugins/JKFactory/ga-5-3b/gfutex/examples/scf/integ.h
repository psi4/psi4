#ifdef __cplusplus
extern "C" {
#endif
  extern double exprjh(double x);
  extern void setfm(void);
  extern void f0(double* value, double t);
  extern void addin(double *g, int *i, int *j, int *k, int *l, double *fock,
		    double *dens, int *iky);
  extern void dfill(int *n, double *val, double *a, int *ia);
#ifdef __cplusplus
}
#endif
