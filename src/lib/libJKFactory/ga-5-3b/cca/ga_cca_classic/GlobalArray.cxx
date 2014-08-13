#include "gacca.h"

#define GA_DATA_TYPE     C_DBL

static int sTmpVar = 0;


/**
 * Constructors and Destructor of GlobalArray          
 */


GA::GlobalArray::GlobalArray(int type, int ndim, int dims[], char *arrayname,
			 int chunk[]) {
  mHandle = NGA_Create(type, ndim, dims, arrayname, chunk);
  if(!mHandle)  GA_Error((char *)" GA creation failed",0);
}

GA::GlobalArray::GlobalArray(int type, int ndim, int dims[], char *arrayname,
                         int block[], int maps[]) {
  mHandle = NGA_Create_irreg(type, ndim, dims, arrayname, block, maps);
  if(!mHandle) GA_Error((char *)" GA creation failed",0);
}

GA::GlobalArray::GlobalArray(int type, int ndim, int dims[], int width[], 
			     char *arrayname, int chunk[], char ghosts) {
  mHandle = NGA_Create_ghosts(type, ndim, dims, width, arrayname, chunk);
  if(!mHandle)  GA_Error((char *)" GA creation failed",0);
}

GA::GlobalArray::GlobalArray(int type, int ndim, int dims[], int width[], 
			     char *arrayname, int block[], int maps[], 
			     char ghosts) {
  mHandle = NGA_Create_ghosts_irreg(type, ndim, dims, width, arrayname, block, 
				    maps);
  if(!mHandle) GA_Error((char *)" GA creation failed",0);
}

GA::GlobalArray::GlobalArray(const GA::GlobalArray &g_a, char *arrayname) {
  mHandle = GA_Duplicate(g_a.mHandle, arrayname);
  if(!mHandle) GA_Error((char *)" GA creation failed",0);
}

GA::GlobalArray::GlobalArray(const GA::GlobalArray &g_a) {
  char temp_name[20];
  
  sprintf(temp_name, "tmpGA%d", sTmpVar++);
  mHandle = GA_Duplicate(g_a.mHandle, temp_name);
  if(!mHandle) GA_Error((char *)" GA creation failed",0);
  GA_Copy(g_a.mHandle, mHandle);
}

GA::GlobalArray::GlobalArray() {
  char temp_name[20];
  int n_dim;
  int *dimensions;
  
  sprintf(temp_name, "tmpGA%d", sTmpVar++);
  n_dim = DEF_NDIM;
  dimensions = new int [n_dim];    
  for(int i=0; i<n_dim; i++) dimensions[i] = DEF_DIMS;
  mHandle = NGA_Create(GA_DATA_TYPE, n_dim, dimensions, temp_name, NULL);
  if(!mHandle) GA_Error((char *)" GA creation failed",0);
  //free(dimensions);
  delete[] dimensions;
}

GA::GlobalArray::~GlobalArray() {
 
}


/*********************************************************************
 *          Member functions of GA::GlobalArray                          *
 *********************************************************************/

void 
GA::GlobalArray::acc(int lo[], int hi[], void *buf, int ld[], void *alpha) const {
  NGA_Acc(mHandle, lo, hi, buf, ld, alpha);
}

void 
GA::GlobalArray::access(int lo[], int hi[], void *ptr, int ld[]) const {
  NGA_Access(mHandle, lo, hi, ptr, ld);
}

void 
GA::GlobalArray::accessGhosts(int dims[], void *ptr, int ld[]) const {
  NGA_Access_ghosts(mHandle, dims, ptr, ld);
}

void 
GA::GlobalArray::accessGhostElement(void *ptr, int subscript[], int ld[]) const {
  NGA_Access_ghost_element(mHandle, ptr, subscript, ld);
}

void 
GA::GlobalArray::add(void *alpha, const GA::GlobalArray * g_a, 
		     void *beta,  const GA::GlobalArray * g_b) const {
  GA_Add(alpha, g_a->mHandle, beta, g_b->mHandle, mHandle);
} 

void  
GA::GlobalArray::addPatch (void *alpha, 
			   const GA::GlobalArray * g_a, int alo[], int ahi[],
			   void *beta,  
			   const GA::GlobalArray * g_b, int blo[], int bhi[],
			   int clo[], int chi[]) const {
  NGA_Add_patch(alpha, g_a->mHandle, alo, ahi, beta, 
		g_b->mHandle, blo, bhi, mHandle, clo, chi);
}

void  
GA::GlobalArray::checkHandle(char* string) const {
  GA_Check_handle(mHandle, string);
}

int 
GA::GlobalArray::compareDistr(const GA::GlobalArray *g_a) const {
  return GA_Compare_distr(mHandle, g_a->mHandle);
}

void 
GA::GlobalArray::copy(const GA::GlobalArray *g_a) const {
  GA_Copy(g_a->mHandle, mHandle);
}

void 
GA::GlobalArray::copyPatch(char trans, const GA::GlobalArray* ga, int alo[], 
			   int ahi[], int blo[], int bhi[]) const {
  NGA_Copy_patch(trans, ga->mHandle, alo, ahi, mHandle, blo, bhi);
}

double 
GA::GlobalArray::ddot(const GA::GlobalArray * g_a) const {
  return GA_Ddot(mHandle, g_a->mHandle);
}

double 
GA::GlobalArray::ddotPatch(char ta, int alo[], int ahi[],
			   const GA::GlobalArray * g_a, char tb, 
			   int blo[], int bhi[]) const {
  return NGA_Ddot_patch(mHandle, ta, alo, ahi, g_a->mHandle, tb, blo, bhi);
}

void 
GA::GlobalArray::destroy() const {
  GA_Destroy(mHandle);
}

void 
GA::GlobalArray::dgemm(char ta, char tb, int m, int n, int k, double alpha,  
		       const GA::GlobalArray *g_a, const GA::GlobalArray *g_b, 
		       double beta) const {
  GA_Dgemm(ta, tb, m, n, k, alpha, g_a->mHandle, g_b->mHandle, 
	   beta, mHandle);
}

void 
GA::GlobalArray::diag(const GA::GlobalArray *g_s, GA::GlobalArray *g_v, 
		      void *eval) const {
  GA_Diag(mHandle, g_s->mHandle, g_v->mHandle, eval);
}

void 
GA::GlobalArray::diagReuse(int control, const GA::GlobalArray *g_s, 
			   GA::GlobalArray *g_v, void *eval) const {
  GA_Diag_reuse(control, mHandle, g_s->mHandle, g_v->mHandle, eval);
}

void 
GA::GlobalArray::diagStd(GlobalArray *g_v, void *eval) const {
  GA_Diag_std(mHandle, g_v->mHandle, eval);
}

void 
GA::GlobalArray::diagSeq(const GA::GlobalArray * g_s, 
			 const GA::GlobalArray * g_v, void *eval) const {
  GA_Diag_seq(mHandle, g_s->mHandle, g_v->mHandle, eval);
}

void 
GA::GlobalArray::diagStdSeq(const GA::GlobalArray * g_v, void *eval) const {
  GA_Diag_std_seq(mHandle, g_v->mHandle, eval);
}

void 
GA::GlobalArray::distribution(int me, int* lo, int* hi) const {
  NGA_Distribution(mHandle, me, lo, hi);
}

float 
GA::GlobalArray::fdot(const GA::GlobalArray * g_a) const {
  return GA_Fdot(mHandle, g_a->mHandle);
}

float 
GA::GlobalArray::fdotPatch(char t_a, int alo[], int ahi[], 
			   const GA::GlobalArray * g_b, char t_b, int blo[], 
			   int bhi[]) const {
  return NGA_Fdot_patch(mHandle, t_a, alo, ahi,
			g_b->mHandle, t_b, blo, bhi);
}

void 
GA::GlobalArray::fill(void *value) const {
  GA_Fill(mHandle, value);
}

void 
GA::GlobalArray::fillPatch (int lo[], int hi[], void *val) const {
  NGA_Fill_patch(mHandle, lo, hi, val);
}

void 
GA::GlobalArray::gather(void *v, int * subsarray[], int n) const {
  NGA_Gather(mHandle, v, subsarray, n);
}

void 
GA::GlobalArray::get(int lo[], int hi[], void *buf, int ld[])  const {
  NGA_Get(mHandle, lo, hi, buf, ld);
}

int 
GA::GlobalArray::hasGhosts()  const {
  return GA_Has_ghosts(mHandle);
}

Integer 
GA::GlobalArray::idot(const GA::GlobalArray * g_a)  const {
  return GA_Idot(mHandle, g_a->mHandle);
}

long 
GA::GlobalArray::idotPatch(char ta, int alo[], int ahi[],
			   const GA::GlobalArray * g_a, char tb, 
			   int blo[], int bhi[])  const {
  return NGA_Idot_patch(mHandle, ta, alo, ahi, g_a->mHandle, tb, blo, bhi);
}

void 
GA::GlobalArray::inquire(int *type, int *ndim, int dims[])  const {
  NGA_Inquire(mHandle, type, ndim, dims);
}

char* 
GA::GlobalArray::inquireName()  const {
  return GA_Inquire_name(mHandle);
}

long 
GA::GlobalArray::ldot(const GA::GlobalArray * g_a)  const {
  return GA_Ldot(mHandle, g_a->mHandle);
}

int 
GA::GlobalArray::lltSolve(const GA::GlobalArray * g_a)  const {
  return GA_Llt_solve(g_a->mHandle, mHandle);
}

int 
GA::GlobalArray::locate(int subscript[])  const {
  return NGA_Locate(mHandle, subscript);
} 

int 
GA::GlobalArray::locateRegion(int lo[], int hi[], int map[], int procs[])  const {
  return NGA_Locate_region(mHandle, lo, hi, map, procs);
}

void 
GA::GlobalArray::luSolve(char trans, const GA::GlobalArray * g_a)  const {
  GA_Lu_solve(trans, g_a->mHandle, mHandle);
}

void 
GA::GlobalArray::matmulPatch(char transa, char transb, void* alpha, 
			     void *beta, const GA::GlobalArray *g_a, 
			     int ailo, int aihi, int ajlo, int ajhi,
			     const GA::GlobalArray *g_b, 
			     int bilo, int bihi, int bjlo, int bjhi,
			     int cilo, int cihi, int cjlo, int cjhi)  const {
  GA_Matmul_patch(transa, transb, alpha, beta, 
		  g_a->mHandle, ailo, aihi, ajlo, ajhi, 
		  g_b->mHandle, bilo, bihi, bjlo, bjhi, 
		  mHandle, cilo, cihi, cjlo, cjhi);
}

void 
GA::GlobalArray::matmulPatch(char transa, char transb, void* alpha, void *beta,
			     const GA::GlobalArray *g_a, int *alo, int *ahi,
			     const GA::GlobalArray *g_b, int *blo, int *bhi,
			     int *clo, int *chi)  const {
  NGA_Matmul_patch(transa, transb, alpha, beta, g_a->mHandle, alo, ahi, 
		  g_b->mHandle, blo, bhi, mHandle, clo, chi);
}

void 
GA::GlobalArray::nblock(int numblock[])  const {
  GA_Nblock(mHandle, numblock);
}

int 
GA::GlobalArray::ndim()  const {
  return GA_Ndim(mHandle);
}

void 
GA::GlobalArray::periodicAcc(int lo[], int hi[], void* buf, int ld[], 
			     void* alpha)  const {
  NGA_Periodic_acc(mHandle, lo, hi, buf, ld, alpha);
}

void 
GA::GlobalArray::periodicGet(int lo[], int hi[], void* buf, int ld[])  const {
  NGA_Periodic_get(mHandle, lo, hi, buf, ld);
}

void 
GA::GlobalArray::periodicPut(int lo[], int hi[], void* buf, int ld[])  const {
  NGA_Periodic_put(mHandle, lo, hi, buf, ld);
}

void 
GA::GlobalArray::print() const {
  GA_Print(mHandle);
}

void 
GA::GlobalArray::printDistribution() const {
  GA_Print_distribution(mHandle);
}

void 
GA::GlobalArray::printFile(FILE *file)  const {
  GA_Print_file(file, mHandle);
}

void 
GA::GlobalArray::printPatch(int* lo, int* hi, int pretty)  const {
  NGA_Print_patch(mHandle, lo, hi, pretty);   
}

void 
GA::GlobalArray::procTopology(int proc, int coord[])  const {
  NGA_Proc_topology(mHandle, proc, coord);
}

void 
GA::GlobalArray::put(int lo[], int hi[], void *buf, int ld[])  const {
  NGA_Put(mHandle, lo, hi, buf, ld);
}

long 
GA::GlobalArray::readInc(int subscript[], long inc)  const {
  return NGA_Read_inc(mHandle, subscript, inc);
}

void 
GA::GlobalArray::release(int lo[], int hi[])  const {
  NGA_Release(mHandle, lo, hi);
}

void 
GA::GlobalArray::releaseUpdate(int lo[], int hi[])  const {
  NGA_Release_update(mHandle, lo, hi);
}

void 
GA::GlobalArray::scale(void *value)  const { 
  GA_Scale(mHandle, value); 
}

void 
GA::GlobalArray::scalePatch (int lo[], int hi[], void *val)  const {
  NGA_Scale_patch(mHandle, lo, hi, val);
}

void 
GA::GlobalArray::scatter(void *v, int *subsarray[], int n)  const {
  NGA_Scatter(mHandle, v, subsarray, n);
}

int 
GA::GlobalArray::solve(const GA::GlobalArray * g_a)  const {
  return GA_Solve(g_a->mHandle, mHandle);
}

int 
GA::GlobalArray::spdInvert()  const {
  return GA_Spd_invert(mHandle);
}

void 
GA::GlobalArray::selectElem(char *op, void* val, int index[])  const {
  NGA_Select_elem(mHandle, op, val, index);
}

void 
GA::GlobalArray::sgemm(char ta, char tb, int m, int n, int k, float alpha,  
		       const GA::GlobalArray *g_a, const GA::GlobalArray *g_b, 
		       float beta)  const {
  GA_Sgemm(ta, tb, m, n, k, alpha, g_a->mHandle, g_b->mHandle, 
	   beta, mHandle);
}

void 
GA::GlobalArray::symmetrize()  const { 
  GA_Symmetrize(mHandle); 
}

void 
GA::GlobalArray::transpose(const GA::GlobalArray * g_a)  const {
  GA_Transpose(mHandle, g_a->mHandle);
}

void 
GA::GlobalArray::updateGhosts()  const {
  GA_Update_ghosts(mHandle);
}

int 
GA::GlobalArray::updateGhostDir(int dimension, int idir, int cflag)  const {
  return NGA_Update_ghost_dir(mHandle, dimension, idir, cflag);
}

DoubleComplex 
GA::GlobalArray::zdot(const GA::GlobalArray * g_a)  const {
  return GA_Zdot(mHandle, g_a->mHandle);
}

DoubleComplex 
GA::GlobalArray::zdotPatch(char ta, int alo[], int ahi[],
			   const GA::GlobalArray * g_a, char tb, 
			   int blo[], int bhi[])  const {
  return NGA_Zdot_patch(mHandle, ta, alo, ahi, g_a->mHandle, tb, blo, bhi);
}

void 
GA::GlobalArray::zero()  const { 
  GA_Zero(mHandle); 
}

void 
GA::GlobalArray::zeroPatch (int lo[], int hi[])  const {
  NGA_Zero_patch(mHandle, lo, hi);
}

void 
GA::GlobalArray::zgemm(char ta, char tb, int m, int n, int k, 
		       DoubleComplex alpha,  
		       const GA::GlobalArray *g_a, const GA::GlobalArray *g_b, 
		       DoubleComplex beta)  const {
  GA_Zgemm(ta, tb, m, n, k, alpha, g_a->mHandle, g_b->mHandle, 
	   beta, mHandle);
}

/* recent additions */

void 
GA::GlobalArray::absValue()  const {
  GA_Abs_value(mHandle);
}

void 
GA::GlobalArray::addConstant(void* alpha)  const {
  GA_Add_constant(mHandle, alpha);
}

void 
GA::GlobalArray::recip()  const {
  GA_Recip(mHandle);
}

void 
GA::GlobalArray::elemMultiply(const GA::GlobalArray * g_a, 
			      const GA::GlobalArray * g_b)  const {
  GA_Elem_multiply(g_a->mHandle, g_b->mHandle, mHandle);
}

void 
GA::GlobalArray::elemDivide(const GA::GlobalArray * g_a, 
			    const GA::GlobalArray * g_b)  const {
  GA_Elem_divide(g_a->mHandle, g_b->mHandle, mHandle);
}


void 
GA::GlobalArray::elemMaximum(const GA::GlobalArray * g_a, 
			     const GA::GlobalArray * g_b)  const {
  GA_Elem_maximum(g_a->mHandle, g_b->mHandle, mHandle);
}


void 
GA::GlobalArray::elemMinimum(const GA::GlobalArray * g_a, 
			     const GA::GlobalArray * g_b)  const {
  GA_Elem_minimum(g_a->mHandle, g_b->mHandle, mHandle);
}

void 
GA::GlobalArray::absValuePatch(int *lo, int *hi)  const {
  GA_Abs_value_patch(mHandle, lo, hi);
}

void 
GA::GlobalArray::addConstantPatch(int *lo,int *hi, void *alpha)  const {
  GA_Add_constant_patch(mHandle, lo, hi, alpha);
}

void 
GA::GlobalArray::recipPatch(int *lo, int *hi)  const {
  GA_Recip_patch(mHandle, lo, hi);
}

void 
GA::GlobalArray::stepMax(const GA::GlobalArray * g_a, double *step)  const {// CHECK all Step Max functions
  GA_Step_max(mHandle, g_a->mHandle, step);
}

void 
GA::GlobalArray::stepMaxPatch(int *alo, int *ahi, 
			      const GA::GlobalArray * g_b, int *blo, int *bhi, 
			      double *step)  const {
  GA_Step_max_patch(mHandle, alo, ahi, g_b->mHandle, blo, bhi, step);
}

void 
GA::GlobalArray::elemMultiplyPatch(const GA::GlobalArray * g_a,
				   int *alo,int *ahi,
				   const GA::GlobalArray * g_b,
				   int *blo,int *bhi,
				   int *clo,int *chi)  const {
  GA_Elem_multiply_patch(g_a->mHandle, alo, ahi, g_b->mHandle, blo, bhi, 
			 mHandle, clo, chi);
}

void
GA::GlobalArray::elemDividePatch(const GA::GlobalArray * g_a,int *alo,int *ahi,
				 const GA::GlobalArray * g_b,int *blo,int *bhi,
				 int *clo,int *chi)  const {
  GA_Elem_divide_patch(g_a->mHandle, alo, ahi, g_b->mHandle, blo, bhi, 
		       mHandle, clo, chi);
}

void 
GA::GlobalArray::elemMaximumPatch(const GA::GlobalArray * g_a,
				  int *alo,int *ahi,
				  const GA::GlobalArray * g_b,
				  int *blo,int *bhi,
				  int *clo,int *chi)  const {
  GA_Elem_maximum_patch(g_a->mHandle, alo, ahi, g_b->mHandle, blo, bhi, 
			mHandle, clo, chi);
}

void 
GA::GlobalArray::elemMinimumPatch(const GA::GlobalArray * g_a,
				  int *alo,int *ahi,
				  const GA::GlobalArray * g_b,
				  int *blo,int *bhi,
				  int *clo,int *chi)  const {
  GA_Elem_minimum_patch(g_a->mHandle, alo, ahi, g_b->mHandle, blo, bhi, 
			mHandle, clo, chi);
}



/*Added by Limin for matrix operations*/

void 
GA::GlobalArray::shiftDiagonal(void *c)  const {
  GA_Shift_diagonal(mHandle, c);
}

void 
GA::GlobalArray::setDiagonal(const GA::GlobalArray * g_v)  const {
  GA_Set_diagonal(mHandle, g_v->mHandle);
}

void 
GA::GlobalArray::zeroDiagonal()  const {
  GA_Zero_diagonal(mHandle);
}

void 
GA::GlobalArray::addDiagonal(const GA::GlobalArray * g_v)  const {
  GA_Add_diagonal(mHandle, g_v->mHandle);
}

void 
GA::GlobalArray::getDiagonal(const GA::GlobalArray * g_a)  const {
  GA_Get_diag(g_a->mHandle, mHandle);
}

void 
GA::GlobalArray::scaleRows(const GA::GlobalArray * g_v)  const {
  GA_Scale_rows(mHandle, g_v->mHandle);
}

void 
GA::GlobalArray::scaleCols(const GA::GlobalArray * g_v)  const {
  GA_Scale_cols(mHandle, g_v->mHandle);
}

void 
GA::GlobalArray::norm1(double *nm)  const {
 GA_Norm1(mHandle, nm);
}

void 
GA::GlobalArray::normInfinity(double *nm)  const {
 GA_Norm_infinity(mHandle, nm);
}

void 
GA::GlobalArray::median(const GA::GlobalArray * g_a, 
			const GA::GlobalArray * g_b, 
			const GA::GlobalArray * g_c)  const {
  GA_Median(g_a->mHandle, g_b->mHandle, g_c->mHandle, mHandle);
}

void 
GA::GlobalArray::medianPatch(const GA::GlobalArray * g_a, int *alo, int *ahi, 
			     const GA::GlobalArray * g_b, int *blo, int *bhi, 
			     const GA::GlobalArray * g_c, int *clo, int *chi, 
			     int *mlo, int *mhi)  const {
  GA_Median_patch(g_a->mHandle, alo, ahi, g_b->mHandle, blo, bhi, 
		  g_c->mHandle, clo, chi, mHandle, mlo, mhi);
}




