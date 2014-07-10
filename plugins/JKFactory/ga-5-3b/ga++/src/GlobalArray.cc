#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include "ga++.h"

#define GA_DATA_TYPE     C_DBL
#define INVALID_HANDLE   -1000

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
                             int chunk[], GA::PGroup * p_handle) {
  mHandle = NGA_Create_config(type, ndim, dims, arrayname, chunk,
                              p_handle->handle());
  if(!mHandle)  GA_Error((char *)" GA creation failed",0);
}

GA::GlobalArray::GlobalArray(int type, int ndim, int64_t dims[],
                             char *arrayname, int64_t chunk[]) {
  mHandle = NGA_Create64(type, ndim, dims, arrayname, chunk);
  if(!mHandle)  GA_Error((char *)" GA creation failed",0);
}

GA::GlobalArray::GlobalArray(int type, int ndim, int64_t dims[],
                             char *arrayname, int64_t chunk[],
                             GA::PGroup * p_handle) {
  mHandle = NGA_Create_config64(type, ndim, dims, arrayname, chunk,
                                p_handle->handle());
  if(!mHandle)  GA_Error((char *)" GA creation failed",0);
}

GA::GlobalArray::GlobalArray(int type, int ndim, int dims[], char *arrayname,
                             int block[], int maps[]) {
  mHandle = NGA_Create_irreg(type, ndim, dims, arrayname, block, maps);
  if(!mHandle) GA_Error((char *)" GA creation failed",0);
}

GA::GlobalArray::GlobalArray(int type, int ndim, int dims[], char *arrayname,
                             int block[], int maps[], GA::PGroup * p_handle) {
  mHandle = NGA_Create_irreg_config(type, ndim, dims, arrayname, block, maps,
                                    p_handle->handle());
  if(!mHandle) GA_Error((char *)" GA creation failed",0);
}

GA::GlobalArray::GlobalArray(int type, int ndim, int64_t dims[],
                             char *arrayname, int64_t block[],
                             int64_t maps[]) {
  mHandle = NGA_Create_irreg64(type, ndim, dims, arrayname, block, maps);
  if(!mHandle) GA_Error((char *)" GA creation failed",0);
}

GA::GlobalArray::GlobalArray(int type, int ndim, int64_t dims[],
                             char *arrayname, int64_t block[], int64_t maps[],
                             GA::PGroup * p_handle) {
  mHandle = NGA_Create_irreg_config64(type, ndim, dims, arrayname, block, maps,
                                      p_handle->handle());
  if(!mHandle) GA_Error((char *)" GA creation failed",0);
}

GA::GlobalArray::GlobalArray(int type, int ndim, int dims[], int width[], 
			     char *arrayname, int chunk[], char ghosts) {
  mHandle = NGA_Create_ghosts(type, ndim, dims, width, arrayname, chunk);
  if(!mHandle)  GA_Error((char *)" GA creation failed",0);
}

GA::GlobalArray::GlobalArray(int type, int ndim, int64_t dims[], int64_t width[], 
			     char *arrayname, int64_t chunk[], char ghosts) {
  mHandle = NGA_Create_ghosts64(type, ndim, dims, width, arrayname, chunk);
  if(!mHandle)  GA_Error((char *)" GA creation failed",0);
}

GA::GlobalArray::GlobalArray(int type, int ndim, int dims[], int width[], 
			     char *arrayname, int chunk[],
                             GA::PGroup * p_handle, char ghosts) {
  mHandle = NGA_Create_ghosts_config(type, ndim, dims, width, arrayname, chunk,
                                     p_handle->handle());
  if(!mHandle)  GA_Error((char *)" GA creation failed",0);
}

GA::GlobalArray::GlobalArray(int type, int ndim, int64_t dims[], int64_t width[], 
			     char *arrayname, int64_t chunk[],
                             GA::PGroup * p_handle, char ghosts) {
  mHandle = NGA_Create_ghosts_config64(type, ndim, dims, width, arrayname, chunk,
                                     p_handle->handle());
  if(!mHandle)  GA_Error((char *)" GA creation failed",0);
}

GA::GlobalArray::GlobalArray(int type, int ndim, int dims[], int width[], 
			     char *arrayname, int block[], int maps[], 
			     char ghosts) {
  mHandle = NGA_Create_ghosts_irreg(type, ndim, dims, width, arrayname, block, 
				    maps);
  if(!mHandle) GA_Error((char *)" GA creation failed",0);
}

GA::GlobalArray::GlobalArray(int type, int ndim, int64_t dims[], int64_t width[], 
			     char *arrayname, int64_t block[], int64_t maps[], 
			     char ghosts) {
  mHandle = NGA_Create_ghosts_irreg64(type, ndim, dims, width, arrayname, block, 
				    maps);
  if(!mHandle) GA_Error((char *)" GA creation failed",0);
}

GA::GlobalArray::GlobalArray(int type, int ndim, int dims[], int width[], 
			     char *arrayname, int block[], int maps[],
                             GA::PGroup * p_handle,
			     char ghosts) {
  mHandle = NGA_Create_ghosts_irreg_config(type, ndim, dims, width, arrayname,
                                           block, maps, p_handle->handle());
  if(!mHandle) GA_Error((char *)" GA creation failed",0);
}

GA::GlobalArray::GlobalArray(int type, int ndim, int64_t dims[], int64_t width[], 
			     char *arrayname, int64_t block[], int64_t maps[],
                             GA::PGroup * p_handle,
			     char ghosts) {
  mHandle = NGA_Create_ghosts_irreg_config64(type, ndim, dims, width, arrayname,
                                           block, maps, p_handle->handle());
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
  mHandle = GA_Create_handle();
  if(!mHandle) GA_Error((char *)" GA creation failed",0);
}

GA::GlobalArray::~GlobalArray() {
  GA_Destroy(mHandle); 
  mHandle = INVALID_HANDLE;
}


/*********************************************************************
 *          Member functions of GA::GlobalArray                          *
 *********************************************************************/

void 
GA::GlobalArray::acc(int lo[], int hi[], void *buf, int ld[], void *alpha) const {
  NGA_Acc(mHandle, lo, hi, buf, ld, alpha);
}

void 
GA::GlobalArray::acc(int64_t lo[], int64_t hi[], void *buf, int64_t ld[], void *alpha) const {
  NGA_Acc64(mHandle, lo, hi, buf, ld, alpha);
}

void 
GA::GlobalArray::access(int lo[], int hi[], void *ptr, int ld[]) const {
  NGA_Access(mHandle, lo, hi, ptr, ld);
}

void 
GA::GlobalArray::access(int64_t lo[], int64_t hi[], void *ptr,
                        int64_t ld[]) const {
  NGA_Access64(mHandle, lo, hi, ptr, ld);
}

void
GA::GlobalArray::accessBlock(int idx, void *ptr, int ld[]) const {
    NGA_Access_block(mHandle, idx, ptr, ld);
}

void
GA::GlobalArray::accessBlock(int64_t idx, void *ptr, int64_t ld[]) const {
    NGA_Access_block64(mHandle, idx, ptr, ld);
}

void
GA::GlobalArray::accessBlockGrid(int index[], void *ptr, int ld[]) const {
    NGA_Access_block_grid(mHandle, index, ptr, ld);
}

void
GA::GlobalArray::accessBlockGrid(int64_t index[], void *ptr, int64_t ld[]) const {
    NGA_Access_block_grid64(mHandle, index, ptr, ld);
}

void
GA::GlobalArray::accessBlockSegment(int proc, void *ptr, int *len) const {
    NGA_Access_block_segment(mHandle, proc, ptr, len);
}

void
GA::GlobalArray::accessBlockSegment(int proc, void *ptr, int64_t *len) const {
    NGA_Access_block_segment64(mHandle, proc, ptr, len);
}

void 
GA::GlobalArray::accessGhosts(int dims[], void *ptr, int ld[]) const {
  NGA_Access_ghosts(mHandle, dims, ptr, ld);
}

void 
GA::GlobalArray::accessGhosts(int64_t dims[], void *ptr, int64_t ld[]) const {
  NGA_Access_ghosts64(mHandle, dims, ptr, ld);
}

void 
GA::GlobalArray::accessGhostElement(void *ptr, int subscript[], int ld[]) const {
  NGA_Access_ghost_element(mHandle, ptr, subscript, ld);
}

void 
GA::GlobalArray::accessGhostElement(void *ptr, int64_t subscript[], int64_t ld[]) const {
  NGA_Access_ghost_element64(mHandle, ptr, subscript, ld);
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
GA::GlobalArray::addPatch (void *alpha, 
			   const GA::GlobalArray * g_a, int64_t alo[], int64_t ahi[],
			   void *beta,  
			   const GA::GlobalArray * g_b, int64_t blo[], int64_t bhi[],
			   int64_t clo[], int64_t chi[]) const {
  NGA_Add_patch64(alpha, g_a->mHandle, alo, ahi, beta, 
		g_b->mHandle, blo, bhi, mHandle, clo, chi);
}

int
GA::GlobalArray::allocate() const {
  return GA_Allocate(mHandle);
}

void  
GA::GlobalArray::allocGatscatBuf(int nelems) const {
  NGA_Alloc_gatscat_buf(nelems);
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

void 
GA::GlobalArray::copyPatch(char trans, const GA::GlobalArray* ga, int64_t alo[], 
			   int64_t ahi[], int64_t blo[], int64_t bhi[]) const {
  NGA_Copy_patch64(trans, ga->mHandle, alo, ahi, mHandle, blo, bhi);
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

double 
GA::GlobalArray::ddotPatch(char ta, int64_t alo[], int64_t ahi[],
			   const GA::GlobalArray * g_a, char tb, 
			   int64_t blo[], int64_t bhi[]) const {
  return NGA_Ddot_patch64(mHandle, ta, alo, ahi, g_a->mHandle, tb, blo, bhi);
}

void 
GA::GlobalArray::destroy() {
  GA_Destroy(mHandle);
  mHandle = INVALID_HANDLE;
}

void 
GA::GlobalArray::dgemm(char ta, char tb, int m, int n, int k, double alpha,  
		       const GA::GlobalArray *g_a, const GA::GlobalArray *g_b, 
		       double beta) const {
  GA_Dgemm(ta, tb, m, n, k, alpha, g_a->mHandle, g_b->mHandle, 
	   beta, mHandle);
}

void 
GA::GlobalArray::dgemm(char ta, char tb, int64_t m, int64_t n, int64_t k, double alpha,  
		       const GA::GlobalArray *g_a, const GA::GlobalArray *g_b, 
		       double beta) const {
  GA_Dgemm64(ta, tb, m, n, k, alpha, g_a->mHandle, g_b->mHandle, 
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

void 
GA::GlobalArray::distribution(int me, int64_t* lo, int64_t* hi) const {
  NGA_Distribution64(mHandle, me, lo, hi);
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

float 
GA::GlobalArray::fdotPatch(char t_a, int64_t alo[], int64_t ahi[], 
			   const GA::GlobalArray * g_b, char t_b, int64_t blo[], 
			   int64_t bhi[]) const {
  return NGA_Fdot_patch64(mHandle, t_a, alo, ahi,
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
GA::GlobalArray::fillPatch (int64_t lo[], int64_t hi[], void *val) const {
  NGA_Fill_patch64(mHandle, lo, hi, val);
}

void  
GA::GlobalArray::freeGatscatBuf() {
  NGA_Free_gatscat_buf();
}

void 
GA::GlobalArray::gather(void *v, int * subsarray[], int n) const {
  NGA_Gather(mHandle, v, subsarray, n);
}

void 
GA::GlobalArray::gather(void *v, int64_t * subsarray[], int64_t n) const {
  NGA_Gather64(mHandle, v, subsarray, n);
}

void 
GA::GlobalArray::get(int lo[], int hi[], void *buf, int ld[])  const {
  NGA_Get(mHandle, lo, hi, buf, ld);
}

void 
GA::GlobalArray::get(int64_t lo[], int64_t hi[], void *buf,
                     int64_t ld[]) const {
  NGA_Get64(mHandle, lo, hi, buf, ld);
}

void
GA::GlobalArray::getBlockInfo(int num_blocks[], int block_dims[]) {
  GA_Get_block_info(mHandle, num_blocks, block_dims);
}

int 
GA::GlobalArray::hasGhosts()  const {
  return GA_Has_ghosts(mHandle);
}

int 
GA::GlobalArray::idot(const GA::GlobalArray * g_a)  const {
  return GA_Idot(mHandle, g_a->mHandle);
}

long 
GA::GlobalArray::idotPatch(char ta, int alo[], int ahi[],
			   const GA::GlobalArray * g_a, char tb, 
			   int blo[], int bhi[])  const {
  return NGA_Idot_patch(mHandle, ta, alo, ahi, g_a->mHandle, tb, blo, bhi);
}

long 
GA::GlobalArray::idotPatch(char ta, int64_t alo[], int64_t ahi[],
			   const GA::GlobalArray * g_a, char tb, 
			   int64_t blo[], int64_t bhi[])  const {
  return NGA_Idot_patch64(mHandle, ta, alo, ahi, g_a->mHandle, tb, blo, bhi);
}

void 
GA::GlobalArray::inquire(int *type, int *ndim, int dims[])  const {
  NGA_Inquire(mHandle, type, ndim, dims);
}

void 
GA::GlobalArray::inquire(int *type, int *ndim, int64_t dims[])  const {
  NGA_Inquire64(mHandle, type, ndim, dims);
}

char* 
GA::GlobalArray::inquireName()  const {
  return GA_Inquire_name(mHandle);
}

long 
GA::GlobalArray::ldot(const GA::GlobalArray * g_a)  const {
  return GA_Ldot(mHandle, g_a->mHandle);
}

long 
GA::GlobalArray::ldotPatch(char ta, int alo[], int ahi[],
			   const GA::GlobalArray * g_a, char tb, 
			   int blo[], int bhi[])  const {
  return NGA_Ldot_patch(mHandle, ta, alo, ahi, g_a->mHandle, tb, blo, bhi);
}

long 
GA::GlobalArray::ldotPatch(char ta, int64_t alo[], int64_t ahi[],
			   const GA::GlobalArray * g_a, char tb, 
			   int64_t blo[], int64_t bhi[])  const {
  return NGA_Ldot_patch64(mHandle, ta, alo, ahi, g_a->mHandle, tb, blo, bhi);
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
GA::GlobalArray::locate(int64_t subscript[])  const {
  return NGA_Locate64(mHandle, subscript);
} 

int 
GA::GlobalArray::locateRegion(int lo[], int hi[], int map[], int procs[])  const {
  return NGA_Locate_region(mHandle, lo, hi, map, procs);
}

int 
GA::GlobalArray::locateRegion(int64_t lo[], int64_t hi[], int64_t map[], int procs[])  const {
  return NGA_Locate_region64(mHandle, lo, hi, map, procs);
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
GA::GlobalArray::matmulPatch(char transa, char transb, void* alpha, 
			     void *beta, const GA::GlobalArray *g_a, 
			     int64_t ailo, int64_t aihi, int64_t ajlo, int64_t ajhi,
			     const GA::GlobalArray *g_b, 
			     int64_t bilo, int64_t bihi, int64_t bjlo, int64_t bjhi,
			     int64_t cilo, int64_t cihi, int64_t cjlo, int64_t cjhi)  const {
  GA_Matmul_patch64(transa, transb, alpha, beta, 
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
GA::GlobalArray::matmulPatch(char transa, char transb, void* alpha, void *beta,
			     const GA::GlobalArray *g_a, int64_t *alo, int64_t *ahi,
			     const GA::GlobalArray *g_b, int64_t *blo, int64_t *bhi,
			     int64_t *clo, int64_t *chi)  const {
  NGA_Matmul_patch64(transa, transb, alpha, beta, g_a->mHandle, alo, ahi, 
		  g_b->mHandle, blo, bhi, mHandle, clo, chi);
}

void
GA::GlobalArray::mergeDistrPatch(int alo[], int ahi[], GlobalArray *g_b,
                                 int blo[], int bhi[]) {
  NGA_Merge_distr_patch(mHandle, alo, ahi, g_b->mHandle, blo, bhi);
}

void
GA::GlobalArray::mergeDistrPatch(int64_t alo[], int64_t ahi[], GlobalArray *g_b,
                                 int64_t blo[], int64_t bhi[]) {
  NGA_Merge_distr_patch64(mHandle, alo, ahi, g_b->mHandle, blo, bhi);
}

int
GA::GlobalArray::isMirrored() {
  return GA_Is_mirrored(mHandle);
}

void
GA::GlobalArray::mergeMirrored() {
  GA_Merge_mirrored(mHandle);
}

void
GA::GlobalArray::nbAcc(int lo[], int hi[], void *buf, int ld[], void *alpha,
                       GANbhdl *nbhandle) {
  NGA_NbAcc(mHandle, lo, hi, buf, ld, alpha, nbhandle);
}

void
GA::GlobalArray::nbAcc(int64_t lo[], int64_t hi[], void *buf, int64_t ld[], void *alpha,
                       GANbhdl *nbhandle) {
  NGA_NbAcc64(mHandle, lo, hi, buf, ld, alpha, nbhandle);
}

void
GA::GlobalArray::nbGet(int lo[], int hi[], void *buf, int ld[], GANbhdl *nbhandle) {
  NGA_NbGet(mHandle, lo, hi, buf, ld, nbhandle);
}

void
GA::GlobalArray::nbGet(int64_t lo[], int64_t hi[], void *buf, int64_t ld[], GANbhdl *nbhandle) {
  NGA_NbGet64(mHandle, lo, hi, buf, ld, nbhandle);
}

void
GA::GlobalArray::nbGetGhostDir(int mask[], GANbhdl *nbhandle)
{
  NGA_NbGet_ghost_dir(mHandle, mask, nbhandle);
}

void
GA::GlobalArray::nbGetGhostDir(int64_t mask[], GANbhdl *nbhandle)
{
  NGA_NbGet_ghost_dir64(mHandle, mask, nbhandle);
}

void 
GA::GlobalArray::nblock(int numblock[])  const {
  GA_Nblock(mHandle, numblock);
}

void
GA::GlobalArray::nbPut(int lo[], int hi[], void *buf, int ld[], GANbhdl *nbhandle) {
  NGA_NbPut(mHandle, lo, hi, buf, ld, nbhandle);
}

void
GA::GlobalArray::nbPut(int64_t lo[], int64_t hi[], void *buf, int64_t ld[], GANbhdl *nbhandle) {
  NGA_NbPut64(mHandle, lo, hi, buf, ld, nbhandle);
}

int 
GA::GlobalArray::ndim()  const {
  return GA_Ndim(mHandle);
}

void
GA::GlobalArray::pack(const GA::GlobalArray *g_dest,
                      const GA::GlobalArray *g_mask,
                      int lo, int hi, int *icount) const {
    GA_Pack(mHandle, g_dest->mHandle, g_mask->mHandle, lo, hi, icount);
}

void
GA::GlobalArray::pack(const GA::GlobalArray *g_dest,
                      const GA::GlobalArray *g_mask,
                      int64_t lo, int64_t hi, int64_t *icount) const {
    GA_Pack64(mHandle, g_dest->mHandle, g_mask->mHandle, lo, hi, icount);
}

void
GA::GlobalArray::patchEnum(int lo, int hi, void *start, void *inc) {
  GA_Patch_enum(mHandle, lo, hi, start, inc);
}

void
GA::GlobalArray::patchEnum(int64_t lo, int64_t hi, void *start, void *inc) {
  GA_Patch_enum64(mHandle, lo, hi, start, inc);
}

void 
GA::GlobalArray::periodicAcc(int lo[], int hi[], void* buf, int ld[], 
			     void* alpha)  const {
  NGA_Periodic_acc(mHandle, lo, hi, buf, ld, alpha);
}

void 
GA::GlobalArray::periodicAcc(int64_t lo[], int64_t hi[], void* buf, int64_t ld[], 
			     void* alpha)  const {
  NGA_Periodic_acc64(mHandle, lo, hi, buf, ld, alpha);
}

void 
GA::GlobalArray::periodicGet(int lo[], int hi[], void* buf, int ld[])  const {
  NGA_Periodic_get(mHandle, lo, hi, buf, ld);
}

void 
GA::GlobalArray::periodicGet(int64_t lo[], int64_t hi[], void* buf, int64_t ld[])  const {
  NGA_Periodic_get64(mHandle, lo, hi, buf, ld);
}

void 
GA::GlobalArray::periodicPut(int lo[], int hi[], void* buf, int ld[])  const {
  NGA_Periodic_put(mHandle, lo, hi, buf, ld);
}

void 
GA::GlobalArray::periodicPut(int64_t lo[], int64_t hi[], void* buf, int64_t ld[])  const {
  NGA_Periodic_put64(mHandle, lo, hi, buf, ld);
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
GA::GlobalArray::printPatch(int64_t* lo, int64_t* hi, int pretty)  const {
  NGA_Print_patch64(mHandle, lo, hi, pretty);   
}

void 
GA::GlobalArray::procTopology(int proc, int coord[])  const {
  NGA_Proc_topology(mHandle, proc, coord);
}

void 
GA::GlobalArray::put(int lo[], int hi[], void *buf, int ld[])  const {
  NGA_Put(mHandle, lo, hi, buf, ld);
}

void 
GA::GlobalArray::put(int64_t lo[], int64_t hi[], void *buf,
                     int64_t ld[])  const {
  NGA_Put64(mHandle, lo, hi, buf, ld);
}

long 
GA::GlobalArray::readInc(int subscript[], long inc)  const {
  return NGA_Read_inc(mHandle, subscript, inc);
}

long 
GA::GlobalArray::readInc(int64_t subscript[], long inc)  const {
  return NGA_Read_inc64(mHandle, subscript, inc);
}

void 
GA::GlobalArray::release(int lo[], int hi[])  const {
  NGA_Release(mHandle, lo, hi);
}

void 
GA::GlobalArray::release(int64_t lo[], int64_t hi[])  const {
  NGA_Release64(mHandle, lo, hi);
}

void
GA::GlobalArray::releaseBlock(int idx) const {    
  NGA_Release_block(mHandle, idx);
}

void
GA::GlobalArray::releaseBlockGrid(int index[]) const {    
  NGA_Release_block_grid(mHandle, index);
}


void
GA::GlobalArray::releaseBlockSegment(int proc) const {    
  NGA_Release_block_segment(mHandle, proc);
}

void
GA::GlobalArray::releaseGhosts() const {    
  NGA_Release_ghosts(mHandle);
}

void
GA::GlobalArray::releaseGhostElement(int subscript[]) const {    
  NGA_Release_ghost_element(mHandle, subscript);
}

void
GA::GlobalArray::releaseGhostElement(int64_t subscript[]) const {    
  NGA_Release_ghost_element64(mHandle, subscript);
}

void 
GA::GlobalArray::releaseUpdate(int lo[], int hi[])  const {
  NGA_Release_update(mHandle, lo, hi);
}

void 
GA::GlobalArray::releaseUpdate(int64_t lo[], int64_t hi[])  const {
  NGA_Release_update64(mHandle, lo, hi);
}

void 
GA::GlobalArray::releaseUpdateBlock(int idx)  const {
  NGA_Release_update_block(mHandle, idx);
}

void 
GA::GlobalArray::releaseUpdateBlockGrid(int index[])  const {
  NGA_Release_update_block_grid(mHandle, index);
}

void 
GA::GlobalArray::releaseUpdateBlockSegment(int idx)  const {
  NGA_Release_update_block_segment(mHandle, idx);
}

void
GA::GlobalArray::releaseUpdateGhosts() const {    
  NGA_Release_update_ghosts(mHandle);
}

void
GA::GlobalArray::releaseUpdateGhostElement(int subscript[]) const {    
  NGA_Release_update_ghost_element(mHandle, subscript);
}

void
GA::GlobalArray::releaseUpdateGhostElement(int64_t subscript[]) const {    
  NGA_Release_update_ghost_element64(mHandle, subscript);
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
GA::GlobalArray::scalePatch (int64_t lo[], int64_t hi[], void *val)  const {
  NGA_Scale_patch64(mHandle, lo, hi, val);
}

void
GA::GlobalArray::scanAdd(const GA::GlobalArray *g_dest,
                         const GlobalArray *g_mask,
                         int lo, int hi, int excl) const {
    GA_Scan_add(mHandle, g_dest->mHandle, g_mask->mHandle, lo, hi, excl);
}

void
GA::GlobalArray::scanAdd(const GA::GlobalArray *g_dest,
                         const GlobalArray *g_mask,
                         int64_t lo, int64_t hi, int excl) const {
    GA_Scan_add64(mHandle, g_dest->mHandle, g_mask->mHandle, lo, hi, excl);
}

void
GA::GlobalArray::scanCopy(const GA::GlobalArray *g_dest,
                          const GA::GlobalArray *g_mask,
                          int lo, int hi) const {
    GA_Scan_copy(mHandle, g_dest->mHandle, g_mask->mHandle, lo, hi);
}

void
GA::GlobalArray::scanCopy(const GA::GlobalArray *g_dest,
                          const GA::GlobalArray *g_mask,
                          int64_t lo, int64_t hi) const {
    GA_Scan_copy64(mHandle, g_dest->mHandle, g_mask->mHandle, lo, hi);
}

void 
GA::GlobalArray::scatter(void *v, int *subsarray[], int n)  const {
  NGA_Scatter(mHandle, v, subsarray, n);
}

void 
GA::GlobalArray::scatter(void *v, int64_t *subsarray[], int64_t n)  const {
  NGA_Scatter64(mHandle, v, subsarray, n);
}

void 
GA::GlobalArray::scatterAcc(void *v, int *subsarray[], int n, void *alpha)  const {
  NGA_Scatter_acc(mHandle, v, subsarray, n, alpha);
}

void 
GA::GlobalArray::scatterAcc(void *v, int64_t *subsarray[], int64_t n, void
    *alpha)  const {
  NGA_Scatter_acc64(mHandle, v, subsarray, n, alpha);
}

void 
GA::GlobalArray::selectElem(char *op, void* val, int index[])  const {
  NGA_Select_elem(mHandle, op, val, index);
}

void 
GA::GlobalArray::selectElem(char *op, void* val, int64_t index[])  const {
  NGA_Select_elem64(mHandle, op, val, index);
}

void
GA::GlobalArray::setArrayName(char *name) const {
    GA_Set_array_name(mHandle, name);
}

void
GA::GlobalArray::setBlockCyclic(int dims[]) const {
    GA_Set_block_cyclic(mHandle, dims);
}

void
GA::GlobalArray::setBlockCyclicProcGrid(int dims[], int proc_grid[]) const{
    GA_Set_block_cyclic_proc_grid(mHandle, dims, proc_grid);    
}

void
GA::GlobalArray::setChunk(int chunk[]) const {
    GA_Set_chunk(mHandle, chunk);    
}

void
GA::GlobalArray::setChunk(int64_t chunk[]) const {
    GA_Set_chunk64(mHandle, chunk);    
}

void
GA::GlobalArray::setData(int ndim, int dims[], int type) const {
    GA_Set_data(mHandle, ndim, dims, type);
}

void
GA::GlobalArray::setData(int ndim, int64_t dims[], int type) const {
    GA_Set_data64(mHandle, ndim, dims, type);
}

void
GA::GlobalArray::setGhosts(int width[]) const {
    GA_Set_ghosts(mHandle, width);
}

void
GA::GlobalArray::setGhosts(int64_t width[]) const {
    GA_Set_ghosts64(mHandle, width);
}

void
GA::GlobalArray::setIrregDistr(int mapc[], int nblock[]) const {
    GA_Set_irreg_distr(mHandle, mapc, nblock);
}

void
GA::GlobalArray::setIrregDistr(int64_t mapc[], int64_t nblock[]) const {
    GA_Set_irreg_distr64(mHandle, mapc, nblock);
}

void
GA::GlobalArray::setRestricted(int list[], int nprocs) const {
    GA_Set_restricted(mHandle, list, nprocs);
}

void
GA::GlobalArray::setRestrictedRange(int lo_proc, int hi_proc) const {
    GA_Set_restricted_range(mHandle, lo_proc, hi_proc);
}

void
GA::GlobalArray::setPGroup(GA::PGroup *pHandle) const {
    GA_Set_pgroup(mHandle, pHandle->handle());
}

void 
GA::GlobalArray::sgemm(char ta, char tb, int m, int n, int k, float alpha,  
		       const GA::GlobalArray *g_a, const GA::GlobalArray *g_b, 
		       float beta)  const {
  GA_Sgemm(ta, tb, m, n, k, alpha, g_a->mHandle, g_b->mHandle, 
	   beta, mHandle);
}

void 
GA::GlobalArray::sgemm(char ta, char tb, int64_t m, int64_t n, int64_t k, float alpha,  
		       const GA::GlobalArray *g_a, const GA::GlobalArray *g_b, 
		       float beta)  const {
  GA_Sgemm64(ta, tb, m, n, k, alpha, g_a->mHandle, g_b->mHandle, 
	   beta, mHandle);
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
GA::GlobalArray::stridedAcc(int lo[], int hi[], int skip[], void*buf,
                            int ld[], void *alpha) const {
    NGA_Strided_acc(mHandle, lo, hi, skip, buf, ld, alpha);
}

void
GA::GlobalArray::stridedAcc(int64_t lo[], int64_t hi[], int64_t skip[], void*buf,
                            int64_t ld[], void *alpha) const {
    NGA_Strided_acc64(mHandle, lo, hi, skip, buf, ld, alpha);
}

void
GA::GlobalArray::stridedGet(int lo[], int hi[], int skip[], void*buf,
                            int ld[]) const {
    NGA_Strided_get(mHandle, lo, hi, skip, buf, ld);
}

void
GA::GlobalArray::stridedGet(int64_t lo[], int64_t hi[], int64_t skip[], void*buf,
                            int64_t ld[]) const {
    NGA_Strided_get64(mHandle, lo, hi, skip, buf, ld);
}

void
GA::GlobalArray::stridedPut(int lo[], int hi[], int skip[], void*buf,
                            int ld[]) const {
    NGA_Strided_put(mHandle, lo, hi, skip, buf, ld);
}

void
GA::GlobalArray::stridedPut(int64_t lo[], int64_t hi[], int64_t skip[], void*buf,
                            int64_t ld[]) const {
    NGA_Strided_put64(mHandle, lo, hi, skip, buf, ld);
}

void
GA::GlobalArray::summarize(int verbose) const {
    GA_Summarize(verbose);
}

void 
GA::GlobalArray::symmetrize()  const { 
  GA_Symmetrize(mHandle); 
}

int
GA::GlobalArray::totalBlocks() const 
{
    return GA_Total_blocks(mHandle);
}

void 
GA::GlobalArray::transpose(const GA::GlobalArray * g_a)  const {
  GA_Transpose(mHandle, g_a->mHandle);
}

void
GA::GlobalArray::unpack(GlobalArray *g_dest, GlobalArray *g_mask,
                        int lo, int hi, int *icount) const {
    GA_Unpack(mHandle, g_dest->mHandle, g_mask->mHandle, lo, hi, icount);
}    

void
GA::GlobalArray::unpack(GlobalArray *g_dest, GlobalArray *g_mask,
                        int64_t lo, int64_t hi, int64_t *icount) const {
    GA_Unpack64(mHandle, g_dest->mHandle, g_mask->mHandle, lo, hi, icount);
}    

void 
GA::GlobalArray::updateGhosts()  const {
  GA_Update_ghosts(mHandle);
}

void
GA::GlobalArray::updateGhostsNb(GANbhdl *nbhandle) const{
  NGA_Update_ghosts_nb(mHandle, nbhandle);
}

int 
GA::GlobalArray::updateGhostDir(int dimension, int idir, int cflag)  const {
  return NGA_Update_ghost_dir(mHandle, dimension, idir, cflag);
}

void
GA::GlobalArray::getGhostBlock(int lo[], int hi[], void *buf, int ld[]) const {
  NGA_Get_ghost_block(mHandle, lo, hi, buf, ld);
}

void
GA::GlobalArray::getGhostBlock(int64_t lo[], int64_t hi[], void *buf, int64_t ld[]) const {
  NGA_Get_ghost_block64(mHandle, lo, hi, buf, ld);
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

DoubleComplex 
GA::GlobalArray::zdotPatch(char ta, int64_t alo[], int64_t ahi[],
			   const GA::GlobalArray * g_a, char tb, 
			   int64_t blo[], int64_t bhi[])  const {
  return NGA_Zdot_patch64(mHandle, ta, alo, ahi, g_a->mHandle, tb, blo, bhi);
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
GA::GlobalArray::zeroPatch (int64_t lo[], int64_t hi[])  const {
  NGA_Zero_patch64(mHandle, lo, hi);
}

void 
GA::GlobalArray::zgemm(char ta, char tb, int m, int n, int k, 
		       DoubleComplex alpha,  
		       const GA::GlobalArray *g_a, const GA::GlobalArray *g_b, 
		       DoubleComplex beta)  const {
  GA_Zgemm(ta, tb, m, n, k, alpha, g_a->mHandle, g_b->mHandle, 
	   beta, mHandle);
}

void 
GA::GlobalArray::zgemm(char ta, char tb, int64_t m, int64_t n, int64_t k, 
		       DoubleComplex alpha,  
		       const GA::GlobalArray *g_a, const GA::GlobalArray *g_b, 
		       DoubleComplex beta)  const {
  GA_Zgemm64(ta, tb, m, n, k, alpha, g_a->mHandle, g_b->mHandle, 
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
GA::GlobalArray::absValuePatch(int64_t *lo, int64_t *hi)  const {
  GA_Abs_value_patch64(mHandle, lo, hi);
}

void 
GA::GlobalArray::addConstantPatch(int *lo,int *hi, void *alpha)  const {
  GA_Add_constant_patch(mHandle, lo, hi, alpha);
}

void 
GA::GlobalArray::addConstantPatch(int64_t *lo,int64_t *hi, void *alpha)  const {
  GA_Add_constant_patch64(mHandle, lo, hi, alpha);
}

void 
GA::GlobalArray::recipPatch(int *lo, int *hi)  const {
  GA_Recip_patch(mHandle, lo, hi);
}

void 
GA::GlobalArray::recipPatch(int64_t *lo, int64_t *hi)  const {
  GA_Recip_patch64(mHandle, lo, hi);
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
GA::GlobalArray::stepMaxPatch(int64_t *alo, int64_t *ahi, 
			      const GA::GlobalArray * g_b, int64_t *blo, int64_t *bhi, 
			      double *step)  const {
  GA_Step_max_patch64(mHandle, alo, ahi, g_b->mHandle, blo, bhi, step);
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
GA::GlobalArray::elemMultiplyPatch(const GA::GlobalArray * g_a,
				   int64_t *alo,int64_t *ahi,
				   const GA::GlobalArray * g_b,
				   int64_t *blo,int64_t *bhi,
				   int64_t *clo,int64_t *chi)  const {
  GA_Elem_multiply_patch64(g_a->mHandle, alo, ahi, g_b->mHandle, blo, bhi, 
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
GA::GlobalArray::elemDividePatch(const GA::GlobalArray * g_a,int64_t *alo,
             int64_t *ahi, const GA::GlobalArray * g_b, int64_t *blo,
             int64_t *bhi, int64_t *clo, int64_t *chi)  const {
  GA_Elem_divide_patch64(g_a->mHandle, alo, ahi, g_b->mHandle, blo, bhi, 
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
GA::GlobalArray::elemMaximumPatch(const GA::GlobalArray * g_a,
				  int64_t *alo, int64_t *ahi,
				  const GA::GlobalArray * g_b,
				  int64_t *blo, int64_t *bhi,
				  int64_t *clo, int64_t *chi)  const {
  GA_Elem_maximum_patch64(g_a->mHandle, alo, ahi, g_b->mHandle, blo, bhi, 
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

void 
GA::GlobalArray::elemMinimumPatch(const GA::GlobalArray * g_a,
				  int64_t *alo, int64_t *ahi,
				  const GA::GlobalArray * g_b,
				  int64_t *blo, int64_t *bhi,
				  int64_t *clo, int64_t *chi)  const {
  GA_Elem_minimum_patch64(g_a->mHandle, alo, ahi, g_b->mHandle, blo, bhi, 
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

void 
GA::GlobalArray::medianPatch(const GA::GlobalArray * g_a, int64_t *alo, int64_t *ahi, 
			     const GA::GlobalArray * g_b, int64_t *blo, int64_t *bhi, 
			     const GA::GlobalArray * g_c, int64_t *clo, int64_t *chi, 
			     int64_t *mlo, int64_t *mhi)  const {
  GA_Median_patch64(g_a->mHandle, alo, ahi, g_b->mHandle, blo, bhi, 
		  g_c->mHandle, clo, chi, mHandle, mlo, mhi);
}




