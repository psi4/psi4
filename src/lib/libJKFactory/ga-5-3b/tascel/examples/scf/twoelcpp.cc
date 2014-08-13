#include "ga.h"
#include "macdecls.h"
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <assert.h>

// #define min(x,y) (((x)<=(y))?(x):(y))

#define maxatom  286
#define maxnbfn  (15*maxatom)
#define mxiter   30
#define maxnnbfn ((maxnbfn)*(maxnbfn+1)/2)
#define pi       3.141592653589793d0)
#define tol      (0.5e-3)
#define tol2e    (1.0e-6)

#define ichunk 20

#ifdef __cplusplus
extern "C" {
#endif
  extern int armci_me;
  void nga_get_(Integer *g_a,
		Integer *lo,
		Integer *hi,
		void    *buf,
		Integer *ld);

  void g_(double *gg, Integer *i, Integer *j, Integer *k, Integer *l);
  double contract_matrices_(Integer *g_a, Integer *g_b);
  void new_nga_acc_(Integer *g_a,
		    Integer *lo,
		    Integer *hi,
		    void    *buf,
		    Integer *ld,
		    void    *alpha);

  void nga_acc_(Integer *g_a,
		Integer *lo,
		Integer *hi,
		void    *buf,
		Integer *ld,
		void    *alpha);

  void ctwoel_(Integer *g_schwarz, Integer *g_dens, Integer *g_fock,
	       double *schwmax, double *etwo,
	       long long *pnbfn, long long *icut1, long long *icut2,
	       long long *icut3, long long *icut4);

#ifdef __cplusplus
}
#endif


static void cpptwoel(Integer g_schwarz, Integer g_dens, Integer g_fock,
		     double schwmax, Integer nbfn,
		     long long *icut1, long long *icut2,
		     long long *icut3, long long *icut4);

void ctwoel_(Integer *g_schwarz, Integer *g_dens, Integer *g_fock,
	     double *schwmax, double *etwo,
	     long long *pnbfn, long long *icut1, long long *icut2,
	     long long *icut3, long long *icut4) {
  const Integer nbfn = (Integer)*pnbfn;
  long long ijcnt,klcnt,ijklcnt;

  ijcnt = *icut1;
  klcnt = *icut2;
  ijklcnt = *icut3;

  cpptwoel(*g_schwarz, *g_dens, *g_fock, *schwmax, nbfn, 
	   icut1, icut2, icut3, icut4);

  *etwo = 0.5*contract_matrices_(g_fock,g_dens);

  ijcnt = *icut1 - ijcnt;
  klcnt = *icut2 - klcnt;
  ijklcnt = (*icut3) - ijklcnt;
  *icut4 = *icut3;
  if (*icut3 <= 0L) {
    /*     no integrals may be calculated if there is no work for  */
    /*     this node (ichunk too big), or, something is wrong */
    printf("no two-electron integrals computed by node %d\n", GA_Nodeid());
    fflush(stdout);
  }
  return;
}

/*************************************************************************
 *
 * the code involving task pool stuff follows
 *
 *************************************************************************/

//#include <UniformTaskCollectionSplit.h>
#include <UniformTaskCollSplitData.h>
#include <DenseArray.h>
#include <algorithm>

using namespace std;
using namespace tascel;

static int next_4chunk(int g_counter, Integer nbfn, 
		       Integer lo[4], Integer hi[4], 
		       Integer *ilo, Integer *jlo,
		       Integer *klo, Integer *llo) {
  const long one=1;
  const int zero=0;
  long imax, itask;
  Integer itmp;
  long long total_ntsks;
  imax = nbfn/ichunk;
  if (nbfn - ichunk*imax > 0) {
    imax = imax + 1;
  }
  total_ntsks = imax*imax*imax*imax;

  itask = NGA_Read_inc(g_counter, (int *)&zero, one);
/*   printf("%d: itask=%d total_ntsks=%lf\n", armci_me, (int)itask, (double)total_ntsks); */
  if (itask < 0) {
    printf("next_4chunk: itask negative: %d imax:%d nbfn:%d ichunk:%d\n",
	   (int) itask, (int)imax, (int)nbfn, (int)ichunk);
    printf("probable GA integer precision problem if %lf > 2^31\n", 
	   pow(imax*1.0, 4));
    fflush(stdout);
    assert(0);
  }
  if (itask < total_ntsks) {
    itmp = itask;
    *ilo = itmp % imax;
    itmp /= imax;
    *jlo = itmp %imax;
    itmp /= imax;
    *klo = itmp % imax;
    itmp /= imax;
    assert(itmp < imax);
    *llo = itmp;

    lo[0] = *ilo*ichunk;
    lo[1] = *jlo*ichunk;
    lo[2] = *klo*ichunk;
    lo[3] = *llo*ichunk;
    hi[0] = min((*ilo+1)*ichunk,nbfn) - 1;
    hi[1] = min((*jlo+1)*ichunk,nbfn) - 1;
    hi[2] = min((*klo+1)*ichunk,nbfn) - 1;
    hi[3] = min((*llo+1)*ichunk,nbfn) - 1;

    return 1;
  }
  else {
    return 0;
  }
}

static void clean_chunk(double *chunk) {
  int i;
  for(i=0; i<ichunk*ichunk; i++) {
    chunk[i] = 0.0;
  }
}

typedef struct {
  int id;
} task_dscr_t;

typedef struct {
  Integer nbfn;
  double schwmax;
  int cnt1, cnt2;
  long long *icut1, *icut2, *icut3, *icut4;
  Integer g_schwarz, g_dens, g_fock;
} task_plo_t;

static void compute_lo_hi(long itask, long nbfn, long lo[4], long hi[4], 
			 long *ilo, long *jlo, long *klo, long *llo) {
  long imax;
  Integer itmp;
  long long total_ntsks;
  imax = nbfn/ichunk;
  if (nbfn - ichunk*imax > 0) {
    imax = imax + 1;
  }
  total_ntsks = imax*imax*imax*imax;

  if (itask < 0) {
    printf("next_4chunk: itask negative: %d imax:%d nbfn:%d ichunk:%d\n",
	   (int) itask, (int)imax, (int)nbfn, (int)ichunk);
    printf("probable GA integer precision problem if %lf > 2^31\n", 
	   pow(imax*1.0, 4));
    fflush(stdout);
    assert(0);
  }
  assert(itask < total_ntsks);

  itmp = itask;
  *ilo = itmp % imax;
  itmp /= imax;
  *jlo = itmp %imax;
  itmp /= imax;
  *klo = itmp % imax;
  itmp /= imax;
  assert(itmp < imax);
  *llo = itmp;
  
  lo[0] = *ilo*ichunk;
  lo[1] = *jlo*ichunk;
  lo[2] = *klo*ichunk;
  lo[3] = *llo*ichunk;
  hi[0] = min((*ilo+1)*ichunk,nbfn) - 1;
  hi[1] = min((*jlo+1)*ichunk,nbfn) - 1;
  hi[2] = min((*klo+1)*ichunk,nbfn) - 1;
  hi[3] = min((*llo+1)*ichunk,nbfn) - 1;  
}

static int compute_owner(long itask, const vector<DataColl*> &colls, int arrid, long nbfn)
{
  const int size = 2*sizeof(int);
  char idx[size];
  long int ilo, jlo, klo, llo, lo[4], hi[4], tmp;
  int *iidx = (int *)idx;

  compute_lo_hi(itask, nbfn, lo, hi, &ilo, &jlo, &klo, &llo);

  switch(arrid) {
  case 0:
    iidx[0] = lo[1]/ichunk;
    iidx[1] = lo[0]/ichunk;
    break;
  case 1:
    iidx[0] = lo[3]/ichunk;
    iidx[1] = lo[2]/ichunk;
    break;
  case 2:
    iidx[0] = lo[3]/ichunk;
    iidx[1] = lo[2]/ichunk;
    break;
  case 3:
    iidx[0] = lo[3]/ichunk;
    iidx[1] = lo[1]/ichunk;
    break;
  case 4:
    iidx[0] = lo[1]/ichunk;
    iidx[1] = lo[0]/ichunk;
    break;
  case 5:
    iidx[0] = lo[2]/ichunk;
    iidx[1] = lo[0]/ichunk;
    break;
  default:
    assert(0);
  } 

  return colls.at(arrid)->getProc(idx,size);
}

static void compute_index(void *dscr, int dscr_len, 
			 void *pldata, int pldata_len, int arrid, void *idx, int idxlen) {
  assert(dscr_len == sizeof(task_dscr_t));
  assert(arrid >=0 && arrid<6);
  assert(idxlen ==  2*sizeof(int));
  const long itask = ((task_dscr_t*)dscr)->id;
  int *iidx = (int *)idx;
  task_plo_t *ptplo = (task_plo_t*)pldata;
  assert(pldata_len == sizeof(task_plo_t));
  const long int nbfn = ptplo->nbfn;
  long int ilo, jlo, klo, llo, lo[4], hi[4], tmp;

  compute_lo_hi(itask, nbfn, lo, hi, &ilo, &jlo, &klo, &llo);

  switch(arrid) {
  case 0:
    iidx[0] = lo[1]/ichunk;
    iidx[1] = lo[0]/ichunk;
    break;
  case 1:
    iidx[0] = lo[3]/ichunk;
    iidx[1] = lo[2]/ichunk;
    break;
  case 2:
    iidx[0] = lo[3]/ichunk;
    iidx[1] = lo[2]/ichunk;
    break;
  case 3:
    iidx[0] = lo[3]/ichunk;
    iidx[1] = lo[1]/ichunk;
    break;
  case 4:
    iidx[0] = lo[1]/ichunk;
    iidx[1] = lo[0]/ichunk;
    break;
  case 5:
    iidx[0] = lo[2]/ichunk;
    iidx[1] = lo[0]/ichunk;
    break;
  default:
    assert(0);
  }  
}


static void twoel_task(UniformTaskCollection *utc, void *bigd, int bigd_len, 
		       void *pldata, int pldata_len, vector<void *> data_bufs);

static void cpptwoel(Integer g_schwarz, Integer g_dens, Integer g_fock,
		     double schwmax, Integer nbfn,
		     long long *icut1, long long *icut2,
		     long long *icut3, long long *icut4) {
  int dotask;
  int g_counter, ione = 1;
  task_dscr_t tdscr;
  task_plo_t tplo;

  {
    tplo.nbfn = nbfn;
    tplo.schwmax = schwmax;
    tplo.icut1 = icut1;
    tplo.icut2 = icut2;
    tplo.icut3 = icut3;
    tplo.icut4 = icut4;
    tplo.g_schwarz = g_schwarz;
    tplo.g_dens = g_dens;
    tplo.g_fock = g_fock;

    TslFuncRegTbl frt;
    TslFunc tf = frt.add(twoel_task);
    vector<DataColl*> colls;
    int block[2] = {ichunk, ichunk};
    DenseArray coll_schwarz(g_schwarz, block, 2);
    DenseArray coll_dens(g_dens, block, 2);
    DenseArray coll_fock(g_fock, block, 2);
    /**FIXME: modes are associated with arrays and not specific data
       movement operation*/
    colls.push_back(&coll_schwarz);
    colls.push_back(&coll_schwarz);
    colls.push_back(&coll_dens);
    colls.push_back(&coll_dens);
    colls.push_back(&coll_fock);
    colls.push_back(&coll_fock);

    vector<AccessMode> modes(4, MODE_RONLY);
    vector<int> idxlens(6, 2*sizeof(int));

    modes.push_back(MODE_ACC);
    modes.push_back(MODE_ACC);

    const long me = GA_Nodeid();
    const long nproc = GA_Nnodes();
    const long imax = nbfn/ichunk + ((nbfn%ichunk) ? 1 : 0);
    const long total_ntasks = imax*imax*imax*imax;
    const long ntasks_per_proc = (long)ceil(1.0*total_ntasks/nproc);
    const long tasklo = ntasks_per_proc * me;
    const long taskhi = min(tasklo+ntasks_per_proc,total_ntasks)-1;

    TaskCollProps props;
    props.functions(tf,frt).taskSize(sizeof(task_dscr_t))
#define EVEN_DISTRIBUTION 0
#if EVEN_DISTRIBUTION
        .maxTasks(ntasks_per_proc)
#else
        .maxTasks(total_ntasks)
#endif
        .localData(&tplo,sizeof(tplo));
    UniformTaskCollSplitData utc(props, colls, modes, idxlens, compute_index);

#if EVEN_DISTRIBUTION
    for(long i=tasklo; i<=taskhi; i++)
#else
    for (long i=0; i<total_ntasks; ++i)
      if (me == compute_owner(i, colls, 0, nbfn))
#endif
      {
        tdscr.id = i;
        utc.addTask(&tdscr, sizeof(tdscr));
      }
    GA_Sync();
    utc.process();
  }
}

typedef struct {
  double buf[ichunk][ichunk];
} tbuf_t;

static void twoel_task(UniformTaskCollection *utc, void *_bigd, int bigd_len, 
		       void *pldata, int pldata_len, vector<void *> data_bufs) {
  assert(_bigd!=NULL);
  assert(bigd_len == sizeof(task_dscr_t));
  assert(pldata!=NULL);
  assert(pldata_len == sizeof(task_plo_t));
  assert(data_bufs.size()==6);

  task_dscr_t *ptdscr = (task_dscr_t*)_bigd;
  task_plo_t *ptplo = (task_plo_t*)pldata;

  tbuf_t *s_ij = (tbuf_t*)data_bufs[0];
  tbuf_t *s_kl = (tbuf_t*)data_bufs[1];
  tbuf_t *d_kl = (tbuf_t*)data_bufs[2];
  tbuf_t *d_jl = (tbuf_t*)data_bufs[3];
  tbuf_t *f_ij = (tbuf_t*)data_bufs[4];
  tbuf_t *f_ik = (tbuf_t*)data_bufs[5];

  int ich;
  Integer lch;
  Integer lo_ik[2],hi_ik[2],lo_jl[2],hi_jl[2];
  Integer i,j,k,l,iloc,jloc,kloc,lloc,ld;
  double gg;
  Integer ione = 1;
  double one = 1.0;
  Integer nbfn;
  long lo[4], hi[4];
  long it, jt, kt, lt;
  long long *icut1, *icut2, *icut3, *icut4;
  Integer g_schwarz, g_dens, g_fock;
  double schwmax;

  long itask = ptdscr->id;
  nbfn = ptplo->nbfn;
  schwmax = ptplo->schwmax;
  icut1 = ptplo->icut1;
  icut2 = ptplo->icut2;
  icut3 = ptplo->icut3;
  icut4 = ptplo->icut4;
  g_schwarz = ptplo->g_schwarz;
  g_dens = ptplo->g_dens;
  g_fock = ptplo->g_fock;  

  ld = maxnbfn;
  ich = lch = ichunk;

  compute_lo_hi(itask, nbfn, lo, hi, &it, &jt, &kt, &lt);

  assert(ich == hi[0]-lo[0]+1);
  assert(ich == hi[2]-lo[2]+1);
  clean_chunk((double*)f_ij->buf);
  clean_chunk((double *)f_ik->buf);
  for(i = lo[0]; i<=hi[0]; i++) {
    iloc = i-lo[0];
    for(j = lo[1]; j<= hi[1]; j++) {
      jloc = j-lo[1];
      if (s_ij->buf[jloc][iloc]*(schwmax) < tol2e) {
	*icut1 = *icut1 + (hi[2]-lo[2]+1)*(hi[3]-lo[3]+1);
      }
      else {
	/* *cnt1+=1;*/
	for(k = lo[2]; k<=hi[2]; k++) {
	  kloc = k-lo[2];
	  for(l = lo[3]; l<=hi[3]; l++) {
	    lloc = l-lo[3];
	    if (s_ij->buf[jloc][iloc] * s_kl->buf[lloc][kloc] < tol2e) {
	      *icut2 = *icut2 + 1;
	    }
	    else  {
	      Integer _i=i+1, _j=j+1, _k=k+1, _l=l+1;
	      g_(&gg, &_i, &_j, &_k, &_l);
	      /* *cnt2 += 1;*/
	      f_ij->buf[jloc][iloc] = f_ij->buf[jloc][iloc] + gg*d_kl->buf[lloc][kloc];
	      f_ik->buf[kloc][iloc] = f_ik->buf[kloc][iloc] - 0.5*gg*d_jl->buf[lloc][jloc];
	      *icut3 = *icut3 + 1;
	    }
	  }
	}
      }
    }
  }
}

