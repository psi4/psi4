/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

/*! \file
    \ingroup TRANSQT
    \brief Enter brief description of file here 
*/
#ifndef _psi3_bin_transqt_yoshimine_h_
#define _psi3_bin_transqt_yoshimine_h_

#include <libmints/typedefs.h>

namespace psi { 

/*
** YOSHIMINE.H 
** Function prototypes for Yoshimine sort object
**
** David Sherrill
** Center for Comptational Quantum Chemistry, UGA
** February 1995
*/

/* need to include iwl.h before including this file */

struct IWLLight {
    int p_;
    int q_;
    int r_;
    int s_;
    double val_;
};

class BaseBucket {
  private:
    long int in_bucket_;
    IWL* IWLBuf_;
    int hi_;
    int lo_;
  // You should never copy a bucket, it's inefficient.
  // Also it contains raw pointers and copying is dangerous.
    BaseBucket(const BaseBucket& input ) {}
    BaseBucket& operator =(const BaseBucket& input) {}

  public:
    virtual ~BaseBucket();
    BaseBucket() : IWLBuf_(NULL) {}
  // Accessor functions to data
    long int& in_bucket()       { return in_bucket_; }
    IWL* iwlbuf()               { return IWLBuf_; }
    int& hi()                   { return hi_; }
    int& lo()                   { return lo_; }

    // Setting up files
    virtual void set_file(IWL* input) { IWLBuf_ = input; }
    virtual void set_file(std::string fileID) {}
    // Operation to close files
    virtual void close_file(int keep);
    // Allocating bucket-owned storage
    virtual void alloc(unsigned long int size) {};
    // Deleting all bucket-owned storage
    virtual void dealloc();
    // Filling the bucket with integrals and labels
    virtual void fill(int pin, int qin, int rin, int sin, double valin) {};
    // Flushing the bucket to disk.
    virtual void flush(const int lastbuf) {};

    // Opening bucket for read
    virtual FILE* open_bucket(unsigned long int &nints) {}
    // Closing bucket after read
    virtual void close_bucket(FILE* handle, int erase) {}
};

// Usual structure of a bucket, containing integrals and labels
// for a certain range of indices
class Bucket : public BaseBucket {
  private:
    int* p_;
    int* q_;
    int* r_;
    int* s_;
    double* val_;
  public:
  // Constructor and destructor
  Bucket() : p_(NULL), q_(NULL), r_(NULL), s_(NULL), val_(NULL) {}
  virtual ~Bucket();
  // Accessor functions to data
    int* p()                    { return p_; }
    int* q()                    { return q_; }
    int* r()                    { return r_; }
    int* s()                    { return s_; }
    double* val()               { return val_; }

    // Allocating bucket-owned storage
    virtual void alloc(unsigned long int size);
    // Deleting all bucket-owned storage
    virtual void dealloc();
    // Filling the bucket with integrals and labels
    virtual void fill(int pin, int qin, int rin, int sin, double valin);
    // Flushing the bucket
    virtual void flush(const int lastbuf);
};

// Light Bucket that avoids IWL buffers.
class BucketLight : public BaseBucket {
private:
    IWLLight* ints_;
    unsigned long int nints_;
    std::string filename_;
    FILE* fh_;
public:
    BucketLight();
    virtual ~BucketLight();

    // Accessor functions to data
    IWLLight* ints()                  {return ints_; }
    unsigned long int& nints()        {return nints_; }
    std::string& filename()           {return filename_; }
    FILE* file_pointer()              {return fh_; }

    // Set up the files
    virtual void set_file(std::string fileID);
    // Close the file
    virtual void close_file(int keep);
    // Allocating bucket-owned storage
    virtual void alloc(unsigned long int size);
    // Deleting all bucket-owned storage
    virtual void dealloc();
    // Filling the bucket with integrals and labels
    virtual void fill(int pin, int qin, int rin, int sin, double valin);
    // Flushing the bucket
    virtual void flush(const int lastbuf);
    // Opening bucket for read
    virtual FILE* open_bucket(unsigned long &nints);
    // Closing bucket after read
    virtual void close_bucket(FILE* handle, int erase);
};

// Base class for a Yoshimine object, handling
// Yoshimine sorting of the integrals.
// First, a pre-sorting is performed based on the
// canonical index PQRS into different files.
// Then, each file is read and the integrals sorted in core
// before to be written to the PK file.
class YoshBase {
  private:
    int core_loads_;
    int nbuckets_;
    int* bucket_for_pq_;
    unsigned long int bucketsize_;
    int first_tmp_file_;
    unsigned int bra_indices_;
    double cutoff_;
   // boost::shared_ptr<PSIO> psio_;
    PSIO* psio_;
    std::vector<BaseBucket*> buckets_;

  // Never copy this, it contains raw pointers and file pointers
    YoshBase(const YoshBase& input) {}
    YoshBase& operator = (const YoshBase& input) {}

  protected:
  // Constructor and destructor, only callable through derived classes
    YoshBase(unsigned bra_idx, long maxcor, long maxcord, const int max_buckets,
         unsigned int first_tmp_file, double cutoff, PSIO* psio);
    virtual ~YoshBase();

  public:
  // Accessor functions
    int& core_loads()                    {return core_loads_; }
    int& nbuckets()                      {return nbuckets_; }
    int* bucket_for_pq()                 {return bucket_for_pq_; }
    unsigned long int& bucketsize()      {return bucketsize_; }
    int& first_tmp_file()                {return first_tmp_file_; }
    unsigned int& bra_indices()          {return bra_indices_; }
    double& cutoff()                     {return cutoff_; }
    std::vector<BaseBucket*>& buckets()  {return buckets_; }
    PSIO* get_psio()                     {return psio_; }

  // Only the functions I need for PK are ported here right now

    void init(long maxcord, unsigned int bra_idx);
    void print();
    virtual void init_buckets();
    void close_buckets(int erase);
    static void rdtwo_pk(YoshBase* YBuffJ, YoshBase* YBuffK, int itapERI,
         int del_tei_file, int nirreps, int* so2rel, int* so2sym, int* pksymoff, int printflag);
    virtual void sort_pk(int is_exch, int out_tape, int keep_bins,
         int* so2ind, int* so2sym, int* pksymoff, int print_lvl);
    void done();
    void flush();
//    void buff_put_val(int *ioff, int pq,
//         int p, int q, int r, int s, double value, int prtflg,
//         std::string OutFileRMR);
}; 

// Usual Yoshimine class with usual buckets
class Yosh : public YoshBase {
private:
  // Never copy this, it contains raw pointers and file pointers
//    Yosh(const Yosh& input) {}
//      Yosh& operator = (const Yosh& input) {}
public:
    // Constructor and destructor
    Yosh(unsigned int bra_idx, long maxcor, long maxcord, const int max_buckets,
         unsigned int first_tmp_file, double cutoff, PSIO *psio);
    virtual ~Yosh() {}
};

// New Yoshimine class without IWL
class YoshLight : public YoshBase {
public:
    // Constructor and destructor
    YoshLight(unsigned int bra_idx, long maxcor, long maxcord, const int max_buckets,
         unsigned int first_tmp_file, double cutoff, PSIO *psio);
    virtual ~YoshLight() {}

    // Bucket initialization needs specialization
    virtual void init_buckets();
    // Bucket sorting needs specialization...
    virtual void sort_pk(int is_exch, int out_tape, int keep_bins,
         int* so2ind, int* so2sym, int* pksymoff, int print_lvl);
    void sort_buffer_pk_light(BaseBucket *bucket, int out_tape, int is_exch,
                              double* ints, unsigned int fpq, unsigned int lpq,
                              int *so2ind, int *so2sym, int *pksymoff, int printflg,
                              std::string out);
};

class Yoshopt : public YoshBase {
private:
  // Never copy this, it contains raw pointers and file pointers
//    Yoshopt(const Yoshopt& input) {}
//      Yoshopt& operator = (const Yoshopt& input) {}
public:
    // Constructor and destructor
    Yoshopt(unsigned bra_idx, long maxcor, long maxcord, const int max_buckets,
         unsigned int first_tmp_file, double cutoff);
    ~Yoshopt() {}
};

namespace transqt {

struct bucket {
   long int in_bucket;
   struct iwlbuf IWLBuf;
   int *p;
   int *q;
   int *r;
   int *s;
   double *val;
   int hi;
   int lo;
   };

struct yoshimine {
   int core_loads;
   int nbuckets;
   int *bucket_for_pq;
   unsigned long int bucketsize;
   struct bucket *buckets;
   int first_tmp_file;
   int pq_per_bucket;
   int bra_indices;
   int ket_indices;
   double cutoff;
   };

union psi_buffer {
  double *pki ;
  int *lbli ;
  unsigned char *lbl;
  double *val;
} ;

#ifdef YEXTERN
#undef YEXTERN
#define YEXTERN  
#else
#define YEXTERN extern
#endif

YEXTERN void yosh_init(struct yoshimine *YBuff, unsigned bra_indices, 
      unsigned ket_indices, long maxcor, long maxcord,
const int max_buckets,
      unsigned int first_tmp_file, double cutoff, std::string OutFileRMR);
YEXTERN void yosh_init_pk(struct yoshimine *YBuff, unsigned bra_indices,
      long maxcor, long maxcord, const int max_buckets,
      unsigned int first_tmp_file, double cutoff, std::string OutFileRMR);
YEXTERN void yosh_print(struct yoshimine *YBuff, std::string OutFileRMR);
YEXTERN void yosh_init_buckets(struct yoshimine *YBuff);
YEXTERN void yosh_close_buckets(struct yoshimine *YBuff, int erase);
YEXTERN void yosh_rdtwo(struct yoshimine *YBuff, int itapERI, int del_tei_file, int *num_so,
      int nirreps, int *ioff, int elbert, int fzcflag, double *P,
      double *Hc, int matrix, int printflag, std::string OutFileRMR);
YEXTERN void yosh_rdtwo_pk(struct yoshimine *YBuffJ, struct yoshimine *YBuffK, int itapERI,
      int del_tei_file, int nirreps, int* so2rel, int* so2sym, int* pksymoff, int printflag);
YEXTERN void yosh_rdtwo_uhf(struct yoshimine *YBuff, int itapERI, int del_tei_file, int *num_so,
      int nirreps, int *ioff, int elbert, int fzcflag, double *Pa, double *Pb,
      double *Hca, double *Hcb, int matrix, int printflag, std::string OutFileRMR);
YEXTERN void yosh_rdtwo_backtr(struct yoshimine *YBuff, int tei_file, 
      int *ioff, int symmetrize, int add_ref_pt, int del_tei_file, int prtflg, 
      std::string OutFileRMR);
YEXTERN void yosh_rdtwo_backtr_uhf(std::string, struct yoshimine *YBuff, int tei_file, 
      int *ioff, int swap_bk, int symm_pq, int del_tei_file, int prtflg, std::string OutFileRMR);
YEXTERN void flush_bucket(struct bucket *bptr, int lastbuf);
YEXTERN void yosh_sort(struct yoshimine *YBuff, int out_tape, int keep_bins,
      int *ioff, int *ioff2, int nbfso, int nbstri, 
      int elbert, int intermediate, int no_pq_perm, int qdim,
      int add, int print_lvl, std::string OutFileRMR);
YEXTERN void yosh_sort_pk(struct yoshimine *YBuff, int is_exch, int out_tape, int keep_bins,
      int* so2ind, int* so2sym, int* pksymoff, int print_lvl);
YEXTERN void yosh_done(struct yoshimine *YBuff);
YEXTERN void yosh_flush(struct yoshimine *YBuff);
YEXTERN void yosh_wrt_arr(struct yoshimine *YBuff, int p, int q, int pq, 
   int pqsym, double *arr, int rmax, int *ioff, 
   int *orbsym, int *firsti, int *lasti, int sortby_rs, int printflag, 
   std::string OutFileRMR);
YEXTERN void yosh_wrt_arr2(struct yoshimine *YBuff, int size, double *arr,
   int p, int q, int *rlist, int *slist, int *ioff, int printflag,
   std::string OutFileRMR);
YEXTERN void yosh_wrt_arr_mp2(struct yoshimine *YBuff, int p, int q, int pq,
                      int pqsym, double **arr, int rsym, int *firstr,
                      int *lastr, int *firsts, int *lasts, int sortby_rs,
                      int ndocc, int nvirt, int *occ, int *vir, int *ioff3,
                      int printflag, std::string OutFileRMR);
YEXTERN void add_2pdm_ref_pt(struct yoshimine *YBuff,int *ioff,int prtflg,
                             std::string OutFileRMR);
YEXTERN void yosh_buff_put_val(struct yoshimine *YBuff, int *ioff, int pq,
                       int p, int q, int r, int s, double value, int prtflg,
                       std::string OutFileRMR);
YEXTERN void yosh_wrt_arr_mp2r12a(struct yoshimine *YBuff, int p, int q, int pq,
                          int pqsym, double **arr, int rsym, int *firstr,
                          int *lastr, int *firsts, int *lasts, int sortby_rs,
                          int *occ, int *ioff3,
                          int printflag, std::string OutFileRMR);
}} // end namespace psi::transqt
#endif // header guard
