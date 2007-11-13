#ifndef _psi3_libiwl_iwl_hpp_
#define _psi3_libiwl_iwl_hpp_

#include <libiwl/config.h>

namespace psi {
  class PSIO;
  
  /**
		IWL is an instance of the libiwl library. Multiple instances are 
		supported by having multiple instances of libpsio available.
		
		Each instance is configured by the constructor:
		IWL new_instance(&libpsio_instance, PSIO_OPEN_OLD);
	*/
  class IWL
  {
    PSIO* psio;
    struct iwlbuf Buf;
    bool keep;
    
  public:
    // Equivalent to iwl_buf_init
    IWL(PSIO* psio_obj, int itape, double cutoff, int oldfile, int readflag);
    // Equivalent to iwl_buf_close
    virtual ~IWL();
    
    // Set flag for keeping the file after closing. Default is true.
    void set_keep(bool k = true) { keep = k; }
    
    // Fetch the next buffer
    void fetch();
    
    // Flush the buffer
    void flush(int lastbuf);
    
    // Put the buffer to disk
    void put();
    
    // Read buffer
    int rd(int target_pq, double *ints, int *ioff_lt, int *ioff_rt, int mp2, int printflg, 
           FILE *outfile);
    int rd_all(double *ints, int *ioff_lt, int *ioff_rt, int no_pq_perm, int *ioff,
           int printflg, FILE *outfile);
    int rd_all2(double **ints, int *ioff_lt, int *ioff_rt, int no_pq_perm, int *ioff, 
           int printflg, FILE *outfile);
    int rd_all_act(double *ints, int *ioff_lt, int *ioff_rt, int no_pq_perm, int *ioff,
           int fstact, int lstact, int printflg, FILE *outfile);
    int rd_all_mp2r12a(double *ints, int *ioff_lt, int *ioff_rt, int bra_ket_symm, 
           int *ioff, int printflg, FILE *outfile);
    int rd_arr(int target_pq, double *ints, int *rlist, int *slist, int *size, int *ioff, 
           int printflg, FILE *outfile);
    int rd_arr2(double *ints, int *plist, int *qlist, int *rlist, int *slist, int *size, 
           int *ioff, int printflg, FILE *outfile);
           
    bool to_end();
    
    // Write buffer
    void wrt(int p, int q, int pq, int pqsym, double *arr, int rmax, int *active, 
           int *ioff, int *orbsym, int *firsti, int *lasti, int sortby_rs, int printflag, 
           FILE *outfile);
    void wrt_all(int nbfso, double *ints, int *ioff, int printflg, FILE *outfile);
    void wrt_arr(double *arr, int *p, int *q, int *r, int *s, long int size);
    void wrt_arr2(double *arr, int p, int q, int *rlist, int *slist, int size, 
           int printflag, FILE *outfile);
    void wrt_arr_SI(double *arr, short int *p, 
           short int *q, short int *r, short int *s, int size);
    void wrt_arr_SI_nocut(double *arr, short int *p, short int *q, short int *r, 
           short int *s, int size);
    void wrt_mat(int ptr, int qtr, double **mat, int rfirst, int rlast, int sfirst, 
           int slast, int *reorder, int reorder_offset, int printflag, int *ioff,
           FILE *outfile);    
    void wrt_mat2(int ptr, int qtr, double **mat, int rfirst, int rlast, int sfirst, 
           int slast, int *reorder, int reorder_offset, int printflag, int *ioff,
           FILE *outfile);
    void wrt_mp2(int p, int q, int pq, int pqsym, double **arr, int rsym, int *firstr, 
           int *lastr, int *firsts, int *lasts, int *occ, int *vir, int *ioff, 
           int printflag, FILE *outfile);
    void wrt_mp2r12a(int p, int q, int pq, int pqsym, double **arr, int rsym, 
           int *firstr, int *lastr, int *firsts, int *lasts, int *occ, int bra_ket_symm, 
           int *ioff, int printflag, FILE *outfile);
    void wrt_val(int p, int q, int r, int s, double value, int printflag, FILE *outfile, int dirac);
    void wrt_val_SI(short int p, short int q, short int r, short int s, double value, int printflag, 
           FILE *outfile, int dirac);
    
    static int rdone(PSIO *psio, int itap, char *label, double *ints, int ntri, int erase, 
           int printflg, FILE *outfile);
    static void rdtwo(PSIO *psio, int itap, double *ints, int *ioff, int norbs, 
           int nfzc, int nfzv, int printflg, FILE *outfile);
           
    static void sortbuf(IWL *Inbuf, IWL *Outbuf,
           double *ints, int fpq, int lpq, int *ioff, int *ioff2, 
           int nbfso, int elbert, int intermediate, int no_pq_perm, 
           int qdim, int add, int printflg, FILE *outfile);
    
    static void wrtone(PSIO* psio, int itap, char *label, int ntri, double *onel_ints);
    static void wrttwo(PSIO* psio, int itap, int nbfso, double *ints, int *ioff, double toler, 
           int printflg, FILE *outfile);
  };
};

#endif // _psi3_libiwl_iwl_hpp_h

