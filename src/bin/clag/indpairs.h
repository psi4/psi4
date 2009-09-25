#ifndef _psi_src_bin_detcas_indpairs_h
#define _psi_src_bin_detcas_indpairs_h

namespace psi { namespace clag {

/*
** INDPAIRS.H
** 
** Contains code pertaining to the "independent pairs" of orbitals for which
** the energy is not invariant.  Only these pairs need to be considered in
** the MCSCF orbital rotation procedure.
**
** C. David Sherrill
** University of California, Berkeley
** May 1998
*/

class IndepPairs {

  protected:

    // The following items are overall...contains all irreps
    int nirreps;           // number of irreducible representations in pt grp
    int npairs;            // number of independent orbital pairs
    int *p;                // array of first elements in pairs 
    int *q;                // array of second elements in pairs
    int *map_pair_ir;      // map a pair to what irrep it is
    int *map_pair_rel;     // map a pair to the relative pair within an irrep

    // The following items are per-irrep
    int *ir_npairs;        // number of pairs for an irrep
    int **ir_p;            // first index of pairs for a given irrep
    int **ir_q;            // second index of pairs for a given irrep
    int **ir_p_rel;        // relative first index of pairs for irrep
    int **ir_q_rel;        // relative second index of pairs for irrep
    int **ir_map_pair;     // absolute number of pair for pair within irrep

  public:
    IndepPairs();
    IndepPairs(int nirr, int num_ras, int **ras_opi, int ***ras_orbs,
      int *fzc, int **fzc_orbs, int *cor, int **cor_orbs,
      int *vir, int **vir_orbs, int *fzv, int **fzv_orbs,
      int *ci2relpitz, int ignore_ras_ras, int ignore_fz);
   ~IndepPairs();
    void set(int nirr, int num_ras, int **ras_opi, int ***ras_orbs,
      int *fzc, int **fzc_orbs, int *cor, int **cor_orbs,
      int *vir, int **vir_orbs, int *fzv, int **fzv_orbs,
      int *ci2relpitz, int ignore_ras_ras, int ignore_fz);
    void set_part(int &count, int *ir_cnt, int *num_orbs_i, 
      int *num_orbs_j, int **orbs_i, int **orbs_j, int *ci2relpitz);
    void print(FILE *outfile);
    int get_num_pairs(void) { return npairs; }
    int * get_p_ptr(void) { return p; }
    int * get_q_ptr(void) { return q; }
    void print_vec(double *arr, const char *label, FILE *outfile);
    int get_ir_num_pairs(int h) { return ir_npairs[h]; }
    int * get_ir_p_ptr(int h) { return ir_p[h]; }
    int * get_ir_q_ptr(int h) { return ir_q[h]; }
    int * get_ir_prel_ptr(int h) { return ir_p_rel[h]; }
    int * get_ir_qrel_ptr(int h) { return ir_q_rel[h]; }
    double * get_irrep_vec(int h, double *arr);
    void put_irrep_vec(int irrep, double *ir_vec, double *tot_vec);

  private:
    void print_selected(int num, int *parr, int *qarr, FILE *outfile);
    void print_selected(int num, int *parr, int *qarr,
                        int *prel, int *qrel, FILE *outfile);

};

}} // end namespace psi::clag

#endif // header guard
