/*! \file
 \ingroup CINTS
 \brief Compute unique te_deriv ints in symmetry block order, contract with correlated tpdm from disk
 */
#include <vector>
#include <sstream>
#include <algorithm>
#include<cmath>
#include<cstring>
#include<cstdio>
#include<memory.h>
#include<cstdlib>

#include<psitypes.h>
#include<libipv1/ip_lib.h>
#include<libiwl/iwl.h>
#include<libciomr/libciomr.h>

#include <libint/libint.h>
#include<libderiv/libderiv.h>
#include"defines.h"
#include"small_fns.h"
#define EXTERN
#include"global.h"
#include <stdexcept>
#include"schwartz.h"
#include"deriv1_quartet_data.h"
#include"iwl_tebuf.h"
#include"norm_quartet.h"
#ifdef USE_TAYLOR_FM
#include"taylor_fm_eval.h"
#else
#include"int_fjt.h"
#endif

using std::vector;
using std::min;

namespace psi {
namespace cints {
// For a set of SALCs, we need to be able to select coefficients of all nonzero derivatives
// of a given shell quartet. It is best done when the relevant subset of the (sparse) SALC coefficient
// matrix is pre-arranged in an array
typedef struct {
    double coef; // contribution of this cartesian derivative to this SALC
    int cart_der; // cartesian derivative (for this quartet) [0,3*num_unique_atoms)
    int salc; // SALC index [0,num_atoms*3)
} cdsalc_elem;

// Class to facilitate iterating through shells
// that contribute to particular irreps of symmetry orbital pairs
class symblock_shell_order {
  private:
    vector<int> shells_[8]; // indices of those shells (one shell can appear in multiple irreps)

  public:
    /// use data in Symmetry to fill std::vector<int> shells_[nirrep]
    void init() {
      int aso = 0;
      for (int g = 0; g < Symmetry.nirreps; g++) {
        shells_[g].clear();
        shells_[g].push_back(Symmetry.uso2shell[aso]);
        aso++;
        for (int so = 1; so < Symmetry.sopi[g]; so++, aso++) {
          if (Symmetry.uso2shell[aso - 1] != Symmetry.uso2shell[aso]) {
            shells_[g].push_back(Symmetry.uso2shell[aso]);
          }
        }
      }
    }

    size_t n(int g) {
      return shells_[g].size();
    }

    // This could benefit from using std::vector<int>::iterator objects
    // instead of keeping track of 'n' externally
    int s(int g, int n) {
      return shells_[g][n];
    }

};

void te_deriv1_corr_symm() {
  FILE *debugfile; // jtf
  symblock_shell_order sso;

  // grabbed from te_deriv1_corr
  struct iwlbuf TPDM; /* IWL buffer for two-pdm matrix elements */
  const double toler = UserOptions.cutoff;
  const int num_coords = 3 * Molecule.num_atoms; // total number of coordinates, symmetry-adapted or otherwise
  // max number of nonzero derivatives for a SO shell quartet = 4 (centers) * 3 (x,y,z) * maximum multiplicity of a nucleus
  const int max_num_cdsalcs_per_quartet = 12 * Symmetry.max_stab_index;

  /*--- ASCII file to print integrals ---*/
  FILE *d1eriout;

  /*--- Various data structures ---*/
  struct iwlbuf* D1ERIOUT; /* IWL buffer for target integrals */
  struct tebuf** tot_data; /* accum. for contracted integrals */
  vector<int> num_of_ints_in_totdata; /* counts how many integrals are in each tot_data */

  struct shell_pair *sp_ij, *sp_kl;
  struct unique_shell_pair *usp_ij, *usp_kl;

  Libderiv_t Libderiv;
#ifndef USE_TAYLOR_FM
  double_array_t fjt_table;
#endif

  std::vector<int> unique_center;
  unique_center.reserve(4);
  std::vector<int> nuc(4);

  PSI_INT_LEAST64 total_te_count = 0;
  // should no longer need this
  // vector<PSI_INT_LEAST64> te_count_per_coord(num_coords, 0); // keeps track of how many integrals were written out

  int ij, kl, ik, jl, ijkl;
  int ioffset, joffset, koffset, loffset;
  int count;
  int dum;
  int total_am, am;
  int orig_am[4];
  int pkblock_end_index = -1;
  int i, j, k, l, m, ii, jj, kk, ll;
  int si, sj, sk, sl;
  int sii, sjj, skk, sll, slll;
  bool switch_ij, switch_kl, switch_ijkl;
  int pi, pj, pk, pl;
  int max_pj, max_pl;
  int upk, num_unique_pk;
  int usi_arr[3], usj_arr[3], usk_arr[3], usl_arr[3];
  int *sj_arr, *sk_arr, *sl_arr;
  int *sj_fbf_arr, *sk_fbf_arr, *sl_fbf_arr;
  int usi, usj, usk, usl;
  int stab_i, stab_j, stab_k, stab_l, stab_ij, stab_kl;
  int *R_list, *S_list, *T_list;
  int R, S, T;
  int dcr_ij, dcr_kl, dcr_ijkl;
  int lambda_T = 1;
  int num_unique_quartets;
  int plquartet;
  int max_num_unique_quartets;
  int max_num_prim_comb;

  int size, class_size;
  int max_cart_class_size;

  int bf_i, bf_j, bf_k, bf_l, so_i, so_j, so_k, so_l, s;
  int np_i, np_j, np_k, np_l;
  int ni, nj, nk, nl;

  int index;
  int iimax, jjmax, kkmax, llmax;
  int irrep, npi_ij, npi_kl, npi_ik, npi_jl, ind_offset;

  int num_prim_comb, p;

  double AB2, CD2;
  double *FourInd;
  double **grad_te_local;
  double pfac;
  double temp;
  double alpha, beta;
  double **dens_i, **dens_j;
  double ddax, dday, ddaz, ddbx, ddby, ddbz, ddcx, ddcy, ddcz, dddx, dddy, dddz;

  double pkblock_end_value = 0.0;

  vector<double> so_int;
  {
    const int max_num_cdsalc_per_quartet = min(max_num_cdsalcs_per_quartet, CDSALCs.nsalcs);
    so_int.reserve(max_num_cdsalc_per_quartet);
  }

  /*---------------
   Initialization
   ---------------*/
  int nirrep = Symmetry.nirreps;
  const int nsalcs = CDSALCs.nsalcs;
  vector<PSI_INT_LEAST64> te_count_per_buffer(nirrep * nirrep * nirrep, 0);
  // initialize symblock shell order
  sso.init();
  printf("Printing shell ordering within irrep \n");
  for (int g = 0; g < nirrep; g++) {
    printf("irrep = %d, num shells in irrep = %d\n", g, (int) sso.n(g));
    for (int s = 0; s < sso.n(g); s++) {
      printf("%d\t", sso.s(g, s));
    }
    printf("\n");
  }

  debugfile = fopen("te_deriv1_debug", "w"); //jtf
  fflush(stdout);
  // may also need init_molecule, init_symm, init_sp, init_uniq_sp or whatever

#if PRINT
  eriout = fopen("d1eriout.dat","w");
#endif
  // will no longer be writing iwl buffers to disk
  //D1ERIOUT = new struct iwlbuf[nsalcs * nsalcs * nsalcs];
  //for (int c = 0; c < num_coords; c++)
  //  iwl_buf_init(&D1ERIOUT[c], IOUnits.itapD1ERI_SO + c, toler, 0, 0);
#ifdef USE_TAYLOR_FM
  init_Taylor_Fm_Eval(BasisSet.max_am*4-4+DERIV_LVL,UserOptions.cutoff);
#else
  init_fjt(BasisSet.max_am * 4 + DERIV_LVL);
  init_fjt_table(&fjt_table);
#endif
  init_libderiv_base();
  printf("after init_libderiv_base()\n");
  fflush(stdout);
  /*-------------------------
   Allocate data structures
   -------------------------*/

  // print uso2ao matrix
  print_mat(Symmetry.usotao, Symmetry.num_so, BasisSet.num_ao, stdout);

  grad_te_local = block_matrix(Molecule.num_atoms, 3);

  // jtf - OK
  max_cart_class_size = (ioff[BasisSet.max_am]) * (ioff[BasisSet.max_am]) * (ioff[BasisSet.max_am])
      * (ioff[BasisSet.max_am]);
  // jtf - OK
  max_num_unique_quartets = Symmetry.max_stab_index * Symmetry.max_stab_index
      * Symmetry.max_stab_index;

  int max_so_per_irrep = 0;
  for (i = 0; i < Symmetry.nirreps; ++i) {
    if (Symmetry.sopi[i] > max_so_per_irrep)
      max_so_per_irrep = Symmetry.sopi[i];
  }

  // jtf - No Good.  Need to come up with new scheme for saving completed integrals
  // this vision was for keeping an entire block in memory
  // 1) that's stupid, it would require N^4 storage
  // 2) we want them direct
  // Final integrals are stored here; size of buffers, too
  //tot_data = new struct tebuf*[Symmetry.nirreps];
  //num_of_ints_in_totdata.resize(Symmetry.nirreps);
  //for (int i = 0; i < Symmetry.nirreps; i++) {
  //  tot_data[i] = new struct tebuf[ipow(max_so_per_irrep, 4)];
  //  num_of_ints_in_totdata[i] = 0;
  //}

  // jtf - don't know if this is OK yet
  // I /think/ it is, b/c the ao integrals will still need to be accumulated
  // into so integrals for a shell quartet in exactly the same way as before
  // These arrays are used to hold cartesian AO derivative integrals
  double* cart_ints[12];
  for (int i = 0; i < 12; i++)
    cart_ints[i] = new double[max_cart_class_size];
  // plist_ints holds all symmetry-unique AO shell quartets for a given SO shell quartet
  // plist_ints[i][j] is the pointer to j-th *total* derivative of i-th AO shell quartet which contributes to the current SO quartet
  // note that only total derivatives are stored, i.e. partial derivatives produced by libderiv are combined
  // to yield total derivatives
  double*** plist_ints = new double**[max_num_unique_quartets];
  for (int i = 0; i < max_num_unique_quartets; i++) {
    plist_ints[i] = new double*[12];
    for (int j = 0; j < 12; j++) {
      plist_ints[i][j] = new double[max_cart_class_size];
    }
  }

  // jtf - unchanged, deals with individual shell quartet things.

  // For every SO shell quartet we need the following information:
  // 1) irreps of all total derivative operators which give nonzero when applied to the SO quartet
  // 2) given an irrep, which total derivative operators contribute
  // For every AO shell quartet in a petit list for a given SO shell quartet we need the following information:
  // 1) number of nonzero total derivatives (same as number of distinct centers, n, times 3, 1 <= n <= 4)
  // 2) irreps of all total derivative operators which give nonzero when applied to the AO quartet
  // 3) given an irrep, which total derivative operators contribute

  // For a set of SALCs, we need to be able to select coefficients of all nonzero derivatives
  // of a given shell quartet. It is best done when the relevant subset of the (sparse) SALC coefficient
  // matrix is pre-arranged in an array
  // plist_salcs[s][g][i] is an i-th nonzero contribution to a SALC of irrep g
  // which results in nonzero when applied to member s of a petite list

  struct cdsalc_elem_vec {
      int nelems;
      cdsalc_elem* elems;
  };
  cdsalc_elem_vec** plist_salcs = new cdsalc_elem_vec*[max_num_unique_quartets];
  for (int i = 0; i < max_num_unique_quartets; i++) {
    plist_salcs[i] = new cdsalc_elem_vec[Symmetry.nirreps];
    for (int j = 0; j < Symmetry.nirreps; j++) {
      plist_salcs[i][j].nelems = 0;
      plist_salcs[i][j].elems = new cdsalc_elem[max_num_cdsalcs_per_quartet];
    }
  }

  int* salc_all2thisquartet = new int[num_coords]; // maps absolute SALC index to the SALC index for this quartet
  int* salc_thisquartet2all = new int[max_num_cdsalcs_per_quartet]; // the reverse of the above

  sj_arr = (int *) malloc(sizeof(int) * max_num_unique_quartets);
  sk_arr = (int *) malloc(sizeof(int) * max_num_unique_quartets);
  sl_arr = (int *) malloc(sizeof(int) * max_num_unique_quartets);
  if (Symmetry.nirreps > 1) {
    sj_fbf_arr = (int *) malloc(sizeof(int) * max_num_unique_quartets);
    sk_fbf_arr = (int *) malloc(sizeof(int) * max_num_unique_quartets);
    sl_fbf_arr = (int *) malloc(sizeof(int) * max_num_unique_quartets);
  }

  max_num_prim_comb = (BasisSet.max_num_prims * BasisSet.max_num_prims) * (BasisSet.max_num_prims
      * BasisSet.max_num_prims);

  init_libderiv1(&Libderiv, BasisSet.max_am - 1, max_num_prim_comb, max_cart_class_size);

  // jtf - start new organizational structure to multiply generate integral derivatives
  // in order to have access to an ordered <ij||kl>*12 set, where all <ij| belong to
  // one irrep and all |kl> belong to another (in this case, same)
  int G_target = 0; // set irrep of target quantity
  for (int G_ij = 0; G_ij < Symmetry.nirreps; ++G_ij) { // loop over all possible irreps of <ij|
    int G_kl = G_ij ^ G_target;
    fprintf(outfile, "Beginning Block %d, %d\n", G_ij, G_kl);
    // For each G_ij / G_kl symmetry block, need to seek through all
    // unique shell quartets, and compute those that will contribute.
    // Contribution depends on 'span' of each shell pair, precomputed
    // in init_shell_pairs.  Not sure what this will do when shells
    // are resorted into 'ascending' AM.  Probably nothing, if permutational
    // symmetry is handled correctly.
#if 0
    // six open braces to generate usii, usjj, uskk, usll to test for
    // conicidence of centers
    int G_i, G_j, G_k, G_l;
    for (G_i = 0; G_i < Symmetry.nirreps; G_i++) {
      G_j = G_i ^ G_ij;
      for (G_k = 0; G_k < Symmetry.nirreps; G_k++) {
        G_l = G_k ^ G_kl;
        for (int i_sso = 0; i_sso < sso.n(G_i); i_sso++) {
          int usii = sso.s(G_i, i_sso);
          for (int j_sso = 0; j_sso < sso.n(G_j) && sso.s(G_j, j_sso) <= usii; j_sso++) {
            int usjj = sso.s(G_j, j_sso);
            for (int k_sso = 0; k_sso < sso.n(G_k); k_sso++) {
              int uskk = sso.s(G_k, k_sso);
              int usll_max = (usii == uskk ? usjj : uskk);
              for (int l_sso = 0; l_sso < sso.n(G_k) && sso.s(G_l, l_sso) <= usll_max; l_sso++) {
                int usll = sso.s(G_l, l_sso);


#endif
#if 1
    // six open braces to generate usii, usjj, uskk, usll to test for
    // conicidence of centers
                { { // G_i -> G_k blanks
    for (int usii = 0; usii < Symmetry.num_unique_shells; usii++) {
      for (int usjj = 0; usjj <= usii; usjj++) {
        sp_ij = &(BasisSet.shell_pairs[Symmetry.us2s[usii]][Symmetry.us2s[usjj]]);
        if (!(sp_ij->span[G_ij])) {
          //printf("<%d %d| rejected, != %d\n", usii, usjj, G_ij);
          continue;
        }
        for (int uskk = 0; uskk < Symmetry.num_unique_shells; uskk++) {
          const int usll_max = (usii == uskk ? usjj : uskk);
          for (int usll = 0; usll <= usll_max; usll++) {
            sp_kl = &(BasisSet.shell_pairs[Symmetry.us2s[uskk]][Symmetry.us2s[usll]]);
            if (!(sp_kl->span[G_kl])) {
              //printf("|%d %d> rejected, != %d\n", uskk, usll, G_kl);
              continue;
            }
#endif

            // if all shells are on same center, shell quartet contributes zero to gradient
            int center_i = BasisSet.shells[Symmetry.us2s[usii]].center - 1;
            int center_j = BasisSet.shells[Symmetry.us2s[usjj]].center - 1;
            int center_k = BasisSet.shells[Symmetry.us2s[uskk]].center - 1;
            int center_l = BasisSet.shells[Symmetry.us2s[usll]].center - 1;
            if (center_i == center_j && center_j == center_k && center_k == center_l) {
              //printf("%d %d %d %d rejected, no unique centers\n", usii, usjj, uskk, usll);
              continue;
            }

            // done eliminating shell quartets from full list that will not contribute to
            // the G_ij / G_kl symmetry block


            int usi = usii;
            int usj = usjj;
            int usk = uskk;
            int usl = usll;

            /* place in "ascending" angular mom-
             my simple way of optimizing PHG recursion (VRR) */
            /* these first two are good for the HRR */
            switch_ij = switch_kl = switch_ijkl = false;
            if (BasisSet.shells[Symmetry.us2s[usi]].am < BasisSet.shells[Symmetry.us2s[usj]].am) {
              dum = usi;
              usi = usj;
              usj = dum;
              switch_ij = true;
            }
            if (BasisSet.shells[Symmetry.us2s[usk]].am < BasisSet.shells[Symmetry.us2s[usl]].am) {
              dum = usk;
              usk = usl;
              usl = dum;
              switch_kl = true;
            }
            /* this should be /good/ for the VRR */
            if (BasisSet.shells[Symmetry.us2s[usi]].am + BasisSet.shells[Symmetry.us2s[usj]].am
                > BasisSet.shells[Symmetry.us2s[usk]].am + BasisSet.shells[Symmetry.us2s[usl]].am) {
              dum = usi;
              usi = usk;
              usk = dum;
              dum = usj;
              usj = usl;
              usl = dum;
              switch_ijkl = true;
            }

            si = Symmetry.us2s[usi];
            sjj = Symmetry.us2s[usj];
            skk = Symmetry.us2s[usk];
            sll = Symmetry.us2s[usl];

            printf("To compute:\n");
            printf("G_ij = %d, G_kl = %d\n", G_ij, G_kl);
            printf("usii = %d, usjj = %d, uskk = %d, usll = %d\n", usii, usjj, uskk, usll);
            printf("si =   %d, sj =   %d, sk =   %d, sl =   %d\n", si, sjj, skk, sll);
            printf("span of usii: %d %d %d %d\n", BasisSet.shells[usii].span[0],
                BasisSet.shells[usii].span[1], BasisSet.shells[usii].span[2],
                BasisSet.shells[usii].span[3]);
            printf("span of usjj: %d %d %d %d\n", BasisSet.shells[usjj].span[0],
                BasisSet.shells[usjj].span[1], BasisSet.shells[usjj].span[2],
                BasisSet.shells[usjj].span[3]);
            printf("span of uskk: %d %d %d %d\n", BasisSet.shells[uskk].span[0],
                BasisSet.shells[uskk].span[1], BasisSet.shells[uskk].span[2],
                BasisSet.shells[uskk].span[3]);
            printf("span of usll: %d %d %d %d\n", BasisSet.shells[usll].span[0],
                BasisSet.shells[usll].span[1], BasisSet.shells[usll].span[2],
                BasisSet.shells[usll].span[3]);
            printf("centers: %d %d %d %d\n", center_i, center_j, center_k, center_l);
            fprintf(
                debugfile,
                "compute unique shell <%d %d||%d %d> for sym blk %d, %d using shell <%d %d||%d %d>\n",
                usii, usjj, uskk, usll, G_ij, G_kl, si, sjj, skk, sll);
            fflush(debugfile);
            fflush(stdout);

            center_i = BasisSet.shells[si].center - 1;
            center_j = BasisSet.shells[sjj].center - 1;
            center_k = BasisSet.shells[skk].center - 1;
            center_l = BasisSet.shells[sll].center - 1;

            // unchanged code block:
            /*--- Generate the petite list of shell quadruplets using DCD approach of Davidson ---*/
            usp_ij = &(Symmetry.us_pairs[usi][usj]);
            usp_kl = &(Symmetry.us_pairs[usk][usl]);
            stab_i = Symmetry.atom_positions[center_i];
            stab_j = Symmetry.atom_positions[center_j];
            stab_k = Symmetry.atom_positions[center_k];
            stab_l = Symmetry.atom_positions[center_l];
            stab_ij = Symmetry.GnG[stab_i][stab_j];
            stab_kl = Symmetry.GnG[stab_k][stab_l];
            R_list = Symmetry.dcr[stab_i][stab_j];
            S_list = Symmetry.dcr[stab_k][stab_l];
            T_list = Symmetry.dcr[stab_ij][stab_kl];
            lambda_T = Symmetry.nirreps / Symmetry.dcr_deg[stab_ij][stab_kl];
            ni = (BasisSet.puream ? 2 * BasisSet.shells[si].am - 1 : ioff[BasisSet.shells[si].am]);
            nj
                = (BasisSet.puream ? 2 * BasisSet.shells[sjj].am - 1
                    : ioff[BasisSet.shells[sjj].am]);
            nk
                = (BasisSet.puream ? 2 * BasisSet.shells[skk].am - 1
                    : ioff[BasisSet.shells[skk].am]);
            nl
                = (BasisSet.puream ? 2 * BasisSet.shells[sll].am - 1
                    : ioff[BasisSet.shells[sll].am]);

            memset(sj_arr, 0, sizeof(int) * max_num_unique_quartets);
            memset(sk_arr, 0, sizeof(int) * max_num_unique_quartets);
            memset(sl_arr, 0, sizeof(int) * max_num_unique_quartets);
            memset(sj_fbf_arr, 0, sizeof(int) * max_num_unique_quartets);
            memset(sk_fbf_arr, 0, sizeof(int) * max_num_unique_quartets);
            memset(sl_fbf_arr, 0, sizeof(int) * max_num_unique_quartets);
            count = 0;
            for (dcr_ij = 0; dcr_ij < Symmetry.dcr_dim[stab_i][stab_j]; dcr_ij++) {
              R = R_list[dcr_ij];
              sj = BasisSet.shells[sjj].trans_vec[R] - 1;
              for (dcr_ijkl = 0; dcr_ijkl < Symmetry.dcr_dim[stab_ij][stab_kl]; dcr_ijkl++) {
                T = T_list[dcr_ijkl];
                sk = BasisSet.shells[skk].trans_vec[T] - 1;
                slll = BasisSet.shells[sll].trans_vec[T] - 1;
                for (dcr_kl = 0; dcr_kl < Symmetry.dcr_dim[stab_k][stab_l]; dcr_kl++) {
                  S = S_list[dcr_kl];
                  sl = BasisSet.shells[slll].trans_vec[S] - 1;

                  total_am = BasisSet.shells[si].am + BasisSet.shells[sj].am
                      + BasisSet.shells[sk].am + BasisSet.shells[sl].am;
                  /*-------------------------------------------------------------
                   Obviously redundant or zero cases should be eliminated here!
                   Right now only zero case is eliminated. Redundancies arising
                   in DCD approach when usi == usj etc. may be eliminated too
                   but lambda_T will have to be replaced by an array (it won't
                   the same for every shell quartet in petite list anymore).
                   -------------------------------------------------------------*/
                  if (!(total_am % 2) || (BasisSet.shells[si].center != BasisSet.shells[sj].center)
                      || (BasisSet.shells[sj].center != BasisSet.shells[sk].center)
                      || (BasisSet.shells[sk].center != BasisSet.shells[sl].center)) {
                    sj_arr[count] = sj;
                    sk_arr[count] = sk;
                    sl_arr[count] = sl;
                    sj_fbf_arr[count] = BasisSet.shells[sj].fbf - 1;
                    sk_fbf_arr[count] = BasisSet.shells[sk].fbf - 1;
                    sl_fbf_arr[count] = BasisSet.shells[sl].fbf - 1;

                    count++;
                  }
                }
              }
            } /* petite list is ready to be used */
            num_unique_quartets = count;

            class_size = ni * nj * nk * nl;

            /* Determine irreps of all derivative operators whose application to this SO quartet is not identically zero */
            /* assume Abelian groups only, i.e. at most 8 irreps */
            int irreps_of_allowed_derivatives = CDSALCs.atom_irreps[center_i]
                | CDSALCs.atom_irreps[center_j] | CDSALCs.atom_irreps[center_k]
                | CDSALCs.atom_irreps[center_l];
            /* Determine the number of all SALC derivative operators which give nonzero when applied to the SO quartet */
            nuc[0] = center_i;
            nuc[1] = center_j;
            nuc[2] = center_k;
            nuc[3] = center_l;
            int salc_count = 0;
            for (int cd = 0; cd < num_coords; cd++)
              salc_all2thisquartet[cd] = -1;
            for (int c = 0; c < 4; c++) {
              int cart_der = 3 * nuc[c];
              for (int xyz = 0; xyz < 3; xyz++, cart_der++) {
                const CDSALC_t::cd2salc_map_t& cd2salc_map = CDSALCs.cd2salc_map[cart_der];
                const int nsalcs = cd2salc_map.nsalcs;
                for (int s = 0; s < nsalcs; s++) {
                  const int salc = cd2salc_map.salcs[s];
                  if (salc_all2thisquartet[salc] == -1) {
                    salc_all2thisquartet[salc] = salc_count;
                    salc_thisquartet2all[salc_count] = salc;
                    ++salc_count;
                  }
                }
              }
            }
            const int nsalcders_for_this_quartet = salc_count;

            np_i = BasisSet.shells[si].n_prims;
            np_j = BasisSet.shells[sjj].n_prims;
            np_k = BasisSet.shells[skk].n_prims;
            np_l = BasisSet.shells[sll].n_prims;

            orig_am[0] = BasisSet.shells[si].am - 1;
            orig_am[1] = BasisSet.shells[sjj].am - 1;
            orig_am[2] = BasisSet.shells[skk].am - 1;
            orig_am[3] = BasisSet.shells[sll].am - 1;
            am = orig_am[0] + orig_am[1] + orig_am[2] + orig_am[3];

            /*----------------------------------
             Compute the nonredundant quartets
             ----------------------------------*/
            for (plquartet = 0; plquartet < num_unique_quartets; plquartet++) {
              sj = sj_arr[plquartet];
              sk = sk_arr[plquartet];
              sl = sl_arr[plquartet];

              int center_j = BasisSet.shells[sj].center - 1;
              int center_k = BasisSet.shells[sk].center - 1;
              int center_l = BasisSet.shells[sl].center - 1;

              sp_ij = &(BasisSet.shell_pairs[si][sj]);
              sp_kl = &(BasisSet.shell_pairs[sk][sl]);

              Libderiv.AB[0] = sp_ij->AB[0];
              Libderiv.AB[1] = sp_ij->AB[1];
              Libderiv.AB[2] = sp_ij->AB[2];
              Libderiv.CD[0] = sp_kl->AB[0];
              Libderiv.CD[1] = sp_kl->AB[1];
              Libderiv.CD[2] = sp_kl->AB[2];

              AB2 = Libderiv.AB[0] * Libderiv.AB[0] + Libderiv.AB[1] * Libderiv.AB[1]
                  + Libderiv.AB[2] * Libderiv.AB[2];
              CD2 = Libderiv.CD[0] * Libderiv.CD[0] + Libderiv.CD[1] * Libderiv.CD[1]
                  + Libderiv.CD[2] * Libderiv.CD[2];

              /*--- Compute data for primitive quartets here ---*/
              num_prim_comb = 0;
              const double pfac = lambda_T;

              for (pi = 0; pi < np_i; pi++) {
                max_pj = (si == sj) ? pi + 1 : np_j;
                for (pj = 0; pj < max_pj; pj++) {
                  m = (1 + (si == sj && pi != pj));
                  for (pk = 0; pk < np_k; pk++) {
                    max_pl = (sk == sl) ? pk + 1 : np_l;
                    for (pl = 0; pl < max_pl; pl++) {
                      const int n = m * (1 + (sk == sl && pk != pl));

#ifdef USE_TAYLOR_FM
                      deriv1_quartet_data(&(Libderiv.PrimQuartet[num_prim_comb++]),
                          NULL, AB2, CD2,
                          sp_ij, sp_kl, am, pi, pj, pk, pl, n*pfac);
#else
                      deriv1_quartet_data(&(Libderiv.PrimQuartet[num_prim_comb++]), &fjt_table,
                          AB2, CD2, sp_ij, sp_kl, am, pi, pj, pk, pl, n * pfac);
#endif
                    }
                  }
                }
              }

              /*--- Compute the derivative integrals ---*/
              size = ioff[BasisSet.shells[si].am] * ioff[BasisSet.shells[sj].am]
                  * ioff[BasisSet.shells[sk].am] * ioff[BasisSet.shells[sl].am];
              build_deriv1_eri[orig_am[0]][orig_am[1]][orig_am[2]][orig_am[3]](&Libderiv,
                  num_prim_comb);

              // OK - Now have real data to manipulate
              // Libderiv.ABCD[x] is a pointer to the integrals.
              // The first 12 (x=0->11) point to first derivatives w/rt x1, y1, z1, x2, y2, z2, ... z4
              // The following 144 (x=12->143) would point to second
              // derivatives in normal order (x1x1, x1y1, x1z1, y1x1, y1y1, y1z1, etc.)
              // Because we've only build_deriv1_eri, only first 12 computed and used
              for (int c = 0; c < 3; c++) {
#ifdef NONDOUBLE_INTS
                for(int j=0;j<size;j++)
                cart_ints[c][j] = (double) Libderiv.ABCD[c][j];
#else
                cart_ints[c] = Libderiv.ABCD[c];
#endif
              }
              for (int c = 6; c < 12; c++) {
#ifdef NONDOUBLE_INTS
                for(int j=0;j<size;j++)
                cart_ints[c][j] = (double) Libderiv.ABCD[c][j];
#else
                cart_ints[c] = Libderiv.ABCD[c];
#endif
              }
              // cart_ints[c] now points to the integrals for a particular nuclear perturbation (x1, etc)

              // reconstruct integrals using translational invariance condition
              for (int j = 0; j < size; j++) {
                cart_ints[3][j] = -(cart_ints[0][j] + cart_ints[6][j] + cart_ints[9][j]);
                cart_ints[4][j] = -(cart_ints[1][j] + cart_ints[7][j] + cart_ints[10][j]);
                cart_ints[5][j] = -(cart_ints[2][j] + cart_ints[8][j] + cart_ints[11][j]);
              }
              // determine all unique centers (0..3)
              unique_center.resize(0);
              nuc[0] = center_i;
              nuc[1] = center_j;
              nuc[2] = center_k;
              nuc[3] = center_l;
              unique_center.push_back(0);
              if (nuc[1] != nuc[0]) {
                unique_center.push_back(1);
              }
              if (nuc[2] != nuc[0] && nuc[2] != nuc[1]) {
                unique_center.push_back(2);
              }
              if (nuc[3] != nuc[0] && nuc[3] != nuc[1] && nuc[3] != nuc[2]) {
                unique_center.push_back(3);
              }
              const int num_unique_centers = unique_center.size();

              // If two shells share a center, sum their derivatives.
              // Store them in the vector corresponding to the first
              // shell on that center, so that if, for example nuc[0] == nuc[2],
              // cart_ints[0..2] winds up = cart_ints[0..2] + cart_ints[6..8]
              // so that later cart_ints[6..8] are no longer needed.  This presumes
              // that the only use of these will be to sum them onto centers...

              // compute total derivatives from partial derivatives
              for (int i = 0; i < num_unique_centers; i++) {
                const int c = unique_center[i];
                for (int j = c + 1; j < 4; j++) {
                  if (nuc[c] == nuc[j]) {
                    int totderindex = 3 * c;
                    int parderindex = 3 * j;
                    for (int xyz = 0; xyz < 3; xyz++, totderindex++, parderindex++) {
                      double* tot_der = cart_ints[totderindex];
                      const double* part_der = cart_ints[parderindex];
                      for (int k = 0; k < size; k++)
                        tot_der[k] += part_der[k];
                    }
                  }
                }
              }
              const int num_tot_der = num_unique_centers * 3;

              // copy over total derivatives to plist_ints
              for (int uc = 0; uc < num_unique_centers; uc++) {
                const int c = unique_center[uc]; // loop over each unique center in shell quartet (0-4)
                int uniquetotderindex = 3 * uc; // where in plist_ints to place ints for that center dx
                int totderindex = 3 * c; // where in cart_ints to find ints for that center dx
                for (int xyz = 0; xyz < 3; xyz++, totderindex++, uniquetotderindex++) { // dx, dy, dz sequentially
                  double* target_ints = plist_ints[plquartet][uniquetotderindex]; // direct pointer to location in petite list to store values
                  double* source_ints = cart_ints[totderindex]; // direct pointer to cart_ints stored values
                  if (am) { // if it's not an <ss||ss> block
                    const double* data = norm_quartet(source_ints, target_ints, orig_am,
                        BasisSet.puream); // normalize a vector of values (one dx worth) and transform to puream if requested
                    if (data != target_ints)
                      memcpy(target_ints, data, sizeof(double) * class_size); // some variants return transformed ints in place, some in target_ints
                  } else { // if it /is/ an <ss||ss> block
                    target_ints[0] = source_ints[0];
                  }
                }
              }

              // It is /right here/ that we should have what we need to contract with TPDM elements.
              // Read a shell quartet of TPDM, manipulate
              // Loop over unique centers
              // Loop down integrals in quartet
              // Contract and sum into gradient matrix; complicated.  Needs to account for AO2SO transformation
              for (int uc = 0; uc < num_unique_centers; uc++) {
                const int c = unique_center[uc]; // loop over each unique center in shell quartet (0-4)
                int uniquetotderindex = 3 * uc; // where in plist_ints to place ints for that center dx

                int totderindex = 3 * c; // where in cart_ints to find ints for that center dx
                for (int xyz = 0; xyz < 3; xyz++, totderindex++, uniquetotderindex++) { // dx, dy, dz sequentially
                  double* source_ints = plist_ints[plquartet][uniquetotderindex]; // direct pointer to location in petite list with values

                  if (am) { // if it's not an <ss||ss> block
                    for (int k = 0; k < size; k++) {

                    }
                  } else { // if it /is/ an <ss||ss> block


                  }
                }
              }

#if 0
              //
              // compute nonzero contributions to SALCs generated by this quartet of unique shells
              //
              // not needed if adding cartesian components to gradient vector

              cdsalc_elem_vec* plist_salc = plist_salcs[plquartet];
              for (int i = 0; i < Symmetry.nirreps; i++)
              plist_salc[i].nelems = 0;
              for (int uc = 0; uc < num_unique_centers; uc++) {
                const int atom = nuc[unique_center[uc]];
                int abs_cart_der = atom * 3;
                int cart_der = uc * 3;
                for (int xyz = 0; xyz < 3; xyz++, abs_cart_der++, cart_der++) {
                  cdsalc_elem elem;
                  elem.cart_der = cart_der;

                  const CDSALC_t::cd2salc_map_t& cd2salc_map = CDSALCs.cd2salc_map[abs_cart_der];
                  const int nsalcs = cd2salc_map.nsalcs;
                  for (int salc = 0; salc < nsalcs; salc++) {
                    elem.salc = cd2salc_map.salcs[salc];
                    elem.coef = Symmetry.cdsalc2cd[abs_cart_der][elem.salc];
                    const int irrep = CDSALCs.salc2irrep[elem.salc];
                    plist_salc[irrep].elems[plist_salc[irrep].nelems++] = elem;
                  }
                }
              }
#endif

            } /* end of computing "petit" list */

#if 0
            const int bf_i_offset = BasisSet.shells[si].fbf - 1;
            //
            // Compute SO integrals and Contract with Density matrix
            //
            bool have_nonzero_integrals = false;

            /*---
             npi_ij - number of pairs of SOs arising from the ij pair of unique shells
             whose direct product transforms as G_ij
             npi_kl - number of pairs of SOs arising from the ij pair of unique shells
             whose direct product transforms as G_kl
             ---*/
            npi_ij = usp_ij->SOpair_npi[G_ij];
            npi_kl = usp_kl->SOpair_npi[G_kl];

            /* product of G_ij, G_kl, and irrep of the derivative operator must be totally symmetric */
            int irrep_deriv = G_ij ^ G_kl;

            /* is there a derivative operator which transforms as irrep_deriv? */
            int exists_nonzero_deriv = irreps_of_allowed_derivatives & (1 << irrep_deriv);
            if (!exists_nonzero_deriv) {
              printf("We should never have arrived here, with no derivatives that count!\n");
              continue;
            }

            for (ij = 0; ij < npi_ij; ij++) { /*--- Loop over SO pairs from usij ---*/
              i = usp_ij->SOpair_bf_i[G_ij][ij]; /*--- BF index ---*/
              j = usp_ij->SOpair_bf_j[G_ij][ij];
              so_i = usp_ij->SOpair_so_i[G_ij][ij]; /*--- Absolute index of this SO from usi ---*/
              so_j = usp_ij->SOpair_so_j[G_ij][ij];
              ind_offset = (i * nj + j) * nk * nl;

              for (kl = 0; kl < npi_kl; kl++) { /*--- Loop over SO pairs from uskl ---*/
                k = usp_kl->SOpair_bf_i[G_kl][kl];
                l = usp_kl->SOpair_bf_j[G_kl][kl];
                so_k = usp_kl->SOpair_so_i[G_kl][kl];
                so_l = usp_kl->SOpair_so_j[G_kl][kl];

                index = ind_offset + k * nl + l; /* position of this integral in quartet */

                // zero out target SO integrals
                for (int salc = 0; salc < nsalcders_for_this_quartet; salc++)
                so_int[salc] = 0.0;

                // compute target SO integrals
                for (s = 0; s < num_unique_quartets; s++) { /*--- Sum over petite list quartets ---*/
                  const int bf_i = bf_i_offset + i; /*--- Absolute basis function index ---*/
                  const int bf_j = sj_fbf_arr[s] + j;
                  const int bf_k = sk_fbf_arr[s] + k;
                  const int bf_l = sl_fbf_arr[s] + l;
                  double so2ao_pfac = Symmetry.usotao[so_i][bf_i] * Symmetry.usotao[so_j][bf_j]
                  * Symmetry.usotao[so_k][bf_k] * Symmetry.usotao[so_l][bf_l];

                  /* loop over all contributions to derivative operator SALCs within this irrep */
                  const cdsalc_elem_vec& salc_elems = plist_salcs[s][irrep_deriv];
                  const int ncontribs = salc_elems.nelems;
                  for (int derop = 0; derop < ncontribs; derop++) {
                    const cdsalc_elem& elem = salc_elems.elems[derop];

                    //if memset is fast then use absolute salc index, otherwise relative, i.e. [0,12), index
                    const int salc = salc_all2thisquartet[elem.salc];
                    const double ao_integral = plist_ints[s][elem.cart_der][index];
                    const double so_contribution = so2ao_pfac * elem.coef * ao_integral;
                    so_int[salc] += so_contribution;

                    // Print out the contribution if needed
                    //if (UserOptions.print_lvl > PRINT_DEBUG) {
                    if (0) {
                      fprintf(outfile, "  -AO integral contribution:\n");
                      fprintf(
                          outfile,
                          "   AO integral    -- deriv wrt coord %d ( %d %d | %d %d ) = %20.10lf\n",
                          elem.cart_der, bf_i, bf_j, bf_k, bf_l, ao_integral);
                      fprintf(outfile, "   SO pfac = %lf  deriv SALC pfac = %lf\n", so2ao_pfac,
                          elem.coef);
                      fprintf(
                          outfile,
                          "   SO integral -- deriv wrt coord %d ( %d %d | %d %d ) += %20.10lf = %20.10lf\n\n",
                          salc_thisquartet2all[salc], so_i, so_j, so_k, so_l, so_contribution,
                          so_int[salc]);
                    }
                  }
                }
                // for symmetry blocked integrals, need to store
                // buffers of integrals indexed by three symmetry irreps: G[ij], G[kl], G[salc]
                // For D2h, then need 8x8x8 tensor of buffers if storing
                // add nonzero target integrals to the appropriate buffer

                for (int salc = 0; salc < nsalcders_for_this_quartet; salc++) {
                  double value = so_int[salc];
                  if (fabs(value) > toler) {
                    const int abs_salc = salc_thisquartet2all[salc];
                    const int curr_integral = num_of_ints_in_totdata[abs_salc];
                    struct tebuf& tbuf = tot_data[abs_salc][curr_integral];
                    tbuf.val = value;
                    tbuf.i = (short int) so_i;
                    tbuf.j = (short int) so_j;
                    tbuf.k = (short int) so_k;
                    tbuf.l = (short int) so_l;
                    ++num_of_ints_in_totdata[abs_salc];
                    have_nonzero_integrals = true;
                    fprintf(outfile, "d/dr%d (%d %d || %d %d) = %lf\n", abs_salc, so_i, so_j, so_k, so_l, value);
                  }
                }

                //
                // Write out the integrals
                //
                if (have_nonzero_integrals) { /* Let's see if we need to write out something */

                  // for each integral buffer which could have contributions from this SO quartet
                  for (int salc = 0; salc < nsalcders_for_this_quartet; salc++) {
                    const int abs_salc = salc_thisquartet2all[salc];
                    const int nints = num_of_ints_in_totdata[abs_salc];
                    if (nints > 0) {
                      te_count_per_coord[abs_salc] += nints;
                      total_te_count += nints;
                      iwl_buf_wrt_struct_nocut(&(D1ERIOUT[abs_salc]), tot_data[abs_salc], nints);

                      // print out the integrals, if needed
                      if (UserOptions.print_lvl >= PRINT_DEBUG) {
                        for (int i = 0; i < nints; i++)
                        fprintf(outfile, "%5d%5d%5d%5d%5d%20.10lf\n", abs_salc,
                            tot_data[abs_salc][i].i, tot_data[abs_salc][i].j,
                            tot_data[abs_salc][i].k, tot_data[abs_salc][i].l,
                            tot_data[abs_salc][i].val);
                      }

                      num_of_ints_in_totdata[abs_salc] = 0;
                    }

                  }
                }
              }
            }/* end getting unique shell combination */

            for (int c = 0; c < num_coords; c++) {
              iwl_buf_flush(&D1ERIOUT[c], 1);
              iwl_buf_close(&D1ERIOUT[c], 1);
            }
#endif
            //{
            // Safest way to print 64-bit integers is to use std::ostringstream
            std::ostringstream oss;
            oss << total_te_count;
            fprintf(outfile, "\n    Wrote %s first-derivative ERIs to IWL files [%d, %d]\n\n",
                oss.str().c_str(), IOUnits.itapD1ERI_SO, IOUnits.itapD1ERI_SO + num_coords - 1);
            //}


          } // usll ...
        } // uskk ...
      } // usjj (unique shell j index)
    } // usii (unique shell i index)
                } } // G_i, G_k (G_j = G_i ^ G_ij, etc)
  } // G_ij (G_kl = G_ij ^ G_target)


  /*---------
   Clean-up
   ---------*/
  free_libderiv(&Libderiv);
  free(sj_arr);
  free(sk_arr);
  free(sl_arr);
  if (Symmetry.nirreps > 1) {
    free(sj_fbf_arr);
    free(sk_fbf_arr);
    free(sl_fbf_arr);
  }
#ifdef USE_TAYLOR_FM
  free_Taylor_Fm_Eval();
#else
  free_fjt_table(&fjt_table);
  free_fjt();
#endif
  fclose(debugfile); // jtf
  exit(1);

  return;
} // end te_deriv1_corr_symm

} // end namespace cints
} // end namespace psi
