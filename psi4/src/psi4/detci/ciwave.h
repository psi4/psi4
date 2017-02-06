/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

#ifndef CIWAVE_H
#define CIWAVE_H

 #include "psi4/pragma.h"
 PRAGMA_WARNING_PUSH
 PRAGMA_WARNING_IGNORE_DEPRECATED_DECLARATIONS
 #include <memory>
 PRAGMA_WARNING_POP
#include "psi4/libmints/dimension.h"
#include "psi4/libmints/wavefunction.h"

// From psi4
namespace psi {
class Options;
class JK;
class DFERI;
class IntegralTransform;
class MOSpace;
typedef std::shared_ptr<Matrix> SharedMatrix;
typedef std::shared_ptr<Matrix> SharedMatrix;
class SOMCSCF;

// Well this is not ideal
struct _SlaterDetSet;
typedef _SlaterDetSet SlaterDetSet;
}

// From the detci module
namespace psi { namespace detci {
class CIvect;
class SlaterDeterminant;
struct calcinfo;
struct params;
struct stringwr;
struct ci_blks;
struct olsen_graph;
struct graph_set;
struct H_zero_block;
struct detci_timings;
struct mcscf_params;
typedef std::shared_ptr<psi::detci::CIvect> SharedCIVector;
}}

namespace psi { namespace detci {


class CIWavefunction : public Wavefunction
{

public:
    CIWavefunction(std::shared_ptr<Wavefunction> reference_wavefunction);
    CIWavefunction(std::shared_ptr<Wavefunction> reference_wavefunction, Options &options);
    virtual ~CIWavefunction();


    double compute_energy();

    /// Simple accessors
    size_t ndet();

    /**!
     * Similar to wavefunction.Ca_subset(); however, this version knows about all of the CI
     * subspaces in the SO basis. We stick to the MCSCF definitions for now.
     * @param  orbital_name fzc, drc, docc, act, ras1, ras2, ras3, ras4, pop, vir, fzv, drv, or all
     * @return C            Returns the appropriate orbitals in the SO basis.
     */
    SharedMatrix get_orbitals(const std::string &orbital_name);

    /**!
     * Similar to wavefunction.Ca_subset(); however, this version knows about all of the CI
     * subspaces in the SO basis. We stick to the MCSCF definitions for now.
     * @param  orbital_name FZC, DRC, DOCC, ACT, RAS1, RAS2, RAS3, RAS4, POP, VIR, FZV, DRV, or ALL
     * @param  orbitals     SharedMatrix to set
     */
    void set_orbitals(const std::string& orbital_name, SharedMatrix orbitals);

    /**!
     * Gets the dimension of the desired subspace.
     * @param  orbital_name FZC, DRC, DOCC, ACT, RAS1, RAS2, RAS3, RAS4, POP, VIR, FZV, DRV, or ALL
     * @return dim          Dimension object
     */
    Dimension get_dimension(const std::string& orbital_name);

    /**!
     * Transform the one and two electron integrals for a CI computation.
     */
    void transform_ci_integrals(void);

    /**!
     * Transform the one and two electron integrals for a MCSCF computation.
     */
    void transform_mcscf_integrals(bool approx_only);

    /**
     * Rotates the integrals for MCSCF computations, places the integrals in CI order.
     * @param K         Rotation matrix (active x N)
     * @param onel_out  one-electron output vector
     * @param twoel_out two-electron output vector
     */
    void rotate_mcscf_integrals(SharedMatrix K, SharedVector onel_out,
                                SharedVector twoel_out);

    /**
     * Takes a pitzer order act x act Matrix and converts it into a ci-ordered vector that
     * CI sigma vectors use.
     * @param src  Source Matrix - can have irreps
     * @param dest Destination Vector - should be of size, ncitri = (nact * (nact + 1)) / 2
     */
    void pitzer_to_ci_order_onel(SharedMatrix src, SharedVector dest);

    /**
     * Takes a pitzer order act^2 x act^2 Matrix and converts it into a ci-ordered vector that
     * CI sigma vectors use.
     * @param src  Source Matrix - Cannot have irreps
     * @param dest Destination Vector - should be of size, (ncitri * (nctri + 1)) / 2
     */
    void pitzer_to_ci_order_twoel(SharedMatrix src, SharedVector dest);

    /**!
     * Obtains the OPDM <Iroot| Epq |Jroot> from the ciwave object. If Jroot is
     * negative then Iroot == Jroot, if both roots are -1 then the "special" CI
     * OPDM is returned.
     * @param Iroot      Left root
     * @param Jroot      Right root
     * @param spin       Selects which spin to return: A, B, or SUM
     * @param full_space If false return only the active OPDM else return full OPDM
     * @return OPDM or TDM shared matrix
     **/
    SharedMatrix get_opdm(int Iroot=-1, int Jroot=-1, const std::string& spin="SUM",
                           bool full_space=false);

    /**!
     * Compute the one-particle density matrix between two CIVectors
     * @param Ivec       The SharedCIVector for I
     * @param Jvec       The SharedCIVector for J
     * @param Jroot      Jroot to use
     * @param Iroot      Iroot to use
     * @return           AA, BB, and summed OPDM's
    **/
    std::vector<SharedMatrix> opdm(SharedCIVector Ivec, SharedCIVector Jvec,
                                    int Iroot, int Jroot);

    /**!
     * Obtains the "special" TPDM, other TPDM roots are not held here.
     * The AA (BB) order: \Gamma_{pqrs} & = \langle 0 | a_{p \alpha}^\dagger a_{r \alpha}^{\dagger} a_{s \alpha} a_{q \alpha} | 0 \rangle
     * The AB order: \Gamma_{pqrs} & = \langle 0 | a_{p \alpha}^\dagger a_{r \beta}^{\dagger} a_{s \beta} a_{q \alpha} | 0 \rangle
     * @param spin       Selects which spin to return AA, AB, BB, or SUM
     * @param symmetrize Symmetrize the TPDM, only works for SUM currently
     * @return           The request 4D active TPDM
     **/
    SharedMatrix get_tpdm(const std::string& spin = "SUM", bool symmetrize=true);

    /**!
     * Compute the two-particle density matrix between two CIVectors
     * @param Ivec       The SharedCIVector for I
     * @param Jvec       The SharedCIVector for J
     * @param Jroot      Jroot to use
     * @param Iroot      Iroot to use
     * @return           AA, AB, BB, and summed TPDM's
    **/
    std::vector<SharedMatrix> tpdm(SharedCIVector Ivec, SharedCIVector Jvec, int Iroot,
                                                   int Jroot);

    /**!
     Builds and returns a new CIvect object.
     * @param maxnvect  Maximum number of vectors in the CIvector
     * @param filenum   File number to use for the on-disk data
     * @param use_disk  If false CIvect disables read/write
     * @param buf_init  If true initializes the buffers internally
     * @return          SharedCIVector object
     **/
    SharedCIVector new_civector(int maxnvect, int filenum, bool use_disk=true,
                                bool buf_init=true);
    /**
     * Returns the "D" vector that contains the current reference CI Wavefunction
     * @return The "D" vector
     */
    SharedCIVector D_vector();

    /**
     * Builds a CIVector that is the diagonal of the Hamiltonian
     * @param hd_type The type of diagonal to be taken, if -1 defaults to current parameters
     *                 0 = HD_EXACT
     *                 1 = HD_KAVE
     *                 2 = ORB_ENER
     *                 3 = EVANGELISTI
     *                 4 = LEININGER
     *                 5 = Z_HD_KAVE
     * @return The Hamiltonian diagonal
     */
    SharedCIVector Hd_vector(int hd_type = -1);

    /**
     * Compute a sigma vector
     * @param C    Input vector
     * @param S    Output vector
     * @param cvec Which vector number to use for the C vec
     * @param svec Which vector number to use for the S vec
     */
    void sigma(SharedCIVector C, SharedCIVector S, int cvec, int svec);

    /**
     * Compute a sigma vector
     * @param C    Input vector
     * @param S    Output vector
     * @param cvec Which vector number to use for the C vec
     * @param svec Which vector number to use for the S vec
     * @param oei  One-electron integrals to use
     * @param tei  Two-electron integrals to use
     */
    void sigma(SharedCIVector C, SharedCIVector S, int cvec, int svec,
               SharedVector oei, SharedVector tei);

    /**
     * Prints the large values in a CIVector
     * @param vec  CIVector to print
     * @param root Which root?
     */
    void print_vector(SharedCIVector vec, int root);

    /**
     * Compute the state-transfer operator. Currently in construction and not to be used.
     */
    void compute_state_transfer(SharedCIVector ref, int ref_vec,
                                SharedMatrix prop, SharedCIVector ret);

    // Compute functions
    void compute_cc();
    int diag_h(double econv = -1.0, double rconv = -1.0);
    void compute_mpn();

    // Build CI quantities

    /**
     * Forms the "special" OPDM of the current Dvec.
     */
    void form_opdm();

    /**
     * Forms the "special" TPDM of the current Dvec.
     */
    void form_tpdm();

    /**
     * Form CI Natural Orbitals (Warning! must build the OPDM beforehand)
     */
    void ci_nat_orbs();

    // Cleanup Functions

    /**
     * Cleans up the strings, H0Block, and CIVectors associated with this Wavefunction.
     */
    void cleanup_ci();

    /**
     * Cleans up the DPD library and associated integrals associated with this Wavefunction.
     */
    void cleanup_dpd();

    /**
     * Sets the diag_h guess option. !Expert option.
     * @param guess CI Guess: (UNIT, H0_BLOCK, or DFLIE)
     */
    void set_ci_guess(std::string guess);

    // Returns a new SOMCSCF object
    std::shared_ptr<SOMCSCF> mcscf_object();

    /// Functions below this line should be used for debug use only

    /**
     *  Builds the full CI hamiltonian for debugging purposes. Currently limits itself to a matrix
     * of 40% of Psi4's memory in size.
     * @param  hsize Compute the first hsize elements of the CI Hamiltonian
     * @return       The requested CI Hamiltonian.
     */
    SharedMatrix hamiltonian(size_t hsize = 0);

private:

    /// => General Helper Functions <= ///

    /// Paramater and CalcInfo setters
    void get_mo_info();
    void get_parameters(Options &options);
    void print_parameters();
    void set_ras_parameters();
    void print_ras_parameters();

    // General setup
    void title(bool is_mcscf);
    void form_strings();
    void set_ciblks();
    void convergence_death();

    /// Sets the ciwavefunction object
    void common_init();
    bool cleaned_up_ci_;

    /// Find out which orbitals belong hwere
    void orbital_locations(const std::string& orbital_name, int* start, int* end);

    /// Symmetry block a matrix
    SharedMatrix symm_block(SharedMatrix x, Dimension dim1, Dimension dim2);

    /// => Integrals <= ///
    bool ints_init_;
    bool df_ints_init_;
    bool mcscf_object_init_;
    bool fzc_fock_computed_;
    std::shared_ptr<IntegralTransform> ints_; // Non-DF
    std::shared_ptr<MOSpace> rot_space_;
    std::shared_ptr<MOSpace> act_space_;
    std::shared_ptr<DFERI> dferi_; // DF
    std::shared_ptr<JK> jk_;
    std::shared_ptr<SOMCSCF> somcscf_;

    /// General transforms
    void init_mcscf_object();
    void tf_onel_ints(SharedVector onel, SharedVector twoel, SharedVector output);
    void form_gmat(SharedVector onel, SharedVector twoel, SharedVector output);
    void onel_ints_from_jk();
    double get_twoel(int i, int j, int k, int l);
    double get_onel(int i, int j);

    /// Non-DF integral functions
    void setup_mcscf_ints();
    void transform_mcscf_ints(bool approx_only = false);
    void read_dpd_ci_ints();
    void rotate_mcscf_twoel_ints(SharedMatrix K, SharedVector twoel_out);

    /// DF integral functions
    void setup_dfmcscf_ints();
    void transform_dfmcscf_ints(bool approx_only = false);
    void rotate_dfmcscf_twoel_ints(SharedMatrix K, SharedVector twoel_out);

    /// MO transformation through TeraCHEM-like
    /// Added in by Kevin Patrick Hannon
    void transform_mcscf_ints_ao(bool approx_only = false);
    void setup_mcscf_ints_ao();
    SharedMatrix tei_raaa_;
    SharedMatrix tei_aaaa_;

    /// => Globals <= //
    struct stringwr **alplist_;
    struct stringwr **betlist_;
    struct calcinfo *CalcInfo_;
    struct params *Parameters_;
    struct ci_blks *CIblks_;
    struct olsen_graph *AlphaG_;
    struct olsen_graph *BetaG_;
    struct H_zero_block *H0block_;
    int **s1_contrib_, **s2_contrib_, **s3_contrib_;
    int ***OV_;
    unsigned char ***Occs_;

    /// => H0block functions <= //
    void H0block_init(unsigned int size);
    void H0block_free(void);
    void H0block_print(void);
    int  H0block_calc(double E);
    void H0block_gather(double **mat, int al, int bl, int cscode, int mscode, int phase);
    void H0block_xy(double *x, double *y, double E);
    void H0block_setup(int num_blocks, int *Ia_code, int *Ib_code);
    void H0block_pairup(int guess);
    void H0block_spin_cpl_chk(void);
    void H0block_filter_setup(void);
    void H0block_fill(void);
    void H0block_coupling_calc(double E);
    std::string print_config(int nbf, int num_alp_el, int num_bet_el,
	  struct stringwr *stralp, struct stringwr *strbet, int num_drc_orbs);

    /// => Slater Matrix Elements <= //
    double matrix_element(SlaterDeterminant* I, SlaterDeterminant* J);
    int sme_first_call_;
    int *same_alpha_, *same_beta_;
    int *common_docc_, *common_alp_socc_, *common_bet_socc_;
    int init_nalp_, init_nbet_;
    int **I_diff_, **J_diff_;

    /// => CI Iterators <= //
    void mitrush_iter(CIvect &Hd, struct stringwr **alplist, struct stringwr
          **betlist, int nroots, double *evals, double conv_rms, double conv_e,
          double enuc, double edrc, int maxiter, int maxnvect);
    void olsen_update(CIvect &C, CIvect &S, CIvect &Hd, double E, double E_est,
          double *norm, double *c1norm, double *ovrlap, double *buffer1,
          double *buffer2,
          int curr, int next, std::string out, int iter, struct stringwr **alplist,
          struct stringwr **betlist);
    void olsen_iter_xy(CIvect &C, CIvect &S, CIvect &Hd, double *x, double *y,
          double *buffer1, double *buffer2, double E, int curvect, int L,
          double **alpha, struct stringwr **alplist, struct stringwr **betlist);
    void mitrush_update(CIvect &C, CIvect &S, double norm, double acur,
       double alast, double *buffer1, double *buffer2, int curr, int next);

    void sem_iter(CIvect &Hd, struct stringwr **alplist, struct stringwr
      **betlist, double *evals, double conv_e,
      double conv_rms, double enuc, double edrc,
      int nroots, int maxiter, int maxnvect);
    void parse_import_vector(SlaterDetSet *sdset, int *ialplist, int *ialpidx,
      int *ibetlist, int *ibetidx, int *blknums);
    void sem_test(double **A, int N, int M, int L, double **evecs, double *evals,
          double **b, double conv_e, double conv_rms, int maxiter, double offst,
          int *vu, int maxnvect);

    /// => Sigma Calculations <= //
    struct sigma_data *SigmaData_;
    void sigma_init(CIvect& C, CIvect &S);
    void sigma_free(void);
    void sigma(CIvect& C, CIvect& S, double *oei, double *tei, int ivec);

    void sigma_a(struct stringwr **alplist, struct stringwr **betlist,
          CIvect& C, CIvect& S, double *oei, double *tei, int fci, int ivec);
    void sigma_b(struct stringwr **alplist, struct stringwr **betlist,
          CIvect& C, CIvect& S, double *oei, double *tei, int fci, int ivec);
    void sigma_c(struct stringwr **alplist, struct stringwr **betlist,
          CIvect& C, CIvect& S, double *oei, double *tei, int fci, int ivec);

    void sigma_block(struct stringwr **alplist, struct stringwr **betlist,
          double **cmat, double **smat, double *oei, double *tei, int fci,
          int cblock, int sblock, int nas, int nbs, int sac, int sbc,
          int cac, int cbc, int cnas, int cnbs, int cnac, int cnbc,
          int sbirr, int cbirr, int Ms0);
    void sigma_get_contrib(struct stringwr **alplist, struct stringwr **betlist,
          CIvect &C, CIvect &S, int **s1_contrib, int **s2_contrib,
          int **s3_contrib);
        void form_ov();
    void sigma_get_contrib_rotf(CIvect &C, CIvect &S,
          int **s1_contrib, int **s2_contrib, int **s3_contrib,
          int *Cnt[2], int **Ij[2], int **Oij[2], int **Ridx[2],
          signed char **Sgn[2], unsigned char **Toccs);

    void print_vec(unsigned int nprint, int *Ialist, int *Iblist,
          int *Iaidx, int *Ibidx, double *coeff);

    /// => MCSCF helpers <= //


    /// => MPn helpers <= //
    void mpn_generator(CIvect &Hd);

    /// => Density Matrix helpers <= //
    std::vector<std::vector<SharedMatrix> > opdm(SharedCIVector Ivec, SharedCIVector Jvec,
                                                std::vector<std::tuple<int, int> > states_vec);
    SharedMatrix opdm_add_inactive(SharedMatrix opdm, double value, bool virt=false);
    void opdm_block(struct stringwr **alplist, struct stringwr **betlist,
            double **onepdm_a, double **onepdm_b, double **CJ, double **CI, int Ja_list,
            int Jb_list, int Jnas, int Jnbs, int Ia_list, int Ib_list,
            int Inas, int Inbs);

    // OPDM holders, opdm_map holds lots of active-active opdms
    // opdm_, opdm_a_, etc are for "the" current OPDM
    bool opdm_called_;
    std::map<std::string, SharedMatrix> opdm_map_;
    SharedMatrix opdm_;
    SharedMatrix opdm_a_;
    SharedMatrix opdm_b_;

    std::vector<SharedMatrix> tpdm(SharedCIVector Ivec, SharedCIVector Jvec,
                                   std::vector<std::tuple<int, int, double> > states_vec);
    void tpdm_block(struct stringwr **alplist, struct stringwr **betlist,
            int nbf, int nalplists, int nbetlists,
            double *twopdm_aa, double *twopdm_bb, double *twopdm_ab, double **CJ, double **CI, int Ja_list,
            int Jb_list, int Jnas, int Jnbs, int Ia_list, int Ib_list,
            int Inas, int Inbs, double weight);

    bool tpdm_called_;
    SharedMatrix tpdm_;
    SharedMatrix tpdm_aa_;
    SharedMatrix tpdm_ab_;
    SharedMatrix tpdm_bb_;

}; // End CIWavefunction

}}

#endif // CIWAVE_H
