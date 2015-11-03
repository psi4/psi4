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

#ifndef CIWAVE_H
#define CIWAVE_H

// Forward declarations
namespace boost {
template<class T> class shared_ptr;
}

// From psi4
namespace psi {
class Wavefunction;
class Options;
class JK;
class DFERI;
class IntegralTransform;
class MOSpace;
typedef boost::shared_ptr<Matrix> SharedMatrix;
}

// From the detci module
namespace psi { namespace detci {
class CIvect;
struct calcinfo;
struct params;
struct stringwr;
struct ci_blks;
struct olsen_graph;
struct graph_set;
struct H_zero_block;
struct detci_timings;
struct mcscf_params;
}}

namespace psi { namespace detci {

// Need to nuke this eventually
extern struct stringwr **alplist;
extern struct stringwr **betlist;

class CIWavefunction : public Wavefunction
{

public:
    CIWavefunction(boost::shared_ptr<Wavefunction> reference_wavefunction);
    CIWavefunction(boost::shared_ptr<Wavefunction> reference_wavefunction, Options &options);
    virtual ~CIWavefunction();

    double compute_energy();

    /**!
     * Similar to wavefunction.Ca_subset(); however, this version knows about all of the CI
     * subspaces in the SO basis. We stick to the MCSCF definitions for now.
     * @param  orbital_name FZC, DOCC, ACT, VIR, FZV
     * @return C            Returns the appropriate orbitals in the SO basis.
     */
    SharedMatrix get_orbitals(const std::string& orbital_name);

    /**!
     * Similar to wavefunction.Ca_subset(); however, this version knows about all of the CI
     * subspaces in the SO basis. We stick to the MCSCF definitions for now.
     * @param  orbital_name FZC, DOCC, ACT, VIR, FZV
     * @param  orbitals     SharedMatrix to set
     * @return C            Returns the appropriate orbitals in the SO basis.
     */
    void set_orbitals(const std::string& orbital_name, SharedMatrix orbitals);

    /**!
     * Gets the dimension of the desired subspace.
     * @param orbital_name FZC, DOCC, ACT, VIR, FZV
     * @return dim         Dimension object
     */
    Dimension get_dimension(const std::string& orbital_name);
    /**
     * Transform the one and two electron integrals.
     * The
     */
    void transform_ci_integrals(void);

    /**!
     * Obtains the OPDM <root| Epq |root> from disk
     * @param root       Root to obtain
     * @param spin       0 returns alpha, 1 returns beta, 2 sums alpha and beta
     * @param transden   Obtain <0 | Epq | root> transition density matrix
     * @return OPDM or TDM shared matrix
     **/
    SharedMatrix get_opdm(int root=0, int spin=2, bool transden=false);

    /**!
     Returns the symmetry block active OPDM
     **/
    SharedMatrix get_active_opdm();

    /**!
     * Loads the OPDM to the wavefunction Da_ and Db_
     * @param use_old_d - If no new OPDM is calculated use the reference density matrices
    **/
    void set_opdm(bool use_old_d=false);

    /**!
     * Obtains the TPDM from disk
     * @param symmetrized Symmetrize the TPDM or not
     * @return TPDM SharedVector
     **/
    SharedVector get_tpdm(bool symmetrized=true, const std::string& = "SUM");

    /**!
     Returns the a full 4D active TPDM
     **/
    SharedMatrix get_active_tpdm(const std::string& tpdm_type = "SUM");

    // Compute functions
    void compute_mcscf();
    void compute_cc();
    void diag_h();
    void compute_mpn();

    // Build CI quantities
    void form_opdm();
    void form_tpdm();

    // Extraneous
    void cleanup();

private:

    /// => General Helper Functions <= ///

    /// Paramater and CalcInfo setters
    void get_mo_info();
    void get_parameters(Options &options);
    void get_mcscf_parameters();
    void print_parameters();
    void set_ras_parameters();
    void print_ras_parameters();

    // General setup
    void title();
    void init_ioff();
    void form_strings();
    void set_ciblks();

    /// Sets the ciwavefunction object
    void common_init();

    /// Find out which orbitals belong hwere
    void orbital_locations(const std::string& orbital_name, int* start, int* end);

    /// Symmetry block a matrix
    SharedMatrix symm_block(SharedMatrix x, Dimension dim1, Dimension dim2);

    /// => Integrals <= ///
    bool ints_init_;
    bool df_ints_init_;
    boost::shared_ptr<IntegralTransform> ints_; // Non-DF
    boost::shared_ptr<MOSpace> rot_space_;
    boost::shared_ptr<MOSpace> act_space_;
    boost::shared_ptr<DFERI> dferi_; // DF
    boost::shared_ptr<JK> jk_;

    /// General transforms
    void tf_onel_ints();
    void form_gmat();
    void onel_ints_from_jk();

    /// Non-DF integral functions
    void setup_mcscf_ints();
    void transform_mcscf_ints(bool approx_only = false);
    void read_dpd_ci_ints();

    /// DF integral functions
    void setup_dfmcscf_ints();
    void transform_dfmcscf_ints(bool approx_only = false);


    /// => Globals <= //
    struct stringwr **alplist_;
    struct stringwr **betlist_;
    struct mcscf_params *MCSCF_Parameters;
    struct calcinfo *CalcInfo_;
    struct params *Parameters_;
    struct ci_blks *CIblks_;
    struct olsen_graph *AlphaG_;
    struct olsen_graph *BetaG_;
    struct H_zero_block *H0block_;
    int **s1_contrib_, **s2_contrib_, **s3_contrib_;
    int ***OV_;

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
    void H0block_fill();
    void H0block_coupling_calc(double E);

    /// => CI Iterators <= //
    //void mitrush_iter(CIvect &Hd, int nroots, double *evals, double conv_rms, double conv_e,
    //  double enuc, double edrc, int maxiter, int maxnvect, std::string out,
    //  int print_lvl);
    void mitrush_iter(CIvect &Hd, struct stringwr **alplist, struct stringwr
          **betlist, int nroots, double *evals, double conv_rms, double conv_e,
          double enuc, double edrc, int maxiter, int maxnvect, std::string out,
          int print_lvl);
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
      int nroots, int maxiter, int maxnvect, std::string out, int print_lvl);

    void sigma_init(CIvect& C, CIvect &S, struct stringwr **alplist,
          struct stringwr **betlist);
    void sigma(struct stringwr **alplist, struct stringwr **betlist,
          CIvect& C, CIvect& S, double *oei, double *tei, int fci, int ivec);
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



    /// => MPn helpers <= //
    void mpn_generator(CIvect &Hd); 

    /// => Density Matrix helpers <= //
    void opdm(struct stringwr **alplist, struct stringwr **betlist,
              int transdens, int dipmom,
              int Inroots, int Iroot, int Inunits, int Ifirstunit,
          int Jnroots, int Jroot, int Jnunits, int Jfirstunit,
          int targetfile, int writeflag, int printflag);
    void opdm_block(struct stringwr **alplist, struct stringwr **betlist,
            double **onepdm_a, double **onepdm_b, double **CJ, double **CI, int Ja_list,
            int Jb_list, int Jnas, int Jnbs, int Ia_list, int Ib_list,
            int Inas, int Inbs);
    void opdm_ave(int targetfile);
    void opdm_ke(double **onepdm);
    void tpdm(struct stringwr **alplist, struct stringwr **betlist,
          int Inroots, int Inunits, int Ifirstunit,
          int Jnroots, int Jnunits, int Jfirstunit,
          int targetfile, int writeflag, int printflag);
    void tpdm_block(struct stringwr **alplist, struct stringwr **betlist,
            int nbf, int nalplists, int nbetlists,
            double *twopdm_aa, double *twopdm_bb, double *twopdm_ab, double **CJ, double **CI, int Ja_list,
            int Jb_list, int Jnas, int Jnbs, int Ia_list, int Ib_list,
            int Inas, int Inbs, double weight);

}; // End CIWavefunction

}}

#endif // CIWAVE_H

