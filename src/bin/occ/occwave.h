#ifndef occwave_h
#define occwave_h

#include <libmints/wavefunction.h>
#include <libdiis/diismanager.h>
#include <libdpd/dpd.h>
#include "arrays.h"

using namespace std;

namespace psi{ 
  
class IntegralTransform;
  
namespace occwave {

class OCCWave : public Wavefunction
{
 
    void common_init();
    
public:
    OCCWave(boost::shared_ptr<Wavefunction> reference_wavefunction, Options &options);

    virtual ~OCCWave();
    virtual bool same_a_b_orbs() const { return reference_wavefunction_->same_a_b_orbs(); }
    virtual bool same_a_b_dens() const { return reference_wavefunction_->same_a_b_dens(); }

    virtual double compute_energy();

protected:
    // General
    void mem_release();
    void mograd();
    void update_mo();
    void ccl_energy();
    void nbo();
    void get_moinfo();
    void title();
    void semi_canonic();
    void ref_energy();
    void fock_alpha();
    void fock_beta();
    void idp();
    void diis(int dimvec, Array2d *vecs, Array2d *errvecs, Array1d *vec_new);
    void kappa_msd();
    void kappa_orb_resp();
    void kappa_orb_resp_iter();
    void orb_resp_pcg_rhf();
    void orb_resp_pcg_uhf();
    void dump_ints();
    void denominators_rhf();
    void denominators_uhf();
    void gfock();
    void trans_ints_rhf();
    void trans_ints_uhf();
    void tpdm_ref();
    void tpdm_corr_opdm();
    void tpdm_oovv();
    void tpdm_oooo();
    void tpdm_ovov();
    void tpdm_vovo();
    void tpdm_ovvo();
    void gfock_diag();
    void coord_grad();
    void dump_pdms();
    void occ_iterations();

    // OMP2
    void omp2_g_int();
    void omp2_response_pdms();
    void omp2_t2_1st_sc(); 
    void omp2_t2_1st_general();
    void omp2_tpdm_oovv();
    void omp2_manager();
    void omp2_mp2_energy();

    // OMP3
    void omp3_manager();
    void omp3_response_pdms();
    void omp3_t2_1st_sc(); 
    void omp3_t2_1st_general();
    void omp3_tpdm_vvvv();
    void omp3_g_int();
    void omp3_mp2_energy();
    void w_1st_order();
    void v_2nd_order();
    void t2_2nd_sc(); 
    void t2_2nd_general();
    void mp3_energy();

    // OCEPA
    void ocepa_manager();
    void cepa_manager();
    void cepa_iterations();
    void ocepa_mp2_energy();
    void ocepa_t2_1st_sc(); 
    void ocepa_tpdm_vvvv();
    void ocepa_g_int();
    void ocepa_response_pdms();
    void t2_amps(); 
    void w_int();
    void v_int();
    void cepa_energy();
    
     class IntegralTransform *ints;
     DIISManager *t2DiisManager;
     //class DIISManager t2DiisManager;

     int nmo;		// Number of MOs
     int nao;		// Number of AOs
     int nso;		// Number of SOs
     int nooA;		// Number of alpha occupied orbitals
     int nooB;		// Number of beta occupied orbitals
     int nvoA;		// Number of alpha virtual orbitals
     int nvoB;		// Number of beta virtual orbitals
     int nacooA;	// Number of active alpha occupied orbitals
     int nacooB;	// Number of active beta occupied orbitals
     int nacso;		// Number of active SOs
     int nacvoA;	// Number of active alpha virtual orbitals
     int nacvoB;	// Number of active beta virtual orbitals
     int nirreps;	// Number of irreducible representations
     int nshell;	// Number of shells
     int nfrzc; 	// Number of frozen cores
     int nfrzv; 	// Number of frozen virtuals
     int npop; 		// Number of populated orbitals: npop=nmo-nfrzv
     int dimtei;	// dimension of tei in pitzer order for all integrals 
     int ntri; 		// square matrix dimension (nmo) -> pitzer order
     int ntri_so;	// square matrix dimension (nso) -> pitzer order
     int cc_maxiter;
     int mo_maxiter;
     int exp_tol_Eod;
     int exp_tol_t2;
     int exp_tol_grad;
     int exp_idp_cutoff;
     int exp_mograd_max;
     int exp_cutoff;
     int itr_occ;
     int nhessroot;
     int idp_return;
     int idp_returnA;
     int idp_returnB;
     int num_vecs; 		// Number of vectors used in diis (diis order)
     int nvar; 			// nvar = num_vecs +1;
     int memory; 
     int multp; 
     int charge;
     int print_;
     int cachelev; 
     int nidp;
     int nidp_tot;             // nidpA + nidpB
     int nidpA;
     int nidpB;
     int conver;
     int mo_optimized; 		// if 0 MOs are NOT optimized, if 1 MOs are optimized.
     int itr_pcg;
     int idp_idx;
     int pcg_maxiter;
     int pcg_conver;
     int do_diis_; 
     int itr_diis;
     int time4grad;             // If 0 it is not the time for grad, if 1 it is the time for grad
     int cc_maxdiis_; 		// MAX Number of vectors used in CC diis
     int cc_mindiis_; 		// MIN Number of vectors used in CC diis

     // Common
     double Enuc;
     double sum;
     double Etotal;
     double Eelec;
     double Escf;
     double Eref;
     double Emp2;
     double Emp2BB;
     double Emp2AA;
     double Emp2AB;
     double Emp2L;
     double Emp2L_old;
     double Ecorr;
     double EcorrL;
     double Ecc_rdm;
     double Escsmp2;
     double Escsmp2BB;
     double Escsmp2AA;
     double Escsmp2AB;
     double Esosmp2AB; 
     double Esosmp2;
     double Escsnmp2;
     double Escsnmp2BB;
     double Escsnmp2AA;
     double Escsmimp2;
     double Escsmimp2BB;
     double Escsmimp2AA;
     double Escsmimp2AB;
     double Escsmp2vdw;
     double Escsmp2vdwBB;
     double Escsmp2vdwAA;
     double Escsmp2vdwAB;
     double Esospimp2AB; 
     double Esospimp2;
     double Eopdm;
     double DE;
     double tol_Eod;
     double tol_grad;
     double idp_cutoff;
     double rms_kappa;
     double rms_kappaA;
     double rms_kappaB;
     double rms_wog;
     double rms_wogA;
     double rms_wogB;
     double step_max;
     double mograd_max;
     double biggest_mograd;
     double biggest_mogradA;
     double biggest_mogradB;
     double biggest_kappa;
     double biggest_kappaA;
     double biggest_kappaB;
     double tol_t2;
     double rms_t2;
     double rms_t2AA;
     double rms_t2AB;
     double rms_t2BB;
     double rms_l2;
     double mu_ls;
     double sc_ls;
     double lshift_parameter;
     double cutoff;
     double os_scale;
     double ss_scale;
     double sos_scale;
     double sos_scale2;
     double a_pcgA;
     double a_pcgB;
     double b_pcgA;
     double b_pcgB;
     double rms_pcgA;
     double rms_pcgB;
     double rms_pcg;
     double tol_pcg;

     // OMP3
     double e3_scale;
     double Emp3;
     double Emp3BB;
     double Emp3AA;
     double Emp3AB;
     double Emp3L;
     double Emp3L_old;
     double Escsmp3BB;
     double Escsmp3AA;
     double Escsmp3AB;
     double Escsmp3;
     double Esosmp3AB; 
     double Esosmp3;
     double Escsnmp3;
     double Escsmimp3;
     double Escsmp3vdw;
     double Esospimp3;

     // OCEPA 
     double Ecepa;
     double Ecepa_old;
     double EcepaAA;
     double EcepaBB;
     double EcepaAB;
     double EcepaL;
     double EcepaL_old;
     double EscscepaBB;
     double EscscepaAA;
     double EscscepaAB;
     double Escscepa;
     double EsoscepaAB; 
     double Esoscepa;
     double Escsncepa;
     double Escsmicepa;
     double Escscepavdw;
     double Esospicepa;
     double cepa_os_scale_;
     double cepa_ss_scale_;
     double cepa_sos_scale_;
     double sos_scale_ocepa;
     
     string wfn;
     string reference;
     string reference_;
     string jobtype;
     string dertype;
     string basis;
     string level_shift;
     string lineq;
     string orth_type;
     string natorb;
     string semicanonic; 
     string opt_method; 
     string hess_type;
     string occ_orb_energy;
     string do_scs;		// Spin-Component-Scaling
     string do_sos;		// Spin-Opposite-Scaling
     string write_mo_coeff;	// Write CmoA to CmoA.psi and CmoB to CmoB.psi
     string read_mo_coeff;	// Read CmoA from CmoA.psi and CmoB from CmoB.psi
     string scs_type_;		
     string sos_type_;		
     string pcg_beta_type_;		
     string compute_mp3l;	// Do compute mp3l energy during iterations?
     string compute_cepal;	// Do compute cepal energy during iterations?
     string twopdm_abcd_type;	// How to handle G_abcd
     string wfn_type_;
     string compute_ccl;
     string orb_resp_solver_;
     string bypass_contract442_;


     int *mopi; 		/* number of all MOs per irrep */
     int *sopi; 		/* number of all SOs per irrep */
     int *occpi;
     int *doccpi; 		/* number of doubly occupied MOs per irrep */
     int *occpiA;		/* number of alpha occupied MOs per irrep */ 
     int *occpiB;		/* number of beta occupied MOs per irrep */ 
     int *soccpi; 		/* number of all singly occupied MOs per irrep */  
     int *virtpiA;		/* number of alpha virtual MOs per irrep */ 
     int *virtpiB;		/* number of beta virtual MOs per irrep */ 
     int *frzcpi; 		/* number of frozen occupied MOs per irrep */
     int *frzvpi; 		/* number of frozen virtual MOs per irrep */  
     int *adoccpi; 		/* number of active doubly occupied MOs per irrep */
     int *aoccpiA; 		/* number of active alpha occupied MOs per irrep */ 
     int *aoccpiB; 		/* number of active beta occupied MOs per irrep */ 
     int *avirtpiA; 		/* number of active alpha virtual MOs per irrep */
     int *avirtpiB; 		/* number of active beta virtual MOs per irrep */
     int *mosym; 		/* symmetry of all MOs in pitzer order */
     int *sosym; 		/* symmetry of all SOs in pitzer order */
     int *mosym_c1; 		/* symmetry of all MOs in energy order */
     int *PitzerOffset; 	/* block offset */
     int *pitzer2symblk;	// convert Pitzer index to sym block index
     int *pitzer2symirrep;	// Return irrep of given SO in Pitzer order
     int *occ2symblkA;		// convert OCC index to sym block index
     int *occ2symblkB;		// convert OCC index to sym block index
     int *virt2symblkA;		// convert VIR index to sym block index
     int *virt2symblkB;		// convert VIR index to sym block index
     int *c1topitzerA; 		// Convert the index of given alpha orbital in C1 subgroup to symmetric subgroup (pitzer ordered)
     int *c1topitzerB; 		// Convert the index of given beta orbital in C1 subgroup to symmetric subgroup (pitzer ordered)
     int *pitzer2c1A;		// Convert the index of given alpha orbital in symmetric subgroup to C1 subgroup (energy ordered)  
     int *pitzer2c1B;		// Convert the index of given alpha orbital in symmetric subgroup to C1 subgroup (energy ordered)  
     int *c1toqtA; 		// Convert the index of given orbital in C1 subgroup to QT order
     int *c1toqtB; 		// Convert the index of given orbital in C1 subgroup to QT order
     int *qt2c1A; 		// Convert the index of given orbital in QT order to C1 subgroup (energy ordered) 
     int *qt2c1B; 		// Convert the index of given orbital in QT order to C1 subgroup (energy ordered) 
     int *qt2pitzerA; 		// Convert the index of given orbital in QT order to symmetric subgroup (pitzer ordered)
     int *qt2pitzerB; 		// Convert the index of given orbital in QT order to symmetric subgroup (pitzer ordered)
     int *pitzer2qtA;  		// Convert the index of given orbital in symmetric subgroup to QT order    
     int *pitzer2qtB;  		// Convert the index of given orbital in symmetric subgroup to QT order    
     int *occ_offA;		/* Alpha OCC block offset */
     int *occ_offB;		/* Beta OCC block offset */
     int *vir_offA;		/* Alpha VIR block offset */
     int *vir_offB;		/* Alpha VIR block offset */
     int *idprowA;
     int *idprowB;
     int *idpcolA;
     int *idpcolB;
     int *idpirrA;
     int *idpirrB;

     double *evalsA; 
     double *evalsB; 
     double *evals_c1A; 
     double *evals_c1B; 
     
     Array1d *wogA; 
     Array1d *wogB; 
     Array1d *kappaA; 
     Array1d *kappaB; 
     Array1d *kappa;          // where kappa = kappaA + kappaB
     Array1d *kappa_barA; 
     Array1d *kappa_barB;   
     Array1d *kappa_newA; 
     Array1d *kappa_newB; 
     Array1d *r_pcgA; 
     Array1d *r_pcgB; 
     Array1d *z_pcgA; 
     Array1d *z_pcgB; 
     Array1d *p_pcgA; 
     Array1d *p_pcgB; 
     Array1d *sigma_pcgA; 
     Array1d *sigma_pcgB; 
     Array1d *Minv_pcgA; 
     Array1d *Minv_pcgB; 
     Array1d *r_pcg_newA; 
     Array1d *r_pcg_newB; 
     Array1d *z_pcg_newA; 
     Array1d *z_pcg_newB; 
     Array1d *p_pcg_newA; 
     Array1d *p_pcg_newB; 
     Array1d *dr_pcgA; 
     Array1d *dr_pcgB; 
     
     Array2d *vecsA;
     Array2d *vecsB;
     Array2d *errvecsA;
     Array2d *errvecsB;
     Array2d *AorbAA;
     Array2d *AorbBB;
     Array2d *AorbAB;
     Array2d *Aorb;


     double **C_pitzerA;     
     double **C_pitzerB;     
     
     SharedMatrix Ca_new;	// New Alpha MO coeff. 
     SharedMatrix Cb_new;	// New Beta MO coeff. 
     SharedMatrix Ca_ref;	// Ref Alpha MO coeff. 
     SharedMatrix Cb_ref;	// Ref Beta MO coeff. 
     SharedMatrix Tso;
     SharedMatrix Vso;
     SharedMatrix Hso;
     SharedMatrix HmoA;
     SharedMatrix HmoB;
     SharedMatrix FsoA;
     SharedMatrix FsoB;
     SharedMatrix FockA;
     SharedMatrix FockB;
     SharedMatrix GFock;
     SharedMatrix GFockA;
     SharedMatrix GFockB;
     SharedMatrix UorbA;
     SharedMatrix UorbB;
     SharedMatrix UorbrotA;
     SharedMatrix UorbrotB;
     SharedMatrix KorbA;
     SharedMatrix KorbB;
     SharedMatrix KsqrA;
     SharedMatrix KsqrB;
     SharedMatrix HG1;
     SharedMatrix HG1A;
     SharedMatrix HG1B;
     SharedMatrix gamma1corr;
     SharedMatrix gamma1corrA;
     SharedMatrix gamma1corrB;
     SharedMatrix g1symm;
     SharedMatrix g1symmA;
     SharedMatrix g1symmB;
     SharedMatrix Worb;
     SharedMatrix WorbA;
     SharedMatrix WorbB;
     SharedMatrix HCA;
     SharedMatrix HCB;
     SharedMatrix FCA;
     SharedMatrix FCB;
     SharedMatrix GooA;
     SharedMatrix GooB;
     SharedMatrix GvvA;
     SharedMatrix GvvB;
    
};

} }

#endif // occwave_h

