#ifndef omp2wave_h
#define omp2wave_h

#include <libmints/wavefunction.h>
#include <libdpd/dpd.h>

using namespace std;

namespace psi{ 
  
class IntegralTransform;
  
namespace omp2wave {

class OMP2Wave : public Wavefunction
{
 
    void common_init();
    
public:
    OMP2Wave(boost::shared_ptr<Wavefunction> reference_wavefunction, Options &options);

    virtual ~OMP2Wave();
    virtual bool same_a_b_orbs() const { return reference_wavefunction_->same_a_b_orbs(); }
    virtual bool same_a_b_dens() const { return reference_wavefunction_->same_a_b_dens(); }

    virtual double compute_energy();

protected:
    void mem_release();
    void response_pdms();
    void GFockmo();
    void mograd();
    //void orbital_hess();
    void update_mo();
    //void korbrot_nr();
    void korbrot_sd();
    //void korbrot_aughess();
    void occ_iterations();
    void t2_1st_sc(); 
    void t2_1st_general();
    void mp2_norm();
    void mp2_energy();
    void twopdm_corr_opdm();
    void twopdm_ref();
    void twopdm_oovv();
    void nbo();
    void mp2l_energy();
    //void idp4stable();
    //void stable();
    void G_int();
    void get_moinfo();
    void title();
    void semi_canonic();
    void trans_ints();
    void denominators();
    void ref_energy();
    void Fockmo_alpha();
    void Fockmo_beta();
    void idp();
    double dot_product(int dim, double *x, double *y);
    void diis(int dimvec, double **vecs, double **errvecs, double *vec_new);
    
     class IntegralTransform *ints;
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
     int dimvvAA; 
     int dimvvBB;
     int dimooAA; 
     int dimooBB; 
     int dimvoAA; 
     int dimvoBB; 
     int dimvoAB; 
     int dimvoBA; 
     int dimss; 
     int dimosA;	// O-all 
     int dimosB;	// o-all 
     int dimvsA; 
     int dimvsB; 
     int dimacsacs;
     int dimacoacoAA;
     int dimacoacoBB;
     int dimacvacvAA;
     int dimacvacvBB;
     int dimacvacoAA;
     int dimacvacoAB;
     int dimacvacoBA;
     int dimacvacoBB;
     int dimacofoAA;
     int dimacofoBB;
     int dimfvacvAA;
     int dimfvacvBB;
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
     int nidpA;
     int nidpB;
     int conver;
     //int stability_flag;
     int mo_optimized; 		// if 0 MOs are NOT optimized, if 1 MOs are optimized.
     int swap_mo_i;
     int swap_mo_j;

     
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
     double Ecorrnew;
     double Emp2_adm;
     double Escsmp2;
     double Escsmp2BB;
     double Escsmp2AA;
     double Escsmp2AB;
     double Esosmp2AB; 
     double Esosmp2;
     double Eopdm;
     double DE;
     double tol_Eod;
     double tol_Emp2;
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
     double cc_norm;
     double lshift_parameter;
     double cutoff;
     double os_scale;
     double ss_scale;
     double sos_scale;
     double sos_scale2;
     double a_cg;
     double b_cg;
     
     string wfn;
     string reference;
     string jobtype;
     string dertype;
     string subgroup;
     string volume;
     string basis;
     string level_shift;
     string lineq;
     string orth_type;
     string stability;
     string natorb;
     string rotation_blocks;
     string semicanonic; 
     string wfn_norm;
     string opt_method; 
     string hess_type;
     string omp2_orb_energy;
     string do_scs;		// Spin-Component-Scaling
     string do_sos;		// Spin-Opposite-Scaling
     string write_mo_coeff;	// Write CmoA to CmoA.psi and CmoB to CmoB.psi
     string read_mo_coeff;	// Read CmoA from CmoA.psi and CmoB from CmoB.psi
     string swap_mo;		// Swap phi_i with phi_j


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

     
     /*
     int **idprowA;
     int **idprowB;
     int **idpcolA;
     int **idpcolB;
     */
     
     /*
     SharedIntVector nidpA;
     SharedIntVector nidpB;
     SharedIntVector idprowA;
     SharedIntVector idpcolA;
     SharedIntVector idprowB;
     SharedIntVector idpcolB;
     SharedVector wogA;
     SharedVector wogB;
     SharedVector kappaA;
     SharedVector kappaB;
     */

     double *evalsA; 
     double *evalsB; 
     double *evals_c1A; 
     double *evals_c1B; 
     double *wogA; 
     double *wogB; 
     double *kappaA; 
     double *kappaB; 
     double *kappa_barA; 
     double *kappa_barB;      

     double **C_pitzerA;     
     double **C_pitzerB;     
     double **vecsA; 
     double **vecsB; 
     double **errvecsA;
     double **errvecsB;
     char **irreplabels; 
     
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
     SharedMatrix HG1A;
     SharedMatrix HG1B;
     SharedMatrix gamma1corrA;
     SharedMatrix gamma1corrB;
     SharedMatrix g1symmA;
     SharedMatrix g1symmB;
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
     SharedMatrix BmatA;
     SharedMatrix BmatB;
     //SharedSimpleMatrix AorbAA;
     //SharedSimpleMatrix AorbAB;
     //SharedSimpleMatrix AorbBB;

    
};

} }

#endif // omp2wave_h

