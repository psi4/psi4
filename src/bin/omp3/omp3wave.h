#ifndef omp3wave_h
#define omp3wave_h

#include <libmints/wavefunction.h>
#include <libdpd/dpd.h>
#include "arrays.h"

using namespace std;

namespace psi{ 
  
class IntegralTransform;
  
namespace omp3wave {

class OMP3Wave : public Wavefunction
{
 
    void common_init();
    
public:
    OMP3Wave(boost::shared_ptr<Wavefunction> reference_wavefunction, Options &options);

    virtual ~OMP3Wave();
    virtual bool same_a_b_orbs() const { return reference_wavefunction_->same_a_b_orbs(); }
    virtual bool same_a_b_dens() const { return reference_wavefunction_->same_a_b_dens(); }

    virtual double compute_energy();

protected:
    void mem_release();
    void response_pdms();
    void GFockmo();
    void mograd();
    void update_mo();
    void korbrot_sd();
    void occ_iterations();
    void t2_1st_sc(); 
    void t2_1st_general();
    void t2_2nd_sc(); 
    void t2_2nd_general();
    void mp2_energy();
    void mp3_energy();
    void twopdm_corr_opdm();
    void twopdm_ref();
    void twopdm_oovv();
    void twopdm_oooo();
    void twopdm_vvvv();
    void twopdm_ovov();
    void twopdm_vovo();
    void twopdm_ovvo();
    void nbo();
    void mp3l_energy();
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
    void diis(int dimvec, Array2d *vecs, Array2d *errvecs, Array1d *vec_new);
    void W_1st_order();
    void V_2nd_order();
    
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
     double Emp3;
     double Emp2BB;
     double Emp2AA;
     double Emp2AB;
     double Emp3BB;
     double Emp3AA;
     double Emp3AB;
     double Emp3L;
     double Emp3L_old;
     double Ecorr;
     double EcorrL;
     double Ecorrnew;
     double Emp3_rdm;
     double Escsmp2;
     double Escsmp2BB;
     double Escsmp2AA;
     double Escsmp2AB;
     double Esosmp2AB; 
     double Esosmp2;
     double Escsmp3BB;
     double Escsmp3AA;
     double Escsmp3AB;
     double Escsmp3;
     double Esosmp3AB; 
     double Esosmp3;
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
     double sos_scale_omp3;
     double e3_scale;
     
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
     string omp3_orb_energy;
     string do_scs;		// Spin-Component-Scaling
     string do_sos;		// Spin-Opposite-Scaling
     string compute_mp3l;	// Do compute mp3l energy during iterations?
     string twopdm_abcd_type;	// How to handle G_abcd
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

     double *evalsA; 
     double *evalsB; 
     double *evals_c1A; 
     double *evals_c1B; 

     Array1d *wogA; 
     Array1d *wogB; 
     Array1d *kappaA; 
     Array1d *kappaB; 
     Array1d *kappa_barA; 
     Array1d *kappa_barB;   
     
     Array2d *vecsA;
     Array2d *vecsB;
     Array2d *errvecsA;
     Array2d *errvecsB;

     double **C_pitzerA;     
     double **C_pitzerB;     
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
    
};

} }

#endif // omp3wave_h

