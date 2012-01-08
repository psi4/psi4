#ifndef ADC_H
#define ADC_H

#include <libmints/matrix.h>
#include <libmints/vector.h>
#include <libmints/wavefunction.h>
#include <libdpd/dpd.h>

#define ID(x) _ints->DPD_ID(x)
#define PSIF_ADC       260
#define PSIF_ADC_SEM   261

#define DEBUG_         false
#define ANGL_TOL_      50.0
#define CUTOFF_DENS_   1e-6

namespace boost {
    template<class T> class shared_ptr;
}

namespace psi{ 

class Options;
class PSIO;
class Chkpt;
class Matrix;
class Vector;
class IntegralTransform;
typedef boost::shared_ptr<Matrix> SharedMatrix;
typedef boost::shared_ptr<Vector> SharedVector;

namespace adc{
    
struct pole{
    int    iter;              // Iterated time 
    double iter_value;        // Converged value of the excitation energy
    double ps_value;          // Pseudo-pertirbative value of the excitation energy
    double osc_strength;      // Oscillator strength
    double renorm_factor;     // Residue of the propagator, which is identical to the squared norm of the singly excited vector
    double rot_angle;         // Rotation angle from corresponding CIS vector
};
        
class ADC: public Wavefunction
{
public:
    ADC();
    ~ADC();
    double compute_energy();
    virtual bool same_a_b_orbs() const { return reference_wavefunction_->same_a_b_orbs(); }
    virtual bool same_a_b_dens() const { return reference_wavefunction_->same_a_b_dens(); }

protected:
    void init();
    void release_mem();
    void rhf_prepare_tensors();
    void amps_write(dpdfile2 *B, int length, FILE *outfile);
    void onestack_insert(struct onestack *stack, double value, int i, int a, int level, int staclen);
    double rhf_init_tensors();
    double rhf_differentiate_omega(int irrep, int root);
    void rhf_diagonalize(int irrep, int num_root, bool first, double omega_in, double *eps);
    void rhf_construct_sigma(int irrep, int root);
    void shift_denom2(int root, int irrep, double omega);
    void shift_denom4(int irrep, double omega);

    // Number of the singly excited configurations
    int nxs_;
    // Number of singly occupied orbitals
    int nopen_;
    // Convergence criteria in Newton-Raphson procedure
    int conv_;
    // Maximum iteration number in Newton-Raphson procedure
    int pole_max_;
    // MAximum iteration number in simultaneous expansion method
    int sem_max_;
    // Norm tolerance for the residual vector 
    int norm_tol_;
    // Number of components of transition amplitudes printed in outfile
    int num_amps_;
    // Number of alpha active occupied MOs per irrep
    int *aoccpi_;
    // Number of alpha active virtual MOs per irrep
    int *avirpi_;
    // Number of beta active occupied MOs per irrep
    int *boccpi_;
    // Number of beta active virtual MOs per irrep
    int *bvirpi_;
    // Number of doubly occupied MOs per irrep
    int *clsdpi_;
    // Roots sought per irrep
    int *rpi_;
    // Number of sngly excited configurations per irrep
    int *nxspi_;
    // Irreps for X, Y and Z
    int *irrep_axis_;
    // Ground state energy
    double gs_energy_; 
    // An array containing alpha occupied orbital energies in DPD order
    double *aocce_;
    // An array containing alpha virtual orbital energies in DPD order
    double *avire_;
    // An array containing beta occupied orbital energies in DPD order
    double *bocce_;
    // An array containing beta virtual orbital energies in DPD order
    double *bvire_;
    // Array of struct that contains all the information on the pole of the propagator
    struct pole **poles_;
    // Integral transformation object
    IntegralTransform *_ints;
    // Guesses for the correlated excitation energies, which are given as CIS/ADC(1) energies
    SharedVector omega_guess_;
};
    
}}

#endif

