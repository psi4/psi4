#ifndef DFGF_H
#define DFGF_H

#include <libmints/matrix.h>
#include <libmints/vector.h>
#include <libmints/wavefunction.h>

#define DEBUG_         false 
#define ANGL_TOL_      50.0

namespace boost {
    template<class T> class shared_ptr;
}

namespace psi{ 

class Options;
class PSIO;
class Chkpt;
class Matrix;
typedef boost::shared_ptr<Matrix> SharedMatrix;

namespace plugin_dfadc{
/*    
struct pole{
    int    iter;              // Iterated time 
    double iter_value;        // Converged value of the excitation energy
    double ps_value;          // Pseudo-pertirbative value of the excitation energy
    double osc_strength;      // Oscillator strength
    double renorm_factor;     // Residue of the propagator, which is identical to the squared norm of the singly excited vector
    double rot_angle;         // Rotation angle from corresponding CIS vector
};
*/    
enum Order {Irow, Icol};
        
class DFADC: public Wavefunction
{
public:
    DFADC();
    ~DFADC();
    double compute_energy();
    virtual bool same_a_b_orbs() const { return true; }
    virtual bool same_a_b_dens() const { return true; }

protected:
    // Initialization of parameters 
    void init();
    // Release all the intermediates and etc.
    void release_mem();
    // Diagonalize response matrix
    void diagonalize(double *&eps, double **&V, int nroot, bool first, bool CISguess);
    // Construct sigma tensor
    void sigma_tensor(double **B, int dim, double **&S, bool do_ADC);
    // Construct DF tensor
    // Order species whether the tensor is prepared as, Q^I_{pq}, or Q^{pq}_I
    // where I and p,q stand for the DF index and MO indices.
    void formDFtensor(SharedMatrix Cl, SharedMatrix Cr, double **J_mhalf, enum Order order, double **&Bout);
    // Construct J^{-1/2} tensor
    void formInvSqrtJ(double **&J_mhalf);
    // Returns the differentiated energy
    double accelerate(double *V);
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
    // Number of alpha active occupied states
    int naocc_;
    // Number of alpha active virtual states
    int navir_;
    // Number of Roots
    int num_roots_;
    // Cut off criterion for amplitudes
    int cutoff_;
    // Initial dimension of the Ritz space
    int init_ritz_;
    // Maximum dimension of the Ritz space
    int max_ritz_;
    // Current energy for the pole searching process
    double omega_;
    // Pseudo-perturbative energy
    double omega_ps_;
    // An array containing occupied orbital energies in DPD order
    double *occe_;
    // An array containing virtual orbital energies in DPD order
    double *vire_;
    // An array that contains diagonal zeroth order term in CIS matrix
    double *diag_;
    // CIS excitation energy
    double *Ecis_;
    // CIS wavefunction
    double **Bcis_;
    // An intermediate tensor apears as controbution from a 3h-3p diagram
    double **Aij_;
    // An intermediate tensor apears as controbution from a 3h-3p diagram
    double **Aab_;    
    // MP1 amplitude
    double **Kiajb_;
    // Q^I_{ij} tensor
    double **Qpij_;
    // Q^I_{ia} tensor
    double **Qpia_;
    // Q^I_{ab} tensor
    double **Qpab_;
    // Information about pole structure
//    struct pole *poles_;
    // SO2occpiedMO transformation coefficient
    SharedMatrix occCa_;
    // SO2virtualMO transformation coefficient
    SharedMatrix virCa_;
    // Dipole integrals in MO basis
    std::vector<SharedMatrix> dipole_ints_;
    // RI Basis
    boost::shared_ptr<BasisSet> ribasis_;
};
    
}}

#endif

