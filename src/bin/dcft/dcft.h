#ifndef _PSI4_SRC_BIN_DCFT_DCFT_H_
#define _PSI4_SRC_BIN_DCFT_DCFT_H_
#include "psi4-dec.h"
#include <libmints/mints.h>
#include <libchkpt/chkpt.hpp>
#include <libpsio/psio.hpp>
#include <libscf_solver/uhf.h>
#include <libtrans/integraltransform.h>
#include <libdpd/dpd.h>
#include <libdiis/diismanager.h>
#include <libmints/matrix.h>


namespace psi{
    using namespace scf;
 namespace dcft{

class DCFTSolver
{
  public:
    DCFTSolver(Options &options);
    ~DCFTSolver();

    void compute();

  protected:
    Options           _options;
    shared_ptr<Chkpt> _chkpt;
    shared_ptr<PSIO>  _psio;
    IntegralTransform *_ints;

    void perform_scf();
    void transform_integrals();
    void build_lambda();
    void init_moinfo();
    void free_moinfo();
    void read_checkpoint();
    void compute_energy();
    void update_lambda_from_residual();
    void compute_scf_energy();
    void mp2_guess();
    void build_tau();
    void print_opdm();
    void write_orbitals_to_checkpoint();
    void check_n_representability();
    void print_orbital_energies();
    void find_occupation(Vector &evals, bool forcePrint = false);
    void build_intermediates();
    void process_so_ints();
    void build_G();
    void build_tensors();
    void build_denominators();
    void dump_density();
    void dpd_buf4_add(dpdbuf4 *A, dpdbuf4 *B, double alpha);
    void scf_guess();
    void half_transform(dpdbuf4 *A, dpdbuf4 *B, Matrix C1, Matrix C2,
            int *mospi_left, int *mospi_right, int **so_row, int **mo_row,
            bool backwards, double alpha, double beta);
    void file2_transform(dpdfile2 *A, dpdfile2 *B, Matrix *C, bool backwards);
    void AO_contribute(dpdbuf4 *tau1_AO, dpdbuf4 *tau2_AO, int p, int q,
            int r, int s, double value, dpdfile2* = NULL, dpdfile2* = NULL, dpdfile2* = NULL);
    //void AO_contribute(dpdfile2 *tau1_AO, dpdfile2 *tau2_AO, int p, int q,
    //        int r, int s, double value);
    bool correct_mo_phases(bool dieOnError = true);
    double compute_lambda_residual();
    double compute_scf_error_vector();
    double update_scf_density();
    /// Whether the DOCC array was provided by the user
    bool _inputDocc;
    /// Whether the SOCC array was provided by the user
    bool _inputSocc;
    /// The maximum number of lambda iterations per update
    int _lambdaMaxIter;
    /// The maximum number of SCF iterations per update
    int _scfMaxIter;
    /// The number of irreps in the current point group
    int _nIrreps;
    /// The amount of information to print
    int _print;
    /// The number of ymmetrized atomic orbitals
    int _nSo;
    /// The number of molecular orbitals
    int _nMo;
    /// The number of unique pairs of symmetrized atomic orbitals
    int _nTriSo;
    /// The number of occupied alpha orbitals
    int _nAOcc;
    /// The number of occupied beta orbitals
    int _nBOcc;
    /// The number of virtual alpha orbitals
    int _nAVir;
    /// The number of virtual beta orbitals
    int _nBVir;
    /// The number of iterations needed to reach lambda convergence
    int _nLambdaIterations;
    /// The number of iterations needed to reach scf convergence
    int _nScfIterations;
    /// The maximum size of the DIIS subspace
    int _maxDiis;
    /// The number of DIIS vectors needed for extrapolation
    int _minDiisVecs;
    /// The maximum number of iterations
    int _maxNumIterations;
    /// The number of unpaired electrons per irrep
    int *_openPI;
    /// The number of doubly-occupied orbitals per irrep
    int *_clsdPI;
    /// The number of frozen core orbitals per irrep
    int *_frzcPI;
    /// The number of frozen virtual orbitals per irrep
    int *_frzvPI;
    /// The number of molecular orbitals per irrep
    int *_moPI;
    /// The number of symmetrized atomic orbitals per irrep
    int *_soPI;
    /// The number of occupied alpha orbitals per irrep
    int *_nAOccPI;
    /// The number of occupied beta orbitals per irrep
    int *_nBOccPI;
    /// The number of virtual alpha orbitals per irrep
    int *_nAVirPI;
    /// The number of virtual beta orbitals per irrep
    int *_nBVirPI;
    /// The nuclear repulsion energy in Hartree
    double _eNuc;
    /// The cutoff below which and integral is assumed to be zero
    double _intTolerance;
    /// The RMS value of the error vector after the SCF iterations
    double _scfConvergence;
    /// The RMS value of the change in lambda after the lambda iterations
    double _lambdaConvergence;
    /// The convergence criterion for the lambda iterations
    double _lambdaThreshold;
    /// The convergence criterion for the scf iterations
    double _scfThreshold;
    /// The convergence that must be achieved before DIIS extrapolation starts
    double _diisStartThresh;
    /// The SCF component of the energy
    double _scfEnergy;
    /// The previous total energy
    double _oldTotalEnergy;
    /// The updated total energy
    double _newTotalEnergy;
    /// The Tikhonow regularizer used to remove singularities (c.f. Taube and Bartlett, JCP, 2009)
    double _regularizer;
    /// The alpha occupied eigenvectors, per irrep
    Matrix _aOccC;
    /// The beta occupied eigenvectors, per irrep
    Matrix _bOccC;
    /// The alpha virtual eigenvectors, per irrep
    Matrix _aVirC;
    /// The beta virtual eigenvectors, per irrep
    Matrix _bVirC;
    /// The Tau matrix in the AO basis, stored by irrep, to perturb the alpha Fock matrix
    double ***_aTau;
    /// The Tau matrix in the AO basis, stored by irrep, to perturb the beta Fock matrix
    double ***_bTau;
    /// The overlap matrix in the AO basis
    Matrix _aoS;
    /// The one-electron integrals in the SO basis
    Matrix _soH;
    /// The alpha Fock matrix in the SO basis
    Matrix _Fa;
    /// The beta Fock matrix in the SO basis
    Matrix _Fb;
    /// The inverse square root overlap matrix in the SO basis
    Matrix _sHalfInv;
    /// The full alpha MO coefficients
    Matrix _Ca;
    /// The full beta MO coefficients
    Matrix _Cb;
    /// The old full alpha MO coefficients
    Matrix _oldCa;
    /// The old full beta MO coefficients
    Matrix _oldCb;
    /// The alpha kappa matrix in the SO basis
    Matrix _aKappa;
    /// The beta kappa matrix in the SO basis
    Matrix _bKappa;
    /// The alpha external potential in the SO basis
    Matrix _aGTau;
    /// The beta external potential in the SO basis
    Matrix _bGTau;
    /// The alpha SCF error vector
    Matrix _aScfError;
    /// The beta SCF error vector
    Matrix _bScfError;
    /// The full alpha MO energies
    Vector _aEvals;
    /// The full beta MO energies
    Vector _bEvals;
    /// Used to align things in the output
    std::string indent;
};

}} // Namespaces

#endif // Header guard
