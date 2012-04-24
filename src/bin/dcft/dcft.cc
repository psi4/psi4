#include "dcft.h"
#include "defines.h"
#include <vector>
#include <cmath>
#include <liboptions/liboptions.h>
#include <libpsio/psio.hpp>
#include <libtrans/integraltransform.h>
#include <libtrans/mospace.h>
#include <libdpd/dpd.h>
#include <libdiis/diismanager.h>
#include <libiwl/iwl.hpp>

using namespace boost;

namespace psi{ namespace dcft{

DCFTSolver::DCFTSolver(boost::shared_ptr<Wavefunction> reference_wavefunction, Options &options):
        Wavefunction(options, _default_psio_lib_)
{
    reference_wavefunction_ = reference_wavefunction;
    scfmaxiter_       = options.get_int("SCF_MAXITER");
    lambdamaxiter_    = options.get_int("LAMBDA_MAXITER");
    maxiter_          = options.get_int("MAXITER");
    print_            = options.get_int("PRINT");
    maxdiis_          = options.get_int("DIIS_MAX_VECS");
    mindiisvecs_      = options.get_int("DIIS_MIN_VECS");
    regularizer_      = options.get_double("TIKHONOW_OMEGA");
    diis_start_thresh_ = options.get_double("DIIS_START_CONVERGENCE");
    scf_threshold_     = options.get_double("R_CONVERGENCE");
    lambda_threshold_  = options.get_double("R_CONVERGENCE");
    int_tolerance_     = options.get_double("INTS_TOLERANCE");
    lock_occupation_   = options.get_bool("LOCK_OCC");
    psio_->open(PSIF_DCFT_DPD, PSIO_OPEN_OLD);

    if(options.get_str("REFERENCE") != "UHF") throw PSIEXCEPTION("You must have reference=UHF in the input file");

    // Sets up the memory, and orbital info
    init();

    energy_tau_squared_ = 0.0;
}

/**
 * Computes A = A + alpha * B, writing the result back to A
 */
void
DCFTSolver::dpd_buf4_add(dpdbuf4 *A, dpdbuf4 *B, double alpha)
{
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(A, h);
        dpd_buf4_mat_irrep_init(B, h);
        dpd_buf4_mat_irrep_rd(A, h);
        dpd_buf4_mat_irrep_rd(B, h);

        #pragma omp parallel for
        for(int row = 0; row < A->params->rowtot[h]; ++row){
            for(int col = 0; col < A->params->coltot[h]; ++col){
                A->matrix[h][row][col] += alpha * B->matrix[h][row][col];
            }
        }
        dpd_buf4_mat_irrep_wrt(A, h);
        dpd_buf4_mat_irrep_close(A, h);
        dpd_buf4_mat_irrep_close(B, h);
    }
}

DCFTSolver::~DCFTSolver()
{
}

}} // Namespaces
