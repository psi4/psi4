#include "dcft.h"
#include "defines.h"
#include <vector>
#include <liboptions/liboptions.h>
#include <libpsio/psio.hpp>
#include <libtrans/integraltransform.h>
#include <libtrans/mospace.h>
#include <libdpd/dpd.h>
#include <libdiis/diismanager.h>
#include <libiwl/iwl.hpp>

namespace psi{ namespace dcft{

DCFTSolver::DCFTSolver(Options &options):
    _psio(_default_psio_lib_),
    _chkpt(_default_chkpt_lib_)
{
    _options          = options;
    _scfMaxIter       = options.get_int("SCF_MAXITER");
    _lambdaMaxIter    = options.get_int("LAMBDA_MAXITER");
    _maxNumIterations = options.get_int("MAXITER");
    _print            = options.get_int("PRINT");
    _maxDiis          = options.get_int("MAX_DIIS");
    _minDiisVecs      = options.get_int("DIIS_NUM_VECS");
    _regularizer      = options.get_double("REGULARIZER");
    _diisStartThresh  = pow(10.0, -options.get_int("DIIS_START"));
    _scfThreshold     = pow(10.0, -options.get_int("CONVERGENCE"));
    _lambdaThreshold  = pow(10.0, -options.get_int("CONVERGENCE"));
    _intTolerance     = pow(10.0, -options.get_int("INT_THRESH"));
    _psio->open(PSIF_DCFT_DPD, PSIO_OPEN_OLD);
    read_checkpoint();
}

/**
 * Computes A = A + alpha * B, writing the result back to A
 */
void
DCFTSolver::dpd_buf4_add(dpdbuf4 *A, dpdbuf4 *B, double alpha)
{
    for(int h = 0; h < _nIrreps; ++h){
        dpd_buf4_mat_irrep_init(A, h);
        dpd_buf4_mat_irrep_init(B, h);
        dpd_buf4_mat_irrep_rd(A, h);
        dpd_buf4_mat_irrep_rd(B, h);
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
    free_moinfo();
    delete _ints;
    _psio->close(PSIF_DCFT_DPD, 1);
}

}} // Namespaces
