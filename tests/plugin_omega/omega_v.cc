#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <utility>
#include <string>
#include <cstring>

#include <psifiles.h>
#include <physconst.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>

#include <libmints/mints.h>
#include <libfunctional/superfunctional.h>
#include <libscf_solver/ks.h>
#include <libscf_solver/integralfunctors.h>
#include <libscf_solver/omegafunctors.h>

#include "omega.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace psi;
using namespace psi::functional;
using namespace boost;

namespace psi{ namespace scf {

OmegaV::OmegaV(boost::shared_ptr<PSIO> psio, boost::shared_ptr<BasisSet> primary,
        boost::shared_ptr<SuperFunctional> functional, boost::shared_ptr<Integrator> integrator,
        boost::shared_ptr<Properties> props) :
    psio_(psio), primary_(primary), functional_(functional), integrator_(integrator),
    properties_(props)
{
    common_init();
}
OmegaV::~OmegaV()
{
}
void OmegaV::common_init()
{
    int nso = primary_->nbf();
    
    Va_ = boost::shared_ptr<Matrix>(new Matrix("Va", nso, nso));
    Vb_ = boost::shared_ptr<Matrix>(new Matrix("Vb", nso, nso));
}
void OmegaV::set_omega(double omega)
{
    functional_->setOmega(omega);
}
void OmegaV::form_V(boost::shared_ptr<Matrix> Da, boost::shared_ptr<Matrix> Ca, int na,
                    boost::shared_ptr<Matrix> Db, boost::shared_ptr<Matrix> Cb, int nb)
{
    // Zero the V_ matrix
    Va_->zero();
    Vb_->zero();

    int nso = primary_->nbf();

    // Some temporary buffers
    boost::shared_ptr<Matrix> Vt1a = boost::shared_ptr<Matrix>(new Matrix("Va", nso, nso));
    boost::shared_ptr<Matrix> Vt2a = boost::shared_ptr<Matrix>(new Matrix("Va", nso, nso));
    double** V_temp1a = Vt1a->pointer();
    double** V_temp2a = Vt2a->pointer();
    boost::shared_ptr<Matrix> Vt1b = boost::shared_ptr<Matrix>(new Matrix("Va", nso, nso));
    boost::shared_ptr<Matrix> Vt2b = boost::shared_ptr<Matrix>(new Matrix("Va", nso, nso));
    double** V_temp1b = Vt1b->pointer();
    double** V_temp2b = Vt2b->pointer();

    // Zero the registers
    double functional_E = 0.0;
    double densityCheckA = 0.0;
    double dipoleCheckXA = 0.0;
    double dipoleCheckYA = 0.0;
    double dipoleCheckZA = 0.0;
    double densityCheckB = 0.0;
    double dipoleCheckXB = 0.0;
    double dipoleCheckYB = 0.0;
    double dipoleCheckZB = 0.0;

    // Grab the integrator's grid block
    int nblocks = integrator_->getNBlocks();
    double *x;
    double *y;
    double *z;
    double *w;

    // Grab the properties references
    double *rho_a       = properties_->getRhoA();
    double *rho_b       = properties_->getRhoB();
    double *gamma_aa    = properties_->getGammaAA();
    double *gamma_ab    = properties_->getGammaAB();
    double *gamma_bb    = properties_->getGammaBB();
    double *tau_a       = properties_->getTauA();
    double *tau_b       = properties_->getTauB();

    double *rho_ax       = properties_->getRhoAX();
    double *rho_ay       = properties_->getRhoAY();
    double *rho_az       = properties_->getRhoAZ();
    double *rho_bx       = properties_->getRhoBX();
    double *rho_by       = properties_->getRhoBY();
    double *rho_bz       = properties_->getRhoBZ();
    // (These will not be used unless needed)

    // Grab sparseness information and scratch array
    int* rel2abs = properties_->rel2absFunctions();
    double **scratch    = properties_->getScratch(); //Already allocated

    // Grab the basis points references
    double **bas = properties_->getPoints();
    double **bas_x = properties_->getGradientsX();
    double **bas_y = properties_->getGradientsY();
    double **bas_z = properties_->getGradientsZ();
    // (These will not be used unless needed)

    // Grab the functional references
    double *zk          = functional_->getFunctionalValue();
    double *v_rho_a     = functional_->getV_RhoA();
    double *v_rho_b     = functional_->getV_RhoB();
    double *v_gamma_aa  = functional_->getV_GammaAA();
    double *v_gamma_ab  = functional_->getV_GammaAB();
    double *v_gamma_bb  = functional_->getV_GammaBB();
    double *v_tau_a     = functional_->getV_TauA();
    double *v_tau_b     = functional_->getV_TauB();
    // (These will not be used unless needed)

    // GGA? Meta?
    bool GGA = functional_->isGGA();
    bool Meta = functional_->isMeta();

    // Some indexing
    int ntrue, nsigf, nbf, index, offset, h, mu, nu, m ,n;
    double contribution;

    nbf = nso;

    // Traverse grid blocks
    for (int N = 0; N < nblocks; N++) {

        // Compute integration points
        boost::shared_ptr<GridBlock> block = integrator_->getBlock(N);
        x = block->getX();
        y = block->getY();
        z = block->getZ();
        w = block->getWeights();
        ntrue = block->getTruePoints();

        // Compute properties and basis points
        properties_->computeUKSProperties(block, Da, Db, Ca, Cb, &na, &nb);
        nsigf = properties_->nSignificantFunctions();

        // Compute functional values and partials
        functional_->computeUKSFunctional(properties_);

        for (index = 0; index < ntrue; index++) {
            zk[index] *= w[index];
            v_rho_a[index] *= w[index];
            v_rho_b[index] *= w[index];
        }
        if (GGA) {
            for (index = 0; index < ntrue; index++) {
                v_gamma_aa[index] *= w[index];
                v_gamma_ab[index] *= w[index];
                v_gamma_bb[index] *= w[index];
            }
        }
        if (Meta) {
            for (index = 0; index < ntrue; index++) {
                v_tau_a[index] *= w[index];
                v_tau_b[index] *= w[index];
            }
        }

        // Compute functional energy
        // And monopole/dipole checks
        for (index = 0; index < ntrue; index++) {
            //printf(" Block %d, Point %d, rho %14.10E\n", N, index, rho_a[index]);

            functional_E += zk[index];
            densityCheckA += w[index]*rho_a[index];
            dipoleCheckXA += w[index]*x[index]*rho_a[index];
            dipoleCheckYA += w[index]*y[index]*rho_a[index];
            dipoleCheckZA += w[index]*z[index]*rho_a[index];
            densityCheckB += w[index]*rho_b[index];
            dipoleCheckXB += w[index]*x[index]*rho_b[index];
            dipoleCheckYB += w[index]*y[index]*rho_b[index];
            dipoleCheckZB += w[index]*z[index]*rho_b[index];
        }

        // LSDA contribution to potential:
        // Alpha
        // V_mn^a += v_rho_a * phi_m * phi_n
        for (index = 0; index < ntrue; index++) {
            memcpy((void*) &scratch[index][0], (void*) &bas[index][0], nsigf*sizeof(double));
            C_DSCAL(nsigf, v_rho_a[index], &scratch[index][0], 1);
        }
        C_DGEMM('T','N', nsigf, nsigf, ntrue, 1.0, &scratch[0][0], nbf, \
            &bas[0][0], nbf, 0.0, &V_temp1a[0][0], nbf);
        // Beta
        // V_mn^b += v_rho_b * phi_m * phi_n
        for (index = 0; index < ntrue; index++) {
            memcpy((void*) &scratch[index][0], (void*) &bas[index][0], nsigf*sizeof(double));
            C_DSCAL(nsigf, v_rho_b[index], &scratch[index][0], 1);
        }
        C_DGEMM('T','N', nsigf, nsigf, ntrue, 1.0, &scratch[0][0], nbf, \
            &bas[0][0], nbf, 0.0, &V_temp1b[0][0], nbf);

        // GGA contribution to potential:
        //V_mn^a += [2 * v_gamma_aa * \nabla rho_a + v_gamma_ab * \nabla rho_b] \dot \nabla (phi_m * phi_n)
        if (GGA) {
            // Alpha
            // Gradient x
            for (index = 0; index < ntrue; index++) {
                memcpy((void*) &scratch[index][0], (void*) &bas[index][0], nsigf*sizeof(double));
                C_DSCAL(nsigf, 2.0*v_gamma_aa[index]*rho_ax[index] + v_gamma_ab[index]*rho_bx[index], &scratch[index][0], 1);
            }
            C_DGEMM('T','N', nsigf, nsigf, ntrue, 1.0, &scratch[0][0], nbf, \
                &bas_x[0][0], nbf, 0.0, &V_temp2a[0][0], nbf);

            // Gradient y
            for (index = 0; index < ntrue; index++) {
                memcpy((void*) &scratch[index][0], (void*) &bas[index][0], nsigf*sizeof(double));
                C_DSCAL(nsigf, 2.0*v_gamma_aa[index]*rho_ay[index] + v_gamma_ab[index]*rho_by[index], &scratch[index][0], 1);
            }
            C_DGEMM('T','N', nsigf, nsigf, ntrue, 1.0, &scratch[0][0], nbf, \
                &bas_y[0][0], nbf, 1.0, &V_temp2a[0][0], nbf);

            // Gradient z
            for (index = 0; index < ntrue; index++) {
                memcpy((void*) &scratch[index][0], (void*) &bas[index][0], nsigf*sizeof(double));
                C_DSCAL(nsigf, 2.0*v_gamma_aa[index]*rho_az[index] + v_gamma_ab[index]*rho_bz[index], &scratch[index][0], 1);
            }
            C_DGEMM('T','N', nsigf, nsigf, ntrue, 1.0, &scratch[0][0], nbf, \
                &bas_z[0][0], nbf, 1.0, &V_temp2a[0][0], nbf);

            // I only did the left term, so V_xc_GGA = V_temp2_ + v_temp2_':
            for (int m = 0; m < nsigf; m++)
                for (int n = 0; n < nsigf; n++)
                    V_temp1a[m][n] += V_temp2a[m][n] + V_temp2a[n][m];

            // Beta
            // Gradient x
            for (index = 0; index < ntrue; index++) {
                memcpy((void*) &scratch[index][0], (void*) &bas[index][0], nsigf*sizeof(double));
                C_DSCAL(nsigf, 2.0*v_gamma_bb[index]*rho_bx[index] + v_gamma_ab[index]*rho_ax[index], &scratch[index][0], 1);
            }
            C_DGEMM('T','N', nsigf, nsigf, ntrue, 1.0, &scratch[0][0], nbf, \
                &bas_x[0][0], nbf, 0.0, &V_temp2b[0][0], nbf);

            // Gradient y
            for (index = 0; index < ntrue; index++) {
                memcpy((void*) &scratch[index][0], (void*) &bas[index][0], nsigf*sizeof(double));
                C_DSCAL(nsigf, 2.0*v_gamma_bb[index]*rho_by[index] + v_gamma_ab[index]*rho_ay[index], &scratch[index][0], 1);
            }
            C_DGEMM('T','N', nsigf, nsigf, ntrue, 1.0, &scratch[0][0], nbf, \
                &bas_y[0][0], nbf, 1.0, &V_temp2b[0][0], nbf);

            // Gradient z
            for (index = 0; index < ntrue; index++) {
                memcpy((void*) &scratch[index][0], (void*) &bas[index][0], nsigf*sizeof(double));
                C_DSCAL(nsigf, 2.0*v_gamma_bb[index]*rho_bz[index] + v_gamma_ab[index]*rho_az[index], &scratch[index][0], 1);
            }
            C_DGEMM('T','N', nsigf, nsigf, ntrue, 1.0, &scratch[0][0], nbf, \
                &bas_z[0][0], nbf, 1.0, &V_temp2b[0][0], nbf);

            // I only did the left term, so V_xc_GGA = V_temp2_ + v_temp2_':
            for (int m = 0; m < nsigf; m++)
                for (int n = 0; n < nsigf; n++)
                    V_temp1b[m][n] += V_temp2b[m][n] + V_temp2b[n][m];

        }

        // TODO determine if the 2.0 is still here for meta-UKS
        // Meta contribution to potential
        //V_mn += 2 * v_tau_a * \nabla \phi_m \dot \nabla \phi_n
        if (Meta) {
            // Alpha
            // phi_x
            for (index = 0; index < ntrue; index++) {
                memcpy((void*) &scratch[index][0], (void*) &bas_x[index][0], nsigf*sizeof(double));
                C_DSCAL(nsigf, 2.0*v_tau_a[index], &scratch[index][0], 1);
            }
            C_DGEMM('T','N', nsigf, nsigf, ntrue, 1.0, &scratch[0][0], nbf, \
                &bas_x[0][0], nbf, 1.0, &V_temp1a[0][0], nbf);

            // phi_y
            for (index = 0; index < ntrue; index++) {
                memcpy((void*) &scratch[index][0], (void*) &bas_y[index][0], nsigf*sizeof(double));
                C_DSCAL(nsigf, 2.0*v_tau_a[index], &scratch[index][0], 1);
            }
            C_DGEMM('T','N', nsigf, nsigf, ntrue, 1.0, &scratch[0][0], nbf, \
                &bas_y[0][0], nbf, 1.0, &V_temp1a[0][0], nbf);

            // phi_z
            for (index = 0; index < ntrue; index++) {
                memcpy((void*) &scratch[index][0], (void*) &bas_z[index][0], nsigf*sizeof(double));
                C_DSCAL(nsigf, 2.0*v_tau_a[index], &scratch[index][0], 1);
            }
            C_DGEMM('T','N', nsigf, nsigf, ntrue, 1.0, &scratch[0][0], nbf, \
                &bas_z[0][0], nbf, 1.0, &V_temp1a[0][0], nbf);

            // Beta
            // phi_x
            for (index = 0; index < ntrue; index++) {
                memcpy((void*) &scratch[index][0], (void*) &bas_x[index][0], nsigf*sizeof(double));
                C_DSCAL(nsigf, 2.0*v_tau_b[index], &scratch[index][0], 1);
            }
            C_DGEMM('T','N', nsigf, nsigf, ntrue, 1.0, &scratch[0][0], nbf, \
                &bas_x[0][0], nbf, 1.0, &V_temp1b[0][0], nbf);

            // phi_y
            for (index = 0; index < ntrue; index++) {
                memcpy((void*) &scratch[index][0], (void*) &bas_y[index][0], nsigf*sizeof(double));
                C_DSCAL(nsigf, 2.0*v_tau_b[index], &scratch[index][0], 1);
            }
            C_DGEMM('T','N', nsigf, nsigf, ntrue, 1.0, &scratch[0][0], nbf, \
                &bas_y[0][0], nbf, 1.0, &V_temp1b[0][0], nbf);

            // phi_z
            for (index = 0; index < ntrue; index++) {
                memcpy((void*) &scratch[index][0], (void*) &bas_z[index][0], nsigf*sizeof(double));
                C_DSCAL(nsigf, 2.0*v_tau_b[index], &scratch[index][0], 1);
            }
            C_DGEMM('T','N', nsigf, nsigf, ntrue, 1.0, &scratch[0][0], nbf, \
                &bas_z[0][0], nbf, 1.0, &V_temp1b[0][0], nbf);
        }

        // Add into the global Vxc matrix
        for (int m = 0; m < nsigf; m++)
            for (int n = 0; n <= m; n++) {
                Va_->add(0, rel2abs[m], rel2abs[n], V_temp1a[m][n]);
                if (m!=n)
                    Va_->add(0, rel2abs[n], rel2abs[m], V_temp1a[n][m]);
            }
        for (int m = 0; m < nsigf; m++)
            for (int n = 0; n <= m; n++) {
                Vb_->add(0, rel2abs[m], rel2abs[n], V_temp1b[m][n]);
                if (m!=n)
                    Vb_->add(0, rel2abs[n], rel2abs[m], V_temp1b[n][m]);
            }

    } // End traverse over grid blocks

    quad_values_["E_xc"] = functional_E;
    quad_values_["<rho_a>"] = densityCheckA;
    quad_values_["<rho_a*x>"] = dipoleCheckXA;
    quad_values_["<rho_a*y>"] = dipoleCheckYA;
    quad_values_["<rho_a*z>"] = dipoleCheckZA;
    quad_values_["<rho_b>"] = densityCheckB;
    quad_values_["<rho_b*x>"] = dipoleCheckXB;
    quad_values_["<rho_b*y>"] = dipoleCheckYB;
    quad_values_["<rho_b*z>"] = dipoleCheckZB;
}
double OmegaV::Exc()
{
    return quad_values_["E_xc"];
}
double OmegaV::variable(const std::string& key)
{
    return quad_values_[key];
}

}} // End Namespaces
