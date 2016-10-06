#include "psi4-dec.h"
#include <libmints/mints.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>

namespace psi{ namespace dfmp2{

void formInvSqrtJ(double **&J_mhalf, shared_ptr<BasisSet> basis,
                       shared_ptr<BasisSet> ribasis, shared_ptr<BasisSet> zero)
{
    // Create integral factories for the RI basis
    shared_ptr<IntegralFactory>
            rifactory_J(new IntegralFactory(ribasis, zero, ribasis, zero));
    shared_ptr<TwoBodyAOInt> Jint(rifactory_J->eri());

    double **J = block_matrix(ribasis->nbf(), ribasis->nbf());
    J_mhalf = block_matrix(ribasis->nbf(), ribasis->nbf());
    const double *Jbuffer = Jint->buffer();

#ifdef TIME_DF_MP2
    timer_on("Form J");
#endif

    int index = 0;
    for (int MU=0; MU < ribasis->nshell(); ++MU) {
        int nummu = ribasis->shell(MU)->nfunction();
        for (int NU=0; NU < ribasis->nshell(); ++NU) {
            int numnu = ribasis->shell(NU)->nfunction();
            Jint->compute_shell(MU, 0, NU, 0);
            index = 0;
            for (int mu=0; mu < nummu; ++mu) {
                int omu = ribasis->shell(MU)->function_index() + mu;
                for (int nu=0; nu < numnu; ++nu, ++index) {
                    int onu = ribasis->shell(NU)->function_index() + nu;
                    J[omu][onu] = Jbuffer[index];
                }
            }
        }
    }

    // First, diagonalize J
    // the C_DSYEV call replaces the original matrix J with its eigenvectors
    int lwork = ribasis->nbf() * 3;
    double* eigval = init_array(ribasis->nbf());
    double* work = init_array(lwork);
    int status = C_DSYEV('v', 'u', ribasis->nbf(), J[0],
                        ribasis->nbf(), eigval, work, lwork);
    if(status){
        throw PsiException("Diagonalization of J failed", __FILE__, __LINE__);
    }
    free(work);

    // Now J contains the eigenvectors of the original J
    // Copy J to J_copy
    double **J_copy = block_matrix(ribasis->nbf(), ribasis->nbf());
    C_DCOPY(ribasis->nbf()*ribasis->nbf(), J[0], 1, J_copy[0], 1);

    // Now form J^{-1/2} = U(T)*j^{-1/2}*U,
    // where j^{-1/2} is the diagonal matrix of the inverse square roots
    // of the eigenvalues, and U is the matrix of eigenvectors of J
    for(int i=0; i<ribasis->nbf(); ++i){
        eigval[i] = (eigval[i] < 1.0E-10) ? 0.0 : 1.0 / sqrt(eigval[i]);
        // scale one set of eigenvectors by the diagonal elements j^{-1/2}
        C_DSCAL(ribasis->nbf(), eigval[i], J[i], 1);
    }
    free(eigval);

    // J_mhalf = J_copy(T) * J
    C_DGEMM('t','n',ribasis->nbf(),ribasis->nbf(),ribasis->nbf(),1.0,
    J_copy[0],ribasis->nbf(),J[0],ribasis->nbf(),0.0,J_mhalf[0],ribasis->nbf());
    free_block(J);
    free_block(J_copy);

#ifdef TIME_DF_MP2
    timer_off("Form J");
#endif
}

}} // Namespaces
