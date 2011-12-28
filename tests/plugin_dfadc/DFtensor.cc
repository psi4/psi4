#include "psi4-dec.h"
#include <libmints/mints.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "dfadc.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace psi{ namespace plugin_dfadc {

//
// At first, I considered to use lib3index instead of these. But I didn't like to write
// any intermediate tensors into disk, so that I coded the functions below to generate
// fitting metrix and DF tensors on-the-fly fashion. Sorry, Dr. Parrish and Dr. Hohenstein.
// All two were developed by referebce to Dr. Justin's code.
//
// I wonder formInvSqrtJ works faster if all the OpenMP macros are removed in many cases.
//
// Masaaki (2011/12/25)
//
    
void 
DFADC::formInvSqrtJ(double **&J_mhalf)
{
    
    int nthread = 1;
    #ifdef _OPENMP
    if(!omp_in_parallel()){
        nthread = omp_get_max_threads();
    }
    #endif

    boost::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
    boost::shared_ptr<IntegralFactory> rifactory_J(new IntegralFactory(ribasis_, zero, ribasis_, zero));
    const double **Jbuffer = new const double *[nthread];
    boost::shared_ptr<TwoBodyAOInt> *Jint = new boost::shared_ptr<TwoBodyAOInt>[nthread];
    for(int I = 0;I < nthread;I++){
        Jint[I] = boost::shared_ptr<TwoBodyAOInt>(rifactory_J->eri());
        Jbuffer[I] = Jint[I]->buffer();
    }
    
    double **J = block_matrix(ribasis_->nbf(), ribasis_->nbf());
    J_mhalf    = block_matrix(ribasis_->nbf(), ribasis_->nbf());

    #pragma omp parallel for schedule (dynamic) num_threads(nthread)
    for(int MU = 0;MU < ribasis_->nshell();MU++){
        int numMU = ribasis_->shell(MU)->nfunction();
        
        int thread = 0;
        #ifdef _OPENMP
        thread = omp_get_thread_num();
        #endif
        
        for(int NU = 0;NU <=MU;NU++){
            int numNU = ribasis_->shell(NU)->nfunction();
            Jint[thread]->compute_shell(MU, 0, NU, 0);
            int index = 0;
            for(int mu = 0;mu < numMU;mu++){
                int omu = ribasis_->shell(MU)->function_index() + mu;
                for(int nu = 0;nu < numNU;nu++, index++){
                    int onu = ribasis_->shell(NU)->function_index() + nu;
                    J[omu][onu] = Jbuffer[thread][index];
                    J[onu][omu] = Jbuffer[thread][index];
                }
            }
        }
    }

    int lwork = ribasis_->nbf() * 3;
    double *eigval = init_array(ribasis_->nbf());
    double *work   = init_array(lwork);
    int status = C_DSYEV('V', 'U', ribasis_->nbf(), J[0], ribasis_->nbf(), eigval, work, lwork);
    
    if(status){
        throw PsiException("Diagonalization of J failed!", __FILE__, __LINE__);
    }
    free(work);
    
    double **J_copy = block_matrix(ribasis_->nbf(), ribasis_->nbf());
    C_DCOPY(ribasis_->nbf()*ribasis_->nbf(), J[0], 1, J_copy[0], 1);
    for(int i = 0;i < ribasis_->nbf();i++){
        eigval[i] = (eigval[i] < 1.0e-10) ? 0.0 : 1.0 / sqrt(eigval[i]);
        C_DSCAL(ribasis_->nbf(), eigval[i], J[i], 1);
    }
    free(eigval);
    C_DGEMM('T', 'N', ribasis_->nbf(), ribasis_->nbf(), ribasis_->nbf(), 1.0, J_copy[0], ribasis_->nbf(), J[0], ribasis_->nbf(), 0.0, J_mhalf[0], ribasis_->nbf());
    
    delete [] Jbuffer;
    delete [] Jint;
    free_block(J);
    free_block(J_copy);
}

void
DFADC::formDFtensor(SharedMatrix Cl, SharedMatrix Cr, double **J_mhalf, enum Order order, double **&Bout)
{    
    int maxPshell = 0;
    for(int Pshell = 0;Pshell < ribasis_->nshell();Pshell++){
        int numPshell = ribasis_->shell(Pshell)->nfunction();
        maxPshell = numPshell > maxPshell ? numPshell : maxPshell;
    }
    double ***temp = new double**[maxPshell];
    for(int P = 0;P < maxPshell;P++) temp[P] = block_matrix(ribasis_->nbf(), ribasis_->nbf());
    
    boost::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
    boost::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(ribasis_, zero, basisset_, basisset_));
    boost::shared_ptr<TwoBodyAOInt> eri(rifactory->eri());
    const double *buffer = eri->buffer();
    
    int Lcol = Cl->ncol();
    int Rrow = Cr->ncol(); 
    double **mo_P_pq = block_matrix(ribasis_->nbf(), Lcol*Rrow);
    
    double **Cleft  = Cl->pointer();
    double **Cright = Cr->pointer();
    double **half   = block_matrix(Lcol, ribasis_->nbf());
    
    for(int Pshell = 0;Pshell < ribasis_->nshell();Pshell++){
        int numPshell = ribasis_->shell(Pshell)->nfunction();
        for(int P = 0;P < numPshell;P++){
            zero_mat(temp[P], ribasis_->nbf(), ribasis_->nbf());
        }
        for(int MU = 0;MU < basisset_->nshell();MU++){
            int numMU = basisset_->shell(MU)->nfunction();
            for(int NU = 0;NU <= MU;NU++){
                int numNU = basisset_->shell(NU)->nfunction();
                eri->compute_shell(Pshell, 0, MU, NU);
                for(int P = 0, index = 0;P < numPshell;P++){
                    for(int mu = 0;mu < numMU;mu++){
                        int omu = basisset_->shell(MU)->function_index() + mu;
                        for(int nu = 0;nu < numNU;nu++, index++){
                            int onu = basisset_->shell(NU)->function_index() + nu;
                            temp[P][omu][onu] = buffer[index];
                            temp[P][onu][omu] = buffer[index];
                        }
                    }
                }
            }
        }

        for(int P = 0, index = 0;P < numPshell;P++){
            int oP = ribasis_->shell(Pshell)->function_index() + P;
            C_DGEMM('T', 'N', Lcol, ribasis_->nbf(), ribasis_->nbf(), 1.0, Cleft[0], Lcol, temp[P][0], ribasis_->nbf(), 0.0, half[0],     ribasis_->nbf());
            C_DGEMM('N', 'N', Lcol, Rrow, ribasis_->nbf(), 1.0, half[0],  ribasis_->nbf(), Cright[0],  Rrow, 0.0, mo_P_pq[oP], Rrow);
        }
    }
    free_block(half);
    
    for(int P = 0;P < maxPshell;P++) free_block(temp[P]);

    if(order == Irow){
        // Tiling B_P^{pq} 
        Bout = block_matrix(ribasis_->nbf(), Lcol*Rrow);
        // B^P_{ia} <-- \sum_{Q} J^{-1/2}_{PQ} (Q|ia)
        C_DGEMM('T', 'N', ribasis_->nbf(), Lcol*Rrow, ribasis_->nbf(), 1.0, J_mhalf[0], ribasis_->nbf(), mo_P_pq[0], Lcol*Rrow, 0.0, Bout[0], Lcol*Rrow);
    }
    else if(order == Icol){
        // Tile B_{pq}^P
        Bout = block_matrix(Lcol*Rrow, ribasis_->nbf());
        // B_{ia}^P <-- \sum_{Q} (ia|Q) J^{-1/2}_{QP}
        C_DGEMM('T', 'N', Lcol*Rrow, ribasis_->nbf(), ribasis_->nbf(), 1.0, mo_P_pq[0], Lcol*Rrow, J_mhalf[0], ribasis_->nbf(), 0.0, Bout[0], ribasis_->nbf());
    }
    free_block(mo_P_pq);
    
}
        
}} // End Namespaces

