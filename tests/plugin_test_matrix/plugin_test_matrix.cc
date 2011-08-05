#include "psi4-dec.h"
#include <liboptions/liboptions.h>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include <libplugin/plugin.h>
#include "libchkpt/chkpt.h"

INIT_PLUGIN

namespace psi{ namespace plugin_test_matrix{

extern "C" int
read_options(std::string name, Options &options){
    if(name == "PLUGIN_TEST_MATRIX") {
        /*- The amount of information printed
            to the output file -*/
        options.add_int("PRINT", 1);
        /*- How to cache quantities within the DPD library -*/
        options.add_int("CACHELEV", 2);
        /*- The amount of memory available (in Mb) -*/
        options.add_int("MEMORY", 2000);

        /*- The size of the matrices -*/
        options.add_int("MATRIX_SIZE", 1000);

        /*- The size of the tiles -*/
        options.add_int("TILE_SIZE", 10);

    }
}


extern "C" PsiReturnType
plugin_test_matrix(Options &options)
{
    int print = options.get_int("PRINT");
    // This will print out all of the user-provided options for this module
    if (print > 2) options.print();
   
    fprintf(outfile, "\n\n***** Lets test out some different matrix multiplication schemes *****\n\n");



    int mat_size = options.get_int("MATRIX_SIZE");
    int T = options.get_int("TILE_SIZE");

    double mb = ( mat_size * mat_size * sizeof(double) * 3) / 1048576;
    double gb = ( mat_size * mat_size * sizeof(double) * 3) / 1073741824;

    double tile_kb = ( 2*T*T*sizeof(double) ) * 0.0009765625;
    std::cout << "memory = " << mb << " MB" << std::endl;
    std::cout << "memory = " << gb << " GB" << std::endl;
    std::cout << "tile_size = " << tile_kb << " KB" << std::endl;



    double **A, **B, **C;
    A = block_matrix(mat_size,mat_size);
    B = block_matrix(mat_size,mat_size);
    C = block_matrix(mat_size,mat_size);

    zero_mat(A,mat_size,mat_size);
    zero_mat(B,mat_size,mat_size);
    zero_mat(C,mat_size,mat_size);

    for (int i=0; i < mat_size; i++) {
        A[i][i] = 1.0;
        B[i][i] = 1.0;
    }

/*    std::cout << "Started Matrix I " << std::endl;
    timer_on("Matrix I");

    for (int i=0; i < mat_size; i++) {
        for (int j=0; j < mat_size; j++) {
            double sum = 0.0;
            for (int k=0; k < mat_size; k++) {
                sum += A[i][k] * B[k][j];
            }
            C[i][j] = sum;
        }
    }

    timer_off("Matrix I");

    std::cout << "Finished Matrix I " << std::endl;

    fprintf(outfile, "Mat A after I");
    print_mat(A, mat_size, mat_size, outfile);
    fprintf(outfile, "Mat B after I");
    print_mat(B, mat_size, mat_size, outfile);
    fprintf(outfile, "Mat C after I");
    print_mat(C, mat_size, mat_size, outfile);
*/


    zero_mat(C,mat_size,mat_size);

    std::cout << "Started Matrix II" << std::endl;
    timer_on("Matrix II");

    for (int i=0; i < mat_size; i++) {
        for (int k=0; k < mat_size; k++) {
            double aik = A[i][k];
            for (int j=0; j < mat_size; j++) {
                C[i][j] += aik * B[k][j];
            }
        }
    }

    timer_off("Matrix II");
    std::cout << "Finished Matrix II " << std::endl;

    zero_mat(C,mat_size,mat_size);


    std::cout << "Started Matrix III" << std::endl;
    timer_on("Matrix III");


    if (mat_size%T != 0) throw PSIEXCEPTION("the block_size must be a multiple of the matix_size");

    for (int ilo=0, ihi=T; ilo < mat_size; ilo+=T, ihi+=T) {
        for (int jlo=0, jhi=T; jlo < mat_size; jlo+=T, jhi+=T) {

            for (int k=0; k < mat_size; k++) {
                for (int i=ilo; i < ihi; i++) {
                    double aik = A[i][k];
                    for (int j=jlo; j < jhi; j++) {
                        C[i][j] += aik * B[k][j];
                    }
                }
            }

        }
    }

    timer_off("Matrix III");
    std::cout << "Finished Matrix III " << std::endl;


    zero_mat(C,mat_size,mat_size);

    std::cout << "Started Matrix IV" << std::endl;
    timer_on("Matrix IV");


    if (mat_size%T != 0) throw PSIEXCEPTION("the block_size must be a multiple of the matix_size");

    double *R = new double[T];

    for (int ilo=0, ihi=T; ilo < mat_size; ilo+=T, ihi+=T) {        
        for (int jlo=0, jhi=T; jlo < mat_size; jlo+=T, jhi+=T) {

            for (int k=0; k < mat_size; k++) {

                for (int i=ilo; i < ihi; i++) {
                    double aik = A[i][k];
                    ::memcpy(&(R[0]), &(C[i][jlo]), T);
                    for (int j=jlo, r=0; j < jhi; j++, r++) {
                        R[r] += aik * B[k][j];
                    }
                    ::memcpy(&(C[i][jlo]), &(R[0]), T);
                }
            }

        }
    }

    timer_off("Matrix IV");
    std::cout << "Finished Matrix IV " << std::endl;


    zero_mat(C,mat_size,mat_size);

    timer_on("blas");
    C_DGEMM('n','n',mat_size, mat_size, mat_size, 1.0, &(A[0][0]), mat_size, &(B[0][0]), mat_size, 0.0, &(C[0][0]), mat_size);
    timer_off("blas");

    free_block(A);
    free_block(B);
    free_block(C);

    return Success;   
}

}} // End Namespaces
