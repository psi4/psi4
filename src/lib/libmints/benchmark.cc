#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include <libutil/libutil.h>
#include <libpsio/psio.hpp>
#include <libpsio/psio.h>
#include "mints.h"
#include <cmath>
#include <cstdlib>
#include <psi4-dec.h>

#include <map>
#include <string>
#include <vector>

#include <psiconfig.h>
#ifdef HAVE_MKL
#include <mkl.h>
#endif

using namespace psi;

namespace psi {

void benchmark_blas1(int N, double min_time)
{
    fprintf(outfile, "\n");
    fprintf(outfile, "                              ------------------------------- \n");
    fprintf(outfile, "                              ======> BLAS1 BENCHMARKS <===== \n");
    fprintf(outfile, "                              ------------------------------- \n");
    fprintf(outfile, "\n");

    fprintf(outfile, "  Parameters:\n");
    fprintf(outfile, "   -Minimum runtime (per operation, per size): %14.10f [s].\n", min_time);
    fprintf(outfile, "   -Maximum dimension exponent N: %d. Arrays are D x D = 2^N x 2^N doubles in size. The D\n", N);
    fprintf(outfile, "        value is reported below\n");
    fprintf(outfile, "\n");

    fprintf(outfile, "  Notes:\n");
    fprintf(outfile, "   -Access: c = A[i]; (stride 1).\n");
    fprintf(outfile, "   -Assign: A[i] = c; (stride 1).\n");
    fprintf(outfile, "   -Cross:  A[i] = B[i]; (stride 1).\n");
    fprintf(outfile, "   -Strides: (XX) indicates the strides used for the various arrays of BLAS 1 operations. 1\n");
    fprintf(outfile, "        indicates stride 1, N indicates stride D.\n");
    fprintf(outfile, "\n");

    double T;
    unsigned long int rounds;
    double t;
    int dim;
    Timer* qq;

    std::map<std::string, std::vector<double> > timings1;
    std::vector<std::string> ops1;

    ops1.push_back("Malloc");
    ops1.push_back("Memset");
    ops1.push_back("Access");
    ops1.push_back("Assign");
    ops1.push_back("Cross");
    ops1.push_back("DSCAL (1)");
    ops1.push_back("DSCAL (N)");
    ops1.push_back("DCOPY (11)");
    ops1.push_back("DCOPY (N1)");
    ops1.push_back("DCOPY (1N)");
    ops1.push_back("DCOPY (NN)");
    ops1.push_back("DSWAP (11)");
    ops1.push_back("DSWAP (N1)");
    ops1.push_back("DSWAP (1N)");
    ops1.push_back("DSWAP (NN)");
    ops1.push_back("DAXPY (11)");
    ops1.push_back("DAXPY (N1)");
    ops1.push_back("DAXPY (1N)");
    ops1.push_back("DAXPY (NN)");
    ops1.push_back("DDOT (11)");
    ops1.push_back("DDOT (N1)");
    ops1.push_back("DDOT (1N)");
    ops1.push_back("DDOT (NN)");
    ops1.push_back("DROT (11)");
    ops1.push_back("DROT (N1)");
    ops1.push_back("DROT (1N)");
    ops1.push_back("DROT (NN)");
    for (int op = 0; op < ops1.size(); op++)
        timings1[ops1[op]].resize(N);

    // Level 1 routines
    dim = 1;
    for (int k = 0; k < N; k++) {

        dim *= 2;
        unsigned long int full_dim = dim * (unsigned long int) dim;

        double* A = init_array(full_dim);
        double* B = init_array(full_dim);

        double alpha = 1.0;
        double beta = 1.0;
        double sina = sqrt(2.0) / 2.0;
        double cosa = sqrt(2.0) / 2.0;

        // Malloc (and free)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            double* C = (double*) malloc(full_dim * sizeof(double));
            free(C);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings1["Malloc"][k] = t;

        // Memset
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            memset((void*) A, '\0', full_dim*sizeof(double));
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings1["Memset"][k] = t;

        // Access
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            for (unsigned long int Q = 0; Q < full_dim; Q++)
                t = A[Q];
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings1["Access"][k] = t;

        // Assign
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            for (unsigned long int Q = 0; Q < full_dim; Q++)
                A[Q] = t;
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings1["Assign"][k] = t;

        // Cross
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            for (unsigned long int Q = 0; Q < full_dim; Q++)
                A[Q] = B[Q];
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings1["Cross"][k] = t;

        for (unsigned long int Q = 0L; Q < full_dim; Q++) {
            A[Q] = rand() / (double) RAND_MAX;
            B[Q] = rand() / (double) RAND_MAX;
        }

        // DSCAL (1)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            C_DSCAL(full_dim, alpha, A, 1);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings1["DSCAL (1)"][k] = t;

        // DSCAL (N)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            for (int h = 0; h < dim; h++)
                C_DSCAL(dim, alpha, &A[h], dim);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings1["DSCAL (N)"][k] = t;

        // DCOPY (11)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            C_DCOPY(full_dim, A, 1, B, 1);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings1["DCOPY (11)"][k] = t;

        // DCOPY (N1)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            for (int h = 0; h < dim; h++)
                C_DCOPY(dim, &A[h], dim, B, 1);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings1["DCOPY (N1)"][k] = t;

        // DCOPY (1N)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            for (int h = 0; h < dim; h++)
                C_DCOPY(dim, A, 1, &B[h], dim);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings1["DCOPY (1N)"][k] = t;

        // DCOPY (NN)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            for (int h = 0; h < dim; h++)
                C_DCOPY(dim, &A[h], dim, &B[h], dim);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings1["DCOPY (NN)"][k] = t;

        // DSWAP (11)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            C_DSWAP(full_dim, A, 1, B, 1);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings1["DSWAP (11)"][k] = t;

        // DSWAP (N1)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            for (int h = 0; h < dim; h++)
                C_DSWAP(dim, &A[h], dim, B, 1);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings1["DSWAP (N1)"][k] = t;

        // DSWAP (1N)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            for (int h = 0; h < dim; h++)
                C_DSWAP(dim, A, 1, &B[h], dim);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings1["DSWAP (1N)"][k] = t;

        // DSWAP (NN)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            for (int h = 0; h < dim; h++)
                C_DSWAP(dim, &A[h], dim, &B[h], dim);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings1["DSWAP (NN)"][k] = t;

        // DAXPY (11)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            C_DAXPY(full_dim, alpha, A, 1, B, 1);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings1["DAXPY (11)"][k] = t;

        // DAXPY (N1)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            for (int h = 0; h < dim; h++)
                C_DAXPY(dim, alpha, &A[h], dim, B, 1);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings1["DAXPY (N1)"][k] = t;

        // DAXPY (1N)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            for (int h = 0; h < dim; h++)
                C_DAXPY(dim, alpha, A, 1, &B[h], dim);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings1["DAXPY (1N)"][k] = t;

        // DAXPY (NN)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            for (int h = 0; h < dim; h++)
                C_DAXPY(dim, alpha, &A[h], dim, &B[h], dim);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings1["DAXPY (NN)"][k] = t;

        // DDOT (11)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            C_DDOT(full_dim, A, 1, B, 1);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings1["DDOT (11)"][k] = t;

        // DDOT (N1)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            for (int h = 0; h < dim; h++)
                C_DDOT(dim, &A[h], dim, B, 1);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings1["DDOT (N1)"][k] = t;

        // DDOT (1N)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            for (int h = 0; h < dim; h++)
                C_DDOT(dim, A, 1, &B[h], dim);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings1["DDOT (1N)"][k] = t;

        // DDOT (NN)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            for (int h = 0; h < dim; h++)
                C_DDOT(dim, &A[h], dim, &B[h], dim);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings1["DDOT (NN)"][k] = t;

        // DROT (11)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            C_DROT(full_dim, A, 1, B, 1, cosa, sina);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings1["DROT (11)"][k] = t;

        // DROT (N1)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            for (int h = 0; h < dim; h++)
                C_DROT(dim, &A[h], dim, B, 1, cosa, sina);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings1["DROT (N1)"][k] = t;

        // DROT (1N)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            for (int h = 0; h < dim; h++)
                C_DROT(dim, A, 1, &B[h], dim, cosa, sina);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings1["DROT (1N)"][k] = t;

        // DROT (NN)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            for (int h = 0; h < dim; h++)
                C_DROT(dim, &A[h], dim, &B[h], dim, cosa, sina);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings1["DROT (NN)"][k] = t;

        free(A);
        free(B);

    }
    fprintf(outfile, "BLAS 1 Timings [s]:\n\n");
    dim = 1;
    fprintf(outfile, "Operation  ");
    for (int k = 0; k < N; k++) {
        dim *= 2;
        fprintf(outfile, "  %9d", dim);
    }
    fprintf(outfile, "\n");
    for (int s = 0; s < ops1.size(); s++) {
        fprintf(outfile, "%-11s", ops1[s].c_str());
        for (int k = 0; k < N; k++) {
            fprintf(outfile, "  %9.3E", timings1[ops1[s]][k]);
        }
        fprintf(outfile, "\n");
    }
    fprintf(outfile, "\n");

    fprintf(outfile, "BLAS 1 Timings Per Double [s]:\n\n");
    dim = 1;
    fprintf(outfile, "Operation  ");
    for (int k = 0; k < N; k++) {
        dim *= 2;
        fprintf(outfile, "  %9d", dim);
    }
    fprintf(outfile, "\n");
    for (int s = 0; s < ops1.size(); s++) {
        fprintf(outfile, "%-11s", ops1[s].c_str());
        dim = 1;
        for (int k = 0; k < N; k++) {
            dim *= 2;
            unsigned long int full_dim = dim * (unsigned long int) dim;
            fprintf(outfile, "  %9.3E", timings1[ops1[s]][k] / (double) full_dim);
        }
        fprintf(outfile, "\n");
    }
    fprintf(outfile, "\n");

    fprintf(outfile, "BLAS 1 FLOPS [Hz]:\n\n");
    dim = 1;
    fprintf(outfile, "Operation  ");
    for (int k = 0; k < N; k++) {
        dim *= 2;
        fprintf(outfile, "  %9d", dim);
    }
    fprintf(outfile, "\n");
    for (int s = 0; s < ops1.size(); s++) {
        fprintf(outfile, "%-11s", ops1[s].c_str());
        dim = 1;
        for (int k = 0; k < N; k++) {
            dim *= 2;
            unsigned long int full_dim = dim * (unsigned long int) dim;
            fprintf(outfile, "  %9.3E", full_dim / timings1[ops1[s]][k]);
        }
        fprintf(outfile, "\n");
    }
    fprintf(outfile, "\n");
    fflush(outfile);
}
void benchmark_blas2(int N, double min_time)
{
    fprintf(outfile, "\n");
    fprintf(outfile, "                              ------------------------------- \n");
    fprintf(outfile, "                              ======> BLAS2 BENCHMARKS <===== \n");
    fprintf(outfile, "                              ------------------------------- \n");
    fprintf(outfile, "\n");

    fprintf(outfile, "  Parameters:\n");
    fprintf(outfile, "   -Minimum runtime (per operation, per size): %14.10f [s].\n", min_time);
    fprintf(outfile, "   -Maximum dimension exponent N: %d. Arrays are D x D = 2^N x 2^N doubles in size. The D\n", N);
    fprintf(outfile, "        value is reported below.\n");
    fprintf(outfile, "\n");

    fprintf(outfile, "  Notes:\n");
    fprintf(outfile, "   -Operations: (OXX) indicates transpose and stride arguments in the order they appear in \n");
    fprintf(outfile, "        the function call. All lda values are D.\n");
    fprintf(outfile, "\n");

    double T;
    unsigned long int rounds;
    double t;
    int dim;
    Timer* qq;

    std::vector<std::string> ops2;
    ops2.push_back("DGEMV (N11)");
    ops2.push_back("DGEMV (NN1)");
    ops2.push_back("DGEMV (N1N)");
    ops2.push_back("DGEMV (NNN)");
    ops2.push_back("DGEMV (T11)");
    ops2.push_back("DGEMV (TN1)");
    ops2.push_back("DGEMV (T1N)");
    ops2.push_back("DGEMV (TNN)");
    ops2.push_back("DGER (11)");
    ops2.push_back("DGER (N1)");
    ops2.push_back("DGER (1N)");
    ops2.push_back("DGER (NN)");

    std::map<std::string, std::vector<double> > timings2;
    for (int op = 0; op < ops2.size(); op++)
        timings2[ops2[op]].resize(N);

    dim = 1;
    for (int k = 0; k < N; k++) {

        dim *= 2;
        unsigned long int full_dim = dim * (unsigned long int) dim;

        double* A = init_array(full_dim);
        double* B = init_array(full_dim);
        double* C = init_array(full_dim);

        for (unsigned long int Q = 0L; Q < full_dim; Q++) {
            A[Q] = rand() / (double) RAND_MAX;
            B[Q] = rand() / (double) RAND_MAX;
        }

        double alpha = 1.0;
        double beta = 0.0;

        // DGEMV (N11)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            C_DGEMV('N', dim, dim, alpha, A, dim, B, 1, beta, C, 1);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings2["DGEMV (N11)"][k] = t;

        // DGEMV (N1N)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            C_DGEMV('N', dim, dim, alpha, A, dim, B, 1, beta, C, dim);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings2["DGEMV (N1N)"][k] = t;

        // DGEMV (NN1)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            C_DGEMV('N', dim, dim, alpha, A, dim, B, dim, beta, C, 1);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings2["DGEMV (NN1)"][k] = t;

        // DGEMV (NNN)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            C_DGEMV('N', dim, dim, alpha, A, dim, B, dim, beta, C, dim);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings2["DGEMV (NNN)"][k] = t;

        // DGEMV (T11)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            C_DGEMV('T', dim, dim, alpha, A, dim, B, 1, beta, C, 1);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings2["DGEMV (T11)"][k] = t;

        // DGEMV (T1N)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            C_DGEMV('T', dim, dim, alpha, A, dim, B, 1, beta, C, dim);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings2["DGEMV (T1N)"][k] = t;

        // DGEMV (TN1)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            C_DGEMV('T', dim, dim, alpha, A, dim, B, dim, beta, C, 1);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings2["DGEMV (TN1)"][k] = t;

        // DGEMV (TNN)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            C_DGEMV('T', dim, dim, alpha, A, dim, B, dim, beta, C, dim);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings2["DGEMV (TNN)"][k] = t;

        // DGER (11)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            C_DGER(dim, dim, alpha, A, 1, B, 1, C, dim);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings2["DGER (11)"][k] = t;

        // DGER (N1)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            C_DGER(dim, dim, alpha, A, dim, B, 1, C, dim);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings2["DGER (N1)"][k] = t;

        // DGER (1N)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            C_DGER(dim, dim, alpha, A, 1, B, dim, C, dim);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings2["DGER (1N)"][k] = t;

        // DGER (NN)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            C_DGER(dim, dim, alpha, A, dim, B, dim, C, dim);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings2["DGER (NN)"][k] = t;

        free(A);
        free(B);
        free(C);


    }

    fprintf(outfile, "BLAS 2 Timings [s]:\n\n");
    dim = 1;
    fprintf(outfile, "Operation  ");
    for (int k = 0; k < N; k++) {
        dim *= 2;
        fprintf(outfile, "  %9d", dim);
    }
    fprintf(outfile, "\n");
    for (int s = 0; s < ops2.size(); s++) {
        fprintf(outfile, "%-11s", ops2[s].c_str());
        for (int k = 0; k < N; k++) {
            fprintf(outfile, "  %9.3E", timings2[ops2[s]][k]);
        }
        fprintf(outfile, "\n");
    }
    fprintf(outfile, "\n");

    fprintf(outfile, "BLAS 2 FLOPS [Hz] (FLOP: += A * B):\n\n");
    dim = 1;
    fprintf(outfile, "Operation  ");
    for (int k = 0; k < N; k++) {
        dim *= 2;
        fprintf(outfile, "  %9d", dim);
    }
    fprintf(outfile, "\n");
    for (int s = 0; s < ops2.size(); s++) {
        fprintf(outfile, "%-11s", ops2[s].c_str());
        dim = 1;
        for (int k = 0; k < N; k++) {
            dim *= 2;
            unsigned long int full_dim = dim * (unsigned long int) dim;
            fprintf(outfile, "  %9.3E", full_dim / timings2[ops2[s]][k]);
        }
        fprintf(outfile, "\n");
    }
    fprintf(outfile, "\n");
    fflush(outfile);
}
void benchmark_blas3(int N, double min_time, int max_threads)
{
    fprintf(outfile, "\n");
    fprintf(outfile, "                              -------------------------------------- \n");
    fprintf(outfile, "                              ======> BLAS3/LAPACK BENCHMARKS <===== \n");
    fprintf(outfile, "                              -------------------------------------- \n");
    fprintf(outfile, "\n");

    fprintf(outfile, "  Parameters:\n");
    fprintf(outfile, "   -Minimum runtime (per operation, per size): %14.10f [s].\n", min_time);
    fprintf(outfile, "   -Maximum dimension exponent N: %d. Arrays are D x D = 2^N x 2^N doubles in size. The D\n", N);
    fprintf(outfile, "        value is reported below.\n");
    fprintf(outfile, "   -Max threads: %d. Currently only supported with MKL.\n", max_threads);
    fprintf(outfile, "\n");

    fprintf(outfile, "  Notes:\n");
    fprintf(outfile, "   -Operations: (OXX) indicates transpose, side, eigenvector request and stride arguments in\n");
    fprintf(outfile, "        the order they appear in the function call. All lda values are D.\n");
    fprintf(outfile, "\n");
    double T;
    unsigned long int rounds;
    double t;
    int dim;
    Timer* qq;

    int max_thread_count = 1;
    #ifdef HAVE_MKL
        max_thread_count = max_threads;
    #endif

    std::vector<std::string> ops3;
    ops3.push_back("DGEMM (NN)");
    ops3.push_back("DGEMM (NT)");
    ops3.push_back("DGEMM (TN)");
    ops3.push_back("DGEMM (TT)");
    ops3.push_back("DSYMM (LU)");
    ops3.push_back("DSYMM (LL)");
    ops3.push_back("DSYMM (RU)");
    ops3.push_back("DSYMM (RL)");
    ops3.push_back("DPOTRF (U)");
    ops3.push_back("DPOTRF (L)");
    ops3.push_back("DGETRF");
    ops3.push_back("DPOTRS (U)");
    ops3.push_back("DPOTRS (L)");
    ops3.push_back("DGETRS");
    ops3.push_back("DGESV");
    ops3.push_back("DPOTRI (U)");
    ops3.push_back("DPOTRI (L)");
    ops3.push_back("DGETRI");
    ops3.push_back("DSYEV (n)");
    ops3.push_back("DSYEV (v)");

    std::map<int, std::map<std::string, std::vector<double> > > timings;
    // Level 3 and LAPACK routines
    for (int thread = 1; thread <= max_thread_count; thread++) {

        if (thread > 4 && thread % 8 != 0) continue;

        std::map<std::string, std::vector<double> > timings3;
        for (int op = 0; op < ops3.size(); op++)
            timings3[ops3[op]].resize(N);

        #ifdef HAVE_MKL
            mkl_set_num_threads(thread);
        #endif

        dim = 1;
        for (int k = 0; k < N; k++) {

            dim *= 2;
            unsigned long int full_dim = dim * (unsigned long int) dim;

            double* A = init_array(full_dim);
            double* Aback = init_array(full_dim);
            double* B = init_array(full_dim);
            double* C = init_array(full_dim);
            int* ipiv = init_int_array(dim);
            double* work = init_array(dim*3);
            int lwork = dim*3;
            double* eig = init_array(dim);

            for (unsigned long int Q = 0; Q < full_dim; Q++)
                A[Q] = rand() / (double) RAND_MAX;

            for (int i = 0; i < dim; i++)
                for (int j = 0; j < dim; j++)
                    B[i*dim + j] = 0.5*(A[i*dim + j] + A[j*dim + i]);

            // Positive definite A
            C_DGEMM('N','N', dim, dim, dim, 1.0, B, dim, B, dim, 0.0, A, dim);

            double alpha = 1.0;
            double beta = 0.0;

            // DGEMM (NN)
            T = 0.0;
            rounds = 0L;
            qq = new Timer();
            while (T < min_time) {
                C_DGEMM('N','N', dim, dim, dim, alpha, A, dim, B, dim, beta, C, dim);
                T = qq->get();
                rounds++;
            }
            delete qq;
            t = T / (double) rounds;
            timings3["DGEMM (NN)"][k] = t;

            // DGEMM (NT)
            T = 0.0;
            rounds = 0L;
            qq = new Timer();
            while (T < min_time) {
                C_DGEMM('N','T', dim, dim, dim, alpha, A, dim, B, dim, beta, C, dim);
                T = qq->get();
                rounds++;
            }
            delete qq;
            t = T / (double) rounds;
            timings3["DGEMM (NT)"][k] = t;

            // DGEMM (TN)
            T = 0.0;
            rounds = 0L;
            qq = new Timer();
            while (T < min_time) {
                C_DGEMM('T','N', dim, dim, dim, alpha, A, dim, B, dim, beta, C, dim);
                T = qq->get();
                rounds++;
            }
            delete qq;
            t = T / (double) rounds;
            timings3["DGEMM (TN)"][k] = t;

            // DGEMM (TT)
            T = 0.0;
            rounds = 0L;
            qq = new Timer();
            while (T < min_time) {
                C_DGEMM('T','T', dim, dim, dim, alpha, A, dim, B, dim, beta, C, dim);
                T = qq->get();
                rounds++;
            }
            delete qq;
            t = T / (double) rounds;
            timings3["DGEMM (TT)"][k] = t;

            // DSYMM (LU)
            T = 0.0;
            rounds = 0L;
            qq = new Timer();
            while (T < min_time) {
                C_DSYMM('L','U', dim, dim, alpha, A, dim, B, dim, beta, C, dim);
                T = qq->get();
                rounds++;
            }
            delete qq;
            t = T / (double) rounds;
            timings3["DSYMM (LU)"][k] = t;

            // DSYMM (LL)
            T = 0.0;
            rounds = 0L;
            qq = new Timer();
            while (T < min_time) {
                C_DSYMM('L','L', dim, dim, alpha, A, dim, B, dim, beta, C, dim);
                T = qq->get();
                rounds++;
            }
            delete qq;
            t = T / (double) rounds;
            timings3["DSYMM (LL)"][k] = t;

            // DSYMM (RU)
            T = 0.0;
            rounds = 0L;
            qq = new Timer();
            while (T < min_time) {
                C_DSYMM('R','U', dim, dim, alpha, A, dim, B, dim, beta, C, dim);
                T = qq->get();
                rounds++;
            }
            delete qq;
            t = T / (double) rounds;
            timings3["DSYMM (RU)"][k] = t;

            // DSYMM (RL)
            T = 0.0;
            rounds = 0L;
            qq = new Timer();
            while (T < min_time) {
                C_DSYMM('R','L', dim, dim, alpha, A, dim, B, dim, beta, C, dim);
                T = qq->get();
                rounds++;
            }
            delete qq;
            t = T / (double) rounds;
            timings3["DSYMM (RL)"][k] = t;

            // DPOTRF (U)
            C_DCOPY(full_dim, A, 1, C, 1);
            C_DCOPY(full_dim, A, 1, Aback, 1);

            T = 0.0;
            rounds = 0L;
            qq = new Timer();
            while (T < min_time) {
                C_DCOPY(full_dim, C, 1, A, 1);
                C_DPOTRF('U', dim, A, dim);
                T = qq->get();
                rounds++;
            }
            delete qq;
            t = T / (double) rounds;
            timings3["DPOTRF (U)"][k] = t;

            // DPOTRS (U)
            C_DCOPY(full_dim, B, 1, C, 1);
            T = 0.0;
            rounds = 0L;
            qq = new Timer();
            while (T < min_time) {
                C_DCOPY(full_dim, C, 1, B, 1);
                C_DPOTRS('U', dim, dim, A, dim, B, dim);
                T = qq->get();
                rounds++;
            }
            delete qq;
            t = T / (double) rounds;
            timings3["DPOTRS (U)"][k] = t;

            C_DCOPY(full_dim, C, 1, B, 1);

            // DPOTRI (U)
            C_DCOPY(full_dim, A, 1, C, 1);
            T = 0.0;
            rounds = 0L;
            qq = new Timer();
            while (T < min_time) {
                C_DCOPY(full_dim, C, 1, A, 1);
                C_DPOTRI('U', dim, A, dim);
                T = qq->get();
                rounds++;
            }
            delete qq;
            t = T / (double) rounds;
            timings3["DPOTRI (U)"][k] = t;
            C_DCOPY(full_dim, Aback, 1, A, 1);

            // DPOTRF (L)
            C_DCOPY(full_dim, A, 1, C, 1);
            C_DCOPY(full_dim, A, 1, Aback, 1);

            T = 0.0;
            rounds = 0L;
            qq = new Timer();
            while (T < min_time) {
                C_DCOPY(full_dim, C, 1, A, 1);
                C_DPOTRF('L', dim, A, dim);
                T = qq->get();
                rounds++;
            }
            delete qq;
            t = T / (double) rounds;
            timings3["DPOTRF (L)"][k] = t;

            // DPOTRS (L)
            C_DCOPY(full_dim, B, 1, C, 1);
            T = 0.0;
            rounds = 0L;
            qq = new Timer();
            while (T < min_time) {
                C_DCOPY(full_dim, C, 1, B, 1);
                C_DPOTRS('L', dim, dim, A, dim, B, dim);
                T = qq->get();
                rounds++;
            }
            delete qq;
            t = T / (double) rounds;
            timings3["DPOTRS (L)"][k] = t;

            C_DCOPY(full_dim, C, 1, B, 1);

            // DPOTRI (L)
            C_DCOPY(full_dim, A, 1, C, 1);
            T = 0.0;
            rounds = 0L;
            qq = new Timer();
            while (T < min_time) {
                C_DCOPY(full_dim, C, 1, A, 1);
                C_DPOTRI('L', dim, A, dim);
                T = qq->get();
                rounds++;
            }
            delete qq;
            t = T / (double) rounds;
            timings3["DPOTRI (L)"][k] = t;
            C_DCOPY(full_dim, Aback, 1, A, 1);

            // DGETRF
            C_DCOPY(full_dim, A, 1, C, 1);
            C_DCOPY(full_dim, A, 1, Aback, 1);

            T = 0.0;
            rounds = 0L;
            qq = new Timer();
            while (T < min_time) {
                C_DCOPY(full_dim, C, 1, A, 1);
                C_DGETRF(dim, dim, A, dim, ipiv);
                T = qq->get();
                rounds++;
            }
            delete qq;
            t = T / (double) rounds;
            timings3["DGETRF"][k] = t;

            // DGETRS
            C_DCOPY(full_dim, B, 1, C, 1);
            T = 0.0;
            rounds = 0L;
            qq = new Timer();
            while (T < min_time) {
                C_DCOPY(full_dim, C, 1, B, 1);
                C_DGETRS('N', dim, dim, A, dim, ipiv, B, dim);
                T = qq->get();
                rounds++;
            }
            delete qq;
            t = T / (double) rounds;
            timings3["DGETRS"][k] = t;

            C_DCOPY(full_dim, C, 1, B, 1);

            // DGETRI
            C_DCOPY(full_dim, A, 1, C, 1);
            T = 0.0;
            rounds = 0L;
            qq = new Timer();
            while (T < min_time) {
                C_DCOPY(full_dim, C, 1, A, 1);
                C_DGETRI(dim, A, dim, ipiv, work, lwork);
                T = qq->get();
                rounds++;
            }
            delete qq;
            t = T / (double) rounds;
            timings3["DGETRI"][k] = t;
            C_DCOPY(full_dim, Aback, 1, A, 1);

            // DGESV
            C_DCOPY(full_dim, B, 1, C, 1);
            C_DCOPY(full_dim, A, 1, Aback, 1);
            T = 0.0;
            rounds = 0L;
            qq = new Timer();
            while (T < min_time) {
                C_DCOPY(full_dim, C, 1, B, 1);
                C_DCOPY(full_dim, Aback, 1, A, 1);
                C_DGESV(dim, dim, A, dim, ipiv, B, dim);
                T = qq->get();
                rounds++;
            }
            delete qq;
            t = T / (double) rounds;
            timings3["DGESV"][k] = t;

            C_DCOPY(full_dim, C, 1, B, 1);
            C_DCOPY(full_dim, Aback, 1, A, 1);

            // DSYEV (n)
            C_DCOPY(full_dim, A, 1, Aback, 1);
            T = 0.0;
            rounds = 0L;
            qq = new Timer();
            while (T < min_time) {
                C_DCOPY(full_dim, Aback, 1, A, 1);
                C_DSYEV('n', 'u', dim, A, dim, eig, work, lwork);
                T = qq->get();
                rounds++;
            }
            delete qq;
            t = T / (double) rounds;
            timings3["DSYEV (n)"][k] = t;

            C_DCOPY(full_dim, Aback, 1, A, 1);

            // DSYEV (v)
            C_DCOPY(full_dim, A, 1, Aback, 1);
            T = 0.0;
            rounds = 0L;
            qq = new Timer();
            while (T < min_time) {
                C_DCOPY(full_dim, Aback, 1, A, 1);
                C_DSYEV('v', 'u', dim, A, dim, eig, work, lwork);
                T = qq->get();
                rounds++;
            }
            delete qq;
            t = T / (double) rounds;
            timings3["DSYEV (v)"][k] = t;

            C_DCOPY(full_dim, Aback, 1, A, 1);

            free(eig);
            free(work);
            free(ipiv);
            free(A);
            free(B);
            free(C);
            free(Aback);

        }
        timings[thread] = timings3;

        fprintf(outfile, "BLAS 3 Timings [s], Threads = %d:\n\n", thread);
        dim = 1;
        fprintf(outfile, "Operation  ");
        for (int k = 0; k < N; k++) {
            dim *= 2;
            fprintf(outfile, "  %9d", dim);
        }
        fprintf(outfile, "\n");
        for (int s = 0; s < ops3.size(); s++) {
            fprintf(outfile, "%-11s", ops3[s].c_str());
            for (int k = 0; k < N; k++) {
                fprintf(outfile, "  %9.3E", timings3[ops3[s]][k]);
            }
            fprintf(outfile, "\n");
        }
        fprintf(outfile, "\n");

        fprintf(outfile, "BLAS 3 Effective DGEMMs [-] (NN):\n\n");
        dim = 1;
        fprintf(outfile, "Operation  ");
        for (int k = 0; k < N; k++) {
            dim *= 2;
            fprintf(outfile, "  %9d", dim);
        }
        fprintf(outfile, "\n");
        for (int s = 0; s < ops3.size(); s++) {
            fprintf(outfile, "%-11s", ops3[s].c_str());
            dim = 1;
            for (int k = 0; k < N; k++) {
                dim *= 2;
                unsigned long int full_dim = dim * (unsigned long int) dim;
                fprintf(outfile, "  %9.3E", timings3[ops3[s]][k] / timings3["DGEMM (NN)"][k]);
            }
            fprintf(outfile, "\n");
        }
        fprintf(outfile, "\n");
    }
    if (max_threads > 1) {
        fprintf(outfile, "BLAS 3 Speedups [-]:\n\n");
        dim = 1;
        fprintf(outfile, "Operation  Threads  ");
        for (int k = 0; k < N; k++) {
            dim *= 2;
            fprintf(outfile, "  %9d", dim);
        }
        fprintf(outfile, "\n");
        for (int s = 0; s < ops3.size(); s++) {
            for (int thread = 1; thread <= max_threads; thread++) {
                if (thread > 4 && thread % 8 != 0) continue;
                fprintf(outfile, "%-11s  %-7d", ops3[s].c_str(), thread);
                dim = 1;
                for (int k = 0; k < N; k++) {
                    dim *= 2;
                    unsigned long int full_dim = dim * (unsigned long int) dim;
                    fprintf(outfile, "  %9.3E", timings[1][ops3[s]][k] / timings[thread][ops3[s]][k]);
                }
                fprintf(outfile, "\n");
            }
        }
        fprintf(outfile, "\n");

        fprintf(outfile, "BLAS 3 Parallel Efficiency [-]:\n\n");
        dim = 1;
        fprintf(outfile, "Operation  Threads  ");
        for (int k = 0; k < N; k++) {
            dim *= 2;
            fprintf(outfile, "  %9d", dim);
        }
        fprintf(outfile, "\n");
        for (int s = 0; s < ops3.size(); s++) {
            for (int thread = 1; thread <= max_threads; thread++) {
                if (thread > 4 && thread % 8 != 0) continue;
                fprintf(outfile, "%-11s  %-7d", ops3[s].c_str(), thread);
                dim = 1;
                for (int k = 0; k < N; k++) {
                    dim *= 2;
                    unsigned long int full_dim = dim * (unsigned long int) dim;
                    fprintf(outfile, "  %9.3E", timings[1][ops3[s]][k] / (timings[thread][ops3[s]][k] * (double) thread));
                }
                fprintf(outfile, "\n");
            }
        }
        fprintf(outfile, "\n");
    }
    fflush(outfile);

}
void benchmark_disk(int N, double min_time)
{
    fprintf(outfile, "\n");
    fprintf(outfile, "                              ------------------------------ \n");
    fprintf(outfile, "                              ======> PSIO BENCHMARKS <===== \n");
    fprintf(outfile, "                              ------------------------------ \n");
    fprintf(outfile, "\n");

    fprintf(outfile, "  Parameters:\n");
    fprintf(outfile, "   -Minimum runtime (per operation, per size): %14.10f [s].\n", min_time);
    fprintf(outfile, "   -Maximum dimension exponent N: %d. Arrays are D x D = 2^N x 2^N doubles in size. The D\n", N);
    fprintf(outfile, "        value is reported below\n");
    fprintf(outfile, "\n");

    fprintf(outfile, "  Operations:\n");
    fprintf(outfile, "   -OPEN/CLOSE: Open and close a file repeatedly witout discard (Data rates are meaningless).\n");
    fprintf(outfile, "   -ZERO: Write the first pass of data, expanding the file. Performed in one op. Timing may\n");
    fprintf(outfile, "        be inaccurate, as only one pass is performed.\n");
    fprintf(outfile, "   -READ (Continuous): Repeatedly read the entire array in one operation of dimension N x N.\n");
    fprintf(outfile, "   -READ (Blocked): Repeatedly read the entire array in N operations of dimension N.\n");
    fprintf(outfile, "   -READ (Transposed): Repeatedly read the entire array in N operations of dimension N. Arrays\n");
    fprintf(outfile, "        are staggered to simulate reading the N/2 blocked transpose of the array.\n");
    fprintf(outfile, "   -WRITE (Continuous): Repeatedly write the entire array in one operation of dimension N x N.\n");
    fprintf(outfile, "   -WRITE (Blocked): Repeatedly write the entire array in N operations of dimension N.\n");
    fprintf(outfile, "   -WRITE (Transposed): Repeatedly write the entire array in N operations of dimension N. Arrays\n");
    fprintf(outfile, "        are staggered to simulate writing the N/2 blocked transpose of the array.\n");
    fprintf(outfile, "\n");

    double T;
    unsigned long int rounds;
    double t;
    int dim;
    Timer* qq;

    std::map<std::string, std::vector<double> > timings;
    std::vector<std::string> ops;

    ops.push_back("OPEN/CLOSE (Reuse)");
    ops.push_back("ZERO (First Touch)");
    ops.push_back("READ (Continuous)");
    ops.push_back("READ (Blocked)");
    ops.push_back("READ (Transposed)");
    ops.push_back("WRITE (Continuous)");
    ops.push_back("WRITE (Blocked)");
    ops.push_back("WRITE (Transposed)");
    for (int op = 0; op < ops.size(); op++)
        timings[ops[op]].resize(N);

    boost::shared_ptr<PSIO> psio_ = PSIO::shared_object();
    psio_address psiadd;
    dim = 1;
    for (int k = 0; k < N; k++) {

        dim *= 2;
        unsigned long int full_dim = dim * (unsigned long int) dim;

        double* A = init_array(full_dim);
        psio_->open(0, PSIO_OPEN_NEW);

        // ZERO
        psiadd = PSIO_ZERO;
        T = 0.0;
        qq = new Timer();
        t = qq->get();
        psio_->write(0,"BENCH_DATA", (char*) &A[0], full_dim * sizeof(double), psiadd, &psiadd);
        delete qq;
        timings["ZERO (First Touch)"][k] = t;

        // Write (Continuous)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            psiadd = PSIO_ZERO;
            psio_->write(0,"BENCH_DATA", (char*) &A[0], full_dim * sizeof(double), psiadd, &psiadd);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings["WRITE (Continuous)"][k] = t;

        // Write (Blocked)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            psiadd = PSIO_ZERO;
            for (int Q = 0; Q < dim; Q++)
                psio_->write(0,"BENCH_DATA", (char*) &A[Q*(unsigned long int) dim], dim * sizeof(double), psiadd, &psiadd);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings["WRITE (Blocked)"][k] = t;

        // Write (Transposed)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            psiadd = PSIO_ZERO;
            for (int Q = 0; Q < dim; Q++) {
                if (Q < dim/2)
                    psiadd = psio_get_address(PSIO_ZERO, 2*Q*dim*sizeof(double));
                else
                    psiadd = psio_get_address(PSIO_ZERO, (Q - dim/2 + 1)*dim*sizeof(double));
                psio_->write(0,"BENCH_DATA", (char*) &A[Q*(unsigned long int) dim], dim * sizeof(double), psiadd, &psiadd);
            }
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings["WRITE (Transposed)"][k] = t;

        // Read (Continuous)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            psiadd = PSIO_ZERO;
            psio_->read(0,"BENCH_DATA", (char*) &A[0], full_dim * sizeof(double), psiadd, &psiadd);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings["READ (Continuous)"][k] = t;

        // Read (Blocked)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            psiadd = PSIO_ZERO;
            for (int Q = 0; Q < dim; Q++)
                psio_->read(0,"BENCH_DATA", (char*) &A[Q*(unsigned long int) dim], dim * sizeof(double), psiadd, &psiadd);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings["READ (Blocked)"][k] = t;

        // Read (Transposed)
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            psiadd = PSIO_ZERO;
            for (int Q = 0; Q < dim; Q++) {
                if (Q < dim/2)
                    psiadd = psio_get_address(PSIO_ZERO, 2*Q*dim*sizeof(double));
                else
                    psiadd = psio_get_address(PSIO_ZERO, (Q - dim/2 + 1)*dim*sizeof(double));
                psio_->read(0,"BENCH_DATA", (char*) &A[Q*(unsigned long int) dim], dim * sizeof(double), psiadd, &psiadd);
            }
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings["READ (Transposed)"][k] = t;

        // Open/Close
        T = 0.0;
        rounds = 0L;
        qq = new Timer();
        while (T < min_time) {
            psio_->close(0, 1);
            psio_->open(0, PSIO_OPEN_OLD);
            T = qq->get();
            rounds++;
        }
        delete qq;
        t = T / (double) rounds;
        timings["OPEN/CLOSE (Reuse)"][k] = t;

        psio_->close(0, 0);
        free(A);
    }
    fprintf(outfile, "PSIO Timings [s]\n\n");
    dim = 1;
    fprintf(outfile, "Operation           ");
    for (int k = 0; k < N; k++) {
        dim *= 2;
        fprintf(outfile, "  %9d", dim);
    }
    fprintf(outfile, "\n");
    for (int s = 0; s < ops.size(); s++) {
        fprintf(outfile, "%-20s", ops[s].c_str());
        for (int k = 0; k < N; k++) {
            fprintf(outfile, "  %9.3E", timings[ops[s]][k]);
        }
        fprintf(outfile, "\n");
    }
    fprintf(outfile, "\n");

    fprintf(outfile, "PSIO Performance [doubles/s]\n\n");
    dim = 1;
    fprintf(outfile, "Operation           ");
    for (int k = 0; k < N; k++) {
        dim *= 2;
        fprintf(outfile, "  %9d", dim);
    }
    fprintf(outfile, "\n");
    for (int s = 0; s < ops.size(); s++) {
        fprintf(outfile, "%-20s", ops[s].c_str());
        dim = 1;
        for (int k = 0; k < N; k++) {
            dim *= 2;
            unsigned long int full_dim = dim * (unsigned long int) dim;
            fprintf(outfile, "  %9.3E", full_dim / timings[ops[s]][k]);
        }
        fprintf(outfile, "\n");
    }
    fprintf(outfile, "\n");

    fprintf(outfile, "PSIO Performance [GiB/s]\n\n");
    dim = 1;
    fprintf(outfile, "Operation           ");
    for (int k = 0; k < N; k++) {
        dim *= 2;
        fprintf(outfile, "  %9d", dim);
    }
    fprintf(outfile, "\n");
    for (int s = 0; s < ops.size(); s++) {
        fprintf(outfile, "%-20s", ops[s].c_str());
        dim = 1;
        for (int k = 0; k < N; k++) {
            dim *= 2;
            unsigned long int full_dim = dim * (unsigned long int) dim;
            fprintf(outfile, "  %9.3E", 8.0E-9 *full_dim / timings[ops[s]][k]);
        }
        fprintf(outfile, "\n");
    }
    fprintf(outfile, "\n");
    fflush(outfile);

}
void benchmark_math(double min_time)
{
    double T;
    unsigned long int rounds;
    double t;
    int dim;
    Timer* qq;

    double a, b, c;

    std::map<std::string, double> timings;
    std::vector<std::string> ops;

    ops.push_back("+");
    ops.push_back("*");
    ops.push_back("-");
    ops.push_back("/");
    ops.push_back("+=");
    ops.push_back("*=");
    ops.push_back("-=");
    ops.push_back("/=");
    ops.push_back("sin");
    ops.push_back("cos");
    ops.push_back("tan");
    ops.push_back("asin");
    ops.push_back("acos");
    ops.push_back("atan");
    ops.push_back("atan2");
    ops.push_back("cosh");
    ops.push_back("sinh");
    ops.push_back("tanh");
    ops.push_back("exp");
    ops.push_back("pow");
    ops.push_back("log");
    ops.push_back("floor");
    ops.push_back("ceil");
    ops.push_back("fabs");

    // In case the compiler gets awesome
    FILE* fh = fopen("dump.dat", "w");

    #define LOOP_SIZE 10000
    #define UNROLL_SIZE 10

    T = 0.0;
    rounds = 0L;
    qq = new Timer();
    while (T < min_time) {
        c = 0.0;
        a = 0.0000001;
        for (int Q = 0; Q < LOOP_SIZE; Q++) {
            c = a + c;
            c = a + c;
            c = a + c;
            c = a + c;
            c = a + c;
            c = a + c;
            c = a + c;
            c = a + c;
            c = a + c;
            c = a + c;
        }

        T = qq->get();
        rounds++;
    }
    fprintf(fh, "%14.10f\n", c);
    delete qq;
    t = T / (double) (rounds * LOOP_SIZE * (unsigned long int)UNROLL_SIZE);
    timings["+"] = t;

    T = 0.0;
    rounds = 0L;
    qq = new Timer();
    while (T < min_time) {
        c = 0.00010;
        a = 1.0000001;
        for (int Q = 0; Q < LOOP_SIZE; Q++) {
            c = a - c;
            c = a - c;
            c = a - c;
            c = a - c;
            c = a - c;
            c = a - c;
            c = a - c;
            c = a - c;
            c = a - c;
            c = a - c;
        }

        T = qq->get();
        rounds++;
    }
    fprintf(fh, "%14.10f\n", c);
    delete qq;
    t = T / (double) (rounds * LOOP_SIZE * (unsigned long int)UNROLL_SIZE);
    timings["-"] = t;

    T = 0.0;
    rounds = 0L;
    qq = new Timer();
    while (T < min_time) {
        a = 1.000000000001;
        c = 1.0;
        for (int Q = 0; Q < LOOP_SIZE; Q++) {
            c = a * c;
            c = a * c;
            c = a * c;
            c = a * c;
            c = a * c;
            c = a * c;
            c = a * c;
            c = a * c;
            c = a * c;
            c = a * c;
        }

        T = qq->get();
        rounds++;
    }
    fprintf(fh, "%14.10f\n", c);
    delete qq;
    t = T / (double) (rounds * LOOP_SIZE * (unsigned long int)UNROLL_SIZE);
    timings["*"] = t;

    T = 0.0;
    rounds = 0L;
    qq = new Timer();
    while (T < min_time) {
        a = 1.00001;
        c = 1.0;
        for (int Q = 0; Q < LOOP_SIZE; Q++) {
            c = c / a;
            c = c / a;
            c = c / a;
            c = c / a;
            c = c / a;
            c = c / a;
            c = c / a;
            c = c / a;
            c = c / a;
            c = c / a;
        }

        T = qq->get();
        rounds++;
    }
    fprintf(fh, "%14.10f\n", c);
    delete qq;
    t = T / (double) (rounds * LOOP_SIZE * (unsigned long int)UNROLL_SIZE);
    timings["/"] = t;
    T = 0.0;
    rounds = 0L;
    qq = new Timer();
    while (T < min_time) {
        c = 1.0;
        for (int Q = 0; Q < LOOP_SIZE; Q++) {
            c += b;
            c += b;
            c += b;
            c += b;
            c += b;
            c += b;
            c += b;
            c += b;
            c += b;
            c += b;
        }

        T = qq->get();
        rounds++;
    }
    fprintf(fh, "%14.10f\n", c);
    delete qq;
    t = T / (double) (rounds * LOOP_SIZE * (unsigned long int)UNROLL_SIZE);
    timings["+="] = t;

    T = 0.0;
    rounds = 0L;
    qq = new Timer();
    while (T < min_time) {
        c = 1.032;
        b = 1.023;
        for (int Q = 0; Q < LOOP_SIZE; Q++) {
            c -= b;
            c -= b;
            c -= b;
            c -= b;
            c -= b;
            c -= b;
            c -= b;
            c -= b;
            c -= b;
            c -= b;
        }

        T = qq->get();
        rounds++;
    }
    fprintf(fh, "%14.10f\n", c);
    delete qq;
    t = T / (double) (rounds * LOOP_SIZE * (unsigned long int)UNROLL_SIZE);
    timings["-="] = t;

    T = 0.0;
    rounds = 0L;
    qq = new Timer();
    while (T < min_time) {
        c = 1.0000001;
        b = 1.000000001;
        for (int Q = 0; Q < LOOP_SIZE; Q++) {
            c *= b;
            c *= b;
            c *= b;
            c *= b;
            c *= b;
            c *= b;
            c *= b;
            c *= b;
            c *= b;
            c *= b;
        }

        T = qq->get();
        rounds++;
    }
    fprintf(fh, "%14.10f\n", c);
    delete qq;
    t = T / (double) (rounds * LOOP_SIZE * (unsigned long int)UNROLL_SIZE);
    timings["*="] = t;

    T = 0.0;
    rounds = 0L;
    qq = new Timer();
    while (T < min_time) {
        c = 1.0;
        b = 1.00001;
        for (int Q = 0; Q < LOOP_SIZE; Q++) {
            c /= b;
            c /= b;
            c /= b;
            c /= b;
            c /= b;
            c /= b;
            c /= b;
            c /= b;
            c /= b;
            c /= b;
        }

        T = qq->get();
        rounds++;
    }
    fprintf(fh, "%14.10f\n", c);
    delete qq;
    t = T / (double) (rounds * LOOP_SIZE * (unsigned long int)UNROLL_SIZE);
    timings["/="] = t;

    T = 0.0;
    rounds = 0L;
    qq = new Timer();
    while (T < min_time) {
        c = 1.0;
        for (int Q = 0; Q < LOOP_SIZE; Q++) {
            c = sin(c);
            c = sin(c);
            c = sin(c);
            c = sin(c);
            c = sin(c);
            c = sin(c);
            c = sin(c);
            c = sin(c);
            c = sin(c);
            c = sin(c);
        }

        T = qq->get();
        rounds++;
    }
    fprintf(fh, "%14.10f\n", c);
    delete qq;
    t = T / (double) (rounds * LOOP_SIZE * (unsigned long int)UNROLL_SIZE);
    timings["sin"] = t;

    T = 0.0;
    rounds = 0L;
    qq = new Timer();
    while (T < min_time) {
        c = 1.0;
        for (int Q = 0; Q < LOOP_SIZE; Q++) {
            c = cos(c);
            c = cos(c);
            c = cos(c);
            c = cos(c);
            c = cos(c);
            c = cos(c);
            c = cos(c);
            c = cos(c);
            c = cos(c);
            c = cos(c);
        }

        T = qq->get();
        rounds++;
    }
    fprintf(fh, "%14.10f\n", c);
    delete qq;
    t = T / (double) (rounds * LOOP_SIZE * (unsigned long int)UNROLL_SIZE);
    timings["cos"] = t;

    T = 0.0;
    rounds = 0L;
    qq = new Timer();
    while (T < min_time) {
        c = 2.0;
        for (int Q = 0; Q < LOOP_SIZE; Q++) {
            c = tan(c);
            c = tan(c);
            c = tan(c);
            c = tan(c);
            c = tan(c);
            c = tan(c);
            c = tan(c);
            c = tan(c);
            c = tan(c);
            c = tan(c);
        }

        T = qq->get();
        rounds++;
    }
    fprintf(fh, "%14.10f\n", c);
    delete qq;
    t = T / (double) (rounds * LOOP_SIZE * (unsigned long int)UNROLL_SIZE);
    timings["tan"] = t;

    b = 0.8;
    T = 0.0;
    rounds = 0L;
    qq = new Timer();
    while (T < min_time) {
        c = 0.9;
        for (int Q = 0; Q < LOOP_SIZE; Q++) {
            c = asin(c);
            c = asin(c);
            c = asin(c);
            c = asin(c);
            c = asin(c);
            c = asin(c);
            c = asin(c);
            c = asin(c);
            c = asin(c);
            c = asin(c);
        }

        T = qq->get();
        rounds++;
    }
    fprintf(fh, "%14.10f\n", c);
    delete qq;
    t = T / (double) (rounds * LOOP_SIZE * (unsigned long int)UNROLL_SIZE);
    timings["asin"] = t;

    T = 0.0;
    rounds = 0L;
    qq = new Timer();
    while (T < min_time) {
        c = 0.9;
        for (int Q = 0; Q < LOOP_SIZE; Q++) {
            c = acos(c);
            c = acos(c);
            c = acos(c);
            c = acos(c);
            c = acos(c);
            c = acos(c);
            c = acos(c);
            c = acos(c);
            c = acos(c);
            c = acos(c);
        }

        T = qq->get();
        rounds++;
    }
    fprintf(fh, "%14.10f\n", c);
    delete qq;
    t = T / (double) (rounds * LOOP_SIZE * (unsigned long int)UNROLL_SIZE);
    timings["acos"] = t;

    T = 0.0;
    rounds = 0L;
    qq = new Timer();
    while (T < min_time) {
        c = 4.0;
        for (int Q = 0; Q < LOOP_SIZE; Q++) {
            c = atan(c);
            c = atan(c);
            c = atan(c);
            c = atan(c);
            c = atan(c);
            c = atan(c);
            c = atan(c);
            c = atan(c);
            c = atan(c);
            c = atan(c);
        }

        T = qq->get();
        rounds++;
    }
    fprintf(fh, "%14.10f\n", c);
    delete qq;
    t = T / (double) (rounds * LOOP_SIZE * (unsigned long int)UNROLL_SIZE);
    timings["atan"] = t;

    T = 0.0;
    rounds = 0L;
    qq = new Timer();
    while (T < min_time) {
        a = 1.0;
        c = 4.0;
        for (int Q = 0; Q < LOOP_SIZE; Q++) {
            c = atan2(c,a);
            c = atan2(c,a);
            c = atan2(c,a);
            c = atan2(c,a);
            c = atan2(c,a);
            c = atan2(c,a);
            c = atan2(c,a);
            c = atan2(c,a);
            c = atan2(c,a);
            c = atan2(c,a);
        }

        T = qq->get();
        rounds++;
    }
    fprintf(fh, "%14.10f\n", c);
    delete qq;
    t = T / (double) (rounds * LOOP_SIZE * (unsigned long int)UNROLL_SIZE);
    timings["atan2"] = t;

    T = 0.0;
    rounds = 0L;
    qq = new Timer();
    while (T < min_time) {
        c = 1.0;
        for (int Q = 0; Q < LOOP_SIZE; Q++) {
            c = cosh(c);
            c = cosh(c);
            c = cosh(c);
            c = cosh(c);
            c = cosh(c);
            c = cosh(c);
            c = cosh(c);
            c = cosh(c);
            c = cosh(c);
            c = cosh(c);
        }

        T = qq->get();
        rounds++;
    }
    fprintf(fh, "%14.10f\n", c);
    delete qq;
    t = T / (double) (rounds * LOOP_SIZE * (unsigned long int)UNROLL_SIZE);
    timings["cosh"] = t;

    T = 0.0;
    rounds = 0L;
    qq = new Timer();
    while (T < min_time) {
        c = 1.0;
        for (int Q = 0; Q < LOOP_SIZE; Q++) {
            c = sinh(c);
            c = sinh(c);
            c = sinh(c);
            c = sinh(c);
            c = sinh(c);
            c = sinh(c);
            c = sinh(c);
            c = sinh(c);
            c = sinh(c);
            c = sinh(c);
        }

        T = qq->get();
        rounds++;
    }
    fprintf(fh, "%14.10f\n", c);
    delete qq;
    t = T / (double) (rounds * LOOP_SIZE * (unsigned long int)UNROLL_SIZE);
    timings["sinh"] = t;

    T = 0.0;
    rounds = 0L;
    qq = new Timer();
    while (T < min_time) {
        c = 1.0;
        for (int Q = 0; Q < LOOP_SIZE; Q++) {
            c = tanh(c);
            c = tanh(c);
            c = tanh(c);
            c = tanh(c);
            c = tanh(c);
            c = tanh(c);
            c = tanh(c);
            c = tanh(c);
            c = tanh(c);
            c = tanh(c);
        }

        T = qq->get();
        rounds++;
    }
    fprintf(fh, "%14.10f\n", c);
    delete qq;
    t = T / (double) (rounds * LOOP_SIZE * (unsigned long int)UNROLL_SIZE);
    timings["tanh"] = t;

    double temp[2];
    temp[0] = 2.3;
    temp[1] = 0.0;
    T = 0.0;
    rounds = 0L;
    qq = new Timer();
    while (T < min_time) {
        for (int Q = 0; Q < LOOP_SIZE; Q++) {
            temp[1] = exp(temp[0]);
            temp[1] = exp(temp[0]);
            temp[1] = exp(temp[0]);
            temp[1] = exp(temp[0]);
            temp[1] = exp(temp[0]);
            temp[1] = exp(temp[0]);
            temp[1] = exp(temp[0]);
            temp[1] = exp(temp[0]);
            temp[1] = exp(temp[0]);
            temp[1] = exp(temp[0]);
        }

        T = qq->get();
        rounds++;
    }
    fprintf(fh, "%14.10f\n", c);
    delete qq;
    t = T / (double) (rounds * LOOP_SIZE * (unsigned long int)UNROLL_SIZE);
    timings["exp"] = t;

    T = 0.0;
    rounds = 0L;
    qq = new Timer();
    while (T < min_time) {
        for (int Q = 0; Q < LOOP_SIZE; Q++) {
            temp[1] = pow(a,temp[0]);
            temp[1] = pow(a,temp[0]);
            temp[1] = pow(a,temp[0]);
            temp[1] = pow(a,temp[0]);
            temp[1] = pow(a,temp[0]);
            temp[1] = pow(a,temp[0]);
            temp[1] = pow(a,temp[0]);
            temp[1] = pow(a,temp[0]);
            temp[1] = pow(a,temp[0]);
            temp[1] = pow(a,temp[0]);
        }

        T = qq->get();
        rounds++;
    }
    fprintf(fh, "%14.10f\n", c);
    delete qq;
    t = T / (double) (rounds * LOOP_SIZE * (unsigned long int)UNROLL_SIZE);
    timings["pow"] = t;

    T = 0.0;
    rounds = 0L;
    qq = new Timer();
    while (T < min_time) {
        for (int Q = 0; Q < LOOP_SIZE; Q++) {
            temp[1] = log(temp[0]);
            temp[1] = log(temp[0]);
            temp[1] = log(temp[0]);
            temp[1] = log(temp[0]);
            temp[1] = log(temp[0]);
            temp[1] = log(temp[0]);
            temp[1] = log(temp[0]);
            temp[1] = log(temp[0]);
            temp[1] = log(temp[0]);
            temp[1] = log(temp[0]);
        }

        T = qq->get();
        rounds++;
    }
    fprintf(fh, "%14.10f\n", c);
    delete qq;
    t = T / (double) (rounds * LOOP_SIZE * (unsigned long int)UNROLL_SIZE);
    timings["log"] = t;

    T = 0.0;
    c = 1.2;
    rounds = 0L;
    qq = new Timer();
    while (T < min_time) {
        for (int Q = 0; Q < LOOP_SIZE; Q++) {
            c = floor(c);
            c = floor(c);
            c = floor(c);
            c = floor(c);
            c = floor(c);
            c = floor(c);
            c = floor(c);
            c = floor(c);
            c = floor(c);
            c = floor(c);
        }

        T = qq->get();
        rounds++;
    }
    fprintf(fh, "%14.10f\n", c);
    delete qq;
    t = T / (double) (rounds * LOOP_SIZE * (unsigned long int)UNROLL_SIZE);
    timings["floor"] = t;

    T = 0.0;
    rounds = 0L;
    qq = new Timer();
    while (T < min_time) {
        for (int Q = 0; Q < LOOP_SIZE; Q++) {
            c = ceil(c);
            c = ceil(c);
            c = ceil(c);
            c = ceil(c);
            c = ceil(c);
            c = ceil(c);
            c = ceil(c);
            c = ceil(c);
            c = ceil(c);
            c = ceil(c);
        }

        T = qq->get();
        rounds++;
    }
    fprintf(fh, "%14.10f\n", c);
    delete qq;
    t = T / (double) (rounds * LOOP_SIZE * (unsigned long int)UNROLL_SIZE);
    timings["ceil"] = t;

    T = 0.0;
    rounds = 0L;
    qq = new Timer();
    while (T < min_time) {
        for (int Q = 0; Q < LOOP_SIZE; Q++) {
            c = fabs(c);
            c = fabs(c);
            c = fabs(c);
            c = fabs(c);
            c = fabs(c);
            c = fabs(c);
            c = fabs(c);
            c = fabs(c);
            c = fabs(c);
            c = fabs(c);
        }

        T = qq->get();
        rounds++;
    }
    fprintf(fh, "%14.10f\n", c);
    delete qq;
    t = T / (double) (rounds * LOOP_SIZE * (unsigned long int)UNROLL_SIZE);
    timings["fabs"] = t;

    fclose(fh);

    fprintf(outfile, "\n");
    fprintf(outfile, "           ------------------------------ \n");
    fprintf(outfile, "           ======> MATH BENCHMARKS <===== \n");
    fprintf(outfile, "           ------------------------------ \n");
    fprintf(outfile, "\n");

    fprintf(outfile, "  Parameters:\n");
    fprintf(outfile, "   -Minimum runtime (per operation): %14.10f [s].\n", min_time);
    fprintf(outfile, "\n");
    fprintf(outfile, "  Notes:\n");
    fprintf(outfile, "   -All operations are for doubles, loops unrolled.\n");
    fprintf(outfile, "   -exp, log, pow, and - are probably optimized out, and unreliable.\n");
    fprintf(outfile, "\n");
    fprintf(outfile, "Operation  ");
    fprintf(outfile, "  %9s  %9s  %9s", "Time", "FLOPS", "Adds");
    fprintf(outfile, "\n");

    double add_time = timings["+"];

    for (int s = 0; s < ops.size(); s++) {
        t = timings[ops[s]];
        fprintf(outfile, "%-11s", ops[s].c_str());
        fprintf(outfile, "    %9.3E  %9.3E  %9.3E", t, 1 / t, t / add_time);
        fprintf(outfile, "\n");
    }
    fprintf(outfile, "\n");
    fflush(outfile);
}
void benchmark_integrals(int max_am, double min_time)
{
    double T;
    unsigned long int rounds;
    double t;
    Timer* qq;

    //We'll contract one each for s, p, and d.
    std::pair<std::vector<std::string>, boost::shared_ptr<BasisSet> > bases = BasisSet::test_basis_set(max_am);
    std::vector<std::string> shell_names = bases.first;
    boost::shared_ptr<BasisSet> basis = bases.second;
    int max_shell = basis->nshell() / basis->molecule()->natom();
    boost::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();

    //Factories
    boost::shared_ptr<IntegralFactory> bbbb(new IntegralFactory(basis, basis, basis, basis));
    boost::shared_ptr<IntegralFactory> b0bb(new IntegralFactory(basis, zero, basis, basis));
    boost::shared_ptr<IntegralFactory> b0b0(new IntegralFactory(basis, zero, basis, zero));

    std::vector<std::string> int_types;
    int_types.push_back("2C Overlap");
    int_types.push_back("2C Kinetic");
    int_types.push_back("2C Potential");
    int_types.push_back("2C Dipole");
    int_types.push_back("2C Quadrupole");
    int_types.push_back("2C ERI");
    int_types.push_back("3C ERI");
    int_types.push_back("3C Overlap");
    int_types.push_back("4C ERI");

    std::vector<int> centers;
    for (int k = 0; k < 6; k++)
        centers.push_back(2);
    for (int k = 0; k < 2; k++)
        centers.push_back(3);
    for (int k = 0; k < 1; k++)
        centers.push_back(4);

    std::map<int, int> ncombinations;
    ncombinations[2] = max_shell * max_shell;
    ncombinations[3] = max_shell * max_shell * max_shell;
    ncombinations[4] = max_shell * max_shell * max_shell * max_shell;

    std::map<int, std::vector<std::string> > combinations;
    combinations[2].resize(ncombinations[2]);
    combinations[3].resize(ncombinations[3]);
    combinations[4].resize(ncombinations[4]);

    std::map<int, std::vector<int> > n_per_combination;
    n_per_combination[2].resize(ncombinations[2]);
    n_per_combination[3].resize(ncombinations[3]);
    n_per_combination[4].resize(ncombinations[4]);

    int index;
    for (int P = 0, index = 0; P < max_shell; P++) {
        for (int Q = 0; Q < max_shell; Q++, index++) {
            combinations[2][index] = "(" + shell_names[P] + "|" + shell_names[Q] + ")";
            n_per_combination[2][index] =  basis->shell(P).nfunction() * basis->shell(Q).nfunction();
        }
    }
    for (int P = 0, index = 0; P < max_shell; P++) {
        for (int Q = 0; Q < max_shell; Q++) {
            for (int R = 0; R < max_shell; R++, index++) {
                combinations[3][index] = "(" + shell_names[P] + "|" + shell_names[Q] + shell_names[R] + ")";
                n_per_combination[3][index] =  basis->shell(P).nfunction() * basis->shell(Q).nfunction() * basis->shell(R).nfunction();
            }
        }
    }
    for (int P = 0, index = 0; P < max_shell; P++) {
        for (int Q = 0; Q < max_shell; Q++) {
            for (int R = 0; R < max_shell; R++) {
               for (int S = 0; S < max_shell; S++, index++) {
                   combinations[4][index] = "(" + shell_names[P] + shell_names[Q] + "|" + shell_names[R] + shell_names[S] + ")";
                   n_per_combination[4][index] =  basis->shell(P).nfunction() * basis->shell(Q).nfunction() * basis->shell(R).nfunction() * basis->shell(S).nfunction();
                }
            }
        }
    }


    std::map<std::string, std::vector<double> > timings;
    for (int k = 0; k < int_types.size(); k++) {
        timings[int_types[k]].resize(ncombinations[centers[k]]);
    }

    std::string this_type;
    int this_ncenter;

    // Overlap Ints
    this_type = "2C Overlap";
    this_ncenter = 2;
    boost::shared_ptr<OneBodyAOInt> o2c(bbbb->ao_overlap());
    for (int P = 0, index = 0; P < max_shell; P++) {
        for (int Q = 0; Q < max_shell; Q++, index++) {
            T = 0.0;
            rounds = 0L;
            qq = new Timer();
            while (T < min_time) {
                o2c->compute_shell(P, Q + max_shell);
                T = qq->get();
                rounds++;
            }
            delete qq;
            t = T / (double) (rounds * n_per_combination[this_ncenter][index]);
            timings[this_type][index] = t;
        }
    }

    // Kinetic Ints
    this_type = "2C Kinetic";
    this_ncenter = 2;
    boost::shared_ptr<OneBodyAOInt> k2c(bbbb->ao_kinetic());
    for (int P = 0, index = 0; P < max_shell; P++) {
        for (int Q = 0; Q < max_shell; Q++, index++) {
            T = 0.0;
            rounds = 0L;
            qq = new Timer();
            while (T < min_time) {
                k2c->compute_shell(P, Q + max_shell);
                T = qq->get();
                rounds++;
            }
            delete qq;
            t = T / (double) (rounds * n_per_combination[this_ncenter][index]);
            timings[this_type][index] = t;
        }
    }

    // Potential Ints
    this_type = "2C Potential";
    this_ncenter = 2;
    boost::shared_ptr<OneBodyAOInt> v2c(bbbb->ao_potential());
    for (int P = 0, index = 0; P < max_shell; P++) {
        for (int Q = 0; Q < max_shell; Q++, index++) {
            T = 0.0;
            rounds = 0L;
            qq = new Timer();
            while (T < min_time) {
                v2c->compute_shell(P, Q + max_shell);
                T = qq->get();
                rounds++;
            }
            delete qq;
            t = T / (double) (rounds * n_per_combination[this_ncenter][index]);
            timings[this_type][index] = t;
        }
    }

    // Dipole Ints
    this_type = "2C Dipole";
    this_ncenter = 2;
    boost::shared_ptr<OneBodyAOInt> d2c(bbbb->ao_dipole());
    for (int P = 0, index = 0; P < max_shell; P++) {
        for (int Q = 0; Q < max_shell; Q++, index++) {
            T = 0.0;
            rounds = 0L;
            qq = new Timer();
            while (T < min_time) {
                d2c->compute_shell(P, Q + max_shell);
                T = qq->get();
                rounds++;
            }
            delete qq;
            t = T / (double) (3L * rounds * n_per_combination[this_ncenter][index]);
            timings[this_type][index] = t;
        }
    }

    // Quadrupole Ints
    this_type = "2C Quadrupole";
    this_ncenter = 2;
    boost::shared_ptr<OneBodyAOInt> q2c(bbbb->ao_quadrupole());
    for (int P = 0, index = 0; P < max_shell; P++) {
        for (int Q = 0; Q < max_shell; Q++, index++) {
            T = 0.0;
            rounds = 0L;
            qq = new Timer();
            while (T < min_time) {
                q2c->compute_shell(P, Q + max_shell);
                T = qq->get();
                rounds++;
            }
            delete qq;
            t = T / (double) (6L * rounds * n_per_combination[this_ncenter][index]);
            timings[this_type][index] = t;
        }
    }

    // 2C ERIs
    this_type = "2C ERI";
    this_ncenter = 2;
    boost::shared_ptr<TwoBodyAOInt> e2c(b0b0->eri());
    for (int P = 0, index = 0; P < max_shell; P++) {
        for (int Q = 0; Q < max_shell; Q++, index++) {
            T = 0.0;
            rounds = 0L;
            qq = new Timer();
            while (T < min_time) {
                e2c->compute_shell(P, 0, Q + max_shell, 0);
                T = qq->get();
                rounds++;
            }
            delete qq;
            t = T / (double) (rounds * n_per_combination[this_ncenter][index]);
            timings[this_type][index] = t;
        }
    }
    /**
    // Poisson Ints
    this_type = "2C Poisson";
    this_ncenter = 2;
    boost::shared_ptr<OneBodyAOInt> p2c(bbbb->poisson_overlap());
    for (int P = 0, index = 0; P < max_shell; P++) {
        for (int Q = 0; Q < max_shell; Q++, index++) {
            T = 0.0;
            rounds = 0L;
            qq = new Timer();
            while (T < min_time) {
                p2c->compute_shell(P, Q + max_shell);
                T = qq->get();
                rounds++;
            }
            delete qq;
            t = T / (double) (rounds * n_per_combination[this_ncenter][index]);
            timings[this_type][index] = t;
        }
    }**/

    // 3C ERIs
    this_type = "3C ERI";
    this_ncenter = 3;
    boost::shared_ptr<TwoBodyAOInt> e3c(b0bb->eri());
    for (int P = 0, index = 0; P < max_shell; P++) {
        for (int Q = 0; Q < max_shell; Q++) {
            for (int R = 0; R < max_shell; R++, index++) {
                T = 0.0;
                rounds = 0L;
                qq = new Timer();
                while (T < min_time) {
                    e3c->compute_shell(P, 0, Q + max_shell, R + max_shell*2);
                    T = qq->get();
                    rounds++;
                }
                delete qq;
                t = T / (double) (rounds * n_per_combination[this_ncenter][index]);
                timings[this_type][index] = t;
            }
        }
    }

    // 3C Overlap
    this_type = "3C Overlap";
    this_ncenter = 3;
    boost::shared_ptr<ThreeCenterOverlapInt> o3c(bbbb->overlap_3c());
    for (int P = 0, index = 0; P < max_shell; P++) {
        for (int Q = 0; Q < max_shell; Q++) {
            for (int R = 0; R < max_shell; R++, index++) {
                T = 0.0;
                rounds = 0L;
                qq = new Timer();
                while (T < min_time) {
                    o3c->compute_shell(P, Q + max_shell, R + max_shell*2);
                    T = qq->get();
                    rounds++;
                }
                delete qq;
                t = T / (double) (rounds * n_per_combination[this_ncenter][index]);
                timings[this_type][index] = t;
            }
        }
    }

    // 4C ERI
    this_type = "4C ERI";
    this_ncenter = 4;
    boost::shared_ptr<TwoBodyAOInt> e4c(bbbb->eri());
    for (int P = 0, index = 0; P < max_shell; P++) {
        for (int Q = 0; Q < max_shell; Q++) {
            for (int R = 0; R < max_shell; R++) {
                for (int S = 0; S < max_shell; S++, index++) {
                    T = 0.0;
                    rounds = 0L;
                    qq = new Timer();
                    while (T < min_time) {
                        e4c->compute_shell(P, Q + max_shell, R + max_shell*2, S + max_shell*3);
                        T = qq->get();
                        rounds++;
                    }
                    delete qq;
                    t = T / (double) (rounds * n_per_combination[this_ncenter][index]);
                    timings[this_type][index] = t;
                }
            }
        }
    }


    fprintf(outfile, "\n");
    fprintf(outfile, "                              ----------------------------------- \n");
    fprintf(outfile, "                              ======> INTEGRALS BENCHMARKS <===== \n");
    fprintf(outfile, "                              ----------------------------------- \n");
    fprintf(outfile, "\n");

    fprintf(outfile, "  Parameters:\n");
    fprintf(outfile, "   -Maximum angular momentum %s\n", shell_names[shell_names.size() - 1].c_str());
    fprintf(outfile, "   -Minimum runtime (per integral, per combination): %14.10f [s].\n", min_time);
    fprintf(outfile, "\n");

    fprintf(outfile, "  Notes:\n");
    fprintf(outfile, "    -All integrals are computed from different centers.\n");
    fprintf(outfile, "    -(s,p,d,f,g,h,i) are single-primitive shells.\n");
    fprintf(outfile, "    -(S,P,D) are [10, 6, 2]-primitive shells, respectively.\n");
    fprintf(outfile, "    -Only up to i functions are currently supported.\n");
    fprintf(outfile, "    -All shells use Spherical Harmonics.\n");
    fprintf(outfile, "    -Timings are reported per double produced (function combination and possibly direction).\n");
    fprintf(outfile, "     Therefore, the cost for a (p|p) overlap shell would be 9x the value reported, while the\n");
    fprintf(outfile, "     cost for a (p|p) dipole shell would be 27x the value reported.\n");
    fprintf(outfile, "\n");


    fprintf(outfile, "Test Basis Set:\n");
    basis->print_by_level(outfile, 3);

    for (int op = 0; op < int_types.size(); op++) {
        this_type = int_types[op];
        this_ncenter = centers[op];
        fprintf(outfile, "  Integral Type: %s\n\n", this_type.c_str());

        // Time per element
        fprintf(outfile, "Combination%11s  %11s    %11s\n", "T [s]", "1/T [Hz]", "1/T [GiB/s]");
        for (index = 0; index < combinations[this_ncenter].size(); index++) {
            t = timings[this_type][index];
            fprintf(outfile, "%-10s    %9.3E    %9.3E    %9.3E\n", combinations[this_ncenter][index].c_str(), t, 1.0 / t, 8.0E-9 / t);
        }
        fprintf(outfile,"\n");
    }
    fflush(outfile);

}

}
