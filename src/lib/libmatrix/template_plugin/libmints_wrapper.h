#include <string>
#include <libmints/vector.h>
#include <libmints/matrix.h>
#include <stdexcept>

namespace psi { namespace libmatrix {

    // Refer to libelemental for example of global data and usage.
    struct libmints_globals {
        static std::string interface_name;

        static void initialize(int argc, char** argv) {
            // Nothing to do.
        }
    };

    struct libmints_vector_wrapper {
        libmints_vector_wrapper(const std::string& name, const Dimension& m) :
            vector_(name, m)
        { }

        void print() const {
            vector_.print(stdout);
        }

        friend struct libmints_matrix_wrapper;
    private:
        Vector vector_;
    };

    struct libmints_matrix_wrapper {
        libmints_matrix_wrapper(const std::string& name, const Dimension& m, const Dimension& n) :
            matrix_(name, m, n)
        { }

        void print() const {
            matrix_.print(stdout);
        }

        void fill(double val) {
            matrix_.set(val);
        }
        
        void gemm(bool ta, bool tb, double alpha, 
                  const libmints_matrix_wrapper& A,
                  const libmints_matrix_wrapper& B,
                  double beta)
        {
            matrix_.gemm(ta, tb, alpha, A.matrix_, B.matrix_, beta);
        }

        // Must assume this is destroyed
        void diagonalize(libmints_matrix_wrapper& X, libmints_vector_wrapper& w) {
            matrix_.diagonalize(X.matrix_, w.vector_);
        }

        // Cause problems if the someone tries to use something other than libmints_matrix_wrapper
        template <typename R>
        void add(const R& rhs) {
            throw std::logic_error("Don't know how to add different matrix types together.");
        }

        void add(const libmints_matrix_wrapper& rhs) {
            matrix_.add(rhs.matrix_);
        }

        double operator()(int h, int m, int n) const {
            return matrix_.get(h, m, n);
        }

    private:
        Matrix matrix_;
    };

}} // end namespaces

