namespace psi { namespace libmatrix {

    // Refer to libelemental for example of global data and usage.
    struct libmints_globals {
        static void initialize(int argc, char** argv) {
            // Nothing to do.
        }
    };

    struct libmints_matrix_wrapper {
        libmints_matrix_wrapper(const std::string& name, const Dimension& m, const Dimension& n) :
            matrix_(name, m, n)
        { }

        void print() const {
            matrix_.print();
        }

        void fill(double val) {
            matrix_.set(val);
        }

        // Cause problems if the someone tries to use something other than libmints_matrix_wrapper
        template <typename R>
        void add(const R& rhs) {
            throw std::logic_error("Don't know how to add different matrix types together.");
        }

        void add(const libmints_matrix_wrapper& rhs) {
            matrix_.add(rhs.matrix_);
        }

    private:
        Matrix matrix_;
    };

}} // end namespaces

