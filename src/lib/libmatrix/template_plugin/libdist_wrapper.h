namespace psi { namespace libmatrix {

    struct libdist_globals {
        static void initialize(int argc, char** argv) {
            // Nothing to do. Psi4 already prepares the environment.
        }
    };

#if defined(HAVE_MADNESS)

    struct libdist_matrix_wrapper {
        libdist_matrix_wrapper(const std::string& name, const Dimension& m, const Dimension& n) : matrix_()
        {
            std::vector<int> dims(2);
            dims[0] = 2; dims[1] = 2;
            process_grid<2> pgrid(dims);

            matrix_.initialize(pgrid, m.sum(), n.sum(), 64, name);
        }

        void print() const {
            matrix_.print();
        }

        void fill(double val) {
            matrix_.fill(val);
        }

        // Cause problems if the someone tries to use something other than libmints_matrix_wrapper
        template <typename R>
        void add(const R& rhs) {
            throw std::logic_error("Don't know how to add different matrix types together.");
        }

        void add(const libdist_matrix_wrapper& rhs) {
            matrix_.add(rhs.matrix_);
        }

        // Is it not unreasonable to require all matrix types to know how to work with the serial matrix object.
        void add(const libmints_matrix_wrapper& rhs) {
            // TODO: implement
        }
    private:
        Distributed_Matrix matrix_;
    };
#else
    // We are not configured with Madness, derive from is_not_available to 
    // cause compiler error when libdist_matrix_wrapper is instantiated in
    // the code.
    UNAVAILABLE_MATRIX(libdist_matrix_wrapper);
#endif

}} // end namespaces

