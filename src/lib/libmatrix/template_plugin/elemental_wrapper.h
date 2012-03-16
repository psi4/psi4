namespace psi { namespace libmatrix {

#if defined(HAVE_ELEMENTAL)
    
    #include <elemental.hpp>
    #include <elemental/basic_internal.hpp>

    struct libelemental {
        static bool initialized = false;
        static mpi::Comm mpi_comm;
        static int rank = 0;
        
        static void initialize(int argc, char** argv) {
            elemental::Init(argc, argv);
            mpi_comm = mpi::COMM_WORLD;
            rank = mpi::commRank(mpi_comm);
        }
    }

    // !! I have elemental installed on my laptop but Psi4 doesn't know
    // !! how to configure with it. This code is untested.
    struct libelemental_matrix_wrapper {
        libelemental_matrix_wrapper(const std::string& name, const Dimension& m, const Dimension& n) :
            name_(name), height_(m.sum()), width_(n.sum())
        { 
            elemental::Grid g(libelemental::mpi_comm);
            matrix_ = Distmatrix<double, MC, MR>(height_, width_, g);
        }

        void print() const {
            matrix_.Print(name_);
        }

        void fill(double val) {
            // elemental doesn't provide a generic fill
            for (int m=0; m<height_; ++m)
                for (int n=0; n<width_; ++n)
                    matrix_.Set(m, n, val);
        }

        // Cause problems if the someone tries to use something other than libmints_matrix_wrapper
        template <typename R>
        void add(const R& rhs) {
            throw std::logic_error("Don't know how to add different matrix types together.");
        }

        void add(const libelemental_matrix_wrapper& rhs) {
            elemental::basic::Axpy(1.0, rhs.matrix_, matrix_);
        }

        // Is it not unreasonable to require all matrix types to know how to work with the serial matrix object.
        void add(const libmints_matrix_wrapper& rhs) {
            // TODO: implement
        }
    private:
        int height_, width_;
        std::string name_;
        DistMatrix<double, MC, MR> matrix_;
    };
#else
    // We are not configured with Elemental, derive from  is_not_available to
    // cause compiler errors when libelemental_matrix_wrapper is instantiated in
    // the code.
    UNAVAILABLE_MATRIX(libelemental_matrix_wrapper);
#endif

}} // end namespaces

