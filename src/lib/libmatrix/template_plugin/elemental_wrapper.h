#if defined(HAVE_ELEMENTAL)
#   include <string>
#   include <libmints/dimension.h>
#   include <elemental.hpp>
#   include <elemental/basic.hpp>
#endif

namespace psi { namespace libmatrix {

    struct libmints_matrix_wrapper;
    struct libelemental_vector_wrapper;
    struct libelemental_matrix_wrapper;

#if defined(HAVE_ELEMENTAL)
    
    struct libelemental_globals {
        static std::string interface_name;
        static elem::mpi::Comm mpi_comm;
        static int rank;
        
        static void initialize(int argc, char** argv) {
            elem::Initialize(argc, argv);
            mpi_comm = elem::mpi::COMM_WORLD;
            rank = elem::mpi::CommRank(mpi_comm);
        }
    };

    struct libelemental_vector_wrapper {
        libelemental_vector_wrapper(const std::string& name, const Dimension& m) :
            name_(name), length_(m.sum())
        {
            elem::Grid g(libelemental_globals::mpi_comm);
            vector_ = elem::DistMatrix<double, elem::VR, elem::STAR>(true, length_, g);
        }

        void print() const {
            vector_.Print(name_);
        }

        friend struct libelemental_matrix_wrapper;
    private:
        int length_;
        std::string name_;
        elem::DistMatrix<double, elem::VR, elem::STAR> vector_;
    };

    // !! I have elemental installed on my laptop but Psi4 doesn't know
    // !! how to configure with it. This code is untested.
    struct libelemental_matrix_wrapper {
        libelemental_matrix_wrapper(const std::string& name, const Dimension& m, const Dimension& n) :
            name_(name), height_(m.sum()), width_(n.sum())
        { 
            elem::Grid g(libelemental_globals::mpi_comm);
            matrix_ = elem::DistMatrix<double, elem::MC, elem::MR>(height_, width_, g);
        }

        void print() const {
            matrix_.Print(name_);
        }

        /// Performs distributed matrix fill.
        void fill(double val) {
            const int colShift    = matrix_.ColShift(); // first row we own
            const int rowShift    = matrix_.RowShift(); // first col we own
            const int colStride   = matrix_.ColStride();
            const int rowStride   = matrix_.RowStride();
            const int localHeight = matrix_.LocalHeight();
            const int localWidth  = matrix_.LocalWidth();
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                for( int iLocal=0; iLocal<localHeight; ++iLocal )
                {
                    // Our process owns the rows colShift:colStride:n,
                    //           and the columns rowShift:rowStride:n
                    const int i = colShift + iLocal*colStride;
                    const int j = rowShift + jLocal*rowStride;
                    matrix_.SetLocalEntry( iLocal, jLocal, val);
                }
            }
        }

        void gemm(bool ta, bool tb, double alpha, 
                  const libelemental_matrix_wrapper& A,
                  const libelemental_matrix_wrapper& B,
                  double beta)
        {
            elem::Orientation opA = ta == true ? elem::TRANSPOSE : elem::NORMAL;
            elem::Orientation opB = tb == true ? elem::TRANSPOSE : elem::NORMAL;

            elem::Gemm(opA, opB, alpha, A.matrix_, B.matrix_, beta, matrix_);
        }

        void diagonalize(libelemental_matrix_wrapper& X, libelemental_vector_wrapper& w) {
            // TODO: w needs to be a vector allocated using: DistMatrix<R,VR,STAR> w( g );
            elem::HermitianEig( elem::LOWER, matrix_, w.vector_, X.matrix_);

            // Sort eigenvalues into ascending order (only functionality provided by elemental)
            elem::SortEig( w.vector_, X.matrix_ );
        }

        // Cause problems if the someone tries to use something other than libmints_matrix_wrapper
        template <typename R>
        void add(const R& rhs) {
            throw std::logic_error("Don't know how to add different matrix types together.");
        }

        void add(const libmints_matrix_wrapper& rhs) {
            const int colShift    = matrix_.ColShift(); // first row we own
            const int rowShift    = matrix_.RowShift(); // first col we own
            const int colStride   = matrix_.ColStride();
            const int rowStride   = matrix_.RowStride();
            const int localHeight = matrix_.LocalHeight();
            const int localWidth  = matrix_.LocalWidth();
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                for( int iLocal=0; iLocal<localHeight; ++iLocal )
                {
                    // Our process owns the rows colShift:colStride:n,
                    //           and the columns rowShift:rowStride:n
                    const int i = colShift + iLocal*colStride;
                    const int j = rowShift + jLocal*rowStride;
                    matrix_.SetLocalEntry(iLocal, jLocal, rhs(0, iLocal, jLocal));  // elemental version doesn't understand symmetry yet
                }
            }
        }

        void add(const libelemental_matrix_wrapper& rhs) {
            elem::Axpy(1.0, rhs.matrix_, matrix_);
        }

    private:
        int height_, width_;
        std::string name_;
        elem::DistMatrix<double, elem::MC, elem::MR> matrix_;
    };

#else
    // We are not configured with Elemental, derive from  is_not_available to
    // cause compiler errors when libelemental_matrix_wrapper is instantiated in
    // the code.
    UNAVAILABLE_MATRIX(libelemental_matrix_wrapper);
#endif

}} // end namespaces

