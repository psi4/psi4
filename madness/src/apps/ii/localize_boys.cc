#include <world/world.h>
#include <mra/mra.h>
#include <utility>
#include <ctime>
#include <cmath>
#include <tensor/tensor.h>
#include <ii/systolic.h>

using namespace madness;

/// A MADNESS functor to compute either x, y, or z
class DipoleFunctor : public FunctionFunctorInterface<double,3> {
private:
    const int axis;
public:
    DipoleFunctor(int axis) : axis(axis) {}
    double operator()(const Vector<double,3>& x) const {
        return x[axis];
    }
};

/// This is an adapter function. returns a distributed and localized set
/// \param mo molecules
/// \param set basis set
template <typename T>
DistributedMatrix<T> plocalize_boys( World & world, const std::vector< Function<T,3> > & mo, const std::vector<int> & set, const double vtol,
        const double thresh = 1e-9, const double thetamax = 0.5, const bool randomize = true )
{
    //START_TIMER(world);
    const bool doprint = false; /// do print
    long nmo = mo.size(); /// number of molecule

    /** make big matrix such like ...
        +---+-------+-------+-------+
        | I | dip_x | dip_y | dip_z |
        +---+-------+-------+-------+
    */
    DistributedMatrix<T> dip[3] = {
        column_distributed_matrix<T>(world, nmo, nmo),
        column_distributed_matrix<T>(world, nmo, nmo),
        column_distributed_matrix<T>(world, nmo, nmo)
    };

    /// for all axis, make dipole function
    for(int axis=0;axis < 3;axis++){
        Function<T,3> fdip = FunctionFactory<T,3>(world).functor(std::shared_ptr< FunctionFunctorInterface<T,3> >(new DipoleFunctor(axis))).initial_level(4);
        dip[axis].copyin( matrix_inner(world, mo, mul_sparse(world, fdip, mo, vtol), true) );
    }

    DistributedMatrix<T> id = idMatrix(dip[0]);
    DistributedMatrix<double> M = concatenate_rows(id, dip[0], dip[1], dip[2]);

    LocalizeBoys<T> lb(M, set, nmo, 3335, thresh, thetamax, randomize );
    lb.solve();

    /// image of usage
    return lb.get_U();
}

int main(int argc, char** argv) {
    initialize(argc, argv);

    madness::World world(MPI::COMM_WORLD);

    redirectio(world); /* redirect a results to file "log.<number of processor>" */

    try {
        print("Test of localize_boys.cc\n");
        print("label size time eig_val");
        for (int64_t n=10; n>1; n-=2) {

            /// image of usage
            /*
            Tensor<double>U = plocalize_boys(world, mo, set);
            double t = cpu_time();
            */

            //print("result:", n, cpu_time()-t, M.localized_set());
        }
    }
    catch (const MPI::Exception& e) {
        //        print(e);
        error("caught an MPI exception");
    }
    catch (const madness::MadnessException& e) {
        print(e);
        error("caught a MADNESS exception");
    }
    catch (const madness::TensorException& e) {
        print(e);
        error("caught a Tensor exception");
    }
    catch (const char* s) {
        print(s);
        error("caught a c-string exception");
    }
    catch (char* s) {
        print(s);
        error("caught a c-string exception");
    }
    catch (const std::string& s) {
        print(s);
        error("caught a string (class) exception");
    }
    catch (const std::exception& e) {
        print(e.what());
        error("caught an STL exception");
    }
    catch (...) {
        error("caught unhandled exception");
    }

    finalize();
}

