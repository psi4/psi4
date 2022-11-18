#include "numinthelper.h"
#include "psi4/libfock/cubature.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libfock/points.h"
#include "psi4/libfock/dft_integrators.h"
#include "psi4/libfock/sap.h"
#include "psi4/libpsi4util/process.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace psi {

NumIntHelper::NumIntHelper(std::shared_ptr<DFTGrid> numint_grid) : print_(0), nthread_(1), numint_grid_(numint_grid) {
#ifdef _OPENMP
    nthread_ = Process::environment.get_n_threads();
#endif

    Options& options_ = Process::environment.options;
    print_ = options_.get_int("PRINT");
}

namespace {

struct GridFunction {
    int n_output;
    const std::vector<SharedMatrix>& grid_data;
    GridFunction(const std::vector<SharedMatrix>& grid_data)
        : n_output(grid_data[0]->rowspi()[0]), grid_data(grid_data) {}

    void process_block(size_t iblock, const std::shared_ptr<BlockOPoints>& block,
                       const std::shared_ptr<PointFunctions>& worker, bool sum_over_atoms, SharedMatrix& output) const {
        auto rho_block = worker->point_values()["RHO_A"];
        double** data_block = grid_data[iblock]->pointer();
        size_t n_data = static_cast<size_t>(n_output);

        for (size_t idata = 0; idata < n_data; ++idata) {
            for (int ig = 0; ig < block->npoints(); ++ig) {
                const size_t icomponent = sum_over_atoms ? 0 : block->parent_atom();
                output->pointer()[icomponent][idata] -= data_block[idata][ig] * rho_block->get(ig) * block->w()[ig];
            }
        }
    }
};
}  // namespace

template <typename Integrand>
SharedMatrix evaluate_density_integral(const std::shared_ptr<DFTGrid>& numint_grid, int nthread, bool sum_over_atoms,
                                       Integrand integrand, const SharedMatrix& D) {
    const int max_points = numint_grid->max_points();
    const int max_functions = numint_grid->max_functions();
    const std::shared_ptr<BasisSet> basisset = numint_grid->primary();
    const int n_atoms = basisset->molecule()->natom();
    const int n_components = sum_over_atoms ? 1 : n_atoms;

    std::vector<SharedMatrix> integral_local;
    std::vector<std::shared_ptr<PointFunctions>> numint_point_workers;
    for (size_t i = 0; i < nthread; i++) {
        integral_local.push_back(std::make_shared<Matrix>("Integral Temp", n_components, integrand.n_output));

        // Use RKS function to have total density at grid points available
        auto point_tmp = std::make_shared<RKSFunctions>(basisset, max_points, max_functions);
        point_tmp->set_ansatz(0);    // LDA is enough
        point_tmp->set_pointers(D);  // TODO What is needed for unrestricted?
        numint_point_workers.push_back(point_tmp);
    }

    int rank = 0;
// Traverse the blocks of points
#pragma omp parallel for private(rank) schedule(guided) num_threads(nthread)
    for (size_t Q = 0; Q < numint_grid->blocks().size(); Q++) {
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        // Extract per-thread info
        std::shared_ptr<BlockOPoints> block = numint_grid->blocks()[Q];
        std::shared_ptr<PointFunctions> worker = numint_point_workers[rank];
        // worker->set_pointers(D);   // TODO Check with parallelised version ... does not seem to be needed

        // Compute Rho on grid.
        parallel_timer_on("Density evaluation", rank);
        worker->compute_points(block, false);
        parallel_timer_off("Density evaluation", rank);

        parallel_timer_on("Process block", rank);
        integrand.process_block(Q, block, worker, sum_over_atoms, integral_local[rank]);
        parallel_timer_off("Process block", rank);
    }  // grid blocks

    // Accumulate results in integral_local[0]
    for (size_t rank = 1; rank < nthread; rank++) {
        integral_local[0]->add(integral_local[rank]);
    }
    integral_local[0]->scale(0.5);  // TODO Why is this needed?
    return integral_local[0];
}

SharedVector NumIntHelper::density_integral(const std::vector<SharedMatrix>& grid_data, const SharedMatrix& D) const {
    if (numint_grid_->blocks().size() != grid_data.size()) {
        throw PSIEXCEPTION("Mismatch of grid data size and DFT integration blocks.");
    }
    timer_on("NumIntHelper: density_integral");
    const bool sum_over_atoms = true;
    auto Vout = evaluate_density_integral(numint_grid_, nthread_, sum_over_atoms, GridFunction(grid_data), D);
    timer_off("NumIntHelper: density_integral");
    return Vout->get_row(0, 0);  // get first (and only) row
}

SharedMatrix NumIntHelper::dd_density_integral(const std::vector<SharedMatrix>& grid_data,
                                               const SharedMatrix& D) const {
    if (numint_grid_->blocks().size() != grid_data.size()) {
        throw PSIEXCEPTION("Mismatch of grid data size and DFT integration blocks.");
    }
    timer_on("NumIntHelper: dd_density_integral");
    const bool sum_over_atoms = false;
    auto Vout = evaluate_density_integral(numint_grid_, nthread_, sum_over_atoms, GridFunction(grid_data), D);
    timer_off("NumIntHelper: dd_density_integral");
    return Vout;
}

SharedMatrix NumIntHelper::potential_integral(const std::vector<SharedVector>& grid_data) const {
    if (numint_grid_->blocks().size() != grid_data.size()) {
        throw PSIEXCEPTION("Mismatch of grid data size and DFT integration blocks.");
    }
    timer_on("NumIntHelper: potential_integral");
    const int max_points = numint_grid_->max_points();
    const int max_functions = numint_grid_->max_functions();
    std::shared_ptr<BasisSet> basisset = numint_grid_->primary();

    std::vector<std::shared_ptr<PointFunctions>> numint_point_workers;
    std::vector<SharedMatrix> V_local;
    for (size_t i = 0; i < nthread_; i++) {
        V_local.push_back(std::make_shared<Matrix>("V Temp", max_functions, max_functions));
        auto point_tmp = std::make_shared<SAPFunctions>(basisset, max_points, max_functions);
        point_tmp->set_ansatz(0);
        numint_point_workers.push_back(point_tmp);
    }

    auto V_AO = std::make_shared<Matrix>("V(PC) operator", basisset->nbf(), basisset->nbf());
    double** Vp = V_AO->pointer();

    int rank = 0;
    // Traverse the blocks of points
#pragma omp parallel for private(rank) schedule(guided) num_threads(nthread_)
    for (size_t Q = 0; Q < numint_grid_->blocks().size(); Q++) {
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif

        // Extract per-thread info
        std::shared_ptr<BlockOPoints> block = numint_grid_->blocks()[Q];
        std::shared_ptr<PointFunctions> worker = numint_point_workers[rank];

        // Compute data on grid points
        parallel_timer_on("Properties", rank);
        worker->compute_points(block, false);
        parallel_timer_off("Properties", rank);

        parallel_timer_on("finish operator", rank);
        // PHI * potential * PHI.
        dft_integrators::sap_integrator(block, grid_data[Q], worker, V_local[rank]);

        // => Unpacking <= //
        double** V2p = V_local[rank]->pointer();
        const std::vector<int>& function_map = block->functions_local_to_global();
        int nlocal = function_map.size();

        for (int ml = 0; ml < nlocal; ml++) {
            int mg = function_map[ml];
            for (int nl = 0; nl < ml; nl++) {
                int ng = function_map[nl];
#pragma omp atomic update
                Vp[mg][ng] += V2p[ml][nl];
#pragma omp atomic update
                Vp[ng][mg] += V2p[ml][nl];
            }
#pragma omp atomic update
            Vp[mg][mg] += V2p[ml][ml];
        }
        parallel_timer_off("finish operator", rank);
    }  // grid blocks
    timer_off("NumIntHelper: potential_integral");
    return V_AO;
}
}  // namespace psi
