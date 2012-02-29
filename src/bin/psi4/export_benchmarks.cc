#include <boost/python.hpp>
#include <libmints/benchmark.h>

using namespace boost::python;

void export_benchmarks()
{
    def("benchmark_blas1",     &psi::benchmark_blas1, "docstring");
    def("benchmark_blas2",     &psi::benchmark_blas2, "docstring");
    def("benchmark_blas3",     &psi::benchmark_blas3, "docstring");
    def("benchmark_disk",      &psi::benchmark_disk, "docstring");
    def("benchmark_math",      &psi::benchmark_math, "docstring");
    def("benchmark_integrals", &psi::benchmark_integrals, "docstring");
}
