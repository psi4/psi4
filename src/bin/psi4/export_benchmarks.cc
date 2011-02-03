#include <boost/python.hpp>
#include <libmints/benchmark.h>

using namespace boost::python;

void export_benchmarks()
{
    def("benchmark_blas1",     &psi::benchmark_blas1);
    def("benchmark_blas2",     &psi::benchmark_blas2);
    def("benchmark_blas3",     &psi::benchmark_blas3);
    def("benchmark_disk",      &psi::benchmark_disk);
    def("benchmark_math",      &psi::benchmark_math);
    def("benchmark_integrals", &psi::benchmark_integrals);
}
