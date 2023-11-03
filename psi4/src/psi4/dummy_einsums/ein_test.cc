/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include <cmath>
#include <cstdlib>

#include "einsums.hpp"
#include "range/v3/algorithm/for_each.hpp"
#include "range/v3/view/zip.hpp"

#include "psi4/libmints/wavefunction.h"

namespace psi {
namespace dummy_einsums {

SharedWavefunction dummy_einsums(SharedWavefunction ref_wfn, Options& options) {

    using namespace einsums;
    using namespace einsums::tensor_algebra;
    using namespace einsums::tensor_algebra::index;

    timer::initialize();
    blas::initialize();

    // Disable HDF5 diagnostic reporting.
    H5Eset_auto(0, nullptr, nullptr);

    // Create a file to hold the data from the DiskTensor tests.
    //einsums::state::data = h5::create("Data.h5", H5F_ACC_TRUNC);
    // TODO need to delete file

    auto [_t, _w] = polynomial::laguerre::gauss_laguerre(40);
    println(_t);
    println(_w);

    auto weights = create_tensor_like(_w);
    {
        using namespace element_operations::new_tensor;
        einsum(Indices{i}, &_w, Indices{i}, _w, Indices{i}, exp(_t));
    }
    println(_w);

    ranges::for_each(ranges::views::zip(_t.vector_data(), _w.vector_data()), [](auto &&v) {
        auto t = std::get<0>(v);
        auto w = std::get<1>(v);

        println("test : {:14.10f} {:14.10f}", t, w);
    });

    // BEGIN: FFT tests

    // auto result = fft::fftfreq(8, 0.1);
    // println(result);

    timer::report();
    blas::finalize();
    timer::finalize();

    return ref_wfn;
}
}  // namespace dummy_einsums
}  // namespace psi
