/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#include <boost/python.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <libmints/mints.h>
#include <libmints/twobody.h>
#include <libmints/integralparameters.h>
#include <libmints/orbitalspace.h>
#include <libmints/view.h>
#include <libmints/pybuffer.h>
#include <libmints/local.h>
#include <libmints/vector3.h>
#include <lib3index/3index.h>
#include <libscf_solver/hf.h>
#include <libscf_solver/rhf.h>

#include <libfock/jk.h>

#include <string>

using namespace boost;
using namespace boost::python;
using namespace psi;

dict matrix_array_interface(SharedMatrix mat, int irrep){
    dict rv;
    int rows = mat->rowspi(irrep);
    int cols = mat->colspi(irrep);
    rv["shape"] = boost::python::make_tuple(rows, cols);
    rv["data"] = boost::python::make_tuple((long)mat->get_pointer(irrep), true);
    std::string typestr = is_big_endian() ? ">" : "<";
    {
        std::stringstream sstr;
        sstr << (int)sizeof(double);
        typestr += "f" + sstr.str();
    }
    rv["typestr"] = typestr;
    return rv;
}

dict matrix_array_interface_c1(SharedMatrix mat){
    if(mat->nirrep() != 1){
        throw PSIEXCEPTION("Pointer export of multiple irrep matrices not yet implemented.");
    }
    return matrix_array_interface(mat, 0);
}

dict vector_array_interface(SharedVector vec, int irrep){
    dict rv;
    int elements = vec->dim(irrep);
    rv["shape"] = boost::python::make_tuple(elements);
    rv["data"] = boost::python::make_tuple((long)vec->pointer(irrep), true);
    std::string typestr = is_big_endian() ? ">" : "<";
    {
        std::stringstream sstr;
        sstr << (int)sizeof(double);
        typestr += "f" + sstr.str();
    }
    rv["typestr"] = typestr;
    return rv;
}

dict vector_array_interface_c1(SharedVector vec){
    if(vec->nirrep() != 1){
        throw PSIEXCEPTION("Pointer export of multiple irrep vectorss not yet implemented.");
    }
    return vector_array_interface(vec, 0);
}

boost::shared_ptr<Vector> py_nuclear_dipole(shared_ptr<Molecule> mol)
{
    //SharedMolecule mol = Process::environment.molecule();
    return DipoleInt::nuclear_contribution(mol, Vector3(0, 0, 0));
}

boost::shared_ptr<MatrixFactory> get_matrix_factory()
{
    // We need a valid molecule with a valid point group to create a matrix factory.
    boost::shared_ptr<Molecule> molecule = Process::environment.molecule();
    if (!molecule) {
        outfile->Printf( "  Active molecule not set!");
        throw PSIEXCEPTION("Active molecule not set!");
    }
    if (!molecule->point_group()) {
        outfile->Printf( "  Active molecule does not have point group set!");
        throw PSIEXCEPTION("Active molecule does not have point group set!");
    }

    // Read in the basis set
    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser);
    boost::shared_ptr<BasisSet> basis = BasisSet::construct(parser, molecule, "BASIS");
    boost::shared_ptr<IntegralFactory> fact(new IntegralFactory(basis, basis, basis, basis));
    boost::shared_ptr<SOBasisSet> sobasis(new SOBasisSet(basis, fact));
    const Dimension& dim = sobasis->dimension();

    boost::shared_ptr<MatrixFactory> matfac(new MatrixFactory);
    matfac->init_with(dim, dim);

    return matfac;
}


BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(CanonicalOrthog, Matrix::canonical_orthogonalization, 1, 2);

/* IntegralFactory overloads */
/* Functions that return OneBodyAOInt objects */
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(ao_overlap_overloads, IntegralFactory::ao_overlap, 0, 1);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(so_overlap_overloads, IntegralFactory::so_overlap, 0, 1);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(ao_kinetic_overloads, IntegralFactory::ao_kinetic, 0, 1);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(ao_potential_overloads, IntegralFactory::ao_potential, 0, 1);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(ao_pseudospectral_overloads, IntegralFactory::ao_pseudospectral, 0, 1);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(ao_dipole_overloads, IntegralFactory::ao_dipole, 0, 1);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(ao_nabla_overloads, IntegralFactory::ao_nabla, 0, 1);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(ao_angular_momentum_overloads, IntegralFactory::ao_angular_momentum, 0, 1);
/* Functions that return TwoBodyAOInt objects */
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(eri_overloads, IntegralFactory::eri, 0, 2);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f12_overloads, IntegralFactory::f12, 1, 3);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f12g12_overloads, IntegralFactory::f12g12, 1, 3);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f12_squared_overloads, IntegralFactory::f12_squared, 1, 3);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(f12_double_commutator_overloads, IntegralFactory::f12_double_commutator, 1, 3);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(erf_eri_overloads, IntegralFactory::erf_eri, 1, 3);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(erf_complement_eri_overloads, IntegralFactory::erf_complement_eri, 1, 3);

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Ca_subset_overloads, Wavefunction::Ca_subset, 0, 2);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Cb_subset_overloads, Wavefunction::Cb_subset, 0, 2);

BOOST_PYTHON_FUNCTION_OVERLOADS(pyconstruct_orb_overloads, BasisSet::pyconstruct_orbital, 3, 4);
BOOST_PYTHON_FUNCTION_OVERLOADS(pyconstruct_aux_overloads, BasisSet::pyconstruct_auxiliary, 5, 6);

void export_mints()
{
    def("nuclear_dipole", py_nuclear_dipole, "docstring");

    // This is needed to wrap an STL vector into Boost.Python. Since the vector
    // is going to contain boost::shared_ptr's we MUST set the no_proxy flag to true
    // (as it is) to tell Boost.Python to not create a proxy class to handle
    // the vector's data type.
    class_<std::vector<SharedMatrix > >("matrix_vector", "docstring").
            def(vector_indexing_suite<std::vector<SharedMatrix >, true >());
    // Other vector types
    class_<std::vector<double> >("vector_of_doubles", "docstring").
            def(vector_indexing_suite<std::vector<double>, true >());

    // Use typedefs to explicitly tell Boost.Python which function in the class
    // to use. In most cases, you should not be making Python specific versions
    // of functions.

    // For example in Vector there are 2 versions of set: a (double*) version and a
    // (int, int, double) version. We create a typedef function pointer to tell
    // Boost.Python we only want the (int, int, double) version.
    typedef void (Vector::*vector_setitem_1)(int, double);
    typedef void (Vector::*vector_setitem_2)(int, int, double);
    typedef double (Vector::*vector_getitem_1)(int);
    typedef double (Vector::*vector_getitem_2)(int, int);
    typedef void (Vector::*vector_setitem_n)(const boost::python::tuple&, double);
    typedef double (Vector::*vector_getitem_n)(const boost::python::tuple&);

    class_<Dimension>("Dimension", "docstring").
            def(init<int>()).
            def(init<int, const std::string&>()).
            def("print_out",
                &Dimension::print,
                "docstring").
            def("init",
                &Dimension::init,
                "Re-initializes the dimension object").
            def("n", &Dimension::n,
                return_value_policy<copy_const_reference>(),
                "The order of the dimension").
            add_property("name",
                         make_function(&Dimension::name, return_value_policy<copy_const_reference>()),
                         &Dimension::set_name,
                         "The name of the dimension. Used in printing.").
            def("__getitem__", &Dimension::get, return_value_policy<copy_const_reference>(), "docstring").
            def("__setitem__", &Dimension::set, "docstring");

    class_<Vector, boost::shared_ptr<Vector> >( "Vector", "docstring").
            def(init<int>()).
            def("get", vector_getitem_1(&Vector::get), "docstring").
            def("get", vector_getitem_2(&Vector::get), "docstring").
            def("set", vector_setitem_1(&Vector::set), "docstring").
            def("set", vector_setitem_2(&Vector::set), "docstring").
            def("print_out", &Vector::print_out, "docstring").
            def("scale", &Vector::scale, "docstring").
            def("dim", &Vector::dim, "docstring").
            def("__getitem__", vector_getitem_1(&Vector::pyget), "docstring").
            def("__setitem__", vector_setitem_1(&Vector::pyset), "docstring").
            def("__getitem__", vector_getitem_n(&Vector::pyget), "docstring").
            def("__setitem__", vector_setitem_n(&Vector::pyset), "docstring").
            def("nirrep", &Vector::nirrep, "docstring").
            add_property("__array_interface__", vector_array_interface_c1, "docstring");

    typedef void  (IntVector::*int_vector_set)(int, int, int);
    class_<IntVector, boost::shared_ptr<IntVector> >( "IntVector", "docstring").
            def(init<int>()).
            def("get", &IntVector::get, "docstring").
            def("set", int_vector_set(&IntVector::set), "docstring").
            def("print_out", &IntVector::print_out, "docstring").
            def("dim", &IntVector::dim, "docstring").
            def("nirrep", &IntVector::nirrep, "docstring");

    enum_<diagonalize_order>("DiagonalizeOrder", "docstring")
            .value("Ascending", ascending)
            .value("Descending", descending)
            .export_values();

    enum_<Molecule::GeometryUnits>("GeometryUnits", "docstring")
            .value("Angstrom", Molecule::Angstrom)
            .value("Bohr", Molecule::Bohr)
            .export_values();


    class_<PyBuffer<double>, shared_ptr<PyBuffer<double> > >("DoublePyBuffer", "Buffer interface to NumPy arrays").
            add_property("__array_interface__", &PyBuffer<double>::array_interface, "docstring");

    typedef void   (Matrix::*matrix_multiply)(bool, bool, double, const SharedMatrix&, const SharedMatrix&, double);
    typedef void   (Matrix::*matrix_diagonalize)(SharedMatrix&, boost::shared_ptr<Vector>&, diagonalize_order);
    typedef void   (Matrix::*matrix_one)(const SharedMatrix&);
    typedef double (Matrix::*double_matrix_one)(const SharedMatrix&);
    typedef void   (Matrix::*matrix_two)(const SharedMatrix&, const SharedMatrix&);
    typedef void   (Matrix::*matrix_save)(const std::string&, bool, bool, bool);
    typedef void   (Matrix::*matrix_set4)(int, int, int, double);
    typedef void   (Matrix::*matrix_set3)(int, int, double);
    typedef double (Matrix::*matrix_get3)(const int&, const int&, const int&) const;
    typedef double (Matrix::*matrix_get2)(const int&, const int&) const;
    typedef void   (Matrix::*matrix_load)(const std::string&);
    typedef const Dimension& (Matrix::*matrix_ret_dimension)() const;

    class_<Matrix, SharedMatrix>("Matrix", "docstring").
            def(init<int, int>()).
            def(init<const std::string&, const Dimension&, const Dimension&>()).
            def(init<const std::string&>()).
            def("clone", &Matrix::clone, "docstring").
            def("set_name", &Matrix::set_name, "docstring").
            def("name", &Matrix::name, return_value_policy<copy_const_reference>(), "docstring").
            def("print_out", &Matrix::print_out, "docstring").
            def("rows", &Matrix::rowdim, "docstring").
            def("cols", &Matrix::coldim, "docstring").
            def("rowdim", matrix_ret_dimension(&Matrix::rowspi), return_value_policy<copy_const_reference>(), "docstring").
            def("coldim", matrix_ret_dimension(&Matrix::colspi), return_value_policy<copy_const_reference>(), "docstring").
            def("nirrep", &Matrix::nirrep, return_value_policy<copy_const_reference>(), "docstring").
            def("symmetry", &Matrix::symmetry, return_value_policy<copy_const_reference>(), "docstring").
            def("identity", &Matrix::identity, "docstring").
            def("copy_lower_to_upper", &Matrix::copy_lower_to_upper, "docstring").
            def("copy_upper_to_lower", &Matrix::copy_upper_to_lower, "docstring").
            def("zero_lower", &Matrix::zero_lower, "docstring").
            def("zero_upper", &Matrix::zero_upper, "docstring").
            def("zero", &Matrix::zero, "docstring").
            def("zero_diagonal", &Matrix::zero_diagonal, "docstring").
            def("trace", &Matrix::trace, "docstring").
            //            def("transpose", &Matrix::transpose).
            def("add", matrix_one(&Matrix::add), "docstring").
            def("subtract", matrix_one(&Matrix::subtract), "docstring").
            def("accumulate_product", matrix_two(&Matrix::accumulate_product), "docstring").
            def("scale", &Matrix::scale, "docstring").
            def("sum_of_squares", &Matrix::sum_of_squares, "docstring").
            def("add_and_orthogonalize_row", &Matrix::add_and_orthogonalize_row, "docstring").
            def("rms", &Matrix::rms, "docstring").
            def("scale_row", &Matrix::scale_row, "docstring").
            def("scale_column", &Matrix::scale_column, "docstring").
            def("transform", matrix_one(&Matrix::transform), "docstring").
            def("transform", matrix_two(&Matrix::transform), "docstring").
            def("transform", matrix_one(&Matrix::back_transform), "docstring").
            def("back_transform", matrix_two(&Matrix::back_transform), "docstring").
            def("vector_dot", double_matrix_one(&Matrix::vector_dot), "docstring").
            def("gemm", matrix_multiply(&Matrix::gemm), "docstring").
            def("diagonalize", matrix_diagonalize(&Matrix::diagonalize), "docstring").
            def("cholesky_factorize", &Matrix::cholesky_factorize, "docstring").
            def("partial_cholesky_factorize", &Matrix::partial_cholesky_factorize, "docstring").
            def("canonical_orthogonalization", &Matrix::canonical_orthogonalization, CanonicalOrthog()).
            def("invert", &Matrix::invert, "docstring").
            def("power", &Matrix::power, "docstring").
            def("get", matrix_get3(&Matrix::get), "docstring").
            def("get", matrix_get2(&Matrix::get), "docstring").
            def("set", matrix_set3(&Matrix::set), "docstring").
            def("set", matrix_set4(&Matrix::set), "docstring").
            def("set", &Matrix::set_by_python_list, "docstring").
            def("project_out", &Matrix::project_out, "docstring").
            def("__getitem__", &Matrix::pyget, "docstring").
            def("__setitem__", &Matrix::pyset, "docstring").
            def("save", matrix_save(&Matrix::save), "docstring").
            def("load", matrix_load(&Matrix::load), "docstring").
            def("load_mpqc", matrix_load(&Matrix::load_mpqc), "docstring").
            def("remove_symmetry", &Matrix::remove_symmetry, "docstring").
            add_property("__array_interface__", matrix_array_interface_c1, "docstring");

    class_<View, boost::noncopyable>("View", no_init).
            def(init<SharedMatrix, const Dimension&, const Dimension&>()).
            def(init<SharedMatrix, const Dimension&, const Dimension&, const Dimension&, const Dimension&>()).
            def("__call__", &View::operator(), "docstring");

    typedef SharedMatrix (MatrixFactory::*create_shared_matrix)();
    typedef SharedMatrix (MatrixFactory::*create_shared_matrix_name)(const std::string&);

    class_<MatrixFactory, boost::shared_ptr<MatrixFactory> >("MatrixFactory", "docstring").
            def("shared_object", &get_matrix_factory, "docstring").
            staticmethod("shared_object").
            def("create_matrix", create_shared_matrix(&MatrixFactory::create_shared_matrix), "docstring").
            def("create_matrix", create_shared_matrix_name(&MatrixFactory::create_shared_matrix), "docstring");

    class_<CdSalcList, boost::shared_ptr<CdSalcList>, boost::noncopyable>("CdSalcList", "docstring", no_init).
            def("print_out", &CdSalcList::print, "docstring").
            def("matrix", &CdSalcList::matrix, "docstring");

    class_<GaussianShell, boost::shared_ptr<GaussianShell> >("GaussianShell", "docstring", no_init).
            add_property("nprimitive", &GaussianShell::nprimitive, "docstring").
            add_property("nfunction", &GaussianShell::nfunction, "docstring").
            add_property("ncartesian", &GaussianShell::ncartesian, "docstring").
            add_property("am", &GaussianShell::am, "docstring").
            add_property("amchar", &GaussianShell::amchar, "docstring").
            add_property("AMCHAR", &GaussianShell::AMCHAR, "docstring").
            add_property("ncenter", &GaussianShell::ncenter, "docstring").
            add_property("function_index", &GaussianShell::function_index, &GaussianShell::set_function_index, "Basis function index where this shell starts.").
//            add_property("center", &GaussianShell::center, "A double* representing the center of the GaussianShell.").
//            add_property("exps", &GaussianShell::exps, "The exponents of all the primitives").
//            add_property("coefs", &GaussianShell::coefs, "The coefficients of all the primitives").
            def("is_cartesian", &GaussianShell::is_cartesian, "docstring").
            def("is_pure", &GaussianShell::is_pure, "docstring").
//            def("normalize_shell", &GaussianShell::normalize_shell, "docstring").
            def("exp", &GaussianShell::exp, "Returns the exponent of the given primitive").
            def("coef", &GaussianShell::coef, "docstring");


    class_<OneBodyAOInt, boost::shared_ptr<OneBodyAOInt>, boost::noncopyable>("OneBodyAOInt", "docstring", no_init).
            def("compute_shell", &OneBodyAOInt::compute_shell, "docstring").
            add_property("origin", &OneBodyAOInt::origin, &OneBodyAOInt::set_origin, "The origin about which the one body ints are being computed.").
            add_property("basis", &OneBodyAOInt::basis, "The basis set on center one").
            add_property("basis1", &OneBodyAOInt::basis1, "The basis set on center one").
            add_property("basis2", &OneBodyAOInt::basis2, "The basis set on center two").
            add_property("py_buffer_object", make_function(&OneBodyAOInt::py_buffer_object, return_internal_reference<>()), "docstring").
            def("set_enable_pybuffer", &OneBodyAOInt::set_enable_pybuffer, "docstring").
            add_property("py_buffer", &OneBodyAOInt::py_buffer, "docstring");

    //typedef void (OneBodySOInt::*matrix_version)(SharedMatrix) const;
    //typedef void (OneBodySOInt::*vector_version)(std::vector<SharedMatrix>) const;
    //class_<OneBodySOInt, boost::shared_ptr<OneBodySOInt>, boost::noncopyable>("OneBodySOInt", "docstring", no_init).
    //        def("compute", matrix_version(&OneBodySOInt::compute_shell), "docstring").
    //        def("compute_list", vector_version(&OneBodySOInt::compute), "docstring").
    //        add_property("basis", &OneBodySOInt::basis, "The basis set on center one").
    //        add_property("basis1", &OneBodySOInt::basis1, "The basis set on center one").
    //        add_property("basis2", &OneBodySOInt::basis2, "The basis set on center two");

    class_<OverlapInt, boost::shared_ptr<OverlapInt>, bases<OneBodyAOInt>, boost::noncopyable>("OverlapInt", "docstring", no_init);
    class_<DipoleInt, boost::shared_ptr<DipoleInt>, bases<OneBodyAOInt>, boost::noncopyable>("DipoleInt", "docstring", no_init);
    class_<QuadrupoleInt, boost::shared_ptr<QuadrupoleInt>, bases<OneBodyAOInt>, boost::noncopyable>("QuadrupoleInt", "docstring", no_init);
    class_<MultipoleInt, boost::shared_ptr<MultipoleInt>, bases<OneBodyAOInt>, boost::noncopyable>("MultipoleInt", "docstring", no_init);
    class_<TracelessQuadrupoleInt, boost::shared_ptr<TracelessQuadrupoleInt>, bases<OneBodyAOInt>, boost::noncopyable>("TracelessQuadrupoleInt", "docstring", no_init);
    class_<ElectricFieldInt, boost::shared_ptr<ElectricFieldInt>, bases<OneBodyAOInt>, boost::noncopyable>("ElectricFieldInt", "docstring", no_init);
    class_<KineticInt, boost::shared_ptr<KineticInt>, bases<OneBodyAOInt>, boost::noncopyable>("KineticInt", "docstring", no_init);
    class_<PotentialInt, boost::shared_ptr<PotentialInt>, bases<OneBodyAOInt>, boost::noncopyable>("PotentialInt", "docstring", no_init);
    class_<PseudospectralInt, boost::shared_ptr<PseudospectralInt>, bases<OneBodyAOInt>, boost::noncopyable>("PseudospectralInt", "docstring", no_init);
    class_<ElectrostaticInt, boost::shared_ptr<ElectrostaticInt>, bases<OneBodyAOInt>, boost::noncopyable>("ElectrostaticInt", "docstring", no_init);
    class_<NablaInt, boost::shared_ptr<NablaInt>, bases<OneBodyAOInt>, boost::noncopyable>("NablaInt", "docstring", no_init);
    class_<AngularMomentumInt, boost::shared_ptr<AngularMomentumInt>, bases<OneBodyAOInt>, boost::noncopyable>("AngularMomentumInt", "docstring", no_init);

    typedef size_t (TwoBodyAOInt::*compute_shell_ints)(int, int, int, int);
    class_<TwoBodyAOInt, boost::shared_ptr<TwoBodyAOInt>, boost::noncopyable>("TwoBodyAOInt", "docstring", no_init).
            def("compute_shell", compute_shell_ints(&TwoBodyAOInt::compute_shell), "docstring").
            add_property("py_buffer_object", make_function(&TwoBodyAOInt::py_buffer_object, return_internal_reference<>()), "docstring").
            add_property("py_buffer", &TwoBodyAOInt::py_buffer, "docstring").
            def("set_enable_pybuffer", &TwoBodyAOInt::set_enable_pybuffer, "docstring");

    class_<TwoElectronInt, boost::shared_ptr<TwoElectronInt>, bases<TwoBodyAOInt>, boost::noncopyable>("TwoElectronInt", "docstring", no_init);
            def("compute_shell", compute_shell_ints(&TwoBodyAOInt::compute_shell), "docstring");

    class_<ERI, boost::shared_ptr<ERI>, bases<TwoElectronInt>, boost::noncopyable>("ERI", "docstring", no_init);
    class_<F12, boost::shared_ptr<F12>, bases<TwoElectronInt>, boost::noncopyable>("F12", "docstring", no_init);
    class_<F12G12, boost::shared_ptr<F12G12>, bases<TwoElectronInt>, boost::noncopyable>("F12G12", "docstring", no_init);
    class_<F12Squared, boost::shared_ptr<F12Squared>, bases<TwoElectronInt>, boost::noncopyable>("F12Squared", "docstring", no_init);
    class_<F12DoubleCommutator, boost::shared_ptr<F12DoubleCommutator>, bases<TwoElectronInt>, boost::noncopyable>("F12DoubleCommutator", "docstring", no_init);
    class_<ErfERI, boost::shared_ptr<ErfERI>, bases<TwoElectronInt>, boost::noncopyable>("ErfERI", "docstring", no_init);
    class_<ErfComplementERI, boost::shared_ptr<ErfComplementERI>, bases<TwoElectronInt>,        boost::noncopyable>("ErfComplementERI", "docstring", no_init);

    class_<AOShellCombinationsIterator, boost::shared_ptr<AOShellCombinationsIterator>, boost::noncopyable>("AOShellCombinationsIterator", no_init).
            add_property("p", &AOShellCombinationsIterator::p, "docstring").
            add_property("q", &AOShellCombinationsIterator::q, "docstring").
            add_property("r", &AOShellCombinationsIterator::r, "docstring").
            add_property("s", &AOShellCombinationsIterator::s, "docstring").
            def("first", &AOShellCombinationsIterator::first, "docstring").
            def("next", &AOShellCombinationsIterator::next, "docstring").
            def("is_done", &AOShellCombinationsIterator::is_done, "docstring");

    class_<ThreeCenterOverlapInt, boost::shared_ptr<ThreeCenterOverlapInt>, boost::noncopyable>("ThreeCenterOverlapInt", "docstring", no_init).
            def("compute_shell", &ThreeCenterOverlapInt::compute_shell, "docstring").
            add_property("py_buffer_object", make_function(&ThreeCenterOverlapInt::py_buffer_object, return_internal_reference<>()), "docstring").
            def("set_enable_pybuffer", &ThreeCenterOverlapInt::set_enable_pybuffer, "docstring");

    class_<IntegralFactory, boost::shared_ptr<IntegralFactory>, boost::noncopyable>("IntegralFactory", "docstring", no_init).
            def(init<boost::shared_ptr<BasisSet>, boost::shared_ptr<BasisSet>, boost::shared_ptr<BasisSet>, boost::shared_ptr<BasisSet> >()).
            def(init<boost::shared_ptr<BasisSet>, boost::shared_ptr<BasisSet> >()).
            def(init<boost::shared_ptr<BasisSet> >()).
            def("shells_iterator", &IntegralFactory::shells_iterator_ptr, return_value_policy<manage_new_object>(), "docstring").
            def("eri", &IntegralFactory::eri, eri_overloads("docstring")[return_value_policy<manage_new_object>()]).
            def("f12", &IntegralFactory::f12, f12_overloads("docstring")[return_value_policy<manage_new_object>()]).
            def("f12g12", &IntegralFactory::f12g12, f12g12_overloads("docstring")[return_value_policy<manage_new_object>()]).
            def("f12_double_commutator", &IntegralFactory::f12_double_commutator, f12_double_commutator_overloads("docstring")[return_value_policy<manage_new_object>()]).
            def("f12_squared", &IntegralFactory::f12_squared, f12_squared_overloads("docstring")[return_value_policy<manage_new_object>()]).
            def("erf_eri", &IntegralFactory::erf_eri, erf_eri_overloads("docstring")[return_value_policy<manage_new_object>()]).
            def("erf_complement_eri", &IntegralFactory::erf_complement_eri, erf_complement_eri_overloads("docstring")[return_value_policy<manage_new_object>()]).
            def("ao_overlap", &IntegralFactory::ao_overlap, ao_overlap_overloads("docstring")[return_value_policy<manage_new_object>()]).
            def("so_overlap", &IntegralFactory::so_overlap, so_overlap_overloads("docstring")[return_value_policy<manage_new_object>()]).
            def("ao_dipole", &IntegralFactory::ao_dipole, ao_dipole_overloads("docstring")[return_value_policy<manage_new_object>()]).
            def("ao_kinetic", &IntegralFactory::ao_kinetic, ao_kinetic_overloads("docstring")[return_value_policy<manage_new_object>()]).
            def("ao_potential", &IntegralFactory::ao_potential, ao_potential_overloads("docstring")[return_value_policy<manage_new_object>()]).
            def("ao_pseudospectral", &IntegralFactory::ao_pseudospectral, ao_pseudospectral_overloads("docstring")[return_value_policy<manage_new_object>()]).
            def("ao_nabla", &IntegralFactory::ao_nabla, ao_nabla_overloads("docstring")[return_value_policy<manage_new_object>()]).
            def("ao_angular_momentum", &IntegralFactory::ao_angular_momentum, ao_angular_momentum_overloads("docstring")[return_value_policy<manage_new_object>()]).
            def("ao_quadrupole", &IntegralFactory::ao_quadrupole, return_value_policy<manage_new_object>(), "docstring").
            def("ao_multipoles", &IntegralFactory::ao_multipoles, return_value_policy<manage_new_object>(), "docstring").
            def("so_multipoles", &IntegralFactory::so_multipoles, return_value_policy<manage_new_object>(), "docstring").
            def("ao_traceless_quadrupole", &IntegralFactory::ao_traceless_quadrupole, return_value_policy<manage_new_object>(), "docstring").
            def("electric_field", &IntegralFactory::electric_field, return_value_policy<manage_new_object>(), "docstring").
            def("electrostatic", &IntegralFactory::electrostatic, return_value_policy<manage_new_object>(), "docstring").
            def("overlap_3c", &IntegralFactory::overlap_3c, return_value_policy<manage_new_object>(), "docstring");

    typedef boost::shared_ptr<PetiteList> (MintsHelper::*petite_list_0)() const;
    typedef boost::shared_ptr<PetiteList> (MintsHelper::*petite_list_1)(bool) const;

    typedef SharedMatrix (MintsHelper::*erf)(double, SharedMatrix, SharedMatrix, SharedMatrix, SharedMatrix);
    typedef SharedMatrix (MintsHelper::*eri)(SharedMatrix, SharedMatrix, SharedMatrix, SharedMatrix);

    class_<MintsHelper, boost::shared_ptr<MintsHelper> >("MintsHelper", "docstring").
            def(init<boost::shared_ptr<BasisSet> >()).
            def("integral", &MintsHelper::integral, "docstring").
            def("integrals", &MintsHelper::integrals, "docstring").
            def("integrals_erf", &MintsHelper::integrals_erf, "docstring").
            def("integrals_erfc", &MintsHelper::integrals_erfc, "docstring").
            def("one_electron_integrals", &MintsHelper::one_electron_integrals, "docstring").
            def("basisset", &MintsHelper::basisset, "docstring").
            def("sobasisset", &MintsHelper::sobasisset, "docstring").
            def("factory", &MintsHelper::factory, "docstring").
            def("ao_overlap", &MintsHelper::ao_overlap, "docstring").
            def("ao_kinetic", &MintsHelper::ao_kinetic, "docstring").
            def("ao_potential", &MintsHelper::ao_potential, "docstring").
            def("one_electron_integrals", &MintsHelper::one_electron_integrals, "docstring").
            def("so_overlap", &MintsHelper::so_overlap, "docstring").
            def("so_kinetic", &MintsHelper::so_kinetic, "docstring").
            def("so_potential", &MintsHelper::so_potential, "docstring").
            def("so_dipole", &MintsHelper::so_dipole, "docstring").
            def("so_quadrupole", &MintsHelper::so_quadrupole, "docstring").
            def("so_traceless_quadrupole", &MintsHelper::so_traceless_quadrupole, "docstring").
            def("ao_nabla", &MintsHelper::ao_nabla, "docstring").
            def("so_nabla", &MintsHelper::so_nabla, "docstring").
            def("so_angular_momentum", &MintsHelper::so_angular_momentum, "docstring").
            def("ao_angular_momentum", &MintsHelper::ao_angular_momentum, "docstring").
            def("ao_eri", &MintsHelper::ao_eri, "docstring").
            def("ao_eri_shell", &MintsHelper::ao_eri_shell, "docstring").
            def("ao_erf_eri", &MintsHelper::ao_erf_eri, "docstring").
            def("ao_f12", &MintsHelper::ao_f12, "docstring").
            def("ao_f12_squared", &MintsHelper::ao_f12_squared, "docstring").
            def("ao_f12g12", &MintsHelper::ao_f12g12, "docstring").
            def("ao_f12_double_commutator", &MintsHelper::ao_f12_double_commutator, "docstring").
            def("mo_eri", eri(&MintsHelper::mo_eri), "docstring").
            def("mo_erf_eri", erf(&MintsHelper::mo_erf_eri), "docstring").
            def("mo_f12", &MintsHelper::mo_f12, "docstring").
            def("mo_f12_squared", &MintsHelper::mo_f12_squared, "docstring").
            def("mo_f12g12", &MintsHelper::mo_f12g12, "docstring").
            def("mo_f12_double_commutator", &MintsHelper::mo_f12_double_commutator, "docstring").
            def("cdsalcs", &MintsHelper::cdsalcs, "docstring").
            def("petite_list", petite_list_0(&MintsHelper::petite_list), "docstring").
            def("petite_list1", petite_list_1(&MintsHelper::petite_list), "docstring").
            def("play", &MintsHelper::play, "docstring");

    class_<FittingMetric, boost::shared_ptr<FittingMetric> >("FittingMetric", "docstring").
            def("get_algorithm", &FittingMetric::get_algorithm, "docstring").
            def("is_poisson", &FittingMetric::is_poisson, "docstring").
            def("is_inverted", &FittingMetric::is_inverted, "docstring").
            def("get_metric", &FittingMetric::get_metric, "docstring").
            def("get_pivots", &FittingMetric::get_pivots, "docstring").
            def("get_reverse_pivots", &FittingMetric::get_reverse_pivots, "docstring").
            def("form_fitting_metric", &FittingMetric::form_fitting_metric, "docstring").
            def("form_cholesky_inverse", &FittingMetric::form_cholesky_inverse, "docstring").
            def("form_QR_inverse", &FittingMetric::form_QR_inverse, "docstring").
            def("form_eig_inverse", &FittingMetric::form_eig_inverse, "docstring").
            def("form_full_inverse", &FittingMetric::form_full_inverse, "docstring");

    class_<PseudoTrial, boost::shared_ptr<PseudoTrial> >("PseudoTrial", "docstring").
            def("getI", &PseudoTrial::getI, "docstring").
            def("getIPS", &PseudoTrial::getIPS, "docstring").
            def("getQ", &PseudoTrial::getQ, "docstring").
            def("getR", &PseudoTrial::getR, "docstring").
            def("getA", &PseudoTrial::getA, "docstring");

    class_<Vector3>("Vector3", "Class for vectors of length three, often Cartesian coordinate vectors, and their common operations").
            def(init<double>()).
            def(init<double, double, double>()).
            def(init<const Vector3&>()).
            //      def(self = other<double>()).
            def(self += self).
            def(self -= self).
            def(self *= other<double>()).
            def(self + self).
            def(self - self).
            def(-self).
            def("dot", &Vector3::dot, "Returns dot product of arg1 and arg2").
            def("distance", &Vector3::distance, "Returns distance between two points represented by arg1 and arg2").
            def("normalize", &Vector3::normalize, "Returns vector of unit length and arg1 direction").
            def("norm", &Vector3::norm, "Returns Euclidean norm of arg1").
            def("cross", &Vector3::cross, "Returns cross product of arg1 and arg2").
            def("__str__", &Vector3::to_string, "Returns a string representation of arg1, suitable for printing.").
            def("__getitem__", &Vector3::get, "Returns the arg2-th element of arg1.");

    typedef void (SymmetryOperation::*intFunction)(int);
    typedef void (SymmetryOperation::*doubleFunction)(double);

    class_<SymmetryOperation>("SymmetryOperation", "Class to provide a 3 by 3 matrix representation of a symmetry operation, such as a rotation or reflection.").
            def(init<const SymmetryOperation& >()).
            def("trace", &SymmetryOperation::trace, "Returns trace of transformation matrix").
            def("zero", &SymmetryOperation::zero, "Zero out the symmetry operation").
            def("operate", &SymmetryOperation::operate, "Performs the operation arg2 * arg1").
            def("transform", &SymmetryOperation::transform, "Performs the transform arg2 * arg1 * arg2~").
            def("unit", &SymmetryOperation::unit, "Set equal to a unit matrix").
            def("E", &SymmetryOperation::E, "Set equal to E").
            def("i", &SymmetryOperation::i, "Set equal to an inversion").
            def("sigma_xy", &SymmetryOperation::sigma_xy, "Set equal to reflection in xy plane").
            def("sigma_yz", &SymmetryOperation::sigma_yz, "Set equal to reflection in yz plane").
            def("sigma_xz", &SymmetryOperation::sigma_xz, "Set equal to reflection in xz plane").
            //        def("sigma_yz", &SymmetryOperation::sigma_yz).
            def("rotate_n", intFunction(&SymmetryOperation::rotation), "Set equal to a clockwise rotation by 2pi/n").
            def("rotate_theta", doubleFunction(&SymmetryOperation::rotation), "Set equal to a clockwise rotation by theta").
            def("c2_x", &SymmetryOperation::c2_x, "Set equal to C2 about the x axis").
            def("c2_y", &SymmetryOperation::c2_y, "Set equal to C2 about the y axis").
            def("c2_z", &SymmetryOperation::c2_z, "Set equal to C2 about the z axis").
            def("transpose", &SymmetryOperation::transpose, "Performs transposition of matrix operation");

    class_<OrbitalSpace>("OrbitalSpace", "docstring", no_init).
            def(init<const std::string&, const std::string&, const SharedMatrix&, const SharedVector&, const boost::shared_ptr<BasisSet>&, const boost::shared_ptr<IntegralFactory>& >()).
            def(init<const std::string&, const std::string&, const SharedMatrix&, const boost::shared_ptr<BasisSet>&, const boost::shared_ptr<IntegralFactory>& >()).
            def(init<const std::string&, const std::string&, const boost::shared_ptr<Wavefunction>& >()).
            def("nirrep", &OrbitalSpace::nirrep, "docstring").
            def("id", &OrbitalSpace::id, return_value_policy<copy_const_reference>(), "docstring").
            def("name", &OrbitalSpace::name, return_value_policy<copy_const_reference>(), "docstring").
            def("C", &OrbitalSpace::C, return_value_policy<copy_const_reference>(), "docstring").
            def("evals", &OrbitalSpace::evals, return_value_policy<copy_const_reference>(), "docstring").
            def("basisset", &OrbitalSpace::basisset, return_value_policy<copy_const_reference>(), "docstring").
            def("integral", &OrbitalSpace::integral, return_value_policy<copy_const_reference>(), "docstring").
            def("dim", &OrbitalSpace::dim, return_value_policy<copy_const_reference>(), "docstring").
            def("print_out", &OrbitalSpace::print, "docstring").
            def("build_cabs_space", &OrbitalSpace::build_cabs_space, "docstring").
            staticmethod("build_cabs_space").
            def("build_ri_space", &OrbitalSpace::build_ri_space, "docstring").
            staticmethod("build_ri_space");

    class_<PointGroup, boost::shared_ptr<PointGroup> >("PointGroup", "docstring").
            def(init<const std::string&>()).
            def("symbol", &PointGroup::symbol, "Returns Schoenflies symbol for point group");
            //def("origin", &PointGroup::origin).
//            def("set_symbol", &PointGroup::set_symbol);

    typedef void (Molecule::*matrix_set_geometry)(const Matrix &);

    class_<Molecule, boost::shared_ptr<Molecule> >("Molecule", "Class to store the elements, coordinates, fragmentation pattern, basis sets, charge, multiplicity, etc. of a molecule.").
            def("set_geometry", matrix_set_geometry(&Molecule::set_geometry), "Sets the geometry, given a (Natom X 3) matrix arg2 of coordinates (in Bohr)").
            def("set_name", &Molecule::set_name, "Sets molecule name").
            def("name", &Molecule::name, "Gets molecule name").
            def("reinterpret_coordentry", &Molecule::set_reinterpret_coordentry, "Do reinterpret coordinate entries during update_geometry().").
            def("fix_orientation", &Molecule::set_orientation_fixed, "Fix the orientation at its current frame").
            def("fix_com", &Molecule::set_com_fixed, "Whether to fix the Cartesian position, or to translate to the C.O.M.").
            def("init_with_checkpoint", &Molecule::init_with_chkpt, "Populate arg1 member data with information from checkpoint file arg2").
            def("save_to_checkpoint", &Molecule::save_to_chkpt, "Saves molecule information to checkpoint file arg2 with prefix arg3").
            def("init_with_io", &Molecule::init_with_psio, "Creates a new checkpoint file with information from arg2").
            def("add_atom", &Molecule::add_atom, "Adds to Molecule arg1 an atom with atomic number arg2, Cartesian coordinates in Bohr (arg3, arg4, arg5), atomic symbol arg6, mass arg7, charge arg8 (optional), and lineno arg9 (optional)").
            def("natom", &Molecule::natom, "Number of real atoms").
            def("multiplicity", &Molecule::multiplicity, "Gets the multiplicity (defined as 2Ms + 1)").
            def("nfragments", &Molecule::nfragments, "Gets the number of fragments in the molecule").
            def("print_in_input_format", &Molecule::print_in_input_format, "Prints the molecule as Cartesian or ZMatrix entries, just as inputted.").
            def("create_psi4_string_from_molecule", &Molecule::create_psi4_string_from_molecule, "Gets a string reexpressing in input format the current states of the molecule").
            def("save_xyz_file", &Molecule::save_xyz_file, "Saves an XYZ file to arg2").
            def("save_string_xyz_file", &Molecule::save_string_xyz_file, "Saves an XYZ file to arg2").
            def("save_string_xyz", &Molecule::save_string_xyz, "Saves the string of an XYZ file to arg2").
            def("Z", &Molecule::Z, return_value_policy<copy_const_reference>(), "Nuclear charge of atom").
            def("x", &Molecule::x, "x position of atom").
            def("y", &Molecule::y, "y position of atom").
            def("z", &Molecule::z, "z position of atom").
            //def("xyz", &Molecule::xyz).
            def("center_of_mass", &Molecule::center_of_mass, "Computes center of mass of molecule (does not translate molecule)").
            def("translate", &Molecule::translate, "Translates molecule by arg2").
            def("move_to_com", &Molecule::move_to_com, "Moves molecule to center of mass").
            def("mass", &Molecule::mass, "Gets mass of atom arg2").
            def("symbol", &Molecule::symbol, "Gets the cleaned up label of atom arg2 (C2 => C, H4 = H)").
            def("label", &Molecule::label, "Gets the original label of the atom as given in the input file (C2, H4)").
            def("charge", &Molecule::charge, "Gets charge of atom").
            def("molecular_charge", &Molecule::molecular_charge, "Gets the molecular charge").
            def("extract_subsets", &Molecule::py_extract_subsets_1, "Returns copy of arg1 with arg2 fragments Real and arg3 fragments Ghost").
            def("extract_subsets", &Molecule::py_extract_subsets_2, "Returns copy of arg1 with arg2 fragments Real and arg3 fragment Ghost").
            def("extract_subsets", &Molecule::py_extract_subsets_3, "Returns copy of arg1 with arg2 fragment Real and arg3 fragments Ghost").
            def("extract_subsets", &Molecule::py_extract_subsets_4, "Returns copy of arg1 with arg2 fragment Real and arg3 fragment Ghost").
            def("extract_subsets", &Molecule::py_extract_subsets_5, "Returns copy of arg1 with arg2 fragments Real").
            def("extract_subsets", &Molecule::py_extract_subsets_6, "Returns copy of arg1 with arg2 fragment Real").
            def("activate_all_fragments", &Molecule::activate_all_fragments, "Sets all fragments in the molecule to be active").
            def("deactivate_all_fragments", &Molecule::deactivate_all_fragments, "Sets all fragments in the molecule to be inactive").
            def("set_active_fragments", &Molecule::set_active_fragments, "Sets the specified list arg2 of fragments to be Real").
            def("set_active_fragment", &Molecule::set_active_fragment, "Sets the specified fragment arg2 to be Real").
            def("set_ghost_fragments", &Molecule::set_ghost_fragments, "Sets the specified list arg2 of fragments to be Ghost").
            def("set_ghost_fragment", &Molecule::set_ghost_fragment, "Sets the specified fragment arg2 to be Ghost").
            def("atom_at_position", &Molecule::atom_at_position1, "Tests to see if an atom is at the position arg2 with a given tolerance arg3").
            def("print_out", &Molecule::print, "Prints the molecule in Cartesians in input units").
            def("print_out_in_bohr", &Molecule::print_in_bohr, "Prints the molecule in Cartesians in Bohr").
            def("print_out_in_angstrom", &Molecule::print_in_angstrom, "Prints the molecule in Cartesians in Angstroms").
            def("rotational_constants", &Molecule::rotational_constants, "Prints the rotational constants of the molecule").
            def("nuclear_repulsion_energy", &Molecule::nuclear_repulsion_energy, "Computes nuclear repulsion energy").
            def("find_point_group", &Molecule::find_point_group, "Finds computational molecular point group, user can override this with the symmetry keyword").
            def("reset_point_group", &Molecule::reset_point_group, "Overrides symmetry from outside the molecule string").
            def("set_point_group", &Molecule::set_point_group, "Sets the molecular point group to the point group object arg2").
            def("get_full_point_group", &Molecule::full_point_group, "Gets point group name such as C3v or S8").
            def("point_group", &Molecule::point_group, "Returns the current point group object").
            def("schoenflies_symbol", &Molecule::schoenflies_symbol, "Returns the Schoenflies symbol").
            def("form_symmetry_information", &Molecule::form_symmetry_information, "Uses the point group object obtain by calling point_group()").
            def("symmetrize", &Molecule::symmetrize_to_abelian_group, "Finds the highest point Abelian point group within the specified tolerance, and forces the geometry to have that symmetry.").
            def("create_molecule_from_string", &Molecule::create_molecule_from_string, "Returns a new Molecule with member data from the geometry string arg1 in psi4 format").
            staticmethod("create_molecule_from_string").
                def("is_variable", &Molecule::is_variable, "Checks if variable arg2 is in the list, returns true if it is, and returns false if not").
            def("set_variable", &Molecule::set_variable, "Assigns the value arg3 to the variable arg2 in the list of geometry variables, then calls update_geometry()").
            def("get_variable", &Molecule::get_variable, "Checks if variable arg2 is in the list, sets it to val and returns true if it is, and returns false if not").
            def("update_geometry", &Molecule::update_geometry, "Reevaluates the geometry with current variable values, orientation directives, etc. Must be called after initial Molecule definition by string.").
            def("set_molecular_charge", &Molecule::set_molecular_charge, "Sets the molecular charge").
            def("set_multiplicity", &Molecule::set_multiplicity, "Sets the multiplicity (defined as 2Ms + 1)").
            def("set_basis_all_atoms", &Molecule::set_basis_all_atoms, "Sets basis set arg2 to all atoms").
            def("set_basis_by_symbol", &Molecule::set_basis_by_symbol, "Sets basis set arg3 to all atoms with symbol (e.g., H) arg2").
            def("set_basis_by_label", &Molecule::set_basis_by_label, "Sets basis set arg3 to all atoms with label (e.g., H4) arg2").
            //def("set_basis_by_number", &Molecule::set_basis_by_number, "Sets basis set arg3 to atom number (1-indexed, incl. dummies) arg2").  // dangerous for user use
            add_property("units", &Molecule::units, &Molecule::set_units, "Units (Angstrom or Bohr) used to define the geometry").
            def("clone", &Molecule::clone, "Returns a new Molecule identical to arg1").
            def("geometry", &Molecule::geometry, "Gets the geometry as a (Natom X 3) matrix of coordinates (in Bohr)");

    class_<PetiteList, boost::shared_ptr<PetiteList>, boost::noncopyable>("PetiteList", "docstring", no_init).
            def("aotoso", &PetiteList::aotoso, "docstring").
            def("sotoao", &PetiteList::sotoao, "docstring").
            def("print", &PetiteList::print, "docstring");

    class_<BasisSetParser, boost::shared_ptr<BasisSetParser>, boost::noncopyable>("BasisSetParser", "docstring", no_init);
    class_<Gaussian94BasisSetParser, boost::shared_ptr<Gaussian94BasisSetParser>, bases<BasisSetParser> >("Gaussian94BasisSetParser", "docstring");

    typedef void (BasisSet::*basis_print_out)() const;
    typedef const GaussianShell& (BasisSet::*no_center_version)(int) const;
    typedef const GaussianShell& (BasisSet::*center_version)(int, int) const;
    typedef boost::shared_ptr<BasisSet> (BasisSet::*ptrversion)(const boost::shared_ptr<BasisSet>&) const;
    class_<BasisSet, boost::shared_ptr<BasisSet>, boost::noncopyable>("BasisSet", "docstring", no_init).
            def("print_out", basis_print_out(&BasisSet::print), "docstring").
            def("print_detail_out", basis_print_out(&BasisSet::print_detail), "docstring").
            def("genbas", &BasisSet::print_detail_cfour, "Returns basis set per atom in CFOUR format").
            def("make_filename", &BasisSet::make_filename, "Returns filename for basis name: pluses, stars, parentheses replaced and gbs extension added").
            staticmethod("make_filename").
            def("construct", &BasisSet::construct, "docstring").
            staticmethod("construct").
            def("zero_ao_basis_set", &BasisSet::zero_ao_basis_set, "Returns a BasisSet object that actually has a single s-function at the origin with an exponent of 0.0 and contraction of 1.0.").
            staticmethod("zero_ao_basis_set").
            def("nbf", &BasisSet::nbf, "Returns number of basis functions (Cartesian or spherical depending on has_puream)").
            def("nao", &BasisSet::nao, "Returns number of atomic orbitals (Cartesian)").
            def("nprimitive", &BasisSet::nprimitive, "Returns total number of primitives in all contractions").
            def("nshell", &BasisSet::nshell, "Returns number of shells").
            def("shell", no_center_version(&BasisSet::shell), return_value_policy<copy_const_reference>(), "docstring").
            def("shell", center_version(&BasisSet::shell), return_value_policy<copy_const_reference>(), "docstring").
            def("max_am", &BasisSet::max_am, "Returns maximum angular momentum used").
            def("has_puream", &BasisSet::has_puream, "Spherical harmonics?").
            def("shell_to_basis_function", &BasisSet::shell_to_basis_function, "docstring").
            def("shell_to_center", &BasisSet::shell_to_center, "docstring").
            def("shell_to_ao_function", &BasisSet::shell_to_ao_function, "docstring").
            def("function_to_shell", &BasisSet::function_to_shell, "docstring").
            def("function_to_center", &BasisSet::function_to_center, "Given a function number, return the number of the center it is on.").
            def("nshell_on_center", &BasisSet::nshell_on_center, "docstring").
            def("ao_to_shell", &BasisSet::ao_to_shell, "docstring").
            def("pyconstruct_orbital", &BasisSet::pyconstruct_orbital, pyconstruct_orb_overloads("Returns new BasisSet for Molecule arg1 for target keyword name arg2 and target keyword value arg3. This suffices for orbital basis sets. For auxiliary basis sets, a default fitting role (e.g., RIFIT, JKFIT) arg4 and orbital keyword value arg5 are required. An optional argument to force the puream setting is arg4 for orbital basis sets and arg6 for auxiliary basis sets.")).
            staticmethod("pyconstruct_orbital").
            def("pyconstruct_auxiliary", &BasisSet::pyconstruct_auxiliary, pyconstruct_aux_overloads("Returns new BasisSet for Molecule arg1 for target keyword name arg2 and target keyword value arg3. This suffices for orbital basis sets. For auxiliary basis sets, a default fitting role (e.g., RIFIT, JKFIT) arg4 and orbital keyword value arg5 are required. An optional argument to force the puream setting is arg4 for orbital basis sets and arg6 for auxiliary basis sets.")).
            staticmethod("pyconstruct_auxiliary").
            def("concatenate", ptrversion(&BasisSet::concatenate), "Concatenates two basis sets together into a new basis without reordering anything. Unless you know what you're doing, you should use the '+' operator instead of this method.").
            def("add", ptrversion(&BasisSet::add), "Combine two basis sets to make a new one.").
            //staticmethod("concatinate").
            def(self + self);

    class_<SOBasisSet, boost::shared_ptr<SOBasisSet>, boost::noncopyable>("SOBasisSet", "docstring", no_init).
            def("petite_list", &SOBasisSet::petite_list, "docstring");

    class_<ExternalPotential, boost::shared_ptr<ExternalPotential>, boost::noncopyable>("ExternalPotential", "docstring").
            def("setName", &ExternalPotential::setName, "docstring").
            def("addCharge", &ExternalPotential::addCharge, "docstring").
            def("addBasis", &ExternalPotential::addBasis, "docstring").
            def("clear", &ExternalPotential::clear, "docstring").
            def("computePotentialMatrix", &ExternalPotential::computePotentialMatrix, "docstring").
            def("print_out", &ExternalPotential::py_print, "docstring");

    class_<DFChargeFitter, boost::shared_ptr<DFChargeFitter>, boost::noncopyable>("DFChargeFitter", "docstring").
            def("setPrimary", &DFChargeFitter::setPrimary, "docstring").
            def("setAuxiliary", &DFChargeFitter::setAuxiliary, "docstring").
            def("setD", &DFChargeFitter::setD, "docstring").
            def("d", &DFChargeFitter::d, "docstring").
            def("fit", &DFChargeFitter::fit, "docstring");

    class_<Wavefunction, boost::shared_ptr<Wavefunction>, boost::noncopyable>("Wavefunction", "docstring", no_init).
            def("nso", &Wavefunction::nso, "docstring").
            def("nmo", &Wavefunction::nmo, "docstring").
            def("nirrep", &Wavefunction::nirrep, "docstring").
            def("Ca_subset", &Wavefunction::Ca_subset, "docstring").
            def("Cb_subset", &Wavefunction::Cb_subset, "docstring").
            def("epsilon_a_subset", &Wavefunction::epsilon_a_subset, "docstring").
            def("epsilon_b_subset", &Wavefunction::epsilon_b_subset, "docstring").
            def("Ca", &Wavefunction::Ca, "docstring").
            def("Cb", &Wavefunction::Cb, "docstring").
            def("Fa", &Wavefunction::Fa, "docstring").
            def("Fb", &Wavefunction::Fb, "docstring").
            def("Da", &Wavefunction::Da, "docstring").
            def("Db", &Wavefunction::Db, "docstring").
            def("epsilon_a", &Wavefunction::epsilon_a, "docstring").
            def("epsilon_b", &Wavefunction::epsilon_b, "docstring").
            def("add_preiteration_callback", &Wavefunction::add_preiteration_callback, "docstring").
            def("add_postiteration_callback", &Wavefunction::add_postiteration_callback, "docstring").
            def("basisset", &Wavefunction::basisset, "docstring").
            def("sobasisset", &Wavefunction::sobasisset, "docstring").
            def("energy", &Wavefunction::reference_energy, "docstring").
            def("gradient", &Wavefunction::gradient, "docstring").
            def("frequencies", &Wavefunction::frequencies, "docstring").
            def("alpha_orbital_space", &Wavefunction::alpha_orbital_space, "docstring").
            def("beta_orbital_space", &Wavefunction::beta_orbital_space, "docstring").
            def("molecule", &Wavefunction::molecule, "docstring").
            def("doccpi", &Wavefunction::doccpi, return_value_policy<copy_const_reference>(), "docstring").
            def("soccpi", &Wavefunction::soccpi, return_value_policy<copy_const_reference>(), "docstring").
            def("nsopi", &Wavefunction::nsopi, return_value_policy<copy_const_reference>(), "docstring").
            def("nmopi", &Wavefunction::nmopi, return_value_policy<copy_const_reference>(), "docstring").
            def("nalphapi", &Wavefunction::nalphapi, return_value_policy<copy_const_reference>(), "docstring").
            def("nbetapi", &Wavefunction::nbetapi, return_value_policy<copy_const_reference>(), "docstring").
            def("frzcpi", &Wavefunction::frzcpi, return_value_policy<copy_const_reference>(), "docstring").
            def("frzvpi", &Wavefunction::frzvpi, return_value_policy<copy_const_reference>(), "docstring").
            def("nalpha", &Wavefunction::nalpha, "docstring").
            def("nbeta", &Wavefunction::nbeta, "docstring");

    class_<scf::HF, boost::shared_ptr<scf::HF>, bases<Wavefunction>, boost::noncopyable>("HF", "docstring", no_init);
    class_<scf::RHF, boost::shared_ptr<scf::RHF>, bases<scf::HF, Wavefunction> >("RHF", "docstring", no_init);

    typedef boost::shared_ptr<Localizer> (*localizer_with_type)(const std::string&, boost::shared_ptr<BasisSet>, boost::shared_ptr<Matrix>);

    class_<Localizer, boost::shared_ptr<Localizer>, boost::noncopyable>("Localizer", "docstring", no_init).
            def("build", localizer_with_type(&Localizer::build), "docstring").
            staticmethod("build").
            def("localize", &Localizer::localize, "Perform the localization procedure").
            add_property("L", &Localizer::L, "Localized orbital coefficients").
            add_property("U", &Localizer::U, "Orbital rotation matrix").
            add_property("converged", &Localizer::converged, "Did the localization procedure converge?");

    class_<BoysLocalizer, boost::shared_ptr<BoysLocalizer>, bases<Localizer> >("BoysLocalizer", "docstring", no_init);
    class_<PMLocalizer, boost::shared_ptr<PMLocalizer>, bases<Localizer> >("PMLocalizer", "docstring", no_init);

    class_<MoldenWriter, boost::shared_ptr<MoldenWriter> >("MoldenWriter", "docstring", no_init).
            def(init<boost::shared_ptr<Wavefunction> >()).
            def("write", &MoldenWriter::write, "docstring");

    class_<NBOWriter, boost::shared_ptr<NBOWriter> >("NBOWriter", "docstring", no_init).
            def(init<boost::shared_ptr<Wavefunction> >()).
            def("write", &NBOWriter::write, "docstring");

    class_<OperatorSymmetry, boost::shared_ptr<OperatorSymmetry> >("MultipoleSymmetry", "docstring", no_init).
            def(init<int, const boost::shared_ptr<Molecule>&,
                const boost::shared_ptr<IntegralFactory>&,
                const boost::shared_ptr<MatrixFactory>&>()).
            def("create_matrices", &OperatorSymmetry::create_matrices, "docstring");

    class_<CorrelationFactor, boost::shared_ptr<CorrelationFactor>, boost::noncopyable>("CorrelationFactor", "docstring", no_init).
            def(init<unsigned int>()).
            def(init<boost::shared_ptr<Vector>, boost::shared_ptr<Vector> >()).
            def("set_params", &CorrelationFactor::set_params, "docstring");
    class_<FittedSlaterCorrelationFactor, bases<CorrelationFactor>, boost::noncopyable>("FittedSlaterCorrelationFactor", "docstring", no_init).
            def(init<double>()).
            def("exponent", &FittedSlaterCorrelationFactor::exponent);


    // LIBFOCK wrappers
    class_<JK, boost::shared_ptr<JK>, boost::noncopyable>("JK", "docstring", no_init)
//            .def(init<boost::shared_ptr<BasisSet> >())
            .def("build_JK", &JK::build_JK, "docstring")
            .staticmethod("build_JK")
            .def("initialize", &JK::initialize)
            .def("compute", &JK::compute)
            .def("finalize", &JK::finalize)
            .def("C_left", &JK::C_left, return_internal_reference<>())
            .def("C_right", &JK::C_right, return_internal_reference<>())
            .def("J", &JK::J, return_internal_reference<>())
            .def("K", &JK::K, return_internal_reference<>())
            .def("wK", &JK::wK, return_internal_reference<>())
            .def("D", &JK::D, return_internal_reference<>())
            .def("print_header", &JK::print_header, "docstring")
            ;
}
