#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <libmints/mints.h>
#include <lib3index/3index.h>
#include <libscf_solver/hf.h>
#include <libscf_solver/rhf.h>

#include <string>

using namespace boost;
using namespace boost::python;
using namespace psi;

boost::shared_ptr<Vector> py_nuclear_dipole(shared_ptr<Molecule> mol)
{
    //SharedMolecule mol = Process::environment.molecule();
    return DipoleInt::nuclear_contribution(mol);
}

boost::shared_ptr<MatrixFactory> get_matrix_factory()
{
    // We need a valid molecule with a valid point group to create a matrix factory.
    boost::shared_ptr<Molecule> molecule = Process::environment.molecule();
    if (!molecule) {
        fprintf(outfile, "  Active molecule not set!");
        throw PSIEXCEPTION("Active molecule not set!");
    }
    if (!molecule->point_group()) {
        fprintf(outfile, "  Active molecule does not have point group set!");
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

void export_mints()
{
    def("nuclear_dipole", py_nuclear_dipole);

    // This is needed to wrap an STL vector into Boost.Python. Since the vector
    // is going to contain boost::shared_ptr's we MUST set the no_proxy flag to true
    // (as it is) to tell Boost.Python to not create a proxy class to handle
    // the vector's data type.
    class_<std::vector<SharedMatrix > >("matrix_vector").
            def(vector_indexing_suite<std::vector<SharedMatrix >, true >());

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

    class_<Vector, boost::shared_ptr<Vector> >( "Vector").
            def(init<int>()).
            def("get", vector_getitem_1(&Vector::get)).
            def("get", vector_getitem_2(&Vector::get)).
            def("set", vector_setitem_1(&Vector::set)).
            def("set", vector_setitem_2(&Vector::set)).
            def("print_out", &Vector::print_out).
            def("scale", &Vector::scale).
            def("dim", &Vector::dim).
            def("__getitem__", vector_getitem_1(&Vector::pyget)).
            def("__setitem__", vector_setitem_1(&Vector::pyset)).
            def("__getitem__", vector_getitem_n(&Vector::pyget)).
            def("__setitem__", vector_setitem_n(&Vector::pyset)).
            def("nirrep", &Vector::nirrep);

    typedef void  (IntVector::*int_vector_set)(int, int, int);
    class_<IntVector, boost::shared_ptr<IntVector> >( "IntVector").
            def(init<int>()).
            def("get", &IntVector::get).
            def("set", int_vector_set(&IntVector::set)).
            def("print_out", &IntVector::print_out).
            def("dim", &IntVector::dim).
            def("nirrep", &IntVector::nirrep);

    enum_<Matrix::DiagonalizeOrder>("DiagonalizeOrder")
            .value("Ascending", Matrix::Ascending)
            .value("Descending", Matrix::Descending)
            .export_values();

    typedef void   (Matrix::*matrix_multiply)(bool, bool, double, const SharedMatrix&, const SharedMatrix&, double);
    typedef void   (Matrix::*matrix_diagonalize)(SharedMatrix&, boost::shared_ptr<Vector>&, Matrix::DiagonalizeOrder);
    typedef void   (Matrix::*matrix_one)(const SharedMatrix&);
    typedef double (Matrix::*double_matrix_one)(const SharedMatrix&);
    typedef void   (Matrix::*matrix_two)(const SharedMatrix&, const SharedMatrix&);
    typedef void   (Matrix::*matrix_save)(const std::string&, bool, bool, bool);
    typedef void   (Matrix::*matrix_set4)(int, int, int, double);
    typedef void   (Matrix::*matrix_set3)(int, int, double);
    typedef double (Matrix::*matrix_get3)(const int&, const int&, const int&) const;
    typedef double (Matrix::*matrix_get2)(const int&, const int&) const;
    typedef void   (Matrix::*matrix_load)(const std::string&);

    class_<Matrix, SharedMatrix >("Matrix").
            def(init<int, int>()).
            def("set_name", &Matrix::set_name).
            def("name", &Matrix::name, return_value_policy<copy_const_reference>()).
            def("print_out", &Matrix::print_out).
            def("rows", &Matrix::rowdim).
            def("cols", &Matrix::coldim).
            def("nirrep", &Matrix::nirrep, return_value_policy<copy_const_reference>()).
            def("symmetry", &Matrix::symmetry, return_value_policy<copy_const_reference>()).
            def("identity", &Matrix::identity).
            def("copy_lower_to_upper", &Matrix::copy_lower_to_upper).
            def("copy_upper_to_lower", &Matrix::copy_upper_to_lower).
            def("zero_lower", &Matrix::zero_lower).
            def("zero_upper", &Matrix::zero_upper).
            def("zero", &Matrix::zero).
            def("zero_diagonal", &Matrix::zero_diagonal).
            def("trace", &Matrix::trace).
            //            def("transpose", &Matrix::transpose).
            def("add", matrix_one(&Matrix::add)).
            def("subtract", matrix_one(&Matrix::subtract)).
            def("accumulate_product", matrix_two(&Matrix::accumulate_product)).
            def("scale", &Matrix::scale).
            def("sum_of_squares", &Matrix::sum_of_squares).
            def("rms", &Matrix::rms).
            def("scale_row", &Matrix::scale_row).
            def("scale_column", &Matrix::scale_column).
            def("transform", matrix_one(&Matrix::transform)).
            def("transform", matrix_two(&Matrix::transform)).
            def("transform", matrix_one(&Matrix::back_transform)).
            def("back_transform", matrix_two(&Matrix::back_transform)).
            def("vector_dot", double_matrix_one(&Matrix::vector_dot)).
            def("gemm", matrix_multiply(&Matrix::gemm)).
            def("diagonalize", matrix_diagonalize(&Matrix::diagonalize)).
            def("cholesky_factorize", &Matrix::cholesky_factorize).
            def("partial_cholesky_factorize", &Matrix::partial_cholesky_factorize).
            def("invert", &Matrix::invert).
            def("power", &Matrix::power).
            def("exp", &Matrix::exp).
            def("get", matrix_get3(&Matrix::get)).
            def("get", matrix_get2(&Matrix::get)).
            def("set", matrix_set3(&Matrix::set)).
            def("set", matrix_set4(&Matrix::set)).
            def("set", &Matrix::set_by_python_list).
            def("project_out", &Matrix::project_out).
            def("__getitem__", &Matrix::pyget).
            def("__setitem__", &Matrix::pyset).
            def("save", matrix_save(&Matrix::save)).
            def("load", matrix_load(&Matrix::load)).
            def("remove_symmetry", &Matrix::remove_symmetry);

    typedef SharedMatrix (MatrixFactory::*create_shared_matrix)();
    typedef SharedMatrix (MatrixFactory::*create_shared_matrix_name)(const std::string&);

    class_<MatrixFactory, boost::shared_ptr<MatrixFactory> >("MatrixFactory").
            def("shared_object", &get_matrix_factory).
            staticmethod("shared_object").
            def("create_matrix", create_shared_matrix(&MatrixFactory::create_shared_matrix)).
            def("create_matrix", create_shared_matrix_name(&MatrixFactory::create_shared_matrix));

    class_<CdSalcList, boost::shared_ptr<CdSalcList>, boost::noncopyable>("CdSalcList", no_init).
            def("print_out", &CdSalcList::print).
            def("matrix", &CdSalcList::matrix);

    class_<MintsHelper, boost::shared_ptr<MintsHelper> >("MintsHelper").
            def("integrals", &MintsHelper::integrals).
            def("one_electron_integrals", &MintsHelper::one_electron_integrals).
            def("basisset", &MintsHelper::basisset).
            def("sobasisset", &MintsHelper::sobasisset).
            def("factory", &MintsHelper::factory).
            def("ao_overlap", &MintsHelper::ao_overlap).
            def("ao_kinetic", &MintsHelper::ao_kinetic).
            def("ao_potential", &MintsHelper::ao_potential).
            def("one_electron_integrals", &MintsHelper::one_electron_integrals).
            def("so_overlap", &MintsHelper::so_overlap).
            def("so_kinetic", &MintsHelper::so_kinetic).
            def("so_potential", &MintsHelper::so_potential).
            def("so_dipole", &MintsHelper::so_dipole).
            def("so_quadrupole", &MintsHelper::so_quadrupole).
            def("so_traceless_quadrupole", &MintsHelper::so_traceless_quadrupole).
            def("ao_nabla", &MintsHelper::ao_nabla).
            def("so_nabla", &MintsHelper::so_nabla).
            def("so_angular_momentum", &MintsHelper::so_angular_momentum).
            def("ao_angular_momentum", &MintsHelper::ao_angular_momentum).
            def("ao_eri", &MintsHelper::ao_eri).
            def("ao_erf_eri", &MintsHelper::ao_erf_eri).
            def("cdsalcs", &MintsHelper::cdsalcs).
            def("petite_list", &MintsHelper::petite_list).
            def("play", &MintsHelper::play);

    class_<FittingMetric, boost::shared_ptr<FittingMetric> >("FittingMetric").
            def("get_algorithm", &FittingMetric::get_algorithm).
            def("is_poisson", &FittingMetric::is_poisson).
            def("is_inverted", &FittingMetric::is_inverted).
            def("get_metric", &FittingMetric::get_metric).
            def("get_pivots", &FittingMetric::get_pivots).
            def("get_reverse_pivots", &FittingMetric::get_reverse_pivots).
            def("form_fitting_metric", &FittingMetric::form_fitting_metric).
            def("form_cholesky_inverse", &FittingMetric::form_cholesky_inverse).
            def("form_QR_inverse", &FittingMetric::form_QR_inverse).
            def("form_eig_inverse", &FittingMetric::form_eig_inverse).
            def("form_full_inverse", &FittingMetric::form_full_inverse);

    class_<PseudoTrial, boost::shared_ptr<PseudoTrial> >("PseudoTrial").
            def("getI", &PseudoTrial::getI).
            def("getIPS", &PseudoTrial::getIPS).
            def("getQ", &PseudoTrial::getQ).
            def("getR", &PseudoTrial::getR).
            def("getA", &PseudoTrial::getA);

    class_<Vector3>("Vector3").
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
            def("dot", &Vector3::dot).
            def("distance", &Vector3::distance).
            def("normalize", &Vector3::normalize).
            def("norm", &Vector3::norm).
            def("cross", &Vector3::cross).
            def("__str__", &Vector3::to_string).
            def("__getitem__", &Vector3::get);

    typedef void (SymmetryOperation::*intFunction)(int);
    typedef void (SymmetryOperation::*doubleFunction)(double);

    class_<SymmetryOperation>("SymmetryOperation").
            def(init<const SymmetryOperation& >()).
            def("trace", &SymmetryOperation::trace).
            def("zero", &SymmetryOperation::zero).
            def("operate", &SymmetryOperation::operate).
            def("transform", &SymmetryOperation::transform).
            def("unit", &SymmetryOperation::unit).
            def("E", &SymmetryOperation::E).
            def("i", &SymmetryOperation::i).
            def("sigma_xy", &SymmetryOperation::sigma_xy).
            def("sigma_yz", &SymmetryOperation::sigma_yz).
            def("sigma_xz", &SymmetryOperation::sigma_xz).
            //        def("sigma_yz", &SymmetryOperation::sigma_yz).
            def("rotate_n", intFunction(&SymmetryOperation::rotation)).
            def("rotate_theta", doubleFunction(&SymmetryOperation::rotation)).
            def("c2_x", &SymmetryOperation::c2_x).
            def("c2_y", &SymmetryOperation::c2_y).
            def("transpose", &SymmetryOperation::transpose);

    class_<PointGroup, boost::shared_ptr<PointGroup> >("PointGroup").
            def(init<const std::string&>()).
            def("symbol", &PointGroup::symbol);
            //def("origin", &PointGroup::origin).
//            def("set_symbol", &PointGroup::set_symbol);

    typedef void (Molecule::*matrix_set_geometry)(const Matrix &);

    class_<Molecule, boost::shared_ptr<Molecule> >("Molecule").
            def("set_geometry", matrix_set_geometry(&Molecule::set_geometry)).
            def("set_name", &Molecule::set_name).
            def("name", &Molecule::name).
            def("fix_orientation", &Molecule::set_orientation_fixed).
            def("fix_com", &Molecule::set_com_fixed).
            def("init_with_checkpoint", &Molecule::init_with_chkpt).
            def("save_to_checkpoint", &Molecule::save_to_chkpt).
            def("init_with_io", &Molecule::init_with_psio).
            def("add_atom", &Molecule::add_atom).
            def("natom", &Molecule::natom).
            def("multiplicity", &Molecule::multiplicity).
            def("nfragments", &Molecule::nfragments).
            def("print_out", &Molecule::print).
            def("save_xyz", &Molecule::save_xyz).
            def("save_string_xyz", &Molecule::save_string_xyz).
            def("update_geometry", &Molecule::update_geometry).
            def("Z", &Molecule::Z, return_value_policy<copy_const_reference>()).
            def("x", &Molecule::x).
            def("y", &Molecule::y).
            def("z", &Molecule::z).
            //def("xyz", &Molecule::xyz).
            def("center_of_mass", &Molecule::center_of_mass).
            def("translate", &Molecule::translate).
            def("move_to_com", &Molecule::move_to_com).
            def("mass", &Molecule::mass).
            def("symbol", &Molecule::symbol).
            def("label", &Molecule::label).
            def("charge", &Molecule::charge).
            def("molecular_charge", &Molecule::molecular_charge).
            def("extract_subsets", &Molecule::py_extract_subsets_1).
            def("extract_subsets", &Molecule::py_extract_subsets_2).
            def("extract_subsets", &Molecule::py_extract_subsets_3).
            def("extract_subsets", &Molecule::py_extract_subsets_4).
            def("extract_subsets", &Molecule::py_extract_subsets_5).
            def("extract_subsets", &Molecule::py_extract_subsets_6).
            def("activate_all_fragments", &Molecule::activate_all_fragments).
            def("deactivate_all_fragments", &Molecule::deactivate_all_fragments).
            def("set_active_fragments", &Molecule::set_active_fragments).
            def("set_active_fragment", &Molecule::set_active_fragment).
            def("set_ghost_fragments", &Molecule::set_ghost_fragments).
            def("set_ghost_fragment", &Molecule::set_ghost_fragment).
            def("atom_at_position", &Molecule::atom_at_position1).
            def("print_out", &Molecule::print).
            def("print_out_in_bohr", &Molecule::print_in_bohr).
            def("nuclear_repulsion_energy", &Molecule::nuclear_repulsion_energy).
            def("find_point_group", &Molecule::find_point_group).
            def("reset_point_group", &Molecule::reset_point_group).
            def("set_point_group", &Molecule::set_point_group).
            def("point_group", &Molecule::point_group).
            def("schoenflies_symbol", &Molecule::schoenflies_symbol).
            def("form_symmetry_information", &Molecule::form_symmetry_information).
            def("create_molecule_from_string", &Molecule::create_molecule_from_string).
            staticmethod("create_molecule_from_string").
            def("is_variable", &Molecule::is_variable).
            def("set_variable", &Molecule::set_variable).
            def("get_variable", &Molecule::get_variable).
            def("update_geometry", &Molecule::update_geometry).
            def("set_basis_all_atoms", &Molecule::set_basis_all_atoms).
            def("set_basis_by_symbol", &Molecule::set_basis_by_symbol).
            def("set_basis_by_label", &Molecule::set_basis_by_label).
            def("set_basis_by_number", &Molecule::set_basis_by_number);

    class_<PetiteList, boost::shared_ptr<PetiteList>, boost::noncopyable>("PetiteList", no_init).
            def("aotoso", &PetiteList::aotoso).
            def("sotoao", &PetiteList::sotoao).
            def("print", &PetiteList::print);

    class_<BasisSetParser, boost::shared_ptr<BasisSetParser>, boost::noncopyable>("BasisSetParser", no_init);
    class_<Gaussian94BasisSetParser, boost::shared_ptr<Gaussian94BasisSetParser>, bases<BasisSetParser> >("Gaussian94BasisSetParser");

    typedef void (BasisSet::*basis_print_out)() const;
    class_<BasisSet, boost::shared_ptr<BasisSet>, boost::noncopyable>("BasisSet", no_init).
            def("print_out", basis_print_out(&BasisSet::print)).
            def("print_detail_out", basis_print_out(&BasisSet::print_detail)).
            def("make_filename", &BasisSet::make_filename).
            staticmethod("make_filename").
            def("construct", &BasisSet::construct).
            staticmethod("construct").
            def("nbf", &BasisSet::nbf).
            def("nao", &BasisSet::nao).
            def("nprimitive", &BasisSet::nprimitive).
            def("nshell", &BasisSet::nshell).
            def("max_am", &BasisSet::max_am);

    class_<SOBasisSet, boost::shared_ptr<SOBasisSet>, boost::noncopyable>("SOBasisSet", no_init).
            def("petite_list", &SOBasisSet::petite_list);

    class_<ExternalPotential, boost::shared_ptr<ExternalPotential>, boost::noncopyable>("ExternalPotential").
            def("setName", &ExternalPotential::setName).
            def("addCharge", &ExternalPotential::addCharge).
            def("addBasis", &ExternalPotential::addBasis).
            def("clear", &ExternalPotential::clear).
            def("computePotentialMatrix", &ExternalPotential::computePotentialMatrix).
            def("print_out", &ExternalPotential::py_print);

    class_<DFChargeFitter, boost::shared_ptr<DFChargeFitter>, boost::noncopyable>("DFChargeFitter").
            def("setPrimary", &DFChargeFitter::setPrimary).
            def("setAuxiliary", &DFChargeFitter::setAuxiliary).
            def("setD", &DFChargeFitter::setD).
            def("d", &DFChargeFitter::d).
            def("fit", &DFChargeFitter::fit);

    class_<Wavefunction, boost::shared_ptr<Wavefunction>, boost::noncopyable>("Wavefunction", no_init).
            def("nso", &Wavefunction::nso).
            def("nmo", &Wavefunction::nmo).
            def("nirrep", &Wavefunction::nirrep).
            def("Ca", &Wavefunction::Ca).
            def("Cb", &Wavefunction::Cb).
            def("Fa", &Wavefunction::Fa).
            def("Fb", &Wavefunction::Fb).
            def("Da", &Wavefunction::Da).
            def("Db", &Wavefunction::Db).
            def("epsilon_a", &Wavefunction::epsilon_a).
            def("epsilon_b", &Wavefunction::epsilon_b).
            def("add_preiteration_callback", &Wavefunction::add_preiteration_callback).
            def("add_postiteration_callback", &Wavefunction::add_postiteration_callback).
            def("basisset", &Wavefunction::basisset).
            def("sobasisset", &Wavefunction::sobasisset).
            def("energy", &Wavefunction::reference_energy).
            def("gradient", &Wavefunction::gradient).
            def("frequencies", &Wavefunction::frequencies);

    class_<scf::HF, boost::shared_ptr<scf::HF>, bases<Wavefunction>, boost::noncopyable>("HF", no_init);
    class_<scf::RHF, boost::shared_ptr<scf::RHF>, bases<scf::HF, Wavefunction> >("RHF", no_init);

    class_<MoldenWriter, boost::shared_ptr<MoldenWriter> >("MoldenWriter", no_init).
            def(init<boost::shared_ptr<Wavefunction> >()).
            def("write", &MoldenWriter::write);

    class_<NBOWriter, boost::shared_ptr<NBOWriter> >("NBOWriter", no_init).
            def(init<boost::shared_ptr<Wavefunction> >()).
            def("write", &NBOWriter::write);

    class_<OperatorSymmetry, boost::shared_ptr<OperatorSymmetry> >("MultipoleSymmetry", no_init).
            def(init<int, const boost::shared_ptr<Molecule>&,
                const boost::shared_ptr<IntegralFactory>&,
                const boost::shared_ptr<MatrixFactory>&>()).
            def("create_matrices", &OperatorSymmetry::create_matrices);
}
