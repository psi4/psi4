#include <string>
#include <boost/python.hpp>
#include <libmints/mints.h>
#include <lib3index/3index.h>
#include <libscf_solver/hf.h>
#include <libscf_solver/rhf.h>

using namespace boost;
using namespace boost::python;
using namespace psi;

shared_ptr<MatrixFactory> get_matrix_factory()
{
    // We need a valid molecule with a valid point group to create a matrix factory.
    shared_ptr<Molecule> molecule = Process::environment.molecule();
    if (!molecule) {
        fprintf(outfile, "  Active molecule not set!");
        throw PSIEXCEPTION("Active molecule not set!");
    }
    if (!molecule->point_group()) {
        fprintf(outfile, "  Active molecule does not have point group set!");
        throw PSIEXCEPTION("Active molecule does not have point group set!");
    }

    // Read in the basis set
    shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser(Process::environment.options.get_str("BASIS_PATH")));
    shared_ptr<BasisSet> basis = BasisSet::construct(parser, molecule, Process::environment.options.get_str("BASIS"));
    shared_ptr<IntegralFactory> fact(new IntegralFactory(basis, basis, basis, basis));
    shared_ptr<SOBasisSet> sobasis(new SOBasisSet(basis, fact));
    const Dimension& dim = sobasis->dimension();

    shared_ptr<MatrixFactory> matfac(new MatrixFactory);
    matfac->init_with(dim, dim);

    return matfac;
}

void export_mints()
{
    typedef void (Vector::*vector_set)(int, int, double);
    class_<Vector, shared_ptr<Vector> >( "Vector").
            def(init<int>()).
            def("get", &Vector::get).
            def("set", vector_set(&Vector::set)).
            def("print_out", &Vector::print_out).
            def("dim", &Vector::dim).
            def("nirrep", &Vector::nirrep);

    class_<IntVector, shared_ptr<IntVector> >( "IntVector").
            def(init<int>()).
            def("get", &IntVector::get).
            def("set", &IntVector::set_python).
            def("print_out", &IntVector::print_out).
            def("dim", &IntVector::dim).
            def("nirrep", &IntVector::nirrep);

    typedef void   (Matrix::*matrix_multiply)(bool, bool, double, const boost::shared_ptr<Matrix>&, const boost::shared_ptr<Matrix>&, double);
    typedef void   (Matrix::*matrix_diagonalize)(boost::shared_ptr<Matrix>&, boost::shared_ptr<Vector>&);
    typedef void   (Matrix::*matrix_one)(const boost::shared_ptr<Matrix>&);
    typedef double (Matrix::*double_matrix_one)(const boost::shared_ptr<Matrix>&);
    typedef void   (Matrix::*matrix_two)(const boost::shared_ptr<Matrix>&, const boost::shared_ptr<Matrix>&);
    typedef void   (Matrix::*matrix_save)(const std::string&, bool, bool, bool);
    typedef void   (Matrix::*matrix_set)(int, int, int, double);

    class_<Matrix, shared_ptr<Matrix> >("Matrix").
            def(init<int, int>()).
            def("set_name", &Matrix::set_name).
            def("name", &Matrix::name).
            def("print_out", &Matrix::print_out).
            def("rows", &Matrix::rowdim).
            def("cols", &Matrix::coldim).
            def("nirrep", &Matrix::nirrep).
            def("identity", &Matrix::identity).
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
            def("invert", &Matrix::invert).
            def("get", &Matrix::get).
            def("set", matrix_set(&Matrix::set)).
            def("__getitem__", &Matrix::pyget).
            def("__setitem__", &Matrix::pyset).
            def("save", matrix_save(&Matrix::save));

    typedef shared_ptr<Matrix> (MatrixFactory::*create_shared_matrix)();
    typedef shared_ptr<Matrix> (MatrixFactory::*create_shared_matrix_name)(std::string);

    class_<MatrixFactory, shared_ptr<MatrixFactory> >("MatrixFactory").
            def("shared_object", &get_matrix_factory).
            staticmethod("shared_object").
            def("create_matrix", create_shared_matrix(&MatrixFactory::create_shared_matrix)).
            def("create_matrix", create_shared_matrix_name(&MatrixFactory::create_shared_matrix));

    class_<MintsHelper, shared_ptr<MintsHelper> >("MintsHelper").
            def("ao_overlap", &MintsHelper::ao_overlap).
            def("ao_kinetic", &MintsHelper::ao_kinetic).
            def("ao_potential", &MintsHelper::ao_potential).
            def("so_overlap", &MintsHelper::so_overlap).
            def("so_kinetic", &MintsHelper::so_kinetic).
            def("so_potential", &MintsHelper::so_potential).
            def("ao_eri", &MintsHelper::ao_eri).
            def("ao_erf_eri", &MintsHelper::ao_erf_eri);

    class_<FittingMetric, shared_ptr<FittingMetric> >("FittingMetric").
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
            def("sigma_h", &SymmetryOperation::sigma_h).
            def("sigma_xz", &SymmetryOperation::sigma_xz).
            //        def("sigma_yz", &SymmetryOperation::sigma_yz).
            def("rotate_n", intFunction(&SymmetryOperation::rotation)).
            def("rotate_theta", doubleFunction(&SymmetryOperation::rotation)).
            def("c2_x", &SymmetryOperation::c2_x).
            def("c2_y", &SymmetryOperation::c2_y).
            def("transpose", &SymmetryOperation::transpose);

    class_<PointGroup, shared_ptr<PointGroup> >("PointGroup").
            def(init<const char*>()).
            def("symbol", &PointGroup::symbol).
            //def("origin", &PointGroup::origin).
            def("set_symbol", &PointGroup::set_symbol);

    class_<Molecule, shared_ptr<Molecule> >("Molecule").
            def("set_name", &Molecule::set_name).
            def("name", &Molecule::name).
            def("fix_orientation", &Molecule::set_orientation_fixed).
            def("init_with_checkpoint", &Molecule::init_with_chkpt).
            def("save_to_checkpoint", &Molecule::save_to_chkpt).
            def("init_with_io", &Molecule::init_with_psio).
            def("add_atom", &Molecule::add_atom).
            def("natom", &Molecule::natom).
            def("multiplicity", &Molecule::multiplicity).
            def("nfragments", &Molecule::nfragments).
            def("print_out", &Molecule::print).
            def("update_geometry", &Molecule::update_geometry).
            def("Z", &Molecule::Z).
            def("x", &Molecule::x).
            def("y", &Molecule::y).
            def("z", &Molecule::z).
            //def("xyz", &Molecule::xyz).
            def("center_of_mass", &Molecule::center_of_mass).
            def("translate", &Molecule::translate).
            def("move_to_com", &Molecule::move_to_com).
            def("mass", &Molecule::mass).
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
            def("nuclear_repulsion_energy", &Molecule::nuclear_repulsion_energy).
//            def("reorient", &Molecule::reorient).
            def("find_point_group", &Molecule::find_point_group).
            def("set_point_group", &Molecule::set_point_group).
            def("point_group", &Molecule::point_group).
            def("schoenflies_symbol", &Molecule::schoenflies_symbol).
            def("form_symmetry_information", &Molecule::form_symmetry_information).
            def("create_molecule_from_string", &Molecule::create_molecule_from_string).
            staticmethod("create_molecule_from_string").
            def("is_variable", &Molecule::is_variable).
            def("set_variable", &Molecule::set_variable).
            def("get_variable", &Molecule::get_variable).
            def("update_geometry", &Molecule::update_geometry);

    class_<PetiteList, shared_ptr<PetiteList>, boost::noncopyable>("PetiteList", no_init).
            def("aotoso", &PetiteList::aotoso).
            def("sotoao", &PetiteList::sotoao);

    typedef void (BasisSet::*basis_print_out)() const;
    class_<BasisSet, shared_ptr<BasisSet>, boost::noncopyable>("BasisSet", no_init).
            def("print_out", basis_print_out(&BasisSet::print)).
            def("print_detail_out", basis_print_out(&BasisSet::print_detail));

    class_<SOBasisSet, shared_ptr<SOBasisSet>, boost::noncopyable>("SOBasisSet", no_init).
            def("petite_list", &SOBasisSet::petitelist);

    class_<Wavefunction, shared_ptr<Wavefunction>, boost::noncopyable>("Wavefunction", no_init).
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
            def("sobasisset", &Wavefunction::sobasisset);

    class_<scf::HF, shared_ptr<scf::HF>, bases<Wavefunction>, boost::noncopyable>("HF", no_init);
    class_<scf::RHF, shared_ptr<scf::RHF>, bases<scf::HF> >("RHF", no_init);

    class_<MoldenWriter, shared_ptr<MoldenWriter> >("MoldenWriter", no_init).
            def(init<shared_ptr<Wavefunction> >()).
            def("write", &MoldenWriter::write);
}
