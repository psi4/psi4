#include "cppe_pe_calc_handler.hh"

#include <cppe/utils/pot_manipulation.hh>

using libcppe::CppeState;

namespace psi {


CppePeCalcHandler::CppePeCalcHandler(libcppe::PeOptions options) :  PeCalcHandler(options),
                                                                    cppe_state(CppeState(options)) {
  dm_count = 0;
  iteration = 0;
}

bool CppePeCalcHandler::initialize() {

  // libcppe::Molecule molecule;
  // const double* nuc_charges = aobasis.nucchg.begin();

  m_nbas = 0;

  // read the geometry in a.u.

  // potentials = cppe_interface.read_potfile(m_options.potfile);
  // potentials = libcppe::PotManipulator(potentials, molecule).manipulate(m_options);

  // cppe_state.set_molecule(molecule);
  // cppe_state.set_potentials(potentials);

  return true;
}

void CppePeCalcHandler::fock_contribution(arma::mat& Ptot,
  arma::mat& out, double* energy) {

  // auto start_fock = std::chrono::high_resolution_clock::now();
  P_gs = Ptot;

  if (iteration == 0) {
    
  }
  cppe_state.update_energies(Ptot);


  // Calculate Electric Field Integrals
  size_t n_sitecoords = 3*cppe_state.get_polarizable_site_number();
  arma::Col<double> elec_fields(n_sitecoords, arma::fill::zeros);

  if (n_sitecoords != 0) {
    // std::cout << "calculating field integrals" << std::endl;
  }

  arma::mat corr = cppe_state.es_operator_copy(); //+ v_ind;

  *energy = cppe_state.get_current_energies().get_total_energy();

  // std::cout << "Total PE energy: " << *energy << std::endl;

  // if (out.get_n_arrays() == 1) {
  //   F_PE += corr;
  // }
  // else {
  //   throw std::runtime_error("PE not implemented for unrestricted calculations.");
  // }
  iteration++;
  // TODO: for debugging
  // cppe_state.print_summary();
}

void CppePeCalcHandler::print_energy() {
  cppe_state.print_summary();
}


// void CppePeCalcHandler::perturbative_exc_energy_correction(double* p_exc, double* energy) {
//   arma::mat Ptot_exc(p_exc, m_nbas, m_nbas, false, true);
//   // Calculate Electric Field Integrals
//   size_t n_sitecoords = 3*cppe_state.get_polarizable_site_number();
//   arma::Col<double> elec_fields(n_sitecoords, arma::fill::zeros);
//   array_view<double> av_fields(elec_fields.memptr(), n_sitecoords);
//   array_view<double> av_p(Ptot_exc.memptr(), m_nbas*m_nbas);
// 
//   // std::cout << "calculating field integrals" << std::endl;
//   // auto start_qints_fields = std::chrono::high_resolution_clock::now();
//   pe_integrals<double>().make_field(aobasis.b1, potentials, av_p, av_fields);
//   // auto end_qints_fields = std::chrono::high_resolution_clock::now();
//   // std::chrono::duration<double> elapsed_qints_fields = end_qints_fields - start_qints_fields;
//   // std::cout << "Time spent for integrals: " << elapsed_qints_fields.count() << std::endl;
//   // std::cout << "Electric fields" << std::endl;
//   // std::cout << elec_fields << std::endl;
// 
//   auto start_inds = std::chrono::high_resolution_clock::now();
//   cppe_state.update_induced_moments(elec_fields, iteration, true);
//   arma::vec induced_moments = cppe_state.get_induced_moments();
//   // std::cout << "induced_moments" << std::endl;
//   // std::cout << induced_moments << std::endl;
//   auto end_inds = std::chrono::high_resolution_clock::now();
//   std::chrono::duration<double> elapsed_inds = end_inds - start_inds;
//   // std::cout << "Time spent for induced moments: " << elapsed_inds.count() << std::endl;
// 
//   auto start_ind_integrals = std::chrono::high_resolution_clock::now();
//   arma::mat v_ind(m_nbas, m_nbas, arma::fill::zeros);
//   array_view<double> av_vind(v_ind.memptr(), m_nbas * m_nbas);
//   pe_integrals<double>().make_ind_integrals(aobasis.b1, potentials, av_vind, induced_moments);
//   auto end_ind_integrals = std::chrono::high_resolution_clock::now();
//   std::chrono::duration<double> elapsed_ind_integrals = end_ind_integrals - start_ind_integrals;
//   // std::cout << "Time spent for field integrals: " << elapsed_ind_integrals.count() << std::endl;
// 
//   *energy = cppe_state.get_current_energies().get("Polarization/Electronic");
// 
//   // QFIT output
//   if (m_options.print_level >= 2) {
//     Ptot_exc.save("P_exc_" + std::to_string(dm_count) + ".txt", arma::raw_ascii);
//   }
//   // induced_moments.save("ind_exc_" + std::to_string(dm_count) + ".txt", arma::raw_ascii);
//   // elec_fields.save("f_exc_"  + std::to_string(dm_count) + ".txt", arma::raw_ascii);
//   dm_count++;
// }


}
