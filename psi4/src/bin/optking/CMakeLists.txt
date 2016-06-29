set(sources_list molecule_sd_step.cc geom_gradients_io.cc lindh_guess.cc frag.cc print.cc getIntcoFileName.cc cart.cc linear_algebra.cc interfrag_orient.cc frag_apply_frozen_constraints.cc fb_frag.cc molecule_fragments.cc molecule_print.cc molecule_linesearch_step.cc mem.cc bend.cc frag_print.cc optking.cc molecule_read_coords.cc combo_coordinates.cc frag_H_guess.cc stre.cc molecule_irc_step.cc molecule_prfo_step.cc tors.cc opt_data_io.cc v3d.cc molecule_nr_step.cc molecule_rfo_step.cc molecule_tests.cc atom_data.cc molecule.cc set_params.cc molecule_backstep.cc opt_data.cc interfrag.cc frag_disp.cc frag_natural.cc oofp.cc print.cc)
psi4_add_binary(optking sources_list)
target_link_libraries(optking mints)

