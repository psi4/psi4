/*-
 * Copyright (c) 2012-2013 Ilya Kaliman
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

#ifndef LIBEFP_EFP_H
#define LIBEFP_EFP_H

/** \file efp.h
 * Public libefp interface.
 *
 * A note on units: masses are in AMU, everything else is in atomic units.
 */

#ifdef __cplusplus
extern "C" {
#endif

/** Version string. */
#define LIBEFP_VERSION_STRING "0.9.9-beta"

/** Result of an operation. */
enum efp_result {
	/** Operation was successful. */
	EFP_RESULT_SUCCESS = 0,
	/** Insufficient memory. */
	EFP_RESULT_NO_MEMORY,
	/** Invalid argument to function was specified. */
	EFP_RESULT_INVALID_ARGUMENT,
	/** EFP structure was not properly initialized. */
	EFP_RESULT_NOT_INITIALIZED,
	/** File not found on disk. */
	EFP_RESULT_FILE_NOT_FOUND,
	/** Syntax error in EFP parameters file. */
	EFP_RESULT_SYNTAX_ERROR,
	/** Unknown fragment type. */
	EFP_RESULT_UNKNOWN_FRAGMENT,
	/** EFP parameters contain duplicate fragments. */
	EFP_RESULT_DUPLICATE_PARAMETERS,
	/** Required callback function was not set. */
	EFP_RESULT_CALLBACK_NOT_SET,
	/** Call to callback function failed. */
	EFP_RESULT_CALLBACK_FAILED,
	/** Gradient computation was not requested. */
	EFP_RESULT_GRADIENT_NOT_REQUESTED,
	/**
	 * Periodic boundary conditions are not supported for current selection
	 * of energy terms. */
	EFP_RESULT_PBC_NOT_SUPPORTED,
	/** Interaction cutoff must be enabled. */
	EFP_RESULT_PBC_REQUIRES_CUTOFF,
	/** Switching function cutoff is too small. */
	EFP_RESULT_SWF_CUTOFF_TOO_SMALL,
	/** Periodic simulation box is too small. */
	EFP_RESULT_BOX_TOO_SMALL,
	/** Fragment must contain at least three atoms. */
	EFP_RESULT_NEED_THREE_ATOMS,
	/** Polarization SCF did not converge. */
	EFP_RESULT_POL_NOT_CONVERGED,
	/** Certain EFP parameters are missing. */
	EFP_RESULT_PARAMETERS_MISSING,
	/** Invalid rotation matrix specified. */
	EFP_RESULT_INVALID_ROTATION_MATRIX,
	/** Index is out of range. */
	EFP_RESULT_INDEX_OUT_OF_RANGE,
	/** Wrong array length. */
	EFP_RESULT_INVALID_ARRAY_SIZE,
	/** Unsupported SCREEN group in EFP parameters file. */
	EFP_RESULT_UNSUPPORTED_SCREEN,
	/**
	 * Inconsistent selection of EFP terms.
	 *
	 * This means that AI/EFP terms were selected without selecting their
	 * EFP/EFP counterparts. Enabling polarization without electrostatics
	 * also produces this error. */
	EFP_RESULT_INCONSISTENT_TERMS
};

/** Flags to specify EFP energy terms. */
enum efp_term {
	/** EFP/EFP electrostatics. */
	EFP_TERM_ELEC = 1 << 0,
	/** EFP/EFP polarization. */
	EFP_TERM_POL = 1 << 1,
	/** EFP/EFP dispersion. */
	EFP_TERM_DISP = 1 << 2,
	/** EFP/EFP exchange repulsion. */
	EFP_TERM_XR = 1 << 3,
	/** EFP/EFP charge transfer, reserved for future use. */
	EFP_TERM_CHTR = 1 << 4,
	/** Ab initio/EFP electrostatics. */
	EFP_TERM_AI_ELEC = 1 << 5,
	/** Ab initio/EFP polarization. */
	EFP_TERM_AI_POL = 1 << 6,
	/** Ab initio/EFP dispersion, reserved for future use. */
	EFP_TERM_AI_DISP = 1 << 7,
	/** Ab initio/EFP exchange repulsion, reserved for future use. */
	EFP_TERM_AI_XR = 1 << 8,
	/** Ab initio/EFP charge transfer, reserved for future use. */
	EFP_TERM_AI_CHTR = 1 << 9
};

/** Fragment-fragment dispersion damping type. */
enum efp_disp_damp {
	EFP_DISP_DAMP_OVERLAP = 0, /**< Overlap-based damping (default). */
	EFP_DISP_DAMP_TT,          /**< Tang-Toennies damping. */
	EFP_DISP_DAMP_OFF          /**< No dispersion damping. */
};

/** Fragment-fragment electrostatic damping type. */
enum efp_elec_damp {
	EFP_ELEC_DAMP_SCREEN = 0,  /**< SCREEN-controlled damping (default). */
	EFP_ELEC_DAMP_OVERLAP,     /**< Overlap-based damping. */
	EFP_ELEC_DAMP_OFF          /**< No electrostatic damping. */
};

/** Fragment-fragment polarization damping type. */
enum efp_pol_damp {
	EFP_POL_DAMP_TT = 0,       /**< Tang-Toennies like damping (default). */
	EFP_POL_DAMP_OFF           /**< No polarization damping. */
};

/** Describes the way fragment coordinates are specified. */
enum efp_coord_type {
	/** Coordinates of center of mass of a fragment and Euler angles. */
	EFP_COORD_TYPE_XYZABC = 0,
	/** Coordinates of three points belonging to a fragment. */
	EFP_COORD_TYPE_POINTS,
	/** Coordinates of the center of mass of a fragment and rotation matrix. */
	EFP_COORD_TYPE_ROTMAT
};

/** \struct efp
 * Main EFP opaque structure.
 */
struct efp;

/** Options controlling EFP computation. */
struct efp_opts {
	/** Specifies which energy terms to compute. */
	unsigned terms;

	/** Dispersion damping type (see #efp_disp_damp). */
	enum efp_disp_damp disp_damp;

	/** Electrostatic damping type (see #efp_elec_damp). */
	enum efp_elec_damp elec_damp;

	/** Polarization damping type (see #efp_pol_damp). */
	enum efp_pol_damp pol_damp;

	/** Enable periodic boundary conditions if nonzero. */
	int enable_pbc;

	/** Enable fragment-fragment interaction cutoff if nonzero. */
	int enable_cutoff;

	/** Cutoff distance for fragment-fragment interactions. */
	double swf_cutoff;
};

/** EFP energy terms. */
struct efp_energy {
	/**
	 * EFP/EFP electrostatic energy. */
	double electrostatic;
	/**
	 * Charge penetration energy from overlap-based electrostatic
	 * damping. Zero if overlap-based damping is turned off. */
	double charge_penetration;
	/**
	 * Interaction energy of EFP electrostatics with point charges. */
	double electrostatic_point_charges;
	/**
	 * All polarization energy goes here. Polarization is computed
	 * self-consistently so it can't be separated into EFP/EFP and AI/EFP
	 * parts. */
	double polarization;
	/**
	 * EFP/EFP dispersion energy. */
	double dispersion;
	/**
	 * EFP/EFP exchange-repulsion energy. */
	double exchange_repulsion;
	/**
	 * Sum of all the above energy terms. */
	double total;
};

/** EFP atom info. */
struct efp_atom {
	char label[32];   /**< Atom label. */
	double x;         /**< X coordinate of atom position. */
	double y;         /**< Y coordinate of atom position. */
	double z;         /**< Z coordinate of atom position. */
	double mass;      /**< Atom mass. */
	double znuc;      /**< Nuclear charge. */
};

/**
 * Callback used to compute electric field in the specified points.
 *
 * This can be used to provide electric field from electrons in \a ab \a initio
 * subsystem. This callback is used during polarization calculation.
 *
 * \param[in] n_pt Number of points in \p xyz array.
 *
 * \param[in] xyz Coordinates of points where electric field should be
 * computed. The size of this array must be [3 * \p n_pt] elements.
 *
 * \param[out] field Computed \a x \a y \a z components of electric field. The
 * size of this array must be at least [3 * \p n_pt] elements.
 *
 * \param[in] user_data User data which was specified during initialization.
 *
 * \return ::EFP_RESULT_CALLBACK_FAILED on error or ::EFP_RESULT_SUCCESS
 *         otherwise.
 */
typedef enum efp_result (*efp_electron_density_field_fn)(
				int n_pt,
				const double *xyz,
				double *field,
				void *user_data);

/**
 * Get a human readable banner string with information about the library.
 *
 * \return Banner string, zero-terminated.
 */
const char *efp_banner(void);

/**
 * Create a new efp object.
 *
 * \return A new efp object or NULL on error.
 */
struct efp *efp_create(void);

/**
 * Get default values of simulation options.
 *
 * \param[out] opts Structure to store the defaults. See ::efp_opts.
 */
void efp_opts_default(struct efp_opts *opts);

/**
 * Set computation options.
 *
 * \param[in] efp The efp structure.
 *
 * \param[in] opts New options for EFP computation. See ::efp_opts.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_set_opts(struct efp *efp, const struct efp_opts *opts);

/**
 * Get currently set computation options.
 *
 * \param[in] efp The efp structure.
 *
 * \param[out] opts Current options for EFP computation. See ::efp_opts.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_opts(struct efp *efp, struct efp_opts *opts);

/**
 * Add EFP potential from a file.
 *
 * \param[in] efp The efp structure.
 *
 * \param[in] path Path to the EFP potential file, zero terminated string.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_add_potential(struct efp *efp, const char *path);

/**
 * Add a new fragment to the EFP subsystem.
 *
 * \param[in] efp The efp structure.
 *
 * \param[in] name Fragment name, zero terminated string.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_add_fragment(struct efp *efp, const char *name);

/**
 * Set the callback function which computes electric field from electrons
 * in \a ab \a initio subsystem.
 *
 * \param[in] efp The efp structure.
 *
 * \param[in] fn The callback function. See ::efp_electron_density_field_fn.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_set_electron_density_field_fn(
			struct efp *efp,
			efp_electron_density_field_fn fn);

/**
 * Set user data to be passed to ::efp_electron_density_field_fn.
 *
 * \param[in] efp The efp structure.
 *
 * \param[in] user_data User data which will be passed as a last parameter to
 * ::efp_electron_density_field_fn.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_set_electron_density_field_user_data(
			struct efp *efp,
			void *user_data);

/**
 * Setup arbitrary point charges interacting with EFP subsystem.
 *
 * This can be used to compute contributions from \a ab \a initio nuclei.
 *
 * \param[in] efp The efp structure.
 *
 * \param[in] n_ptc Number of point charges.
 *
 * \param[in] q Array of \p n_ptc elements with charge values.
 *
 * \param[in] xyz Array of [3 * \p n_ptc] elements with \a x \a y \a z
 * coordinates of charge positions.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_set_point_charges(
			struct efp *efp,
			int n_ptc,
			const double *q,
			const double *xyz);

/**
 * Get the number of currently set point charges.
 *
 * \param[in] efp The efp structure.
 * \param[out] n_ptc Number of point charges.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_point_charge_count(
			struct efp *efp,
			int *n_ptc);

/**
 * Get values of currently set point charges.
 *
 * \param[in] efp The efp structure.
 * \param[out] q Array of \p n_ptc elements where charges will be stored.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_point_charge_values(
			struct efp *efp,
			double *q);

/**
 * Get coordinates of currently set point charges.
 *
 * \param[in] efp The efp structure.
 *
 * \param[out] xyz Array of [3 * \p n_ptc] elements where \a x \a y \a z
 * coordinates of point charges will be stored.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_point_charge_coordinates(
			struct efp *efp,
			double *xyz);

/**
 * Get gradient on point charges from EFP subsystem.
 *
 * \param[in] efp The efp structure.
 *
 * \param[out] grad For each point charge \a x \a y \a z components of energy
 * gradient are stored. The size of this array must be at least [3 * \p n_ptc]
 * elements.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_point_charge_gradient(
			struct efp *efp,
			double *grad);

/**
 * Update positions and orientations of effective fragments.
 *
 * \param[in] efp The efp structure.
 *
 * \param[in] coord_type Specifies the type of coordinates in the \a coord
 * array (see #efp_coord_type).
 *
 * \param[in] coord Array of fragment coordinates.
 *
 * If \p coord_type is \a EFP_COORD_TYPE_XYZABC then for each fragment the \p
 * coord array should contain \a x \a y \a z components of the center of mass
 * position and three Euler rotation angles representing orientation of a
 * fragment. The size of the \p coord array must be at least [6 * \a n]
 * elements, where \a n is the number of fragments.
 *
 * If \p coord_type is \a EFP_COORD_TYPE_POINTS then for each fragment the \p
 * coord array should contain the coordinates of 3 points in space. For each
 * fragment point 1 and first atom of fragment are made to coincide. The vector
 * connecting points 1 and 2 is aligned with the corresponding vector
 * connecting fragment atoms. The plane defined by points 1, 2, and 3 is made
 * to coincide with the corresponding fragment plane. The size of the \p coord
 * array must be at least [9 * \a n] elements, where \a n is the number of
 * fragments.
 *
 * If \p coord_type is \a EFP_COORD_TYPE_ROTMAT then for each fragment the \p
 * coord array should contain \a x \a y \a z components of the center of mass
 * position and nine elements of the rotation matrix representing orientation
 * of a fragment. The size of the \p coord array must be at least [12 * \a n]
 * elements, where \a n is the number of fragments.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_set_coordinates(struct efp *efp,
				    enum efp_coord_type coord_type,
				    const double *coord);

/**
 * Update position and orientation of the specified effective fragment.
 *
 * \param[in] efp The efp structure.
 *
 * \param[in] frag_idx Index of a fragment. Must be a value between zero and
 * the total number of fragments minus one.
 *
 * \param[in] coord_type Specifies the type of coordinates in the \a coord
 * array (see #efp_coord_type).
 *
 * \param[in] coord Array of coordinates specifying fragment position and
 * orientation.
 *
 * If \p coord_type is \a EFP_COORD_TYPE_XYZABC then the \p coord array should
 * contain \a x \a y \a z components of the center of mass position and three
 * Euler rotation angles representing orientation of a fragment. The \p coord
 * array must contain a total of 6 elements.
 *
 * If \p coord_type is \a EFP_COORD_TYPE_POINTS then the \p coord array should
 * contain the coordinates of 3 points in space. Point 1 and first atom of
 * fragment are made to coincide. The vector connecting points 1 and 2 is
 * aligned with the corresponding vector connecting fragment atoms. The plane
 * defined by points 1, 2, and 3 is made to coincide with the corresponding
 * fragment plane. The \p coord array must contain a total of 9 elements.
 *
 * If \p coord_type is \a EFP_COORD_TYPE_ROTMAT then the \p coord array should
 * contain \a x \a y \a z components of the center of mass position and nine
 * elements of the rotation matrix representing orientation of a fragment. The
 * \p coord array must contain a total of 12 elements.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_set_frag_coordinates(struct efp *efp,
					 int frag_idx,
					 enum efp_coord_type coord_type,
					 const double *coord);

/**
 * Get center of mass positions and Euler angles of the effective fragments.
 *
 * \param[in] efp The efp structure.
 *
 * \param[in] n_frags Expected number of fragments.
 *
 * \param[out] xyzabc Upon return the coordinates of the center of mass and
 * Euler rotation angles for each fragment will be written to this array. The
 * size of the \p xyzabc array must be at least [6 * \p n_frags] elements.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_coordinates(struct efp *efp, int n_frags, double *xyzabc);

/**
 * Setup periodic box size.
 *
 * \param[in] efp The efp structure.
 * \param[in] x Box size in x dimension.
 * \param[in] y Box size in y dimension.
 * \param[in] z Box size in z dimension.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_set_periodic_box(struct efp *efp, double x, double y, double z);

/**
 * Get the stress tensor.
 *
 * \param[in] efp The efp structure.
 *
 * \param[out] stress Array of 9 elements where the stress tensor will be
 * stored.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_stress_tensor(struct efp *efp, double *stress);

/**
 * Update wave function dependent energy terms.
 *
 * This function must be called during \a ab \a initio SCF.
 *
 * \param[in] efp The efp structure.
 * \param[out] energy Wave function dependent EFP energy.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_wavefunction_dependent_energy(
			struct efp *efp,
			double *energy);

/**
 * Perform the EFP computation.
 *
 * \param[in] efp The efp structure.
 *
 * \param[in] do_gradient If nonzero value is specified in addition to energy
 * compute the gradient.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_compute(struct efp *efp, int do_gradient);

/**
 * Get total charge of a fragment.
 *
 * \param[in] efp The efp structure.
 *
 * \param[in] frag_idx Index of a fragment. Must be a value between zero and
 * the total number of fragments minus one.
 *
 * \param[out] charge Total charge of a fragment.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_frag_charge(struct efp *efp, int frag_idx, double *charge);

/**
 * Get spin multiplicity of a fragment.
 *
 * \param[in] efp The efp structure.
 *
 * \param[in] frag_idx Index of a fragment. Must be a value between zero and
 * the total number of fragments minus one.
 *
 * \param[out] mult Spin multiplicity of a fragment.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_frag_multiplicity(struct efp *efp, int frag_idx, int *mult);

/**
 * Get total number of multipoles from EFP electrostatics.
 *
 * \param[in] efp The efp structure.
 * \param[out] n_mult Number of electrostatics multipoles.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_multipole_count(struct efp *efp, int *n_mult);

/**
 * Get coordinates of electrostatics multipoles.
 *
 * \param[in] efp The efp structure.
 *
 * \param[out] xyz Array where coordinates of EFP electrostatics multipoles
 * will be stored. Size of the \p xyz array must be at least [3 * \p n_mult]
 * elements, where \p n_mult is the value returned by the
 * ::efp_get_multipole_count function.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_multipole_coordinates(
			struct efp *efp,
			double *xyz);

/**
 * Get electrostatics multipoles from EFP fragments.
 *
 * \param[in] efp The efp structure.
 *
 * \param[out] mult Array where charges, dipoles, quadrupoles, and octupoles
 * for each point will be stored.
 *
 * The size of the \p mult array must be at least [(1 + 3 + 6 + 10) * \p
 * n_mult] elements (charges + dipoles + quadrupoles + octupoles), where \p
 * n_mult is the value returned by the ::efp_get_multipole_count function.
 *
 * Quadrupoles are stored in the following order:
 *    \a xx, \a yy, \a zz, \a xy, \a xz, \a yz
 *
 * Octupoles are stored in the following order:
 *    \a xxx, \a yyy, \a zzz, \a xxy, \a xxz,
 *    \a xyy, \a yyz, \a xzz, \a yzz, \a xyz
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_multipole_values(struct efp *efp, double *mult);

/**
 *  Get the number of polarization induced dipoles.
 *
 * \param[in] efp The efp structure.
 * \param[out] n_dip Number of polarization induced dipoles.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_induced_dipole_count(struct efp *efp, int *n_dip);

/**
 * Get coordinates of induced dipoles.
 *
 * \param[in] efp The efp structure.
 *
 * \param[out] xyz Array where the coordinates of polarizable points will be
 * stored. The size of the \p xyz array must be at least [3 * \p n_dip]
 * elements.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_induced_dipole_coordinates(struct efp *efp, double *xyz);

/**
 * Get values of polarization induced dipoles.
 *
 * \param[in] efp The efp structure.
 *
 * \param[out] dip Array where induced dipoles will be stored. The size of the
 * array must be at least [3 * \p n_dip] elements.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_induced_dipole_values(struct efp *efp, double *dip);

/**
 * Get values of polarization conjugated induced dipoles.
 *
 * \param[in] efp The efp structure.
 *
 * \param[out] dip Array where induced dipoles will be stored. The size of the
 * array must be at least [3 * \p n_dip] elements.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_induced_dipole_conj_values(struct efp *efp, double *dip);

/**
 * Get computed energy components.
 *
 * \param[in] efp The efp structure.
 * \param[out] energy Computed EFP energy components (see efp_energy).
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_energy(struct efp *efp, struct efp_energy *energy);

/**
 * Get computed EFP energy gradient.
 *
 * \param[in] efp The efp structure.
 *
 * \param[in] n_frags Expected number of fragments.
 *
 * \param[out] grad For each fragment \a x \a y \a z components of negative
 * force and torque will be written to this array. The size of this array must
 * be at least [6 * \p n_frags] elements.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_gradient(struct efp *efp, int n_frags, double *grad);

/**
 * Get the number of fragments in this computation.
 *
 * \param[in] efp The efp structure.
 * \param[out] n_frag Total number of fragments in this simulation.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_frag_count(struct efp *efp, int *n_frag);

/**
 * Get the name of the specified effective fragment.
 *
 * \param[in] efp The efp structure.
 *
 * \param[in] frag_idx Index of a fragment. Must be a value between zero and
 * the total number of fragments minus one.
 *
 * \param[in] size Size of a \p frag_name buffer.
 *
 * \param[out] frag_name A buffer where name of the fragment will be stored.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_frag_name(struct efp *efp,
				  int frag_idx,
				  int size,
				  char *frag_name);

/**
 * Get total mass of a fragment.
 *
 * \param[in] efp The efp structure.
 *
 * \param[in] frag_idx Index of a fragment. This must be a value between zero
 * and the total number of fragments minus one.
 *
 * \param[out] mass_out Output mass value.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_frag_mass(struct efp *efp,
				  int frag_idx,
				  double *mass_out);

/**
 * Get fragment principal moments of inertia.
 *
 * \param[in] efp The efp structure.
 *
 * \param[in] frag_idx Index of a fragment. Must be a value between zero and
 * the total number of fragments minus one.
 *
 * \param[out] inertia_out Array of 3 elements where principal moments of
 * inertia of a fragment will be stored.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_frag_inertia(struct efp *efp,
				     int frag_idx,
				     double *inertia_out);

/**
 * Get the number of atoms in the specified fragment.
 *
 * \param[in] efp The efp structure.
 *
 * \param[in] frag_idx Index of a fragment. Must be a value between zero and
 * the total number of fragments minus one.
 *
 * \param[out] n_atoms Total number of atoms in the fragment.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_frag_atom_count(struct efp *efp,
					int frag_idx,
					int *n_atoms);

/**
 * Get atoms comprising the specified fragment.
 *
 * \param[in] efp The efp structure.
 *
 * \param[in] frag_idx Index of a fragment. Must be a value between zero and
 * the total number of fragments minus one.
 *
 * \param[in] size Size of the \p atoms array. Must be greater than or equal to
 * the size returned by the ::efp_get_frag_atom_count function.
 *
 * \param[out] atoms Array where atom information will be stored.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_frag_atoms(struct efp *efp,
				   int frag_idx,
				   int size,
				   struct efp_atom *atoms);

/**
 * Convert rigid body torque to derivatives of energy by Euler angles.
 *
 * \param[in] euler Array of 3 elements containing the values of Euler angles
 * for the current orientation of the rigid body.
 *
 * \param[in] torque Array of 3 elements containing torque components.
 *
 * \param[out] deriv Array of 3 elements where derivatives of energy by
 * Euler angles will be stored. This can point to the same memory location as
 * the \p torque argument.
 */
void efp_torque_to_derivative(const double *euler, const double *torque, double *deriv);

/**
 * Release all resources used by this \a efp.
 *
 * \param[in] efp The efp structure to be released.
 */
void efp_shutdown(struct efp *efp);

/**
 * Convert #efp_result to a human readable message.
 *
 * \param res Result value to be converted to string.
 *
 * \return Human readable string with the description of the result,
 *         zero-terminated.
 */
const char *efp_result_to_string(enum efp_result res);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* LIBEFP_EFP_H */
