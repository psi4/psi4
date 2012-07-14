/*-
 * Copyright (c) 2012 Ilya Kaliman
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
 */

/** \mainpage libefp - The Effective Fragment Potential method implementation
 *
 * \section API Public API Documentation
 *
 * The library is still in active development and API is UNSTABLE.
 * Documentation for public API is available <a href="efp_8h.html">here</a>.
 * For additional information see
 * <a href="http://github.com/libefp/libefp/blob/master/README.md">README</a>
 * file.
 *
 * \section Repo Git Repository
 *
 * Latest development version of code can be found in git
 * <a href="http://github.com/libefp/libefp">repository</a>.
 *
 * \copyright Copyright (c) 2012 Ilya Kaliman.
 * Distributed under the terms of BSD 2-clause license. See
 * <a href="http://github.com/libefp/libefp/blob/master/LICENSE">LICENSE</a>
 * file for details.
 */

#ifdef __cplusplus
extern "C" {
#endif

/** Result of an operation. */
enum efp_result {
	/** Operation was successful. */
	EFP_RESULT_SUCCESS = 0,
	/** Insufficient memory. */
	EFP_RESULT_NO_MEMORY,
	/** Operation is not implemented. */
	EFP_RESULT_NOT_IMPLEMENTED,
	/** Unexpected NULL argument to function was specified. */
	EFP_RESULT_ARGUMENT_NULL,
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
	/**
	 * Exchange-repulsion must be turned on.
	 *
	 * Exchange-repulsion must be enabled to compute overlap-based
	 * electrostatic and dispersion damping. Overlap integrals need to be
	 * computed for these damping types. */
	EFP_RESULT_OVERLAP_INTEGRALS_REQUIRED,
	/** Gradient computation was not requested. */
	EFP_RESULT_GRADIENT_NOT_REQUESTED,
	/** Certain EFP parameters are missing. */
	EFP_RESULT_PARAMETERS_MISSING,
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
	EFP_TERM_ELEC = 1 << 0,     /**< EFP/EFP electrostatics */
	EFP_TERM_POL = 1 << 1,      /**< EFP/EFP polarization */
	EFP_TERM_DISP = 1 << 2,     /**< EFP/EFP dispersion */
	EFP_TERM_XR = 1 << 3,       /**< EFP/EFP exchange repulsion */
	EFP_TERM_CHTR = 1 << 4,     /**< EFP/EFP charge transfer */
	EFP_TERM_AI_ELEC = 1 << 5,  /**< Ab initio/EFP electrostatics */
	EFP_TERM_AI_POL = 1 << 6,   /**< Ab initio/EFP polarization */
	EFP_TERM_AI_DISP = 1 << 7,  /**< Ab initio/EFP dispersion */
	EFP_TERM_AI_XR = 1 << 8,    /**< Ab initio/EFP exchange repulsion */
	EFP_TERM_AI_CHTR = 1 << 9   /**< Ab initio/EFP charge transfer */
};

/** Fragment-fragment dispersion damping type. */
enum efp_disp_damp {
	EFP_DISP_DAMP_OVERLAP = 0, /**< Overlap-based damping (default). */
	EFP_DISP_DAMP_TT,          /**< Tang-Toennies damping. */
	EFP_DISP_DAMP_OFF          /**< No dispersion damping. */
};

/** Fragment-fragment electrostatic damping type. */
enum efp_elec_damp {
	EFP_ELEC_DAMP_SCREEN = 0,  /**< SCREEN controlled damping (default). */
	EFP_ELEC_DAMP_OVERLAP,     /**< Overlap-based damping. */
	EFP_ELEC_DAMP_OFF          /**< No electrostatic damping. */
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
};

/** EFP energy terms. */
struct efp_energy {
	/**
	 * Electrostatic energy. */
	double electrostatic;
	/**
	 * Charge penetration energy from overlap-based electrostatic
	 * damping. Zero if overlap-based damping is turned off. */
	double charge_penetration;
	/**
	 * All polarization energy goes here. Polarization is computed
	 * self-consistently so we can't separate it into EFP/EFP and AI/EFP
	 * parts. */
	double polarization;
	/**
	 * Dispersion energy. */
	double dispersion;
	/**
	 * Exchange-repulsion energy. */
	double exchange_repulsion;
	/**
	 * Charge transfer energy. */
	double charge_transfer;
	/**
	 * AI/EFP electrostatic energy. */
	double ai_electrostatic;
	/**
	 * AI/EFP dispersion energy. */
	double ai_dispersion;
	/**
	 * AI/EFP exchange-repulsion energy. */
	double ai_exchange_repulsion;
	/**
	 * AI/EFP charge transfer energy. */
	double ai_charge_transfer;
	/**
	 * Sum of all the above EFP energy terms. */
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
 * Compute electric field form electrons in the specified points.
 *
 * \param[in] n_pt Number of points in \a xyz array.
 * \param[in] xyz Coordinates of points where electric field should be computed.
 * \param[out] field Computed \a x \a y \a z components of electric field.
 * \param[in] user_data User data which was specified during initialization.
 *
 * \return ::EFP_RESULT_CALLBACK_FAILED on error or ::EFP_RESULT_SUCCESS
 *         otherwise.
 */
typedef enum efp_result (*efp_electron_density_field_fn)(int n_pt,
							 const double *xyz,
							 double *field,
							 void *user_data);

/** Callback function information. */
struct efp_callbacks {
	/** Callback function, see ::efp_electron_density_field_fn. */
	efp_electron_density_field_fn get_electron_density_field;

	/** Will be passed as a last parameter to
	 * efp_callbacks::get_electron_density_field. */
	void *get_electron_density_field_user_data;
};

/**
 * Get a human readable banner string with information about the library.
 *
 * \return Banner string, zero-terminated.
 */
const char *efp_banner(void);

/**
 * Get default values of simulation options.
 *
 * \param[out] opts Structure to store the defaults.
 */
void efp_opts_default(struct efp_opts *opts);

/**
 * Initialize the EFP computation.
 *
 * \param[out] out Initialized efp structure.
 * \param[in] opts User defined options controlling the computation.
 * \param[in] callbacks User supplied callback functions (see efp_callbacks).
 * \param[in] potential_file_list NULL-terminated array with paths to the files
 *                                with EFP parameters.
 * \param[in] frag_name_list NULL-terminated array with names of the fragments
 *                           comprising the molecular system.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_init(struct efp **out,
			 const struct efp_opts *opts,
			 const struct efp_callbacks *callbacks,
			 const char **potential_file_list,
			 const char **frag_name_list);

/**
 * Update information about \a ab \a initio region.
 *
 * \param[in] efp The efp structure.
 * \param[in] n_atoms Number of atoms in \a ab \a initio region.
 * \param[in] znuc Array of \a n_atoms elements with atom nuclear charges.
 * \param[in] xyz Array of [3 * \a n_atoms] elements with \a x \a y \a z
 *                coordinates of atom positions.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_set_qm_atoms(struct efp *efp,
				 int n_atoms,
				 const double *znuc,
				 const double *xyz);

/**
 * Update positions and orientations of effective fragments.
 *
 * \param[in] efp The efp structure.
 * \param[in] xyzabc For each fragment specifies \a x \a y \a z components of
 *                   the center of mass position and three Euler rotation
 *                   angles representing orientation of a fragment.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_set_coordinates(struct efp *efp, const double *xyzabc);

/**
 * Update positions and orientations of effective fragments.
 *
 * This is a convenience function. It does the same as ::efp_set_coordinates.
 * However to specify position and orientation of each fragment it takes
 * coordinates of 3 points in space. For each fragment point 1 and first atom
 * of fragment are made to coincide. The vector connecting points 1 and 2 is
 * aligned with the corresponding vector connecting fragment atoms. The plane
 * defined by points 1, 2, and 3 is made to coincide with the corresponding
 * fragment plane.
 *
 * \param[in] efp The efp structure.
 * \param[in] pts Array of points used to determine fragment positions and
 *                orientations. The size of this array must be at least
 *                [9 * \a n] elements, where \a n is the number of fragments.
 *                For each fragment \a x \a y \a z coordinates of 3 fragment
 *                points should be specified to setup the position and
 *                orientation of a fragment.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_set_coordinates_2(struct efp *efp, const double *pts);

/**
 * Get center of mass positions and Euler angles of the effective fragments.
 *
 * \param[in] efp The efp structure.
 * \param[in] size Size of the xyzabc array. Must be at least [6 * \a n]
 *                 elements, where \a n is the number of fragments.
 * \param[out] xyzabc Upon return [6 * \a n] elements will be written to this
 *                    array. The coordinates of the center of mass and Euler
 *                    rotation angles for each fragment will be stored.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_coordinates(struct efp *efp, int size, double *xyzabc);

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
enum efp_result efp_scf_update(struct efp *efp, double *energy);

/**
 * Perform the EFP computation.
 *
 * This call can take long time to complete.
 *
 * \param[in] efp The efp structure.
 * \param[in] do_gradient Also compute energy gradient if nonzero value is
 *                        specified.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_compute(struct efp *efp, int do_gradient);

/**
 * Get total number of EFP multipoles.
 *
 * \param[in] efp The efp structure.
 * \param[out] n_mult Array of 4 integers where the total number of charges,
 *                    dipoles, quadrupoles and octupoles will be stored.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_multipole_count(struct efp *efp, int *n_mult);

/**
 * Get all electrostatic multipoles from all fragments.
 *
 * \param[in] efp The efp structure.
 *
 * \param[out] xyz Array of four pointers to arrays where the corresponding
 *                 multipole point coordinates will be stored. The size of
 *                 each array must be at least [3 * \a n_mult] elements, where
 *                 \a n_mult is the corresponding value in the array returned
 *                 by the ::efp_get_multipole_count function.
 *
 * \param[out] z Array of four pointers to arrays where the charges, dipoles,
 *               quadrupoles, and octupoles will be stored.
 *
 *               The size of the first array (charges) must be at least
 *               [\a n_mult] elements, where \a n_mult is the corresponding
 *               first element of the array returned by the
 *               ::efp_get_multipole_count function.
 *
 *               The size of the second array (dipoles) must be at least
 *               [3 * \a n_mult] elements, where \a n_mult is the corresponding
 *               second element of the array returned by the
 *               ::efp_get_multipole_count function.
 *
 *               The size of the third array (quadrupoles) must be at least
 *               [6 * \a n_mult] elements, where \a n_mult is the corresponding
 *               third element of the array returned by the
 *               ::efp_get_multipole_count function.
 *
 *               Quadrupoles are stored in the following order:
 *                   \a xx, \a yy, \a zz, \a xy, \a xz, \a yz
 *
 *               The size of the fourth array (octupoles) must be at least
 *               [10 * \a n_mult] elements, where \a n_mult is the
 *               corresponding fourth element of the array returned by the
 *               ::efp_get_multipole_count function.
 *
 *               Octupoles are stored in the following order:
 *                   \a xxx, \a yyy, \a zzz, \a xxy, \a xxz,
 *                   \a xyy, \a yyz, \a xzz, \a yzz, \a xyz
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_multipoles(struct efp *efp, double **xyz, double **z);

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
 * \param[in] size Size of the grad array. Must be at least [6 * \a n]
 *                 elements, where \a n is the number of fragments.
 * \param[out] grad For each fragment contains \a x \a y \a z components of
 *                  negative force and torque.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_gradient(struct efp *efp, int size, double *grad);

/**
 * Get the gradient on atoms in \a ab \a initio subsystem due to EFP
 * interaction.
 *
 * This does not include additional correction due to changes in the
 * \a ab \a initio wave function affected by EFP subsystem.
 *
 * \param[in] efp The efp structure.
 * \param[in] size Size of the grad array. Must be at least [3 * \a n]
 *                 elements, where \a n is the number of atoms in the
 *                 \a ab \a initio subsystem.
 * \param[out] grad For each atom \a x \a y \a z components of energy
 *                  gradient are stored.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_qm_gradient(struct efp *efp, int size, double *grad);

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
 * \param[in] frag_idx Index of a fragment. Must be a value between zero and
 *                     the total number of fragments minus one.
 * \param[in] size Size of a \a frag_name buffer.
 * \param[out] frag_name A buffer where name of the fragment will be stored.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_frag_name(struct efp *efp,
				  int frag_idx,
				  int size,
				  char *frag_name);

/**
 * Get the number of atoms in the specified fragment.
 *
 * \param[in] efp The efp structure.
 * \param[in] frag_idx Index of a fragment. Must be a value between zero and
 *                     the total number of fragments minus one.
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
 * \param[in] frag_idx Index of a fragment. Must be a value between zero and
 *                     the total number of fragments minus one.
 * \param[in] size Size of the \a atoms array. Must be greater than or equal to
 *                 the size returned by the ::efp_get_frag_atom_count function.
 * \param[out] atoms Array where atom information will be stored.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_frag_atoms(struct efp *efp,
				   int frag_idx,
				   int size,
				   struct efp_atom *atoms);

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
 * \return Human readable string with description of the result,
 *         zero-terminated.
 */
const char *efp_result_to_string(enum efp_result res);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* LIBEFP_EFP_H */
