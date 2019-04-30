/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
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

#ifndef _PSI_SRC_LIB_LIBTRANS_INTEGRALTRANSFORM_H_
#define _PSI_SRC_LIB_LIBTRANS_INTEGRALTRANSFORM_H_

#include <array>
#include <map>
#include <vector>
#include <string>
#include "psi4/libmints/dimension.h"
#include "psi4/libmints/typedefs.h"
#include "mospace.h"

#ifndef INDEX
#define INDEX(i, j) (((i) > (j)) ? (((i) * ((i) + 1) / 2) + (j)) : (((j) * ((j) + 1) / 2) + (i)))
#endif

namespace psi {

struct dpdbuf4;
class Matrix;
class Dimension;
class Wavefunction;
class PSIO;

typedef std::vector<std::shared_ptr<MOSpace> > SpaceVec;

/**
   The IntegralTransform class transforms one- and two-electron integrals
   within general spaces
 */

class PSI_API IntegralTransform {
    // TODO check usage of restricted, to make sure that it's correct everywhere
   public:
    /**
     * How to handle reuse of intermediates in two-electron integral transformations
     *
     * MakeAndKeep - Compute the half-transformed integrals for this transformation, and keep after use
     * ReadAndKeep - Read the half-transformed integrals for this transformation from disk, and keep after use
     * MakeAndNuke - Compute the half-transformed integrals for this transformation, and delete after use
     * ReadAndNuke - Read the half-transformed integrals for this transformation from disk, and delete after use
     */
    enum class HalfTrans { MakeAndKeep, ReadAndKeep, MakeAndNuke, ReadAndNuke };

    /**
     * Possible Transformations:-
     *
     * Restricted    - Same alpha and beta orbitals
     * Unrestricted  - Different alpha and beta orbitals
     * SemiCanonical - Start from restricted orbitals and diagnonalize alpha
     *               - and beta occ-occ and vir-vir block seperately, leading
     *               - to different alpha and beta orbitals
     */
    enum class TransformationType { Restricted, Unrestricted, SemiCanonical };
    /**
     * Ordering of the resulting integrals:-
     *
     * QTOrder - Ordered by class (frozen docc < docc < socc < virtual < frozen virtual )
     *           then by irrep within these classes
     * Pitzer  - Ordered by irreps, then by orbital energy within irreps
     */
    enum class MOOrdering { QTOrder, PitzerOrder };
    /**
     * Output format for the resulting integrals:-
     *
     * DPDOnly   - Write the integrals to (a) DPD structure(s)
     * IWLOnly   - Write the integrals to (an) IWL structure(s)
     * IWLAndDPD - Write the integrals to an IWL-formatted file in addition to the
     *           - DPD buffer
     */
    enum class OutputType { DPDOnly, IWLOnly, IWLAndDPD };
    /**
     * Which orbitals are to be "frozen" i.e. excluded from the transformation.
     * N.B. Orbitals are only frozen if input detected a request to do so.  That
     * means that requesting a frozen-core transformation will result in an
     * all-electron transformation, unless some orbitals were actually frozen.
     *
     * None      - No orbitals are excluded
     * OccOnly   - Only the frozen occupied orbitals are excluded
     * VirOnly   - Only the frozen virtual orbitals are excluded
     * OccAndVir - The frozen occupied and the frozen core orbitals are excluded
     */
    enum class FrozenOrbitals { None, OccOnly, VirOnly, OccAndVir };
    /**
     * The spin of the electron.  This could, of course, be a boolean but an
     * enum makes things a little more transparent
     *
     * Alpha = spin up
     * Beta  = spin down
     */
    enum class SpinType { Alpha, Beta };
    /**
     * Set up a transformation involving four MO spaces
     *
     * @param wfn                A (shared pointer to a) wavefunction object with the orbital info
     * @param spaces             A vector containing smart pointers to the unique space(s) involved
     *                           in any transformations that this object will perform
     * @param transformationType The type of transformation, described by the
     *                           enum TransformationType
     * @param moOrdering         The ordering convention of the resulting integrals, see
     *                           enum MOOrdering.  This only affects IWL output.
     * @param outputType         The storage format of the transformed integrals, see
     *                           enum OutputType
     * @param frozenOrbitals     Which orbitals are to be excluded from the transformation, see
     *                           enum FrozenOrbitals
     * @param initialize         Whether to initialize during construction or not.  Useful if some
     *                           options need to be tweaked before initialization.
     */
    IntegralTransform(std::shared_ptr<Wavefunction> wfn, SpaceVec spaces,
                      TransformationType transformationType = TransformationType::Restricted,
                      OutputType outputType = OutputType::DPDOnly, MOOrdering moOrdering = MOOrdering::QTOrder,
                      FrozenOrbitals frozenOrbitals = FrozenOrbitals::OccAndVir, bool initialize = true);

    IntegralTransform(SharedMatrix H, SharedMatrix c, SharedMatrix i, SharedMatrix a, SharedMatrix v, SpaceVec spaces,
                      TransformationType transformationType = TransformationType::Restricted,
                      OutputType outputType = OutputType::DPDOnly, MOOrdering moOrdering = MOOrdering::QTOrder,
                      FrozenOrbitals frozenOrbitals = FrozenOrbitals::OccAndVir, bool initialize = true);

    ~IntegralTransform();

    void initialize();
    void presort_so_tei();
    void generate_oei();
    void update_orbitals();
    void transform_T_plus_V(const std::shared_ptr<MOSpace> s1, const std::shared_ptr<MOSpace> s2);
    void transform_oei(const std::shared_ptr<MOSpace> s1, const std::shared_ptr<MOSpace> s2,
                       const std::array<std::string, 4> &labels);
    void transform_tei(const std::shared_ptr<MOSpace> s1, const std::shared_ptr<MOSpace> s2,
                       const std::shared_ptr<MOSpace> s3, const std::shared_ptr<MOSpace> s4,
                       HalfTrans = HalfTrans::MakeAndNuke);
    void transform_tei_first_half(const std::shared_ptr<MOSpace> s1, const std::shared_ptr<MOSpace> s2);
    void transform_tei_second_half(const std::shared_ptr<MOSpace> s1, const std::shared_ptr<MOSpace> s2,
                                   const std::shared_ptr<MOSpace> s3, const std::shared_ptr<MOSpace> s4);
    void backtransform_density();
    void backtransform_tpdm_restricted();
    void backtransform_tpdm_unrestricted();
    void print_dpd_lookup();
    std::vector<SharedMatrix> compute_fock_like_matrices(SharedMatrix Hcore, std::vector<SharedMatrix> Cmats);

    int DPD_ID(const char c);
    int DPD_ID(char *str);
    int DPD_ID(const char *str);
    int DPD_ID(const std::string &str);
    int DPD_ID(const std::shared_ptr<MOSpace> s1, const std::shared_ptr<MOSpace> s2, SpinType spin, bool pack);

    /*===== The set/get accessor functions =====*/

    /// Sets the SO IWL file to read the TEIs from.
    void set_so_tei_file(int so_tei_file) { soIntTEIFile_ = so_tei_file; }
    /// Set whether to write a DPD formatted SO basis TPDM to disk after density transformations
    void set_write_dpd_so_tpdm(bool t_f) { write_dpd_so_tpdm_ = t_f; }
    /// Set the level of printing used during transformations (0 -> 6)
    void set_print(int n) { print_ = n; }
    /// Whether to build the Fock matrix in the MO basis during integral presort
    void set_build_mo_fock(bool t_f) { buildMOFock_ = t_f; }
    /// Sets the orbitals to the given C matrix. This is a hack for MCSCF wavefunctions.
    /// Use with caution.
    void set_orbitals(SharedMatrix C);

    /// The level of printing used during transformations
    int get_print() const { return print_; }

    /// Returns the frozen-core energy
    double get_frozen_core_energy() const { return frozen_core_energy_; }

    /// Set the library to keep or delete the half-transformed integrals in DPD form after processing
    void set_keep_ht_ints(bool val) { keepHtInts_ = val; }
    /// Whether the library will keep or delete the half-transformed integrals in DPD form after processing
    bool get_keep_ht_ints() const { return keepHtInts_; }

    /// Set the library to keep or delete the SO integrals in DPD form after processing
    void set_keep_dpd_so_ints(bool val) { keepDpdSoInts_ = val; }
    /// Whether the library will keep or delete the SO integrals in DPD form after processing
    bool get_keep_dpd_so_ints() const { return keepDpdSoInts_; }

    /// Set the library to keep or delete the SO integrals in IWL form after processing
    void set_keep_iwl_so_ints(bool val) { keepIwlSoInts_ = val; }
    /// Whether the library will keep or delete the SO integrals in IWL form after processing
    bool get_keep_iwl_so_ints() const { return keepIwlSoInts_; }
    /// Whether TPDM has already presorted
    void set_tpdm_already_presorted(bool val) { tpdmAlreadyPresorted_ = val; }

    /// Whether SO intergals are already presorted
    bool get_tei_already_presorted() { return alreadyPresorted_; }
    void set_tei_already_presorted(bool val) { alreadyPresorted_ = val; }

    /// Set the memory (in MB) available to the library
    void set_memory(size_t memory) { memory_ = memory; }
    /// The amount of memory (in MB) available to the library
    size_t get_memory() const { return memory_; }

    /// Set the number of the DPD instance to be used in the transformation
    void set_dpd_id(int n) { myDPDNum_ = n; }
    /// The number of the DPD instance used in the transformation
    int get_dpd_id() const { return myDPDNum_; }

    /// Get the psio object being used by this object
    std::shared_ptr<PSIO> get_psio() const;
    /// Set the psio object to be used.  You must delay initialization in the ctor for this to work.
    void set_psio(std::shared_ptr<PSIO> psio);

    /// The file to output DPD integrals to
    void set_dpd_int_file(int file) { dpdIntFile_ = file; }
    /// Set the name used for the Alpha-Alpha integrals in the DPD file.  Needs to be
    /// called before each transformation in order to override the default name.
    void set_aa_int_name(const std::string &name) { aaIntName_ = name; }
    /// Set the name used for the Alpha-Beta integrals in the DPD file.  Needs to be
    /// called before each transformation in order to override the default name.
    void set_ab_int_name(const std::string &name) { abIntName_ = name; }
    /// Set the name used for the Beta-Beta integrals in the DPD file.  Needs to be
    /// called before each transformation in order to override the default name.
    void set_bb_int_name(const std::string &name) { bbIntName_ = name; }

    /// Get the alpha correlated to Pitzer ordering array, used in backtransforms
    const int *alpha_corr_to_pitzer() const { return aCorrToPitzer_; }
    /// Get the beta correlated to Pitzer ordering array, used in backtransforms
    const int *beta_corr_to_pitzer() const { return bCorrToPitzer_; }

    int nirrep() const { return nirreps_; }

    void reset_so_int() { alreadyPresorted_ = false; }

   protected:
    void check_initialized();
    void common_initialize();

    void process_eigenvectors();
    void process_spaces();
    void presort_mo_tpdm_restricted();
    void presort_mo_tpdm_unrestricted();
    void setup_tpdm_buffer(const dpdbuf4 *D);
    void sort_so_tpdm(const dpdbuf4 *B, int irrep, size_t first_row, size_t num_rows, bool first_run);

    void transform_oei_restricted(const std::shared_ptr<MOSpace> s1, const std::shared_ptr<MOSpace> s2,
                                  const std::vector<double> &soInts, std::string label);
    void transform_oei_unrestricted(const std::shared_ptr<MOSpace> s1, const std::shared_ptr<MOSpace> s2,
                                    const std::vector<double> &soInts, std::string A_label, std::string B_label);
    void trans_one(int m, int n, double *input, double *output, double **C, int soOffset, int *order,
                   bool backtransform = false, double scale = 0.0);

    // Has this instance been initialized yet?
    bool initialized_;

    // The number of SO tpdm elements in each SO shell quartet
    std::vector<size_t> tpdm_buffer_sizes_;
    // The buffer used in sorting the SO basis tpdm
    double *tpdm_buffer_;
    // Frozen core energy
    double frozen_core_energy_;
    // The wavefunction object, containing the orbital infomation
    std::shared_ptr<Wavefunction> wfn_;
    // Pointer to the PSIO object to use for file I/O
    std::shared_ptr<PSIO> psio_;
    // The type of transformation
    TransformationType transformationType_;
    // The unique MO spaces provided to this object's constructor
    SpaceVec uniqueSpaces_;
    // The ordering of the resulting integrals
    MOOrdering moOrdering_;
    // The format of the outputted integrals
    OutputType outputType_;
    // How to handle frozen orbitals
    FrozenOrbitals frozenOrbitals_;
    // The unique orbital spaces involved in this transformation
    std::vector<char> spacesUsed_;
    // A list of the arrays to pass into libDPD
    std::vector<int *> spaceArray_;
    // The alpha orbitals per irrep for each space
    std::map<char, int *> aOrbsPI_;
    // The beta orbitals per irrep for each space
    std::map<char, int *> bOrbsPI_;
    // The alpha MO coefficients for all unique spaces needed
    std::map<char, SharedMatrix> aMOCoefficients_;
    // The beta MO coefficients for all unique spaces needed
    std::map<char, SharedMatrix> bMOCoefficients_;
    // The alpha orbital indexing arrays
    std::map<char, int *> aIndices_;
    // The beta orbital indexing arrays
    std::map<char, int *> bIndices_;
    // The lookup table for DPD indexing
    std::map<std::string, int> dpdLookup_;
    // Whether the SO integrals have already been presorted
    bool alreadyPresorted_;
    // Whether to also write DPD formatted SO TPDMs after density transformations
    bool write_dpd_so_tpdm_;
    // The file to which DPD formatted integrals are written
    int dpdIntFile_;
    // The file from which IWL SO TEIs are read
    int soIntTEIFile_;
    // The file containing alpha half-transformed integrals in DPD format
    int aHtIntFile_;
    // The file containing beta half-transformed integrals in DPD format
    int bHtIntFile_;
    // The file containing alpha-alpha IWL formatted integrals
    int iwlAAIntFile_;
    // The file containing alpha-beta IWL formatted integrals
    int iwlABIntFile_;
    // The file containing beta-beta IWL formatted integrals
    int iwlBBIntFile_;
    // The number of irreps
    int nirreps_;
    // The number of molecular orbitals
    int nmo_;
    // The number of symmetrized atomic orbitals
    int nso_;
    // The number of pairs of symmetrized atomic orbitals
    int nTriSo_;
    // The number of pairs of molecular orbitals
    int nTriMo_;
    // The number of frozen doubly occupied orbitals
    int nfzc_;
    // The number of frozen virtual orbitals
    int nfzv_;
    // The name of the alpha-alpha DPD integral buffer
    std::string aaIntName_;
    // The name of the alpha-beta DPD integral buffer
    std::string abIntName_;
    // The name of the beta-beta DPD integral buffer
    std::string bbIntName_;
    // A string describing the spaces in which the integrals are to be transformed
    char *spaces_;
    // An array containing labels for each irrep
    std::vector<std::string> labels_;
    // The definition of zero
    double tolerance_;
    // The amount of memory, in MB
    size_t memory_;
    // The PSI file number for the alpha-alpha integrals
    int moIntFileAA_;
    // The PSI file number for the alpha-beta integrals
    int moIntFileAB_;
    // The PSI file number for the beta-beta integrals
    int moIntFileBB_;
    // The DPD id to use internally
    int myDPDNum_;
    // The amount of information to print
    int print_;
    // Just an array of zeros! Used in the null MOSpace "transforms"
    int *zeros_;
    // The alpha Pitzer->QT reordering array
    int *aQT_;
    // The alpha Pitzer->QT reordering array
    int *bQT_;
    // The alpha correlated to Pitzer ordering arrays, used in backtransforms
    int *aCorrToPitzer_;
    // The beta correlated to Pitzer ordering arrays, used in backtransforms
    int *bCorrToPitzer_;
    // The number of symmetrized orbitals per irrep
    Dimension sopi_;
    // The symmetry (irrep number) of each symmetrized atomic orbital
    int *sosym_;
    // The symmetry (irrep number) of each molecular orbital
    int *mosym_;
    // The number of molecular orbitals per irrep
    Dimension mopi_;
    // The number of doubly-occupied orbitals per irrep
    Dimension clsdpi_;
    // The number of singly-occupied orbitals per irrep
    Dimension openpi_;
    // The number of frozen doubly occupied orbitals per irrep
    Dimension frzcpi_;
    // The number of frozen virtual orbitals per irrep
    Dimension frzvpi_;
    // The number of alpha electrons per irrep
    Dimension nalphapi_;
    // The number of beta electrons per irrep
    Dimension nbetapi_;
    // The cache files used by libDPD
    int *cacheFiles_, **cacheList_;
    // The alpha MO coefficients for each irrep
    std::shared_ptr<Matrix> Ca_;
    // The alpha MO coefficients for each irrep
    std::shared_ptr<Matrix> Cb_;
    // The one electron Hamiltonian matrix for each irrep
    std::shared_ptr<Matrix> H_;
    // Whether to keep the IWL SO integral file after processing
    bool keepIwlSoInts_;
    // Whether to keep the IWL MO two particle density matrix
    bool keepIwlMoTpdm_;
    // Whether to keep the DPD SO integral file after processing
    bool keepDpdSoInts_;
    // Whether to keep the DPD MO to particle density matrix after processing
    bool keepDpdMoTpdm_;
    // Whether to keep the half-transformed two electron integrals
    bool keepHtInts_;
    // Whether to keep the half-transformed TPDM
    bool keepHtTpdm_;
    // Whether to print the two-electron integrals or not
    bool printTei_;
    // Whether to output the results to an IWL buffer
    bool useIWL_;
    // Whether to output the results to a DPD buffer
    bool useDPD_;
    // Has this object already pre-sorted?
    bool tpdmAlreadyPresorted_;
    // Whether to form the MO basis Fock matrix during TEI presort
    bool buildMOFock_;
    // This keeps track of which labels have been assigned by other spaces
    std::map<char, int> labelsUsed_;
};

}  // namespace psi

#endif  // Header guard
