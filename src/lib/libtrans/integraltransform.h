#ifndef _PSI_SRC_LIB_LIBTRANS_INTEGRALTRANSFORM_H_
#define _PSI_SRC_LIB_LIBTRANS_INTEGRALTRANSFORM_H_

#include <psi4-dec.h>
#include <map>
#include <vector>
#include <string>
#include <libdpd/dpd.h>
#include <libchkpt/chkpt.hpp>
#include <libmints/dimension.h>
#include <psifiles.h>
#include "mospace.h"

#ifndef INDEX
    #define INDEX(i,j) (((i)>(j)) ? (((i)*((i)+1)/2)+(j)) : (((j)*((j)+1)/2)+(i)))
#endif

namespace psi{

class Matrix;
class Dimension;

typedef std::vector<boost::shared_ptr< MOSpace> > SpaceVec;

  /**
     The IntegralTransform class transforms one- and two-electron integrals
     within general spaces
   */

class IntegralTransform{
// TODO check usage of restricted, to make sure that it's correct everywhere
    public:
        /**
         * Possible Transformations:-
         *
         * Restricted    - Same alpha and beta orbitals
         * Unrestricted  - Different alpha and beta orbitals
         * SemiCanonical - Start from restricted orbitals and diagnonalize alpha
         *               - and beta occ-occ and vir-vir block seperately, leading
         *               - to different alpha and beta orbitals
         */
        enum TransformationType {Restricted, Unrestricted, SemiCanonical};
        /**
         * Ordering of the resulting integrals:-
         *
         * QTOrder - Ordered by class (frozen docc < docc < socc < virtual < frozen virtual )
         *           then by irrep within these classes
         * Pitzer  - Ordered by irreps, then by orbital energy within irreps
         */
        enum MOOrdering {QTOrder, PitzerOrder};
        /**
         * Output format for the resulting integrals:-
         *
         * DPDOnly   - Write the integrals to (a) DPD structure(s)
         * IWLOnly   - Write the integrals to (an) IWL structure(s)
         * IWLAndDPD - Write the integrals to an IWL-formatted file in addition to the
         *           - DPD buffer
         */
        enum OutputType {DPDOnly, IWLOnly, IWLAndDPD};
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
        enum FrozenOrbitals {None, OccOnly, VirOnly, OccAndVir};
        /**
         * The spin of the electron.  This could, of course, be a boolean but an
         * enum makes things a little more transparent
         *
         * Alpha = spin up
         * Beta  = spin down
         */
        enum SpinType {Alpha, Beta};
        /**
         * Set up a transformation involving four MO spaces
         *
         * @param chkpt              The checkpoint object to use for reading orbital info
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
        IntegralTransform(boost::shared_ptr<Chkpt> chkpt,
                          SpaceVec spaces,
                          TransformationType transformationType = Restricted,
                          OutputType outputType = DPDOnly,
                          MOOrdering moOrdering = QTOrder,
                          FrozenOrbitals frozenOrbitals = OccAndVir,
                          bool initialize = true);

        IntegralTransform(boost::shared_ptr<Wavefunction> wave,
                          SpaceVec spaces,
                          TransformationType transformationType = Restricted,
                          OutputType outputType = DPDOnly,
                          MOOrdering moOrdering = QTOrder,
                          FrozenOrbitals frozenOrbitals = OccAndVir,
                          bool initialize = true);

        IntegralTransform(SharedMatrix c,
                          SharedMatrix i,
                          SharedMatrix a,
                          SharedMatrix v,
                          SpaceVec spaces,
                          TransformationType transformationType = Restricted,
                          OutputType outputType = DPDOnly,
                          MOOrdering moOrdering = QTOrder,
                          FrozenOrbitals frozenOrbitals = OccAndVir,
                          bool initialize = true);

        ~IntegralTransform();

        void initialize();
        void presort_so_tei();
        void generate_oei();
        void update_orbitals();
        void transform_oei(const boost::shared_ptr<MOSpace> s1, const boost::shared_ptr<MOSpace> s2, const char *label);
        void transform_tei(const boost::shared_ptr<MOSpace> s1, const boost::shared_ptr<MOSpace> s2,
                           const boost::shared_ptr<MOSpace> s3, const boost::shared_ptr<MOSpace> s4);
        void transform_tei_first_half(const boost::shared_ptr<MOSpace> s1, const boost::shared_ptr<MOSpace> s2);
        void transform_tei_second_half(const boost::shared_ptr<MOSpace> s1, const boost::shared_ptr<MOSpace> s2,
                                       const boost::shared_ptr<MOSpace> s3, const boost::shared_ptr<MOSpace> s4);
        void backtransform_density();
        void backtransform_tpdm_restricted();
        void backtransform_tpdm_unrestricted();
        void print_dpd_lookup();

        int DPD_ID(const char c);
        int DPD_ID(char *str);
        int DPD_ID(const char *str);
        int DPD_ID(const std::string &str);
        int DPD_ID(const boost::shared_ptr<MOSpace> s1, const boost::shared_ptr<MOSpace> s2, SpinType spin, bool pack);

        /*===== The set/get accessor functions =====*/

        /// Set the level of printing used during transformations (0 -> 6)
        void set_print(int n) {_print = n;}
        /// The level of printing used during transformations
        int get_print() const {return _print;}

        /// Set the library to keep or delete the half-transformed integrals in DPD form after processing
        void set_keep_ht_ints(bool val) {_keepHtInts = val;}
        /// Whether the library will keep or delete the half-transformed integrals in DPD form after processing
        bool get_keep_ht_ints() const {return _keepHtInts;}

        /// Set the library to keep or delete the SO integrals in DPD form after processing
        void set_keep_dpd_so_ints(bool val) {_keepDpdSoInts = val;}
        /// Whether the library will keep or delete the SO integrals in DPD form after processing
        bool get_keep_dpd_so_ints() const {return _keepDpdSoInts;}

        /// Set the library to keep or delete the SO integrals in IWL form after processing
        void set_keep_iwl_so_ints(bool val) {_keepIwlSoInts = val;}
        /// Whether the library will keep or delete the SO integrals in IWL form after processing
        bool get_keep_iwl_so_ints() const {return _keepIwlSoInts;}

        /// Set the memory (in MB) available to the library
        void set_memory(size_t memory) {_memory = memory;}
        /// The amount of memory (in MB) available to the library
        size_t get_memory() const {return _memory;}

        /// Set the number of the DPD instance to be used in the transformation
        void set_dpd_id(int n) {_myDPDNum = n;}
        /// The number of the DPD instance used in the transformation
        int get_dpd_id() const {return _myDPDNum;}

        /// Get the psio object being used by this object
        boost::shared_ptr<PSIO> get_psio() const {return _psio;}
        /// Set the psio object to be used.  You must delay initialization in the ctor for this to work.
        void set_psio(boost::shared_ptr<PSIO> psio) {_psio = psio;}

        // Get the alpha correlated to Pitzer ordering array, used in backtransforms
        const int *alpha_corr_to_pitzer() const { return _aCorrToPitzer; }
        // Get the beta correlated to Pitzer ordering array, used in backtransforms
        const int *beta_corr_to_pitzer() const { return _bCorrToPitzer; }

    protected:
        void check_initialized();
        void semicanonicalize();
        void common_moinfo_initialize();

        void process_eigenvectors();
        void process_spaces();
        void presort_mo_tpdm_restricted();
        void presort_mo_tpdm_unrestricted();
        void setup_backtrans_reordering();

        void trans_one(int m, int n, double *input, double *output, double **C, int soOffset,
                       int *order, bool backtransform = false, bool accumulate = false);
        void build_fzc_and_fock(int p, int q, int r, int s, double value,
                          double *aFzcD, double *bFzcD, double *aFzcOp, double *bFzcOp,
                          double *aD, double *bD, double *aFock, double *bFock);
        void idx_permute_presort(dpdfile4 *File, int &thisBucket, int **&bucketMap,
                                 int **&bucketOffset, int &p, int &q, int &r,
                                 int &s, double value, bool symmetrize = false);
        void idx_error(const char *message, int p, int q, int r, int s,
                       int pq, int rs, int pq_sym, int rs_sym);

        // Has this instance been initialized yet?
        bool _initialized;

        // Pointer to the PSIO object to use for file I/O
        boost::shared_ptr<PSIO>  _psio;
        // Pointer to the checkpoint object to use
        boost::shared_ptr<Chkpt> _chkpt;
        // The type of transformation
        TransformationType _transformationType;
        // The unique MO spaces provided to this object's constructor
        SpaceVec _uniqueSpaces;
        // The ordering of the resulting integrals
        MOOrdering _moOrdering;
        // The format of the outputted integrals
        OutputType _outputType;
        // How to handle frozen orbitals
        FrozenOrbitals _frozenOrbitals;
        // The unique orbital spaces involved in this transformation
        std::vector <char> _spacesUsed;
        // A list of the arrays to pass into libDPD
        std::vector<int*> _spaceArrays;
        // The alpha orbitals per irrep for each space
        std::map<char, int* > _aOrbsPI;
        // The beta orbitals per irrep for each space
        std::map<char, int* > _bOrbsPI;
        // The alpha MO coefficients for all unique spaces needed
        std::map<char, double*** > _aMOCoefficients;
        // The beta MO coefficients for all unique spaces needed
        std::map<char, double*** > _bMOCoefficients;
        // The alpha orbital indexing arrays
        std::map<char, int* > _aIndices;
        // The beta orbital indexing arrays
        std::map<char, int* > _bIndices;
        // The lookup table for DPD indexing
        std::map<std::string, int> _dpdLookup;
        // Whether the SO integrals have already been presorted
        bool _alreadyPresorted;
        // The file to which DPD formatted integrals are written
        int _dpdIntFile;
        // The file containing alpha half-transformed integrals in DPD format
        int _aHtIntFile;
        // The file containing beta half-transformed integrals in DPD format
        int _bHtIntFile;
        // The file containing alpha-alpha IWL formatted integrals
        int _iwlAAIntFile;
        // The file containing alpha-beta IWL formatted integrals
        int _iwlABIntFile;
        // The file containing beta-beta IWL formatted integrals
        int _iwlBBIntFile;
        // The number of irreps
        int _nirreps;
        // The number of molecular orbitals
        int _nmo;
        // The number of symmetrized atomic orbitals
        int _nso;
        // The number of pairs of symmetrized atomic orbitals
        int _nTriSo;
        // The number of pairs of molecular orbitals
        int _nTriMo;
        // The number of frozen doubly occupied orbitals
        int _nfzc;
        // The number of frozen virtual orbitals
        int _nfzv;
        // A string describing the spaces in which the integrals are to be transformed
        char *_spaces;
        // An array containing labels for each irrep
        char **_labels;
        // The definition of zero
        double _tolerance;
        // The amount of memory, in MB
        size_t _memory;
        // The PSI file number for the alpha-alpha integrals
        int _moIntFileAA;
        // The PSI file number for the alpha-beta integrals
        int _moIntFileAB;
        // The PSI file number for the beta-beta integrals
        int _moIntFileBB;
        // The DPD id to use internally
        int _myDPDNum;
        // The amount of information to print
        int _print;
        // Just an array of zeros! Used in the null MOSpace "transforms"
        int *_zeros;
        // The alpha correlated to Pitzer ordering arrays, used in backtransforms
        int *_aCorrToPitzer;
        // The beta correlated to Pitzer ordering arrays, used in backtransforms
        int *_bCorrToPitzer;
        // The number of symmetrized orbitals per irrep
        Dimension _sopi;
        // The symmetry (irrep number) of each symmetrized atomic orbital
        int *_sosym;
        // The number of molecular orbitals per irrep
        Dimension _mopi;
        // The number of doubly-occupied orbitals per irrep
        Dimension _clsdpi;
        // The number of singly-occupied orbitals per irrep
        Dimension _openpi;
        // The number of frozen doubly occupied orbitals per irrep
        Dimension _frzcpi;
        // The number of frozen virtual orbitals per irrep
        Dimension _frzvpi;
        // The cache files used by libDPD
        int *_cacheFiles, **_cacheList;
        // Matrix objects of Ca and Cb (these are copies of _Ca, _Cb below).
        SharedMatrix _mCa;
        SharedMatrix _mCb;
        // The alpha MO coefficients for each irrep
        double ***_Ca;
        // The alpha MO coefficients for each irrep
        double ***_Cb;
        // Whether to keep the IWL SO integral file after processing
        bool _keepIwlSoInts;
        // Whether to keep the IWL MO two particle density matrix
        bool _keepIwlMoTpdm;
        // Whether to keep the DPD SO integral file after processing
        bool _keepDpdSoInts;
        // Whether to keep the DPD MO to particle density matrix after processing
        bool _keepDpdMoTpdm;
        // Whether to keep the half-transformed two electron integrals
        bool _keepHtInts;
        // Whether to keep the half-transformed TPDM
        bool _keepHtTpdm;
        // Whether to print the two-electron integrals or not
        bool _printTei;
        // Whether to output the results to an IWL buffer
        bool _useIWL;
        // Whether to output the results to a DPD buffer
        bool _useDPD;
        // Has this object already pre-sorted?
        bool _tpdmAlreadyPresorted;
};

} // End namespaces

#endif // Header guard
