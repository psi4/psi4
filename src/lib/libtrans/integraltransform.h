#ifndef _PSI_SRC_LIB_LIBTRANS_INTEGRALTRANSFORM_H_
#define _PSI_SRC_LIB_LIBTRANS_INTEGRALTRANSFORM_H_

#include <map>
#include <vector>
#include <string>
#include <libdpd/dpd.h>
#include <libmints/dimension.h>
#include <libmints/typedefs.h>
#include "mospace.h"

#ifndef INDEX
    #define INDEX(i,j) (((i)>(j)) ? (((i)*((i)+1)/2)+(j)) : (((j)*((j)+1)/2)+(i)))
#endif

namespace boost {
    template <class T>
    class shared_ptr;
}

namespace psi{

class Matrix;
class Dimension;
class Wavefunction;

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
        IntegralTransform(boost::shared_ptr<Wavefunction> wfn,
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
        void set_print(int n) {print_ = n;}
        /// The level of printing used during transformations
        int get_print() const {return print_;}

        /// Set the library to keep or delete the half-transformed integrals in DPD form after processing
        void set_keep_ht_ints(bool val) {keepHtInts_ = val;}
        /// Whether the library will keep or delete the half-transformed integrals in DPD form after processing
        bool get_keep_ht_ints() const {return keepHtInts_;}

        /// Set the library to keep or delete the SO integrals in DPD form after processing
        void set_keep_dpd_so_ints(bool val) {keepDpdSoInts_ = val;}
        /// Whether the library will keep or delete the SO integrals in DPD form after processing
        bool get_keep_dpd_so_ints() const {return keepDpdSoInts_;}

        /// Set the library to keep or delete the SO integrals in IWL form after processing
        void set_keep_iwl_so_ints(bool val) {keepIwlSoInts_ = val;}
        /// Whether the library will keep or delete the SO integrals in IWL form after processing
        bool get_keep_iwl_so_ints() const {return keepIwlSoInts_;}

        /// Set the memory (in MB) available to the library
        void set_memory(size_t memory) {memory_ = memory;}
        /// The amount of memory (in MB) available to the library
        size_t get_memory() const {return memory_;}

        /// Set the number of the DPD instance to be used in the transformation
        void set_dpd_id(int n) {myDPDNum_ = n;}
        /// The number of the DPD instance used in the transformation
        int get_dpd_id() const {return myDPDNum_;}

        /// Get the psio object being used by this object
        boost::shared_ptr<PSIO> get_psio() const;
        /// Set the psio object to be used.  You must delay initialization in the ctor for this to work.
        void set_psio(boost::shared_ptr<PSIO> psio);

        // Get the alpha correlated to Pitzer ordering array, used in backtransforms
        const int *alpha_corr_to_pitzer() const { return aCorrToPitzer_; }
        // Get the beta correlated to Pitzer ordering array, used in backtransforms
        const int *beta_corr_to_pitzer() const { return bCorrToPitzer_; }

    protected:
        void check_initialized();
        void semicanonicalize();
        void common_moinfo_initialize();

        void process_eigenvectors();
        void process_spaces();
        void presort_mo_tpdm_restricted();
        void presort_mo_tpdm_unrestricted();
        void sort_so_tpdm(const dpdbuf4 *B, size_t first_row, size_t num_rows);

        void trans_one(int m, int n, double *input, double *output, double **C, int soOffset,
                       int *order, bool backtransform = false, double scale = 0.0);
        void build_fzc_and_fock(int p, int q, int r, int s, double value,
                          double *aFzcD, double *bFzcD, double *aFzcOp, double *bFzcOp,
                          double *aD, double *bD, double *aFock, double *bFock);
        void idx_permute_presort(dpdfile4 *File, int &thisBucket, int **&bucketMap,
                                 int **&bucketOffset, int &p, int &q, int &r,
                                 int &s, double value, bool symmetrize = false, bool have_bra_ket_sym=true);
        void idx_error(const char *message, int p, int q, int r, int s,
                       int pq, int rs, int pq_sym, int rs_sym);

        // Has this instance been initialized yet?
        bool initialized_;

        // The wavefunction object, containing the orbital infomation
        boost::shared_ptr<Wavefunction> wfn_;
        // Pointer to the PSIO object to use for file I/O
        boost::shared_ptr<PSIO>  psio_;
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
        std::vector <char> spacesUsed_;
        // A list of the arrays to pass into libDPD
        std::vector<int*> spaceArray_;
        // The alpha orbitals per irrep for each space
        std::map<char, int* > aOrbsPI_;
        // The beta orbitals per irrep for each space
        std::map<char, int* > bOrbsPI_;
        // The alpha MO coefficients for all unique spaces needed
        std::map<char, SharedMatrix> aMOCoefficients_;
        // The beta MO coefficients for all unique spaces needed
        std::map<char, SharedMatrix> bMOCoefficients_;
        // The alpha orbital indexing arrays
        std::map<char, int* > aIndices_;
        // The beta orbital indexing arrays
        std::map<char, int* > bIndices_;
        // The lookup table for DPD indexing
        std::map<std::string, int> dpdLookup_;
        // Whether the SO integrals have already been presorted
        bool alreadyPresorted_;
        // The file to which DPD formatted integrals are written
        int dpdIntFile_;
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
        // A string describing the spaces in which the integrals are to be transformed
        char *spaces_;
        // An array containing labels for each irrep
        char **labels_;
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
        // The alpha correlated to Pitzer ordering arrays, used in backtransforms
        int *aCorrToPitzer_;
        // The beta correlated to Pitzer ordering arrays, used in backtransforms
        int *bCorrToPitzer_;
        // The number of symmetrized orbitals per irrep
        Dimension sopi_;
        // The symmetry (irrep number) of each symmetrized atomic orbital
        int *sosym_;
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
        // The cache files used by libDPD
        int *cacheFiles_, **cacheList_;
        // The alpha MO coefficients for each irrep
        boost::shared_ptr<Matrix> Ca_;
        // The alpha MO coefficients for each irrep
        boost::shared_ptr<Matrix> Cb_;
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
};

inline void
IntegralTransform::idx_permute_presort(dpdfile4 *File, int &thisBucket, int **&bucketMap,
                                       int **&bucketOffset, int &p, int &q, int &r,
                                       int &s, double value, bool symmetrize, bool have_bra_ket_sym)
{
    dpdparams4 *Params = File->params;

    if(symmetrize){
        // Symmetrize the quantity (used in density matrix processing)
        if(p!=q) value *= 0.5;
        if(r!=s) value *= 0.5;
    }

    bool bra_ket_different = !(p==r && q==s);

    /* Get the orbital symmetries */
    int p_sym = Params->psym[p];
    int q_sym = Params->qsym[q];
    int r_sym = Params->rsym[r];
    int s_sym = Params->ssym[s];
    int pq_sym = p_sym^q_sym;
    int rs_sym = r_sym^s_sym;

    /* The allowed (Mulliken) permutations are very simple in this case */
    if(bucketMap[p][q] == thisBucket) {
        /* Get the row and column indices and assign the value */
        int pq = Params->rowidx[p][q];
        int rs = Params->colidx[r][s];
        int offset = bucketOffset[thisBucket][pq_sym];
        if((pq-offset >= Params->rowtot[pq_sym]) || (rs >= Params->coltot[rs_sym]))
            idx_error("MP Params_make: pq, rs", p,q,r,s,pq,rs,pq_sym,rs_sym);
        File->matrix[pq_sym][pq-offset][rs] += value;
    }

    /*
     * We also add in the bra-ket transposed value, as a result of the matrix
     * storage, but we need to make sure we don't duplicate "diagonal" values.
     * We don't do this if the quantity does not have bra-ket symmetry, like
     * in the Alpha-Beta TPDM.
     */
    if(bucketMap[r][s] == thisBucket && bra_ket_different && have_bra_ket_sym) {
        int rs = Params->rowidx[r][s];
        int pq = Params->colidx[p][q];
        int offset = bucketOffset[thisBucket][rs_sym];
        if((rs-offset >= Params->rowtot[rs_sym])||(pq >= Params->coltot[pq_sym]))
            idx_error("MP Params_make: rs, pq", p,q,r,s,rs,pq,rs_sym,pq_sym);
        File->matrix[rs_sym][rs-offset][pq] += value;
    }
}

inline void
IntegralTransform::idx_error(const char *message, int p, int q, int r, int s,
                             int pq, int rs, int pq_sym, int rs_sym)
{

    fprintf(outfile, "\n\tDPD Parameter Error in %s\n", message);
    fprintf(outfile,"\t-------------------------------------------------\n");
    fprintf(outfile,"\t    p      q      r      s  [   pq]  [   rs] pq_symm rs_symm\n");
    fprintf(outfile,"\t%5d  %5d  %5d  %5d  [%5d]  [%5d]   %1d   %1d\n", p,q,r,s,
      pq,rs,pq_sym,rs_sym);
    throw PsiException("DPD idx failure.", __FILE__, __LINE__);
}


} // End namespaces

#endif // Header guard
