#ifndef _PSI_SRC_LIB_LIBTRANS_INTEGRALTRANSFORM_H_
#define _PSI_SRC_LIB_LIBTRANS_INTEGRALTRANSFORM_H_

#include <psi4-dec.h>
#include <map>
#include <vector>
#include <libdpd/dpd.h>

#define INDEX(i,j) ((i>j) ? ((i*(i+1)/2)+j) : ((j*(j+1)/2)+i))

using namespace psi;

namespace psi{ namespace libtrans{

class SpaceInfo;
class MOSpace;

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
         * DPDOnly   - Write the integrals to a DPD structure(s)
         * IWLAndDPD - Write the integrals to an IWL-formatted file in addition to the
         *           - DPD buffer
         */
        enum OutputType {DPDOnly, IWLAndDPD};
        /**
         * Which orbitals are to be "frozen" i.e. excluded from the transformation:-
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
         * @param options            An Options object, passed by reference
         * @param s1                 An MOSpace object describing one of the spaces to transform.
         * @param s2                 An MOSpace object describing one of the spaces to transform.
         * @param s3                 An MOSpace object describing one of the spaces to transform.
         * @param s4                 An MOSpace object describing one of the spaces to transform.
         * @param transformationType The type of transformation, described by the
         *                           enum TransformationType
         * @param moOrdering         The ordering convention of the resulting integrals, see
         *                           enum MOOrdering
         * @param outputType         The storage format of the transformed integrals, see
         *                           enum OutputType
         * @param frozenOrbitals     Which orbitals are to be excluded from the transformation, see
         *                           enum FrozenOrbitals
         */
        IntegralTransform(Options &options,
                          shared_ptr<MOSpace> s1,
                          shared_ptr<MOSpace> s2,
                          shared_ptr<MOSpace> s3,
                          shared_ptr<MOSpace> s4,
                          TransformationType transformationType = Restricted,
                          MOOrdering moOrdering = QTOrder,
                          OutputType outputType = IWLAndDPD,
                          FrozenOrbitals frozenOrbitals = None);
        ~IntegralTransform();

        void presort_so_tei();
        void transform_oei(shared_ptr<MOSpace> s1, shared_ptr<MOSpace> s2);
        void transform_tei(shared_ptr<MOSpace> s1, shared_ptr<MOSpace> s2,
                           shared_ptr<MOSpace> s3, shared_ptr<MOSpace> s4);

    protected:
        void process_spaces(std::vector<shared_ptr<MOSpace> >);

        void trans_one(int m, int n, double *input, double *output,
                                double **C, int nc, int *order);
        void frozen_core(int p, int q, int r, int s, double value);
        void idx_permute_presort(dpdfile4 *File, int &thisBucket, int **&bucketMap,
                                 int **&bucketOffset, int &p, int &q, int &r, int &s,
                                 double &value);
        void idx_error(const char *message, int p, int q, int r, int s,
                       int pq, int rs, int pq_sym, int rs_sym);
        int DPD_ID(shared_ptr<MOSpace> s1, shared_ptr<MOSpace> s2, SpinType spin, bool pack);

        // The options object
        Options _options;
        // The type of transformation
        TransformationType _transformationType;
        // The ordering of the resulting integrals
        MOOrdering _moOrdering;
        // The format of the outputted integrals
        OutputType _outputType;
        // How to handle frozen orbitals
        FrozenOrbitals _frozenOrbitals;
        // The unique orbital spaces involved in this transformation
        std::vector <shared_ptr<MOSpace> > _spacesUsed;
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
        // The order in which each alpha space was added
        std::map<char, int> _aSpaceNum;
        // The order in which each beta space was added
        std::map<char, int> _bSpaceNum;
        // The number of irreps
        int _nirreps;
        // The number of molecular orbitals
        int _nmo;
        // The number of symmetrized atomic orbitals
        int _nso;
        // The number of pairs of symmetrized atomic orbitals
        int _nTriSo;
        // The number of atomic orbitals
        int _nao;
        // The number of frozen doubly occupied orbitals
        int _nfzc;
        // The number of frozen virtual orbitals
        int _nfzv;
        // A string describing the spaces in which the integrals are to be transformed
        char *_spaces;
        // An array containing labels for each irrep
        char **_labels;
        // The nuclear repulsion energy
        double _enuc;
        // The SCF energy (from the checkpoint file)
        double _escf;
        // The definition of zero
        double _tolerance;
        // The alpha frozen core operator
        double *_aFzcOp;
        // The beta frozen core operator
        double *_bFzcOp;
        // The alpha frozen core density operator
        double *_aFzcD;
        // The beta frozen core density operator
        double *_bFzcD;
        // The amount of memory, in bytes
        size_t _memory;
        // The PSI file number for the alpha-alpha integrals
        int _moIntFileAA;
        // The PSI file number for the alpha-beta integrals
        int _moIntFileAB;
        // The PSI file number for the beta-beta integrals
        int _moIntFileBB;
        // The amount of information to print
        int _print;
        // The number of symmetrized orbitals per irrep
        int *_sopi;
        // The symmetry (irrep number) of each symmetrized atomic orbital
        int *_sosym;
        // The number of molecular orbitals per irrep
        int *_mopi;
        // The number of doubly-occupied orbitals per irrep
        int *_clsdpi;
        // The number of singly-occupied orbitals per irrep
        int *_openpi;
        // The number of frozen doubly occupied orbitals per irrep
        int *_frzcpi;
        // The number of frozen virtual orbitals per irrep
        int *_frzvpi;
        // The cache files used by libDPD
        int *_cacheFiles, **_cacheList;
        // The alpha frozen core operator
        double *_aFzcOperator;
        // The beta frozen core operator
        double *_bFzcOperator;
        // The full alpha MO coefficients
        double **_fullCa;
        // The full beta MO coefficients
        double **_fullCb;
        // The alpha MO coefficients for each irrep
        double ***_Ca;
        // The alpha MO coefficients for each irrep
        double ***_Cb;
        // Whether the moinfo routine has already been called
        bool _moinfo_initialized;
        // Whether to keep the IWL SO integral file after processing
        bool _deleteIwlSoTei;
        // Whether to print the two-electron integrals or not
        bool _printTei;
};

}} // End namespaces

// This is here so that files including this have clean(er) syntax
using namespace psi::libtrans;

#endif // Header guard
