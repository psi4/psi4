#ifndef __CINT_H__
#define __CINT_H__


#include <stdint.h>


struct ERD;
struct OED;
struct BasisSet;


typedef struct OED *OED_t;
typedef struct ERD *ERD_t;
typedef struct BasisSet *BasisSet_t;


typedef enum
{
    CINT_STATUS_SUCCESS = 0,
    CINT_STATUS_NOT_INITIALIZED = 1,
    CINT_STATUS_ALLOC_FAILED = 2,
    CINT_STATUS_INVALID_VALUE = 3,
    CINT_STATUS_EXECUTION_FAILED = 4,
    CINT_STATUS_INTERNAL_ERROR = 5,
    CINT_STATUS_FILEIO_FAILED = 6,
    CINT_STATUS_OFFLOAD_ERROR = 7
} CIntStatus_t;


#ifdef __INTEL_OFFLOAD
extern __declspec(target(mic)) ERD_t erd_mic;
extern __declspec(target(mic)) BasisSet_t basis_mic;
#endif

#ifdef __cplusplus
extern "C" {
#endif

// basisset
#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

CIntStatus_t CInt_createBasisSet( BasisSet_t *basis );

CIntStatus_t CInt_destroyBasisSet( BasisSet_t basis );

CIntStatus_t CInt_unpackBasisSet( BasisSet_t basis,
                                  void *buf);

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif

CIntStatus_t CInt_importBasisSet(BasisSet_t basis,
                                 int natoms, int *Zs,
                                 double *X, double *Y, double *Z,
                                 int nprims, int nshells, int pure,
                                 int *shells_p_atom,
                                 int *prims_p_shell,
                                 int *L, double *cc, double *alpha);

CIntStatus_t CInt_loadBasisSet( BasisSet_t basis,
                                char *bsfile,
                                char *xyzfile );

CIntStatus_t CInt_packBasisSet( BasisSet_t basis,
                                void **buf,
                                int *bufsize );

int CInt_getNumShells( BasisSet_t basis );

int CInt_getNumFuncs( BasisSet_t basis );

int CInt_getNumAtoms( BasisSet_t basis );

int CInt_getMaxShellDim( BasisSet_t basis );

int CInt_getNumOccOrb( BasisSet_t basis );

int CInt_getFuncEndInd( BasisSet_t basis,
                        int shellid );

int CInt_getAtomStartInd( BasisSet_t basis,
                          int atomid );

// one electron integrals

CIntStatus_t CInt_createOED( BasisSet_t basis,
                             OED_t *oed );

CIntStatus_t CInt_destroyOED( OED_t oed );
    
CIntStatus_t CInt_computePairKin( BasisSet_t basis,
                                  OED_t oed,
                                  int A,
                                  int B,
                                  double **integrals,
                                  int *nints );

CIntStatus_t CInt_computePairOvl( BasisSet_t basis,
                                  OED_t oed,
                                  int A,
                                  int B,
                                  double **integrals,
                                  int *nints );

CIntStatus_t CInt_computePairPot( BasisSet_t basis,
                                  OED_t oed,
                                  int A,
                                  int B,
                                  double **integrals,
                                  int *nints );

CIntStatus_t CInt_computePairCoreH( BasisSet_t basis,
                                    OED_t oed,
                                    int A,
                                    int B,
                                    double **integrals,
                                    int *nints );

void CInt_getShellxyz ( BasisSet_t basis,
                        int shellid,
                        double *x,
                        double *y,
                        double *z );

int CInt_getShellDim( BasisSet_t basis,
                      int shellid );

int CInt_getFuncStartInd( BasisSet_t basis,
                          int shellid );

void CInt_getInitialGuess( BasisSet_t basis,
                           int atomid,
                           double **guess,
                           int *spos,
                           int *epos );

int CInt_getTotalCharge(BasisSet_t basis);

int CInt_getNneutral(BasisSet_t basis);

double CInt_getNucEnergy (BasisSet_t basis);

// two electron integrals
#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

CIntStatus_t CInt_createERD( BasisSet_t basis,
                             ERD_t *erd,
                             int nthreads );

CIntStatus_t CInt_destroyERD( ERD_t erd );


CIntStatus_t CInt_computeShellQuartet( BasisSet_t basis,
                                       ERD_t erd,
                                       int tid,
                                       int A,
                                       int B,
                                       int C,
                                       int D,
                                       double **integrals,
                                       int *nints );

CIntStatus_t CInt_computeShellQuartets(BasisSet_t basis,
                                       ERD_t erd,
                                       uint32_t threadId,
                                       uint32_t shellIndixA,
                                       const uint32_t* shellIndicesB,
                                       uint32_t shellIndixC,
                                       const uint32_t* shellIndicesD,
                                       uint32_t shellIndicesCount,
                                       double **integrals,
                                       int *integralsCount);

void CInt_getMaxMemory( ERD_t erd,
                        double *memsize );

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif


#ifdef __INTEL_OFFLOAD
CIntStatus_t CInt_offload_createBasisSet( BasisSet_t *_basis );

CIntStatus_t CInt_offload_destroyBasisSet( BasisSet_t basis );

CIntStatus_t CInt_offload_pushBasisSet( BasisSet_t basis );

CIntStatus_t CInt_offload_loadBasisSet( BasisSet_t basis,
                                        char *bsfile,
                                        char *molfile );

CIntStatus_t CInt_offload_createERD( BasisSet_t basis,
                                     ERD_t *erd,
                                     int nthreads, int nthreads_mic);

CIntStatus_t CInt_offload_destroyERD( ERD_t erd );

void CInt_offload_getMaxMemory( ERD_t erd,
                                double *mem_mic,
                                double *mem_cpu );
#endif

#ifdef __cplusplus
}
#endif

#define CINT_ASSERT(condition) if (!(condition)) { \
    dprintf(2, "ASSERTION FAILED: %s in %s:%d\n", #condition, __FILE__, __LINE__); \
    fsync(2); \
    abort(); \
}

#endif /* __CINT_H__ */
