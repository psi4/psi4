#ifndef __CINT_DEF_H__
#define __CINT_DEF_H__

#include "cint_type.h"

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
    CINT_STATUS_FILEIO_FAILED = 6
} CIntStatus_t;


// basisset
// TODO: develop basisset parser

CIntStatus_t CInt_createBasisSet( BasisSet_t *basis );

CIntStatus_t CInt_loadBasisSet( BasisSet_t basis,
                                char *bsfile,
                                char *xyzfile );

CIntStatus_t CInt_destroyBasisSet( BasisSet_t basis );

CIntStatus_t CInt_packBasisSet( BasisSet_t basis,
                                void **buf,
                                int *bufsize );

CIntStatus_t CInt_unpackBasisSet( BasisSet_t basis,
                                void *buf);

int CInt_getNumShells( BasisSet_t basis );

int CInt_getNumFuncs( BasisSet_t basis );

int CInt_getNumAtoms( BasisSet_t basis );

int CInt_getShellDim( BasisSet_t basis,
                      int shellid );

int CInt_getMaxShellDim( BasisSet_t basis );

int CInt_getNumOccOrb( BasisSet_t basis );

int CInt_getFuncStartInd( BasisSet_t basis,
                          int shellid );

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


// two electron integrals

CIntStatus_t CInt_createERD( BasisSet_t basis,
                             ERD_t *erd );

CIntStatus_t CInt_destroyERD( ERD_t erd );


CIntStatus_t CInt_computeShellQuartet( BasisSet_t basis,
                                       ERD_t erd, 
                                       int A,
                                       int B,
                                       int C,
                                       int D,
                                       double **integrals,
                                       int *nints );


#endif /* __CINT_DEF_H__ */
