#ifndef __PFOCK_DEF_H__
#define __PFOCK_DEF_H__


#include "pfock_type.h"


typedef struct CoreH *CoreH_t;
typedef struct Ovl *Ovl_t;    
typedef struct PFock *PFock_t;


typedef enum
{
    PFOCK_STATUS_SUCCESS = 0,
    PFOCK_STATUS_INIT_FAILED = 1,
    PFOCK_STATUS_ALLOC_FAILED = 2,
    PFOCK_STATUS_INVALID_VALUE = 3,
    PFOCK_STATUS_EXECUTION_FAILED = 4,
    PFOCK_STATUS_INTERNAL_ERROR = 5
} PFockStatus_t;


typedef enum
{
    PFOCK_MAT_TYPE_D = 0,
    PFOCK_MAT_TYPE_F = 1,
    PFOCK_MAT_TYPE_J = 2,
    PFOCK_MAT_TYPE_K = 3
} PFockMatType_t;


// must be called by all processes "loosely synchronously"
PFockStatus_t PFock_create( BasisSet_t basis,
                            int nprow,      // process topology is nprow by npcol
                            int npcol,
                            int ntasks,     // number of tasks
                            int maxnumdmat,
                            int symm,       // 0 non symmetric, others symmetric
                            double tolscr,  // screening threshold
                            PFock_t *pfock );

PFockStatus_t PFock_destroy( PFock_t pfock );

PFockStatus_t PFockGetProcessCoord( PFock_t pfock,
                                    int *rowid,
                                    int *colid );

PFockStatus_t PFock_setNumDenMat( PFock_t pfock,
                                  int numdmat );

PFockStatus_t PFock_putDenMat( PFock_t pfock,
                               int index, // 0 to num-1, index of density matrix
                               int rowstart,
                               int rowend,
                               int colstart,
                               int colend,
                               double *dmat,
                               int stride );

PFockStatus_t PFock_fillDenMat( PFock_t pfock,
                                int index, // 0 to num-1, index of density matrix
                                double value );

// called loosely synchronously by all processes
//  after all PFockPutDenMat calls for one density matrix
PFockStatus_t PFock_commitDenMats( PFock_t pfock );

PFockStatus_t PFock_getMat( PFock_t pfock,
                            int index,
                            PFockMatType_t type,
                            int rowstart,
                            int rowend,
                            int colstart,
                            int colend,
                            double *mat,
                            int stride );

// return local indices
PFockStatus_t PFock_getLocalMatInds( PFock_t pfock,
                                     int *rowstart,
                                     int *rowend,
                                     int *colstart,
                                     int *colend );

PFockStatus_t PFock_getLocalMatPtr ( PFock_t pfock,
                                     int index,
                                     PFockMatType_t type,
                                     int *rowstart,
                                     int *rowend,
                                     int *colstart,
                                     int *colend,
                                     double **mat,
                                     int *stride );

PFockStatus_t PFock_getMatGAHandle( PFock_t pfock,
                                    int index,
                                    PFockMatType_t type,
                                    int *ga );

// called loosely synchronously by all processes
//  to compute all J and K matrices
PFockStatus_t PFock_computeFock ( PFock_t pfock,
                                  BasisSet_t basis );

PFockStatus_t PFock_createCoreHMat ( PFock_t pfock,
                                     BasisSet_t basis,
                                     CoreH_t *hmat );

PFockStatus_t PFock_destroyCoreHMat ( CoreH_t hmat );

PFockStatus_t PFock_getCoreHMatGAHandle( CoreH_t hmat,
                                         int *ga );

PFockStatus_t PFock_getCoreHMat( CoreH_t hmat,
                                 int rowstart,
                                 int rowend,
                                 int colstart,
                                 int colend,
                                 double *mat,
                                 int stride );

PFockStatus_t PFock_createOvlMat ( PFock_t pfock,
                                   BasisSet_t basis,
                                   Ovl_t *omat );

PFockStatus_t PFock_destroyOvlMat ( Ovl_t omat );

PFockStatus_t PFock_getOvlMatGAHandle( Ovl_t omat,
                                       int *ga );
                                      
PFockStatus_t PFock_getOvlMat( Ovl_t omat,
                               int rowstart,
                               int rowend,
                               int colstart,
                               int colend,
                               double *mat,
                               int stride );

PFockStatus_t PFock_GAInit( int numFuncs,
                            int nprow,
                            int npcol,
                            int numdenmat,
                            int sizeheap,
                            int sizestack );

PFockStatus_t PFock_GAFinalize ( void );


#endif /* #ifndef __PFOCK_DEF_H__ */
