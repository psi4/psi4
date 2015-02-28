#ifndef __BASISSET_H__
#define __BASISSET_H__


#include "cint_def.h"
struct BasisSet;
typedef struct BasisSet* BasisSet_t;

void _maxMomentum (BasisSet_t basis, int *max_momentum);

void _maxPrimid (BasisSet_t basis, int *max_primid);

void _maxnumExp (BasisSet_t basis, int *max_nexp);

CIntStatus_t import_basis (char *file, BasisSet_t basis);

CIntStatus_t import_molecule (char *file, BasisSet_t basis);

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

CIntStatus_t parse_molecule (BasisSet_t basis);

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif


#endif /* __BASISSET_H__ */
