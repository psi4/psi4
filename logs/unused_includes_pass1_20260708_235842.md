# Unused includes verification pass — unused_includes_pass1_20260708_235842

- Files processed: 20
- Files with >=1 include actually removed: 20
- Total includes removed (compile-verified): 163
- Total includes kept (false positive / needed): 26
- Branch: `unused-includes-cleanup` (based on `master` @ 1496b4e98, plus cherry-picked `26ad5540a` timer.cc fix)

## End-of-batch validation

- `ninja core` in `objdir_psidev/psi4-core-prefix/src/psi4-core-build`: succeeded (30/30 build steps, link OK, only pre-existing deprecation warnings unrelated to this change).
- `ninja install` to refresh `stage/lib/psi4/core*.so`: succeeded.
- `ctest` smoke tests, all passed:
  - `scf1` (2.42s)
  - `tu1-h2o-energy` (0.92s)
  - `dlpnomp2-1` (1.49s) — exercises `dlpno/mp2.cc`, `dlpno/triples.cc`
  - `psimrcc-ccsd_t-1` (1.24s) — exercises `psimrcc/*.cc`
  - `scf-guess-read1` (1.02s)
  - `scf2` (1.10s) — exercises `libfock/*JK.cc`, `libscf_solver/rhf.cc`, `scfgrad/scf_grad.cc`

All 6 tests passed with no regressions.

| File | Removed | Kept (needed) | Commit |
|---|---|---|---|
| `/home/tamago/psi4/psi4/src/psi4/psimrcc/transform_read_so.cc` | 9 | 3 | 80aa81a5d |
| `/home/tamago/psi4/psi4/src/psi4/occ/omp3_ip_poles.cc` | 10 | 1 | 252ee937c |
| `/home/tamago/psi4/psi4/src/psi4/occ/omp2_ip_poles.cc` | 10 | 1 | 6402d3372 |
| `/home/tamago/psi4/psi4/src/psi4/libfock/DirectJK.cc` | 11 | 0 | 45b14a513 |
| `/home/tamago/psi4/psi4/src/psi4/scfgrad/scf_grad.cc` | 6 | 4 | fea3e053d |
| `/home/tamago/psi4/psi4/src/psi4/psimrcc/transform.cc` | 8 | 2 | 2adebd32a |
| `/home/tamago/psi4/psi4/src/psi4/libfock/snLinK.cc` | 9 | 1 | 227cc22ef |
| `/home/tamago/psi4/psi4/src/psi4/libfock/apps.cc` | 7 | 3 | d74d889a4 |
| `/home/tamago/psi4/psi4/src/psi4/dlpno/mp2.cc` | 9 | 1 | d5149cfc8 |
| `/home/tamago/psi4/psi4/src/psi4/psimrcc/sort.cc` | 7 | 2 | fb483b97b |
| `/home/tamago/psi4/psi4/src/psi4/libscf_solver/rhf.cc` | 8 | 1 | 27fdf993e |
| `/home/tamago/psi4/psi4/src/psi4/liboptions/print.cc` | 8 | 1 | eba692cc3 |
| `/home/tamago/psi4/psi4/src/psi4/libfock/wrapper.cc` | 9 | 0 | 0d9517bd0 |
| `/home/tamago/psi4/psi4/src/psi4/libfock/MemDFJK.cc` | 8 | 1 | 49c4a6a8a |
| `/home/tamago/psi4/psi4/src/psi4/dlpno/triples.cc` | 7 | 2 | 975d18c8d |
| `/home/tamago/psi4/psi4/src/psi4/psimrcc/transform_mrpt2.cc` | 6 | 2 | 8b0d1c213 |
| `/home/tamago/psi4/psi4/src/psi4/psimrcc/main.cc` | 7 | 1 | 70f09a84e |
| `/home/tamago/psi4/psi4/src/psi4/libfock/jk.cc` | 8 | 0 | 8c04c445e |
| `/home/tamago/psi4/psi4/src/psi4/libfock/PKJK.cc` | 8 | 0 | bf3f3473e |
| `/home/tamago/psi4/psi4/src/psi4/libfock/CDJK.cc` | 8 | 0 | beb2d4341 |

## `/home/tamago/psi4/psi4/src/psi4/psimrcc/transform_read_so.cc`

Removed:
- `"psi4/psifiles.h"` (was line 57)
- `"psi4/libqt/qt.h"` (was line 56)
- `"psi4/libiwl/iwl.h"` (was line 55)
- `"psi4/libpsio/psio.hpp"` (was line 54)
- `"matrix.h"` (was line 40)
- `"blas.h"` (was line 38)
- `"algebra_interface.h"` (was line 37)
- `"psi4/libpsi4util/libpsi4util.h"` (was line 33)
- `<algorithm>` (was line 30)

Kept (compile failed without it):
- `"psi4/libmints/matrix.h"` (line 35)
- `"psi4/libmints/mintshelper.h"` (line 34)
- `"psi4/libmoinfo/libmoinfo.h"` (line 32)

## `/home/tamago/psi4/psi4/src/psi4/occ/omp3_ip_poles.cc`

Removed:
- `"psi4/libtrans/mospace.h"` (was line 47)
- `"psi4/libqt/qt.h"` (was line 46)
- `"psi4/libiwl/iwl.h"` (was line 45)
- `"psi4/libpsio/psio.hpp"` (was line 44)
- `"psi4/libpsio/psio.h"` (was line 43)
- `"psi4/libciomr/libciomr.h"` (was line 42)
- `<iomanip>` (was line 37)
- `<fstream>` (was line 35)
- `<sstream>` (was line 34)
- `<iostream>` (was line 30)

Kept (compile failed without it):
- `"psi4/libmints/factory.h"` (line 52)

## `/home/tamago/psi4/psi4/src/psi4/occ/omp2_ip_poles.cc`

Removed:
- `"psi4/libtrans/mospace.h"` (was line 47)
- `"psi4/libqt/qt.h"` (was line 46)
- `"psi4/libiwl/iwl.h"` (was line 45)
- `"psi4/libpsio/psio.hpp"` (was line 44)
- `"psi4/libpsio/psio.h"` (was line 43)
- `"psi4/libciomr/libciomr.h"` (was line 42)
- `<iomanip>` (was line 37)
- `<fstream>` (was line 35)
- `<sstream>` (was line 34)
- `<iostream>` (was line 30)

Kept (compile failed without it):
- `"psi4/libmints/factory.h"` (line 52)

## `/home/tamago/psi4/psi4/src/psi4/libfock/DirectJK.cc`

Removed:
- `"psi4/libpsi4util/PsiOutStream.h"` (was line 51)
- `<unordered_set>` (was line 50)
- `<sstream>` (was line 49)
- `<limits>` (was line 48)
- `<algorithm>` (was line 47)
- `"psi4/libiwl/iwl.hpp"` (was line 36)
- `"psi4/psifiles.h"` (was line 35)
- `"psi4/libpsio/aiohandler.h"` (was line 32)
- `"psi4/libpsio/psio.h"` (was line 31)
- `"psi4/libpsio/psio.hpp"` (was line 30)
- `"psi4/lib3index/3index.h"` (was line 29)

## `/home/tamago/psi4/psi4/src/psi4/scfgrad/scf_grad.cc`

Removed:
- `"psi4/libpsi4util/PsiOutStream.h"` (was line 57)
- `"psi4/libdisp/dispersion.h"` (was line 53)
- `"psi4/libpsio/psio.hpp"` (was line 39)
- `<omp.h>` (was line 35)
- `<sstream>` (was line 33)
- `<algorithm>` (was line 31)

Kept (compile failed without it):
- `"psi4/libscf_solver/rhf.h"` (line 56)
- `"psi4/libscf_solver/uhf.h"` (line 55)
- `"psi4/libfunctional/superfunctional.h"` (line 52)
- `"psi4/libmints/molecule.h"` (line 44)

## `/home/tamago/psi4/psi4/src/psi4/psimrcc/transform.cc`

Removed:
- `"matrix.h"` (was line 55)
- `"index.h"` (was line 54)
- `"algebra_interface.h"` (was line 52)
- `"psi4/libqt/qt.h"` (was line 49)
- `"psi4/libpsio/psio.hpp"` (was line 47)
- `"psi4/libciomr/libciomr.h"` (was line 46)
- `"psi4/libpsi4util/libpsi4util.h"` (was line 33)
- `<algorithm>` (was line 30)

Kept (compile failed without it):
- `"blas.h"` (line 53)
- `"psi4/libmoinfo/libmoinfo.h"` (line 32)

## `/home/tamago/psi4/psi4/src/psi4/libfock/snLinK.cc`

Removed:
- `<libint2.hpp>` (was line 57)
- `<omp.h>` (was line 54)
- `<variant>` (was line 51)
- `<unordered_set>` (was line 50)
- `<tuple>` (was line 49)
- `<map>` (was line 48)
- `<iostream>` (was line 47)
- `<algorithm>` (was line 46)
- `"psi4/libpsi4util/PsiOutStream.h"` (was line 43)

Kept (compile failed without it):
- `"psi4/pybind11.h"` (line 44)

## `/home/tamago/psi4/psi4/src/psi4/libfock/apps.cc`

Removed:
- `<omp.h>` (was line 58)
- `<tuple>` (was line 54)
- `<functional>` (was line 53)
- `<algorithm>` (was line 52)
- `"psi4/libpsi4util/PsiOutStream.h"` (was line 48)
- `"psi4/libscf_solver/rhf.h"` (was line 40)
- `"psi4/physconst.h"` (was line 38)

Kept (compile failed without it):
- `"psi4/libmints/basisset.h"` (line 44)
- `"psi4/libmints/molecule.h"` (line 43)
- `"psi4/libmints/factory.h"` (line 42)

## `/home/tamago/psi4/psi4/src/psi4/dlpno/mp2.cc`

Removed:
- `<omp.h>` (was line 52)
- `"psi4/libpsi4util/PsiOutStream.h"` (was line 45)
- `"psi4/libmints/twobody.h"` (was line 43)
- `"psi4/libmints/orthog.h"` (was line 42)
- `"psi4/libmints/molecule.h"` (was line 41)
- `"psi4/libmints/mintshelper.h"` (was line 40)
- `"psi4/libmints/integral.h"` (was line 37)
- `"psi4/libfock/points.h"` (was line 35)
- `"psi4/lib3index/3index.h"` (was line 32)

Kept (compile failed without it):
- `"psi4/libmints/basisset.h"` (line 36)

## `/home/tamago/psi4/psi4/src/psi4/psimrcc/sort.cc`

Removed:
- `"psi4/psifiles.h"` (was line 43)
- `"psi4/libqt/qt.h"` (was line 42)
- `"psi4/libiwl/iwl.h"` (was line 41)
- `"psi4/libpsio/psio.hpp"` (was line 40)
- `"psi4/libciomr/libciomr.h"` (was line 39)
- `"psi4/libpsi4util/libpsi4util.h"` (was line 38)
- `<algorithm>` (was line 34)

Kept (compile failed without it):
- `"blas.h"` (line 45)
- `"psi4/libmoinfo/libmoinfo.h"` (line 36)

## `/home/tamago/psi4/psi4/src/psi4/libscf_solver/rhf.cc`

Removed:
- `"psi4/libpsio/psio.h"` (was line 57)
- `"psi4/libpsi4util/PsiOutStream.h"` (was line 55)
- `"psi4/libmints/pointgrp.h"` (was line 53)
- `"psi4/libiwl/iwl.hpp"` (was line 49)
- `"psi4/physconst.h"` (was line 42)
- `<omp.h>` (was line 38)
- `<algorithm>` (was line 30)
- `<any>` (was line 29)

Kept (compile failed without it):
- `"psi4/libmints/factory.h"` (line 50)

## `/home/tamago/psi4/psi4/src/psi4/liboptions/print.cc`

Removed:
- `"psi4/libpsi4util/process.h"` (was line 45)
- `"psi4/libpsi4util/libpsi4util.h"` (was line 43)
- `"psi4/libpsi4util/exception.h"` (was line 42)
- `<memory>` (was line 41)
- `"psi4/pragma.h"` (was line 40)
- `<algorithm>` (was line 38)
- `<stdexcept>` (was line 33)
- `<cstddef>` (was line 32)

Kept (compile failed without it):
- `"psi4/libpsi4util/PsiOutStream.h"` (line 44)

## `/home/tamago/psi4/psi4/src/psi4/libfock/wrapper.cc`

Removed:
- `<sstream>` (was line 47)
- `<iostream>` (was line 43)
- `"psi4/libpsi4util/process.h"` (was line 41)
- `"psi4/libpsi4util/PsiOutStream.h"` (was line 40)
- `"psi4/psi4-dec.h"` (was line 38)
- `"psi4/libiwl/iwl.h"` (was line 35)
- `"psi4/libpsio/psio.hpp"` (was line 34)
- `"psi4/libpsio/psio.h"` (was line 33)
- `"psi4/psifiles.h"` (was line 31)

## `/home/tamago/psi4/psi4/src/psi4/libfock/MemDFJK.cc`

Removed:
- `"psi4/libpsi4util/process.h"` (was line 49)
- `<omp.h>` (was line 48)
- `"psi4/libpsi4util/PsiOutStream.h"` (was line 46)
- `<sstream>` (was line 45)
- `"psi4/psifiles.h"` (was line 34)
- `"psi4/libpsio/aiohandler.h"` (was line 31)
- `"psi4/libpsio/psio.h"` (was line 30)
- `"psi4/libpsio/psio.hpp"` (was line 29)

Kept (compile failed without it):
- `"psi4/libmints/vector.h"` (line 37)

## `/home/tamago/psi4/psi4/src/psi4/dlpno/triples.cc`

Removed:
- `"psi4/libpsi4util/PsiOutStream.h"` (was line 45)
- `"psi4/libmints/twobody.h"` (was line 43)
- `"psi4/libmints/orthog.h"` (was line 42)
- `"psi4/libmints/mintshelper.h"` (was line 40)
- `"psi4/libmints/integral.h"` (was line 37)
- `"psi4/libdiis/diismanager.h"` (was line 33)
- `"psi4/lib3index/3index.h"` (was line 32)

Kept (compile failed without it):
- `"psi4/libmints/molecule.h"` (line 41)
- `"psi4/libmints/basisset.h"` (line 36)

## `/home/tamago/psi4/psi4/src/psi4/psimrcc/transform_mrpt2.cc`

Removed:
- `"psi4/libiwl/iwl.h"` (was line 52)
- `"psi4/libciomr/libciomr.h"` (was line 50)
- `"algebra_interface.h"` (was line 36)
- `"psi4/libpsi4util/libpsi4util.h"` (was line 35)
- `"matrix.h"` (was line 34)
- `<algorithm>` (was line 30)

Kept (compile failed without it):
- `"blas.h"` (line 37)
- `"psi4/libmoinfo/libmoinfo.h"` (line 32)

## `/home/tamago/psi4/psi4/src/psi4/psimrcc/main.cc`

Removed:
- `"transform.h"` (was line 58)
- `"mrcc.h"` (was line 56)
- `"sort.h"` (was line 55)
- `"main.h"` (was line 54)
- `"psi4/libciomr/libciomr.h"` (was line 49)
- `<complex>` (was line 45)
- `<iostream>` (was line 44)

Kept (compile failed without it):
- `"blas.h"` (line 53)

## `/home/tamago/psi4/psi4/src/psi4/libfock/jk.cc`

Removed:
- `<omp.h>` (was line 51)
- `"psi4/libpsi4util/PsiOutStream.h"` (was line 45)
- `"psi4/libiwl/iwl.hpp"` (was line 38)
- `"psi4/psifiles.h"` (was line 37)
- `"psi4/libpsio/aiohandler.h"` (was line 34)
- `"psi4/libpsio/psio.h"` (was line 33)
- `"psi4/libpsio/psio.hpp"` (was line 32)
- `"psi4/lib3index/3index.h"` (was line 31)

## `/home/tamago/psi4/psi4/src/psi4/libfock/PKJK.cc`

Removed:
- `<omp.h>` (was line 46)
- `"psi4/libpsi4util/PsiOutStream.h"` (was line 44)
- `<sstream>` (was line 43)
- `"psi4/libmints/mintshelper.h"` (was line 40)
- `"psi4/libmints/matrix.h"` (was line 38)
- `"psi4/libiwl/iwl.hpp"` (was line 35)
- `"psi4/libpsio/aiohandler.h"` (was line 31)
- `"psi4/libpsio/psio.h"` (was line 30)

## `/home/tamago/psi4/psi4/src/psi4/libfock/CDJK.cc`

Removed:
- `<omp.h>` (was line 47)
- `"psi4/libpsi4util/PsiOutStream.h"` (was line 45)
- `<sstream>` (was line 44)
- `"psi4/libiwl/iwl.hpp"` (was line 35)
- `"psi4/psifiles.h"` (was line 34)
- `"psi4/libpsio/aiohandler.h"` (was line 31)
- `"psi4/libpsio/psio.h"` (was line 30)
- `"psi4/libpsio/psio.hpp"` (was line 29)
