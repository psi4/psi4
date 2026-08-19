list(APPEND sources
  aio_handler.cc
  change_namespace.cc
  close.cc
  compose_err_msg.cc
  done.cc
  decode_errno.cc
  error.cc
  filemanager.cc
  filescfg.cc
  get_address.cc
  get_filename.cc
  get_global_address.cc
  get_numvols.cc
  get_volpath.cc
  getpid.cc
  init.cc
  open.cc
  open_check.cc
  read.cc
  read_entry.cc
  rename_file.cc
  rw.cc
  tocclean.cc
  toclast.cc
  toclen.cc
  tocprint.cc
  tocread.cc
  tocscan.cc
  tocwrite.cc
  volseek.cc
  write.cc
  write_entry.cc
  )
psi4_add_module(lib psio sources)

# The standalone round-trip test links a curated subset of static archives plus a
# small stubs TU. That resolves with GNU/Clang on Linux/macOS but not with
# clang-cl/lld-link, which emits inline and implicit functions that GCC discards
# as unused, so the Windows objects carry undefined references the curated link
# line does not provide: psio.lib wants DPDMOSpace::~DPDMOSpace (instantiated by
# std::vector<DPDMOSpace> in dpd.h, reached via close.cc), and psi4util.lib wants
# Molecule::Molecule(const Molecule&) (from inline Molecule::clone) and
# ExternalPotential::print (from inline py_print). Not a unity or PCH artifact --
# Linux is unity-built too, and Azure Windows builds with PCH off. Restrict the
# test to non-Windows rather than drag in the full libmints/Libint closure.
if(NOT WIN32)
    add_subdirectory(test)
endif()

if (psi4_ENABLE_PRECOMPILE_HEADERS)
    # Clearly from the paths, libpsio is a weird place to be defining the pch. However, if one sensibly
    #   puts it in mints, CMake creates a dependency of, say, fnocc on mints, and since mints is the
    #   slowest to build, enabling pch slows the build. Revisit if https://gitlab.kitware.com/cmake/cmake/-/work_items/25392 solved.
    # Demonstrated to lengthen build (# inclusions): <vector> (138), <string> (141), "../../../include/psi4/pybind11.h",
    #   "../libpsi4util/libpsi4util.h" (83), "../libiwl/iwl.hpp" (58), "../libmints/vector.h" (108),
    #   "../libmints/integral.h" (92), "../libmints/wavefunction.h" (106), "../../../include/psi4/psifiles.h" (228),
    # Optimal pch Apr 2026: 173MB = libciomr, psio, qt, matrix, molecule, basisset
    target_precompile_headers(psio PRIVATE
        "../libciomr/libciomr.h"  # (333)
        "psio.hpp"                # (198)
        "../libqt/qt.h"           # (529)
        "../libmints/matrix.h"    # (272)
        "../libmints/molecule.h"  # (116)
        "../libmints/basisset.h"  # (120)
    )
endif()
