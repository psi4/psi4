# AGENTS.md

Guidance for AI coding agents (and new human contributors) working in the Psi4
repository. This file is intentionally practical: it covers how the project is
laid out, how to build/run/test it, and the conventions and gotchas that are not
obvious from the code. For deep user- and developer-facing documentation, defer
to the manual rather than duplicating it here:

- User manual & programmer's guide: https://psicode.org/psi4manual/master/index.html
- Contributing: [.github/CONTRIBUTING.md](.github/CONTRIBUTING.md)
- Forum (support/discussion): https://forum.psicode.org
- Developer Slack (invite): https://join.slack.com/t/psi4/shared_invite/zt-5s36s4rb-SQH6_AWyfWOqlKYN3cFs4Q

## What Psi4 is

Psi4 is an open-source suite of *ab initio* quantum chemistry programs for
high-accuracy simulation of molecular properties (energies, gradients,
Hessians, and many one- and two-electron properties). It implements methods
across the accuracy/cost spectrum — Hartree–Fock and DFT, many-body
perturbation theory (MP2–MP4), coupled-cluster (through CCSD(T) and beyond),
intermolecular interaction analysis (SAPT, including DFT-SAPT / F-SAPT), and
configuration interaction and multireference methods (CASSCF, FCI, MRCC).
Density fitting (resolution-of-the-identity) is used heavily for speed — it is
the default for MP2.5 and lower-cost methods — though not universally.
Computations with >2500 basis functions on multi-core machines are routine.

Reference publication (cite this):
D. G. A. Smith *et al.*, "Psi4 1.4: Open-Source Software for High-Throughput
Quantum Chemistry", *J. Chem. Phys.* **152**, 184108 (2020).
https://doi.org/10.1063/5.0006002

## Architecture at a glance

Psi4 is a **hybrid C++/Python** package, and needs *both* halves to be a usable
quantum-chemistry program — the C++ compiled alone is not.

- **C++** (`psi4/src/psi4/`) holds the compute-heavy kernels: integrals, SCF,
  correlated methods, tensor machinery. Performance-critical.
- **Pybind11** exports C++ classes/functions into the `psi4.core` Python module.
  Bindings live in `psi4/src/export_*.cc` and `psi4/src/core.cc`. This is the
  C++↔Python boundary.
- **Python driver** (`psi4/driver/`) orchestrates everything: input parsing,
  method dispatch, composite/CBS extrapolation, finite differences, many-body
  expansion, and the user-facing API. Most *logic about what to run* lives here;
  most *number crunching* lives in C++.

There are **three user interfaces**, and tests/samples exist for them:

- **PSIthon** — a Python superset with chemistry sugar (`molecule {...}` blocks,
  `set` keyword syntax). Run with `psi4 input.in` (input/output file extensions
  are conventional, not enforced).
- **PsiAPI** — plain Python: `import psi4`. Run with `python job.py`.
- **QCSchema** — a single-step job driven by a JSON/MessagePack payload:
  `psi4 --qcschema in.json`.

## Repository layout

```
psi4/
  src/psi4/          C++ source, organized into ~40 modules (core libs below)
  src/export_*.cc    Pybind11 bindings (C++ -> psi4.core)
  src/core.cc        Pybind11 module definition
  src/read_options.cc  Central registry of ALL keywords/options and defaults
  driver/            Python driver (orchestration + user API)
    driver.py          energy()/gradient()/optimize()/frequency() entry points
    procrouting/       method dispatch: proc.py, proc_table.py, dft/, sapt/, ...
    p4util/            helper utilities
  share/psi4/         runtime data: basis sets (basis/), integration grids,
                      databases, and Python scripts/data
  include/            public C++ headers, physical constants
tests/               test suite
  pytests/           pytest-based tests (PsiAPI); conftest.py, markers
  <name>/            individual CTest cases (input.dat + expected refs)
doc/sphinxman/       Sphinx manual sources (prog_*.rst = programmer's guide)
cmake/, external/    build system and bundled/fetched dependencies
conda/, devtools/    packaging and CI helper scripts
```

### C++ core libraries (`psi4/src/psi4/`) — orientation

These underpin every method. The individual method modules (SCF, MP2, coupled
cluster, SAPT, CI, ...) are best learned from the manual rather than enumerated
here.

- `libmints` — molecules, basis sets, integrals, `Matrix`/`Vector`,
  `Wavefunction` (the central object passed around).
- `scf_solver`, `scfgrad` — SCF (HF/DFT) energies and gradients; foundational,
  as most correlated methods build on an SCF reference.
- `libfock` — J/K (Coulomb/exchange) build engines.
- `lib3index`, `libdpd`, `libtrans`, `libiwl` — density-fitting 3-index tensors,
  DPD block-tensor algebra, integral transformations, integrals-with-labels I/O.
- `libpsio` — binary scratch-file I/O; `liboptions` — options plumbing;
  `libqt`, `libpsi4util`, `libciomr` — math/utility helpers.
- `libfunctional`, `libdisp` — DFT exchange-correlation functionals and
  empirical dispersion.

## The central object: `Wavefunction`

`Wavefunction` (defined in `libmints`) is the object that ties Psi4 together. A
computation returns one, and it carries the molecule, basis set, reference
energy, orbitals (`Ca`/`Cb`), densities, Fock matrices, orbital energies, and
computed properties — the shared state passed between the driver, methods, and
post-processing. When you extend Psi4, you are almost always reading from or
writing to a `Wavefunction`. It is exported to Python, so most of it is reachable
as `wfn.Ca()`, `wfn.epsilon_a()`, `wfn.energy()`, etc.

A great, less-obvious way to learn this object and the exposed C++ API is the
**Psi4NumPy** project (https://github.com/psi4/psi4numpy), whose tutorial
notebooks reimplement methods in Python against Psi4's building blocks — e.g.
[the `Wavefunction` walkthrough](https://github.com/psi4/psi4numpy/blob/master/Tutorials/01_Psi4NumPy-Basics/1d_wavefunction.ipynb).

## Building

Psi4 builds with **CMake + Ninja** using system tools (compilers, Python,
BLAS/LAPACK) on any architecture, including Linux, macOS (Intel & Apple
Silicon), and Windows. In principle you supply the toolchain and dependencies
yourself.

In practice, developers use a **conda-centric bootstrap** that provisions all
build tools and dependencies, via
[`conda/psi4-path-advisor.py`](conda/psi4-path-advisor.py) (which reads the
dependency spec [codedeps.yaml](codedeps.yaml) and emits copy/paste-able
commands):

```bash
git clone https://github.com/psi4/psi4.git && cd psi4

# (1) create & activate a conda dev env with all deps. Invaluable — do not skip.
eval $(conda/psi4-path-advisor.py env)          # defaults: env "p4dev"

# (2) configure & build against that env. Often skipped in favor of simple `cmake -S. -Bobjdir_p4dev -GNinja`
eval $(conda/psi4-path-advisor.py cmake)        # defaults: build dir "objdir_p4dev"
```

- For build-configuration options, run `conda/psi4-path-advisor.py cmake --help`.
  The full menu of tweakable CMake options is documented in the top-level
  [CMakeLists.txt](CMakeLists.txt).
- Build products stage into `<objdir>/stage/` (a full install tree, see Gotchas).
  The runnable program is `<objdir>/stage/bin/psi4`; the importable package is
  `<objdir>/stage/lib/psi4`.
- Incremental rebuilds are just re-running `cmake --build <objdir>`.

## Running

First, set up your shell environment — this is broadly needed, not just for
PsiAPI. `psi4 --psiapi` emits the PATH/PYTHONPATH exports; without them even
PSIthon only works for the simplest inputs (those that don't dispatch to
`qcengine`, e.g. for dispersion or external gradients):

```bash
eval $(<objdir>/stage/bin/psi4 --psiapi)

psi4 input.in                             # PSIthon; writes input.out alongside
python job.py                             # PsiAPI (import psi4)
psi4 --qcschema in.json                   # QCSchema single-step job
```

Set a scratch location on fast local (non-network) disk:
`export PSI_SCRATCH=/path/to/local/scratch`. Useful `psi4` flags: `-n<N>`
(threads), `-o` (output file), `--scratch <dir>`, `-v` (verbose).

## Testing

Two complementary systems — **both** matter, and new features should add tests
in whichever fits:

**CTest**, run from the build dir, drives the PSIthon test inputs:

```bash
cd <objdir>
ctest -j$(getconf _NPROCESSORS_ONLN)      # all
ctest -R quick                            # common fast subset
ctest -R <regex>                          # by test name
ctest -L <label>                          # by label (e.g. scf, cc, sapt)
```

**pytest** drives the PsiAPI tests (config in [pytest.ini](pytest.ini)):

```bash
pytest -m "quick and not long"            # select by marker
pytest tests/pytests/test_<x>.py -k <name>

cd <objdir> && pytest ../tests/           # runs the FULL suite (CTest cases too)
```

Markers are defined in [pytest.ini](pytest.ini): scope markers like `quick`,
`smoke`, `long`, `slow`, `stdsuite`, plus per-method markers (`scf`, `cc`,
`mp2`, `sapt`, `casscf`, ...). Prefer running a relevant subset locally; the full
suite is long.

**Optional-dependency (addon) tests** are guarded by `using(<addon>)` /
`uusing(<addon>)`. For **pytest** this adapts to the active environment at
runtime (tests skip if the addon is absent). For **CTest** it is fixed at
configure time — a newly installed addon is not picked up until you re-run cmake.

## Style

Full conventions are in flux and will be revised soon, so keep it minimal for now:

- Line length ~120 characters.
- Every source file carries the LGPL-3.0 license header block — preserve it when
  editing and include it in new files.

## Making changes — where things live

- **Add/modify a keyword/option**: options and their defaults are registered in
  [psi4/src/read_options.cc](psi4/src/read_options.cc). Read from C++ via the
  options object; read from Python via `psi4.core.get_option(...)`.
- **Add/route a method**: method dispatch is in
  `psi4/driver/procrouting/` (`proc.py`, `proc_table.py`). The public entry
  points are `energy()`, `gradient()`, `optimize()`, `frequencies()`,
  `properties()` in [psi4/driver/driver.py](psi4/driver/driver.py).
- **Expose new C++ to Python**: add bindings in the relevant
  `psi4/src/export_*.cc` (and `psi4/src/core.cc`).
- **Add basis sets / grids / databases**: under `psi4/share/psi4/`. DFT
  functional definitions live in `psi4/driver/procrouting/dft/`.
- New features **must** include tests and documentation (per CONTRIBUTING.md).
  Deeper how-to guides: the `prog_*` pages of the manual, e.g.
  https://psicode.org/psi4manual/master/prog_newcode.html

## Gotchas

- **`<objdir>/stage/` is a full installation.** There is no need to
  `cmake install`; most developers just work directly against the staged tree.
- **Even Python-only edits require a rebuild.** The staged files in
  `<objdir>/stage/lib/psi4/` are *copies*, not symlinks. After editing anything
  under `psi4/driver/`, run `cmake --build <objdir>` to copy the `.py` into the
  position where the C++/staged package can find it. Editing a staged file
  directly is silently clobbered on the next build — always edit under `psi4/`.
- **Git tags must be even with the upstream repo** for the version to compute
  correctly. A shallow clone or missing tags yields a wrong/undefined version.
- **Symmetry**: molecules carry point-group symmetry; many quantities are stored
  per-irrep (e.g. `Dimension`, symmetry-blocked matrices). Don't assume C1.
- **Density fitting is the default for MP2.5 and lower-cost methods** (but not
  universally); "conventional" (pk/direct) paths are separate code paths and
  separate keywords.
- **Scratch files**: methods read/write large binary scratch via `libpsio`.
  Point `PSI_SCRATCH` at fast local storage; stale scratch can corrupt restarts.
- **Two test systems** (CTest PSIthon + pytest PsiAPI) — a change may need
  updating in both, and reference values live next to each test.
- **Generated/staged artifacts**: never commit anything under `objdir*/`,
  `stage/`, or `*.out` scratch outputs.

## Contributing workflow

Fork → feature branch → PR against `psi4/psi4` `master`. PRs run CI (Azure
Pipelines + GitHub Actions) and require passing tests + core-developer review.
**Rebase onto upstream is preferred over merge commits** to keep history linear.
See [.github/CONTRIBUTING.md](.github/CONTRIBUTING.md) and the
[PR template](.github/PULL_REQUEST_TEMPLATE.md).
Importrant: AI agents that write PR descriptions, code reviews, GitHub comments or responses must identify themselves as such, by ending the text with "Generated by: <LLM model name>".

# Local overrides

Developers may create `AGENTS.local.md` for machine-, workflow-, or preference-specific notes. This file is ignored by git and should not contain instructions required for normal project development.
