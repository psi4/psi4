# Standard Suite ("stdsuite") Guide for AI Agents

This guide is for agents (and humans) extending psi4's standard-suite regression tests: new methods, new
routes through psi4 modules, or new molecule/basis systems. Read it before editing any of the files below.
It records conventions that are not obvious from the code and lessons learned the hard way.

## 1. What the stdsuite tests

For every method x driver (energy/gradient/Hessian) x reference (RHF/UHF/ROHF) x frozen-core (AE/FC) x
algorithm (CONV/DF/CD) x psi4 module (occ, fnocc, ccenergy, dfmp2, detci, scf, ...) combination, a test
checks that psi4

* runs (or raises the *expected* error for not-yet-implemented combinations),
* routes to the expected module (`qc_module`),
* produces the right numbers: total energy/gradient/Hessian and the full set of QCVariables an ideal
  implementation of the method should set (the "contract"), each in both `psi4.core` globals and on the
  returned `Wavefunction`,
* agrees with conventional-integral reference values to a looser tolerance when run with DF/CD,
* (with `--runnonroutine`) gives analytic derivatives that match 5-point finite differences (FD).

Results also feed the capabilities tables in the docs (Sec. 9).

## 2. File map

| file | repo | role |
|---|---|---|
| `tests/pytests/test_standard_suite.py` | psi4 | ~120 `test_<method>_<driver>_<type>` functions; all parametrization lives here |
| `tests/pytests/standard_suite_runner.py` | psi4 | `runner_asserter()`: derives types, picks reference blocks, runs psi4, asserts, records |
| `tests/pytests/standard_suite_ref_local.py` | psi4 | psi4-local *additions to existing* QCEngine blocks; psi4-local contracts (`contractual_mp2f12`, `contractual_current`) |
| `tests/pytests/conftest.py` | psi4 | `--runnonroutine` option; per-test `clean_options()` |
| `pytest.ini` | psi4 | marker registration (`stdsuite`, `nonroutine`, `findif`, `long`, `d2ints`, `ecp`, ...) |
| `qcengine/programs/tests/standard_suite_ref.py` | QCEngine | `std_molecules`, `std_refs`, `_std_generics`, reference constants, `_std_suite` blocks, `compute_derived_qcvars`, `answer_hash` |
| `qcengine/programs/tests/standard_suite_contracts.py` | QCEngine | `contractual_<method>()` generators: which QCVariables are expected, with per-program excusals |
| `psi4/share/psi4/scripts/merge_stdsuite.py` | psi4 | merges `stdsuite_psi4.txt` records for the docs tables |

QCEngine's `standard_suite_runner.py`/`test_standard_suite.py` serve other programs and are not used by psi4.

**Which QCEngine?** psi4 imports whatever `qcengine` is on the path. Check with
`python -c "import qcengine; print(qcengine.__file__)"` before editing; edit *that* checkout.

**Local extension limits.** `standard_suite_ref_local.py` merges `data` into blocks that already exist in
QCEngine (matched by hash). It cannot add new systems/blocks; those must go into QCEngine
`standard_suite_ref.py` (then released, or the psi4 test env pointed at a QCEngine branch).

## 3. Reference data model (QCEngine `standard_suite_ref.py`)

* `std_molecules[name]`: psi4-style molecule strings. Tests build molecules with
  `psi4.core.Molecule.from_string(s, name=name)`; **generate references the same way** so orientation and
  frame match. Use non-equilibrium geometries so gradients are not ~0.
* `std_refs = ["rhf", "uhf", "rohf"]`; each `basis, subjects` param row lists three molecules
  `[rhf_mol, uhf_mol, rohf_mol]` (closed shell, open shell, open shell).
* `_std_generics["<system>_<basis>_<ae|fc>"] = (nbf, nmo, nalpha, nbeta)`. **Required for every block**
  (`compute_derived_qcvars` indexes it unconditionally).
* `_std_suite`: list of `{"meta": {...}, "data": {QCVariable: value}}` blocks. Key:
  `answer_hash(system, basis, scf_type, reference, fcae, corl_type, sdsc)` =
  `"<system>_<basis>_<scf_type>_<reference>_<fcae>_<corl_type>_<sdsc>"`.
  * `scf_type` in `pk|df|cd` (`direct`, `mem_df`, `disk_df` map onto these); `corl_type` in `conv|df|cd`.
  * Sections, marked `# <<<  CONV-AE-CONV  >>>` etc., read SCF-FC-CORL. Mixed ("lopsided") sections exist.
  * `sdsc`: `sd` (standard) or `sc` (semicanonical); only meaningful for ROHF. ROHF tests default to `sc`
    unless the param row sets `"xptd": {"sdsc": "sd"}`.
* Store only primary quantities (e.g., `X CORRELATION ENERGY`, `X SINGLES ENERGY`,
  `X SAME-SPIN CORRELATION ENERGY`, `(T) CORRECTION ENERGY`, gradients, Hessians);
  `compute_derived_qcvars` builds totals, doubles, opposite-spin, SCS, MP2.5, etc. Read it before adding a
  new method so you know which inputs it needs.
* HF/SCF quantities are repeated in every block of a system (via module-level constants like
  `_scf_<sys>_<bas>_<scf>_<ref>`, `_grad_scf_...`, `_hess_scf_...`).
* `_knownmissing`: placeholder for a CONV value that can't be obtained, so DF/CD comparison is skipped.
  Last resort; see the comment block above `_knownmissing`.
* Comment provenance on every number you add (program/module, analytic vs findif, any caveat).

## 4. How one test executes

```
test_<mtd>_<drv>_<type>(inp, dertype, basis, subjects, ...)
  -> _processor(): picks the molecule by reference, resolves per-dertype "error"/"wrong"/"xptd"/"marks",
                   sets keywords.basis and function_kwargs.dertype
  -> runner_asserter():
       * module skips (detci only for cc-pvdz; ECP-basis CD rows; ...)
       * derive corl_type from method + *_type keywords (hard-coded read_options defaults!), scf_type from
         keywords or natural pairing (conv<->pk); DFT and double hybrids special-cased
       * tolerances: atol_e/g/h (1e-6, 1e-6, 2e-6); FD Hessians 1e-5; atol_conv 1.5e-4 vs CONV refs
       * chash (exact-match block) and chash_conv (pk/conv block, same system/ref/fcae/sdsc)
       * set convergence options (see "reference generation" block), call psi4
       * "error" rows: assert exception type + regex, record "error", return
       * routing check: wfn.module() vs qc_module requested/expected
       * _asserter(contractual_<method>): for each expected QCVariable, compare to exact block (atol_*)
         AND to CONV block (atol_conv); for excused ones assert *absent*
       * "wrong" rows: expect an AssertionError containing the named QCVariable, then xfail
       * contractual_current (CURRENT ENERGY/GRADIENT/HESSIAN aliases), return values, generics
       * _recorder() appends a line to stdsuite_psi4.txt for docs
```

## 5. Test types and param rows

See the ADVICE comment at the top of `test_standard_suite.py`. Short version:

* `*_module`: most important. Every module route for the method, all three `basis, subjects` rows,
  analytic and FD dertypes. Each `pytest.param({"keywords": {...}, "error": {dertype: _pNN}, "wrong":
  {dertype: _wNN}, "xptd": {...}}, id="...")`.
* `*_default`: dz only. What happens when the user doesn't set `qc_module` and/or `*_type`.
* `*_scftype`: RHF dz only. Every `scf_type` under a correlated method.
* Error tuples `_pNN = (ExceptionType, regex, human_reason)`; `wrong` tuples `_wNN = (QCVariable, reason)`.
* dertypes: `grd1/grd0`, `hes2/hes1/hes0`. FD dertypes carry `findif` and usually `nonroutine` (skipped unless
  `--runnonroutine`). `hes2` usually carries `d2ints` (needs Libint2 second-derivative ERIs).
* **Param ids must not contain dashes** in the basis-row id; stdsuite tooling splits test ids on dashes
  (hence `id="def2ecp"`, not `"def2-svp"`).

## 6. Running

```bash
cd tests/pytests          # or objdir with ../tests/pytests/...
pytest test_standard_suite.py -k "test_mp2_ and module and def2ecp" -n 8
pytest test_standard_suite.py -k "test_mp2_ and module" --runnonroutine -n 8   # include FD checks
pytest test_standard_suite.py -m "ecp" --runnonroutine -n 12
```
* `--tb=line -rf` gives one line per failure; that's usually all you need, because the assertion messages name
  the QCVariable, the store (`P::e` vs `wfn`), and whether the miss is vs exact ref or `v. CONV`.
* Long FD jobs (CC FD gradients, FD Hessians): run in the background and poll by PID; don't `pgrep -f`
  for a pattern that also appears in your own wait loop's command line.

## 7. Workflows

### 7a. Adding a method or module route
1. Write the `*_module` test first (copy the closest method), then `*_default`, then `*_scftype`.
2. Check `contractual_<method>` exists in QCEngine; add or extend it. Excusals belong in the contract with the
   precise `qc_module`/reference/corl_type/method condition, never in the runner.
3. Add reference values (7c). If DF/CD values exist but no CONV, prefer finding CONV values over `_knownmissing`.
4. Add the method's branch in `runner_asserter.qcvar_assertions` and `corl_natural_values` if new.
5. Iterate with ADVICE: comment out lesser methods' contracts at first to focus, then restore.
6. Doing "good" things (enabling a route, setting a missing QCVariable) can make tests fail. That's a feature.

### 7b. Adding a system/basis (template: the ECP systems `ch2i2`/`ch2i` with `def2-svp`)
1. QCEngine: add molecules to `std_molecules`, `_std_generics` for ae and fc (check nfzc logic for your atoms).
2. QCEngine: add HF constants and *skeleton blocks* for every section you will test, for each
   `(system, reference, sdsc)` (ROHF: both `sd` and `sc`). HF tests need CONV-AE-CONV, CONV-AE-DF (pk scf with the
   dummy df corl type HF uses), DF-AE-DF, and CD-*-DF if CD is in scope. Correlated tests need FC twins.
3. psi4: add the `basis, subjects` row to each `*_module` test you cover, with marks
   `[*using("<addon>"), pytest.mark.<topic>]` (+ `pytest.mark.long` where FD runs are routine and expensive).
4. If the system needs special handling (skips, tolerances, FD options), key it off one flag in the runner
   (e.g., `is_ecp = basis in [...]`) and comment *why* with numbers.
5. Generate refs (7c), insert, run, triage (Sec. 8).

### 7c. Generating reference values
* Use the runner's "reference generation" convergence (`e_convergence 10`, `d_convergence 9`,
  `r_convergence 9`, ...); uncomment that block in the runner while working, and **re-comment before merge**.
* Script generation (one job per system/reference/fcae/algorithm), dump JSON with every relevant QCVariable,
  then insert with a script that locates blocks by full meta and asserts uniqueness. Hand-pasting is error-prone.
* Cross-check every number at least one independent way before storing it:
  * analytic gradient vs 5-point FD gradient (agreement 1e-7 or better; store the analytic one);
  * analytic Hessian vs FD-of-gradients vs FD-of-energies;
  * two modules for the same quantity (ccenergy vs fnocc; dfmp2 vs dfocc; occ vs ccenergy ROHF-MP2);
  * DF/CD vs CONV within `atol_conv`.
* When analytic and FD disagree, **don't** store either until you know which one is right (Sec. 8).
* FD caveats (see 8): default findif displaces in *mass-weighted* coordinates; for heavy atoms use
  `disp_mass_weight false` for gradients-from-energies; FD-of-energies Hessians of *non-variational*
  (correlated) energies on heavy atoms are noise-limited, so prefer FD-of-analytic-gradients.

## 8. Triage: failure message -> likely meaning

| symptom | first suspects |
|---|---|
| `KeyError` on `std_suite[...]` | missing block/generics for that hash (check `scf_type`/`corl_type`/`sdsc` the runner derived) |
| `KeyError` on a QCVariable in a block | ref not stored; or the contract expects it and shouldn't (fix contract) |
| `<var> wrongly present` | contract excuses it but psi4 now sets it (update contract) |
| `... v. CONV` miss only | DF/CD error vs `atol_conv`; quantify for this system before relaxing |
| `P::e` passes, `wfn` fails (or vice versa) | variable set on only one store; set on wfn in C++, copy to globals in proc.py |
| grd1 passes, grd0 fails | FD noise (check step/mass-weighting/convergence) or a real analytic bug (compare a fine scan) |
| hes2 differs from hes1 *and* hes0 (which agree) | analytic Hessian bug |
| hes0 differs from hes1 (analytic grads verified) | FD-of-energies noise (mass-weighted, non-variational energy) |
| one module disagrees with all others | that module's bug or a different (valid) approximation; test with a common setting |
| `Engine::lmax_exceeded` | Libint2 built with too-low derivative AM (`LIBINT_MAX_AM_LIST`); `d2ints`-type issue |
| error only with frozen core spread over several irreps | ordering bug (see fnocc MP2 below) |

Tools that settled real cases: run the same quantity in C1 vs natural symmetry; vary `disp_size`; scan the
energy along one coordinate and fit a polynomial (smoothness residual and slope); compare per-displacement
energies; force a common auxiliary basis across modules.

## 9. Docs / capabilities tables

Each non-`scftype` test appends a dict line to `stdsuite_psi4.txt` (in pytest's basetemp). The record has no
system/basis, so all basis rows (including ECP) merge into the same capability entries. See
`psi4/share/psi4/scripts/merge_stdsuite.py` for regenerating the tables.

## 10. ECP coverage (ch2i2 / ch2i, def2-svp) and lessons

<!-- ECP-STATUS -->
