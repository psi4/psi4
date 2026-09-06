# libpsio on libspill — design record

The I/O core of `libpsio` is a thin shim (`psio_libspill.cc`) over
[libspill](https://github.com/ReLibQC/libspill). `psio.h`, `psio.hpp` and
`config.h` are unchanged, so every consumer sees the identical API — the ~2000
lines of hand-rolled I/O behind those headers are what the shim replaces. Of the
~3500 call sites, 2468 go through the `PSIO` class and 1043 through the free
functions; **none of them changed.**

This file records why the port is safe, what it deliberately left alone, and the
decisions still open. libspill itself arrives as a build-from-source add-on
(`external/upstream/libspill`); the shim is compiled with
`-DPSIO_USE_PSI4_HEADERS` so it binds against Psi4's own `psio.h`.

## What the shim provides

Both halves of libpsio's API, from one shared implementation so they cannot
drift apart:

- the **`PSIO` class** — all 21 declared methods plus `PSIO::PSIO`, `~PSIO`,
  `shared_object()`, the `_default_psio_lib_` singleton, `PSIO_ZERO` and
  `default_namespace_` (all of which lived in the deleted `init.cc`/`done.cc`);
- the **free functions**, which delegate to the default instance exactly as
  Psi4's `init.cc` did.

State is per instance, keyed on `this` (see *Open decisions*).

## What stays Psi4's (untouched)

`error.cc`, `decode_errno.cc`, `compose_err_msg.cc`, `getpid.cc`,
`filemanager.cc` (PSIOManager), `filescfg.cc`, `get_filename.cc`,
`rename_file.cc`, `change_namespace.cc`. These are Psi4 *policy* — scratch paths,
namespaces, retention, error text — not I/O, and libspill has no opinion about
them. They touch only `pid_` and `files_keywords_`, never `psio_unit`/`state_`,
so the shim owning the handles introduces no dual ownership.

`aio_handler.cc` (Psi4's async layer, 27 call sites in nine files) also stays —
see *Stage 2*.

## Why the swap is mechanical — three invariants

Each was checked against the tree (true at `a0e6ba5c4`; re-verify with the greps):

1. **`psio_address` is linear.** `psio_get_address(start, shift)` is exactly
   `start + shift` in (page, offset) coordinates, so `page * PSIO_PAGELEN +
   offset` is a faithful byte offset. `.page` appears **zero** times outside
   `libpsio`.
2. **`psio_tocscan` is only an existence test.** `->sadd`/`->eadd` are never
   dereferenced outside `libpsio`; every consumer writes `if (!psio_tocscan(...))`
   or compares to `nullptr`. The entry the shim hands back need only be valid and
   correctly filled, not a node in a live list.
3. **The on-disk table of contents is nobody's business.** `rd_toclen`/`tocread`/
   `toclen` have no consumers outside `libpsio`; `tocwrite` has one.

## The clearest win

`PSIO::zero_disk(unit, key, rows, cols)` used to loop `rows` times, writing a
zeroed row through the whole stack each time. The shim is one `ls_reserve` — on
Linux, one `fallocate`. ~59 call sites.

## The conformance harness (issue #1)

`test/psio_conformance_test.cc` takes the address of every free function **and**
every `PSIO` class method Psi4 declares, so a declared-but-undefined entry point
is a link error there rather than in the full core link, and then drives a real
round trip that also asserts two `PSIO` instances don't share unit state.

The history it guards against: the first shim implemented only the free
functions, and the first harness checked only the free functions — so it reported
success while the port could not link (the class is ~2400 of the call sites). A
trap if you extend it: comparing a member pointer to `nullptr` is **not** enough
— `&C::f` is never null, the optimiser folds the test, and no reference is
emitted. The addresses must be materialised (the harness memcpys them into a
volatile sink). **If you extend the shim, extend the harness in the same commit.**

## Stage 2, deliberately separate

`aio_handler.cc` (~630 lines) could be replaced by libspill's `ls_aread`/
`ls_awrite`, but a first port should change nothing that could alter results or
timings, so it is untouched. `AIOHandler` covers nine files; the other 443 that
touch psio have no async path today.

## Open decisions

- **Per-instance state.** The shim keeps each `PSIO` instance's units in a side
  table keyed on `this`, because `psio.hpp` is unchanged and there is nowhere to
  put a member. Correct (Psi4 creates several instances — dfmp2, mcscf, DiskJK,
  sortintegrals) but it costs a map lookup per call. Adding one opaque member to
  `psio.hpp` would be cheaper at runtime and change no call site.
- **Scratch-path policy.** libspill's `ls_opts.dir` takes a directory;
  `PSIOManager` already computes one per unit. The shim currently uses a single
  directory set via `psio_set_scratch_dir`. Wiring it to `PSIOManager`'s per-unit
  paths is the natural next step.
- **Memory budget.** Not yet passed to libspill. Feeding it a fraction of Psi4's
  limit unlocks the in-memory tier — the single largest measured win.
- **Error mapping.** libspill's codes are richer than `PSIO_ERROR_*`; the shim
  narrows them and reports the detail through a log callback. Confirm this suits
  Psi4's error conventions.

## Verification that still matters

The signature/semantics level is covered (conformance harness, in-tree compile).
What turns "the API fits" into evidence: a full build, the `ctest` suite, and a
byte-exact DF-CCSD(T) / disk-based MP2 spot-check against the old layer.
