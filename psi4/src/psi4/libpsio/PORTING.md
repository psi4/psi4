# libpsio on libspill — design record

The I/O core of `libpsio` is a thin shim (`psio_libspill.cc`) over
[libspill](https://github.com/ReLibQC/libspill). `psio.h`, `psio.hpp` and
`config.h` are unchanged, so every consumer sees the identical API — the ~2000
lines of hand-rolled I/O behind those headers are what the shim replaces. Of the
~3500 call sites, 2468 go through the `PSIO` class and 1043 through the free
functions; **none of them changed.**

This file records why the port is safe, what it deliberately left alone, and the
decisions still open. libspill itself arrives as a build-from-source add-on
(`external/upstream/libspill`).

The shim is one translation unit, `psio_libspill.cc`, with no header of its own
and no configuration. It carried both for a while: libspill hosted the port
before Psi4 did, and a `PSIO_USE_PSI4_HEADERS` macro chose between Psi4's real
`psio.h` and a standalone copy of the types, so the port could be built without a
Psi4 tree. In Psi4 that second branch is unbuildable — the copy of the types
never came across — and it was a duplicate implementation of the free functions
sitting next to the real one, which is exactly the shape that produced issue #1.
It is gone, along with the three free functions (`psio_done`, `psio_tocdel`,
`psio_zero_disk`) that existed only to serve it: none is declared in `psio.h`,
and none has ever had a caller in Psi4.

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

`PSIO::get_unit_filename()` — the helper #3515 factored out of `open.cc` — moved
into `get_filename.cc` when `open.cc` went, because it is filename composition
rather than I/O. The shim, `rename_file.cc` and `change_namespace.cc` all build
their paths from it, so there is one composition rule and not three.

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

## The PSIOManager contract

The part of this port that is easiest to get wrong and hardest to notice.

`PSIOManager` is Psi4's ledger of scratch files. It decides where a unit lives
(`get_file_path`, which is how `PSI_SCRATCH` and per-file overrides reach the
I/O layer), it records every file that is opened, and `psiclean()` unlinks
**exactly** the paths in that ledger — nothing else. Retention
(`mark_file_for_retention`, `set_specific_retention`) is keyed on the same
strings.

So a shim that reads and writes perfectly can still be wrong in three ways at
once, and the first version of this one was:

- it opened `<TMPDIR>/psi.<unit>.libspill`, so `PSI_SCRATCH` was ignored;
- it never called `open_file`/`close_file`, so the ledger stayed empty and
  `psiclean()` cleaned nothing;
- retention, reading that same empty ledger, was inert.

The fix is `ls_opts.exact_name`, which libspill added for this: it suppresses
libspill's own `<dir>/<name>.libspill` layout and takes the path verbatim, so
the store sits exactly where `PSIOManager` says it does. `PSIO::open` then
registers the file and `PSIO::close` reports it back, as `open.cc`/`close.cc`
did. `PSIO::exists` asks `ls_store_exists` rather than `stat`ing a path the shim
reconstructed itself.

The cost of `exact_name` is the protection it removes: the `.libspill` suffix is
how `ls_open` tells one of its stores from an unrelated file. With an exact name
the namespace is the caller's to manage — which is fine here, because Psi4's
names already carry pid, namespace and unit. One consequence worth stating:
a scratch directory left over from a **pre-port** Psi4 now holds files with the
right names and the wrong format, and `PSIO_OPEN_OLD` on one of those will fail
rather than silently misread it.

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

It also checks the *PSIOManager contract* above, because none of that is
visible in a read/write round trip: after opening a unit it reads back the
ledger `PSIOManager` mirrors to `psi.<pid>.clean`, requires the unit to appear
there, requires the recorded path to be under `PSIOManager`'s directory, and
stats it. A shim that stores the data somewhere of its own choosing passes
every other check in the file and fails these.

## One loose end in psio.h

`psio_volseek` is still declared in `psio.h`, but its definition (`volseek.cc`)
went with the I/O core and it has no callers anywhere in the tree — it is a
declaration for a function that cannot be linked. Removing it is the obvious
tidy-up, but it is a change to the public header this port otherwise leaves
untouched, so it is left as a separate decision rather than folded in here.

## Stage 2, deliberately separate

`aio_handler.cc` (~630 lines) could be replaced by libspill's `ls_aread`/
`ls_awrite`, but a first port should change nothing that could alter results or
timings, so it is untouched. `AIOHandler` covers nine files; the other 443 that
touch psio have no async path today.

## Open decisions

(Scratch-path policy was one of these; it is settled above, under *The
PSIOManager contract*.)

- **Per-instance state.** The shim keeps each `PSIO` instance's units in a side
  table keyed on `this`, because `psio.hpp` is unchanged and there is nowhere to
  put a member. Correct (Psi4 creates several instances — dfmp2, mcscf, DiskJK,
  sortintegrals) but it costs a map lookup per call. Adding one opaque member to
  `psio.hpp` would be cheaper at runtime and change no call site.
- **Memory budget.** Not yet passed to libspill. Feeding it a fraction of Psi4's
  limit unlocks the in-memory tier — the single largest measured win.
- **Error mapping.** libspill's codes are richer than `PSIO_ERROR_*`; the shim
  narrows them and reports the detail through a log callback. Confirm this suits
  Psi4's error conventions.

## Verification that still matters

The signature/semantics level is covered (conformance harness, in-tree compile).
What turns "the API fits" into evidence is a full build, the `ctest` suite, and a
numerical comparison against the old layer — done properly:

**Do not use a bit-identical / "machine precision" criterion.** Psi4 is threaded,
and OpenMP reduction order alone moves the last digits between runs of the *same
unmodified binary*. Establish that envelope first: run the unmodified code on the
chosen input several times and record the spread. The defensible claim is then
that the libspill runs fall **inside the code's own run-to-run spread**, not that
they reproduce a single reference bit for bit — a bit-identical test would "fail"
against the baseline itself and tell you nothing. (Measured elsewhere for
comparison: a water/cc-pVDZ SCF spread of 2e-13 across runs of untouched code.)

Two traps that would otherwise produce a false green:

1. **Pick a deterministic method.** Anything stochastic (stochastic-PT2-style)
   is useless as a comparison case.
2. **Make sure the case actually spills.** This is the important one for a
   *libpsio* port: `DiskDFJK` and the out-of-core CC paths only engage above a
   memory threshold, so a small test case silently exercises the in-core path and
   verifies nothing about the I/O layer that was replaced. Choose an input that
   genuinely goes to disk, or force it by setting the memory limit low enough —
   and confirm from the output that the disk path really ran.
