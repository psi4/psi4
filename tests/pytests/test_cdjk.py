"""
Tests for CDJK's in-core-only handling of SCF_SUBTYPE.

The Cholesky JK algorithm (SCF_TYPE=CD) has no out-of-core sub-algorithm: CDJK
implements neither preiterations() nor iaia(), so both of those run their
DiskDFJK implementations, and both consult is_core() to pick a sub-algorithm.
CDJK::is_core() therefore has to *override* DiskDFJK::is_core(), not merely hide
it, or the DF-flavored decision (a pessimistic memory estimate, or an explicit
SCF_SUBTYPE=OUT_OF_CORE) can route a CDJK object onto a disk path that does not
exist.

Each test below distinguishes an overriding is_core() from a hiding one.
"""

import pytest

import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api]


def _water():
    return psi4.geometry("""
    0 1
    O
    H 1 0.96
    H 1 0.96 2 104.5
    """)


def _build_cdjk(primary, memory):
    """Build an uninitialized CDJK holding `memory` doubles.

    JK.build()'s memory argument only steers the DF -> MEM_DF/DISK_DF choice and
    never reaches the object, so the budget has to be set separately -- same as
    scf_iterator.initialize_jk() does.
    """

    jk = psi4.core.JK.build(primary, jk_type="CD")
    assert jk.name() == "CDJK"
    jk.set_memory(memory)
    return jk


@pytest.mark.parametrize("scf_subtype", ["OUT_OF_CORE", "YOSHIMINE_OUT_OF_CORE", "REORDER_OUT_OF_CORE"])
def test_cdjk_rejects_out_of_core_subtype(scf_subtype):
    """CD is in-core only, so every disk-flavored SCF_SUBTYPE must be rejected by
    CDJK itself, with a message that names CD and points at the alternative.

    If is_core() only hides its parent, DiskDFJK::is_core() runs instead and the
    user gets one of two misleading errors: for OUT_OF_CORE it happily returns
    false and the failure surfaces much later out of CDJK::initialize_JK_disk(),
    and for the two PK-only spellings it blames DiskDFJK for an SCF_TYPE=CD run.
    """

    _water()
    psi4.set_options({
        "basis": "cc-pvdz",
        "scf_type": "cd",
        "scf_subtype": scf_subtype,
    })

    with pytest.raises(RuntimeError) as e:
        psi4.energy("scf")

    assert "Invalid SCF_SUBTYPE option in CDJK" in str(e.value), str(e.value)


@pytest.mark.parametrize("scf_subtype", ["AUTO", "INCORE"])
def test_cdjk_low_memory_reports_a_cd_memory_error(scf_subtype):
    """Starved of memory, CD must fail on a CD memory check, not get diverted by
    the DF-style memory gate in DiskDFJK::is_core().

    CDJK::memory_estimate() is deliberately crude -- it assumes 4*nbf Cholesky
    vectors, because the true count is not known until the decomposition has
    run. DiskDFJK::is_core() compares that guess against the memory budget, so
    with a hiding is_core() this raises 'Disk algorithm for CD JK not
    implemented.' (AUTO) or 'SCF_SUBTYPE=INCORE was specified, but there is not
    enough memory to do in-core!' (INCORE). With an overriding is_core() the
    guess is bypassed and the decomposition is attempted for real.
    """

    mol = _water()
    psi4.set_options({
        "basis": "cc-pvdz",
        "scf_subtype": scf_subtype,
    })
    primary = psi4.core.BasisSet.build(mol, "ORBITAL", "cc-pVDZ")

    # A quarter of CDJK::memory_estimate(), which is 4 * nbf**3 doubles.
    memory = primary.nbf() ** 3
    jk = _build_cdjk(primary, memory)

    # Precondition: this is the regime where the pessimistic estimate used to
    # send CD to disk. If it ever stops holding, the test proves nothing.
    assert jk.memory_estimate() > memory

    with pytest.raises(RuntimeError) as e:
        jk.initialize()

    # Which of the two genuine CD memory guards trips first depends on the
    # number of Cholesky vectors and on how many function pairs survive
    # screening, neither of which is known ahead of time. Both are CD-side
    # errors; neither is reachable while is_core() merely hides its parent.
    msg = str(e.value)
    assert ("Not enough memory for CD." in msg) or ("Cholesky: Memory constraints exceeded" in msg), msg


def test_cdjk_low_memory_loads_cached_integrals():
    """Below the pessimistic memory estimate, CD must still run when it can.

    Reading cached Cholesky vectors back off disk (DF_INTS_IO=LOAD) needs no
    decomposition and no working memory beyond the vectors themselves, so it
    succeeds on a budget far below CDJK::memory_estimate(). That only gets a
    chance to happen if CDJK::is_core() is the one deciding: with a hiding
    is_core(), DiskDFJK::is_core() sees the pessimistic estimate exceed the
    budget and diverts to CDJK::initialize_JK_disk(), which throws.
    """

    mol = _water()
    psi4.set_options({
        "basis": "cc-pvdz",
        "scf_subtype": "AUTO",
    })
    primary = psi4.core.BasisSet.build(mol, "ORBITAL", "cc-pVDZ")

    # Comfortably above the estimate, so this half runs either way; it exists
    # only to write the (Q|mn) integrals out.
    psi4.set_options({"df_ints_io": "SAVE"})
    jk_save = _build_cdjk(primary, int(1e8))
    jk_save.initialize()
    naux = psi4.variable("NAUX (SCF)")
    assert naux > 0

    # Same starved budget as the test above: a quarter of the estimate.
    memory = primary.nbf() ** 3

    psi4.set_options({"df_ints_io": "LOAD"})
    jk_load = _build_cdjk(primary, memory)
    assert jk_load.memory_estimate() > memory

    jk_load.initialize()

    assert psi4.variable("NAUX (SCF)") == naux
