"""
Same answer under memory pressure.

Each case runs one post-HF method twice with independent reference wavefunctions: first with ample memory, then
with memory squeezed so that the module takes its blocked, tiled, or out-of-core paths. The energies must agree,
and the squeezed output must show that the blocked path actually ran (else a change of basis, default, or memory
formula would silently turn the case into a duplicate of the ample run).

Memory is squeezed with psi4.core.set_memory_bytes, which, unlike psi4.set_memory, skips the 250 MiB floor, so small
molecules suffice. Only the post-HF step is squeezed; the SCF runs with ample memory.

To add coverage, add a row to `pressure_cases`: molecule key, method, options, squeezed memory in bytes, and output
regexes that prove the blocked path ran. Find a squeeze by scanning memory downward until the patterns appear and
the module has not yet refused to run.
"""

import math
import re

import psi4
import pytest
from utils import compare_values

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]

_molecules = {
    "ethene": """
C  0.000000  0.000000  0.667480
C  0.000000  0.000000 -0.667480
H  0.000000  0.922832  1.237695
H  0.000000 -0.922832  1.237695
H  0.000000  0.922832 -1.237695
H  0.000000 -0.922832 -1.237695
""",
}
_molecules["ethene_c1"] = "symmetry c1\n" + _molecules["ethene"]
_molecules["h2o"] = "O\nH 1 0.96\nH 1 0.96 2 104.5"
_molecules["h2o+"] = "1 2\n" + _molecules["h2o"]
_molecules["ch2"] = "0 3\nC\nH 1 1.1\nH 1 1.1 2 109.0"
_molecules["n2"] = "N\nN 1 1.1"

_common_options = {
    "basis": "cc-pvdz",
    "freeze_core": True,
    "e_convergence": 10,
    "d_convergence": 9,
    "r_convergence": 9,
}

# output patterns that prove a blocked path ran
fnocc_tiled = r"v\(ab,cd\) diagrams will be evaluated in\s+(?:[2-9]|[1-9]\d+) blocks"
fnocc_multifile_sort = r"Number of \(ab\|cd\) temporary files:\s+(?:[3-9]|[1-9]\d+)\n"
fnocc_t2_on_disk = r"T2 will be stored on disk"
fnocc_highmem_t = r"\(T\) algorithm:.*\(high-memory\)"
dpd_ooc_sort = r"DPD: not enough memory for an in-core buf4_sort"
dpd_cache_evict = r"DPD: memory is full; writing cached files to disk"
libtrans_tei_buckets = r"Sorting File: SO Ints \(nn\|nn\) nbuckets = (?:[2-9]|[1-9]\d+)"
libtrans_tpdm_buckets = r"Sorting File: MO TPDM \(AA\|AA\) nbuckets = (?:[2-9]|[1-9]\d+)"

fnocc_conv = {"qc_module": "fnocc", "scf_type": "pk", "cc_type": "conv", "mp_type": "conv"}
fnocc_df = {"qc_module": "fnocc", "scf_type": "disk_df", "df_ints_io": "save", "cc_type": "df"}  # fnocc DF-CC wants C1
cc_conv = {"qc_module": "ccenergy", "scf_type": "pk", "cc_type": "conv"}

pressure_cases = [
    # fnocc: CCSD tiles v(ab,cd) and v(ab,ci) and sorts integrals into many files; below ~1.6 MB T2 goes to disk
    pytest.param("ethene", "ccsd(t)", fnocc_conv, 3_000_000, [fnocc_tiled, fnocc_multifile_sort], id="fnocc-ccsd_t_-tiled", marks=pytest.mark.fnocc),
    pytest.param("ethene", "ccsd(t)", fnocc_conv, 1_500_000, [fnocc_tiled, fnocc_t2_on_disk], id="fnocc-ccsd_t_-t2disk", marks=pytest.mark.fnocc),
    pytest.param("ethene", "ccsd(t)", {**fnocc_conv, "triples_low_memory": False}, 3_000_000, [fnocc_tiled, fnocc_highmem_t], id="fnocc-ccsd_t_-highmemT-tiled", marks=pytest.mark.fnocc),
    pytest.param("ethene", "qcisd(t)", fnocc_conv, 1_500_000, [fnocc_t2_on_disk], id="fnocc-qcisd_t_-t2disk", marks=pytest.mark.fnocc),
    pytest.param("ethene", "mp4", fnocc_conv, 3_000_000, [fnocc_tiled], id="fnocc-mp4-tiled", marks=pytest.mark.fnocc),
    pytest.param("ethene", "mp4", fnocc_conv, 1_500_000, [fnocc_t2_on_disk], id="fnocc-mp4-t2disk", marks=pytest.mark.fnocc),
    pytest.param("ethene", "cepa(0)", fnocc_conv, 1_500_000, [fnocc_t2_on_disk], id="fnocc-cepa_0_-t2disk", marks=pytest.mark.fnocc),
    pytest.param("ethene", "fno-ccsd(t)", fnocc_conv, 1_500_000, [fnocc_t2_on_disk], id="fnocc-fno-ccsd_t_-t2disk", marks=pytest.mark.fnocc),
    # fnocc DF-CCSD: T2 goes to disk only within one T2's worth (~0.46 MB here) below the in-core requirement
    pytest.param("ethene_c1", "ccsd(t)", fnocc_df, 8_800_000, [fnocc_t2_on_disk], id="fnocc-dfccsd_t_-t2disk", marks=[pytest.mark.fnocc, pytest.mark.df]),
    # ccenergy/cceom: buf4 sorts run out of core, and libtrans transforms the integrals in buckets. At these squeezes a
    # transform bucket count once divided the rows evenly, dropping the last bucket (PSIO_ERROR 17, block start address).
    pytest.param("h2o", "ccsd", {**cc_conv, "basis": "cc-pvtz"}, 1_000_000, [dpd_ooc_sort, libtrans_tei_buckets], id="ccenergy-ccsd-ooc_sort-trans_buckets", marks=pytest.mark.cc),
    pytest.param("h2o+", "eom-ccsd", {**cc_conv, "basis": "cc-pvtz", "reference": "rohf", "roots_per_irrep": [1, 0, 0, 0]}, 3_500_000, [dpd_ooc_sort], id="cceom-rohf-eom_ccsd-ooc_sort-trans_buckets", marks=pytest.mark.cc),
    # cceom: out-of-core sorts, including qpsr
    pytest.param("ch2", "eom-ccsd", {**cc_conv, "basis": "cc-pvtz", "reference": "uhf", "roots_per_irrep": [1, 0, 0, 1]}, 3_000_000, [dpd_ooc_sort], id="cceom-uhf-eom_ccsd-ooc_sort", marks=pytest.mark.cc),
]  # yapf: disable

gradient_cases = [
    # ccenergy CCSD gradient: out-of-core sorts (incl. buf4_sort_axpy) and a TPDM presort in buckets, whose map once
    # missed the (q,p) half of packed pairs (segfault)
    pytest.param("n2", "ccsd", {**cc_conv, "basis": "cc-pvtz", "freeze_core": False}, 3_000_000, [dpd_ooc_sort, libtrans_tpdm_buckets], id="ccenergy-ccsd-ooc_sort-tpdm_buckets", marks=pytest.mark.cc),
    # cceom EOM-CCSD gradient: the DPD cache evicts to disk, where one label once held two differently shaped buffers
    # (PSIO_ERROR 18, write into next entry)
    pytest.param("h2o+", "eom-ccsd", {**cc_conv, "basis": "cc-pvdz", "reference": "rohf", "roots_per_irrep": [1, 0, 0, 0], "freeze_core": False}, 1_250_000, [dpd_cache_evict], id="cceom-rohf-eom_ccsd-cache_evict", marks=pytest.mark.cc),
]  # yapf: disable


@pytest.fixture
def restore_memory_and_threads():
    memory, nthreads = psi4.core.get_memory(), psi4.core.get_num_threads()
    yield
    psi4.core.set_memory_bytes(memory)
    psi4.core.set_num_threads(nthreads)


def _run(method, outfile, memory, ample, driver=psi4.energy):
    psi4.core.set_output_file(str(outfile.with_suffix(".scf.out")), False)
    psi4.core.set_memory_bytes(ample)
    _, ref_wfn = psi4.energy("scf", return_wfn=True)

    psi4.core.set_output_file(str(outfile), False)
    psi4.core.set_memory_bytes(memory)
    result = driver(method, ref_wfn=ref_wfn)
    psi4.core.set_output_file(str(outfile.with_suffix(".after")), False)  # flushes outfile; closing would leave none open
    return result, outfile.read_text()


@pytest.mark.parametrize("molecule,method,options,squeezed,triggers", pressure_cases)
def test_memory_pressure_and_blocking(molecule, method, options, squeezed, triggers, tmp_path, restore_memory_and_threads):
    psi4.set_num_threads(1)  # memory formulas (e.g., high-memory (T)) scale with threads
    psi4.geometry(_molecules[molecule])
    psi4.set_options({**_common_options, **options})

    ample = 2 * 1024**3
    e_ample, out_ample = _run(method, tmp_path / "ample.out", ample, ample)
    e_squeezed, out_squeezed = _run(method, tmp_path / "squeezed.out", squeezed, ample)

    missing = [pat for pat in triggers if not re.search(pat, out_squeezed)]
    assert not missing, f"squeezed run ({squeezed} bytes) did not take the blocked path: {missing}"
    assert any(not re.search(pat, out_ample) for pat in triggers), "ample run already took every blocked path"

    assert compare_values(e_ample, e_squeezed, 9, f"{method} energy: squeezed vs ample memory")


@pytest.mark.fnocc
def test_fnocc_lowmemory_triples_thread_fallback(tmp_path, restore_memory_and_threads):
    """Low-memory (T) that cannot afford its requested threads must fall back to fewer, not more.

    Ethene/STO-3G with frozen core has o = v = 6, so one-thread low-memory (T) needs 8*(2o^2v^2 + o^3v + ov + 5o^3)
    = 40032 bytes and each further thread 8*5o^3 = 8640 bytes. Squeezing to 7 bytes over the one-thread minimum
    leaves room for no extra thread (and puts T2 on disk). The old fallback, leftover / 5 * o^3, asked for 217 threads
    here; that is wrong but cheap, so a regression fails the assertion rather than exhausting the machine.
    """
    psi4.set_num_threads(4)
    psi4.geometry(_molecules["ethene"])
    psi4.set_options({**_common_options, **fnocc_conv, "basis": "sto-3g", "triples_low_memory": True})

    ample = 2 * 1024**3
    e_ample, out_ample = _run("ccsd(t)", tmp_path / "ample.out", ample, ample)
    e_squeezed, out_squeezed = _run("ccsd(t)", tmp_path / "squeezed.out", 40032 + 7, ample)

    assert "Not enough memory for requested threading" not in out_ample
    assert "Not enough memory for requested threading" in out_squeezed
    assert re.search(fnocc_t2_on_disk, out_squeezed)
    threads = re.search(r"Attempting to proceed with (\d+) threads", out_squeezed)
    assert threads and int(threads.group(1)) == 1, f"fallback chose {threads and threads.group(1)} threads, expected 1"

    assert compare_values(e_ample, e_squeezed, 9, "ccsd(t) energy: low-memory (T) thread fallback vs ample memory")


@pytest.mark.parametrize("molecule,method,options,squeezed,triggers", gradient_cases)
def test_memory_pressure_gradient(molecule, method, options, squeezed, triggers, tmp_path, restore_memory_and_threads):
    psi4.set_num_threads(1)
    psi4.geometry(_molecules[molecule])
    psi4.set_options({**_common_options, **options})

    ample = 2 * 1024**3
    g_ample, out_ample = _run(method, tmp_path / "ample.out", ample, ample, driver=psi4.gradient)
    g_squeezed, out_squeezed = _run(method, tmp_path / "squeezed.out", squeezed, ample, driver=psi4.gradient)

    missing = [pat for pat in triggers if not re.search(pat, out_squeezed)]
    assert not missing, f"squeezed run ({squeezed} bytes) did not take the blocked path: {missing}"
    assert any(not re.search(pat, out_ample) for pat in triggers), "ample run already took every blocked path"

    assert compare_values(g_ample, g_squeezed, 8, f"{method} gradient: squeezed vs ample memory")


@pytest.mark.cc
def test_dpd_32_gib_memory_accounting(tmp_path, restore_memory_and_threads):
    """A large configured budget must not overflow DPD's internal word counts on LLP64 platforms."""
    psi4.set_num_threads(2)
    psi4.geometry(
        """
        O
        H 1 1.0
        H 1 1.0 2 90.0
        """
    )
    psi4.set_options(
        {
            "basis": "sto-3g",
            "scf_type": "pk",
            "cc_type": "conv",
            "qc_module": "ccenergy",
            "freeze_core": True,
        }
    )
    outfile = tmp_path / "dpd-large-memory.out"
    psi4.core.set_output_file(str(outfile), False)
    psi4.core.set_memory_bytes(32 * 1024**3)

    energy = psi4.energy("ccsd(t)")
    psi4.core.set_output_file(str(tmp_path / "dpd-large-memory.after"), False)

    assert math.isfinite(energy)
    assert re.search(r"Number of threads for explicit ijk threading:\s+2", outfile.read_text())
