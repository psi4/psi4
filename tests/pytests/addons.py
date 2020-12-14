import pytest

from qcelemental.util import parse_version, which, which_import
from qcengine.testing import using

__all__ = [
    'hardware_nvidia_gpu',
    'using_cppe',
    'using_dftd3',
    'using_dftd3_321',
    'using_gcp',
    'using_mdi',
    'using_mp2d',
    'using_memory_profiler',
    'using_networkx',
    'using_psi4',
    'using_qcdb',
]


def is_psi4_new_enough(version_feature_introduced):
    if not which_import('psi4', return_bool=True):
        return False
    import psi4
    return parse_version(psi4.__version__) >= parse_version(version_feature_introduced)


def is_nvidia_gpu_present():
    try:
        import GPUtil
    except ModuleNotFoundError:
        try:
            import gpu_dfcc
        except ModuleNotFoundError:
            # who knows?
            return False
        else:
            return gpu_dfcc.cudaGetDeviceCount() > 0
    else:
        try:
            ngpu = len(GPUtil.getGPUs())
        except FileNotFoundError:
            # no `nvidia-smi`
            return False
        else:
            return ngpu > 0


hardware_nvidia_gpu = pytest.mark.skipif(True, #is_nvidia_gpu_present() is False,
                                         reason='Psi4 not detecting Nvidia GPU via `nvidia-smi`. Install one')

using_memory_profiler = pytest.mark.skipif(
    which_import('memory_profiler', return_bool=True) is False,
    reason='Not detecting module memory_profiler. Install package if necessary and add to envvar PYTHONPATH')

using_psi4 = pytest.mark.skipif(
    False, reason='Not detecting module psi4. Install package if necessary and add to envvar PYTHONPATH')

using_qcdb = pytest.mark.skipif(
    True, reason='Not detecting common driver. Install package if necessary and add to envvar PYTHONPATH')

using_gcp = pytest.mark.skipif(
    which("gcp", return_bool=True) is False,
    reason="Not detecting executable gcp. Install package if necessary and add to envvar PATH")

using_cppe = pytest.mark.skipif(
    which_import('cppe', return_bool=True) is False,
    reason="Not detecting module cppe. Rebuild with -DENABLE_cppe")

using_networkx = pytest.mark.skipif(
    which_import('networkx', return_bool=True) is False,
    reason='Not detecting module networkx. Install package if necessary and add to envvar PYTHONPATH')

using_mdi = pytest.mark.skipif(
    which_import('mdi', return_bool=True) is False,
    reason="Not detecting module mdi. Install package if necessary and add to envvar PYTHONPATH (or rebuild Psi with -DENABLE_mdi)")

using_dftd3 = using('dftd3')
using_dftd3_321 = using('dftd3')
using_mp2d = using('mp2d')

using_psixas = pytest.mark.skipif(which_import("psixas", return_bool=True) is False, reason="Not detecting plugin psixas")

