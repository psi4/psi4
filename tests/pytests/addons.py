import subprocess

import pytest
from qcelemental.util import parse_version, which, which_import


def is_psi4_new_enough(version_feature_introduced):
    if not which_import('psi4', return_bool=True):
        return False
    import psi4
    return parse_version(psi4.__version__) >= parse_version(version_feature_introduced)


def is_numpy_new_enough(version_feature_introduced):
    if not which_import('numpy', return_bool=True):
        return False
    import numpy
    return parse_version(numpy.version.version) >= parse_version(version_feature_introduced)


def is_dftd3_new_enough(version_feature_introduced):
    if not which('dftd3', return_bool=True):
        return False
    # Note: anything below v3.2.1 will return the help menu here. but that's fine as version compare evals to False.
    command = [which('dftd3'), '-version']
    proc = subprocess.run(command, stdout=subprocess.PIPE)
    candidate_version = proc.stdout.decode('utf-8').strip()

    return parse_version(candidate_version) >= parse_version(version_feature_introduced)


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


hardware_nvidia_gpu = pytest.mark.skipif(is_nvidia_gpu_present() is False,
                                         reason='Psi4 not detecting Nvidia GPU via `nvidia-smi`. Install one')

using_memory_profiler = pytest.mark.skipif(
    which_import('memory_profiler', return_bool=True) is False,
    reason='Not detecting module memory_profiler. Install package if necessary and add to envvar PYTHONPATH')

using_psi4 = pytest.mark.skipif(
    False, reason='Not detecting module psi4. Install package if necessary and add to envvar PYTHONPATH')

using_qcdb = pytest.mark.skipif(
    True, reason='Not detecting common driver. Install package if necessary and add to envvar PYTHONPATH')

using_dftd3 = pytest.mark.skipif(
    which('dftd3', return_bool=True) is False,
    reason='Not detecting executable dftd3. Install package if necessary and add to envvar PATH or PSIPATH')

using_dftd3_321 = pytest.mark.skipif(is_dftd3_new_enough("3.2.1") is False,
                                     reason='DFTD3 does not include 3.2.1 features. Update package and add to PATH')

using_gcp = pytest.mark.skipif(
    which("gcp", return_bool=True) is False,
    reason="Not detecting executable gcp. Install package if necessary and add to envvar PATH")

using_mp2d = pytest.mark.skipif(
    which('mp2d', return_bool=True) is False,
    reason='Not detecting executable mp2d. Install package if necessary and add to envvar PATH')

#using_psi4_libxc = pytest.mark.skipif(is_psi4_new_enough("1.2a1.dev100") is False,
#                                reason="Psi4 does not include DFT rewrite to use Libxc. Update to development head")

using_networkx = pytest.mark.skipif(
    which_import('networkx', return_bool=True) is False,
    reason='Not detecting module networkx. Install package if necessary and add to envvar PYTHONPATH')
