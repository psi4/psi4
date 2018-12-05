import os
import shutil

import pytest


def _which(command, return_bool=False):
    # environment is $PSIPATH:$PATH, less any None values
    lenv = {'PATH': ':'.join([os.path.abspath(x) for x in os.environ.get('PSIPATH', '').split(':') if x != '']) +
                    ':' + os.environ.get('PATH')}
    lenv = {k: v for k, v in lenv.items() if v is not None}

    ans = shutil.which(command, mode=os.F_OK | os.X_OK, path=lenv['PATH'])

    if return_bool:
        return bool(ans)
    else:
        return ans


def _plugin_import(plug):
    import sys
    if sys.version_info >= (3, 4):
        from importlib import util
        plug_spec = util.find_spec(plug)
    else:
        import pkgutil
        plug_spec = pkgutil.find_loader(plug)
    if plug_spec is None:
        return False
    else:
        return True


def is_psi4_new_enough(version_feature_introduced):
    if not _plugin_import('psi4'):
        return False
    import psi4
    from pkg_resources import parse_version
    return parse_version(psi4.__version__) >= parse_version(version_feature_introduced)


def is_numpy_new_enough(version_feature_introduced):
    if not _plugin_import('numpy'):
        return False
    import numpy
    from pkg_resources import parse_version
    return parse_version(numpy.version.version) >= parse_version(version_feature_introduced)


def is_nvidia_gpu_present():
    try:
        import GPUtil
    except ImportError:  # py36 ModuleNotFoundError
        try:
            import gpu_dfcc
        except ImportError:  # py36 ModuleNotFoundError
            # who knows?
            return False
        else:
            return gpu_dfcc.cudaGetDeviceCount() > 0
    else:
        try:
            ngpu = len(GPUtil.getGPUs())
        except OSError:  # py3 FileNotFoundError
            # no `nvidia-smi`
            return False
        else:
            return ngpu > 0


hardware_nvidia_gpu = pytest.mark.skipif(is_nvidia_gpu_present() is False,
                                reason='Psi4 not detecting Nvidia GPU via `nvidia-smi`. Install one')

using_memory_profiler = pytest.mark.skipif(_plugin_import('memory_profiler') is False,
                                reason='Not detecting module memory_profiler. Install package if necessary and add to envvar PYTHONPATH')

using_psi4 = pytest.mark.skipif(False,
                                reason='Not detecting module psi4. Install package if necessary and add to envvar PYTHONPATH')

using_qcdb = pytest.mark.skipif(True,
                                reason='Not detecting common driver. Install package if necessary and add to envvar PYTHONPATH')

using_dftd3 = pytest.mark.skipif(_which('dftd3', return_bool=True) is False,
                                reason='Not detecting executable dftd3. Install package if necessary and add to envvar PATH or PSIPATH')

using_gcp = pytest.mark.skipif(_which("gcp", return_bool=True) is False,
                                reason="Not detecting executable gcp. Install package if necessary and add to envvar PSIPATH or PATH")

#using_scipy = pytest.mark.skipif(_plugin_import('scipy') is False,
#                                reason='Not detecting module scipy. Install package if necessary and add to envvar PYTHONPATH')
#
#using_psi4_libxc = pytest.mark.skipif(is_psi4_new_enough("1.2a1.dev100") is False,
#                                reason="Psi4 does not include DFT rewrite to use Libxc. Update to development head")
#
#using_psi4_efpmints = pytest.mark.skipif(is_psi4_new_enough("1.2a1.dev507") is False,
#                                reason="Psi4 does not include EFP integrals in mints. Update to development head")
#
#using_numpy_113 = pytest.mark.skipif(is_numpy_new_enough("1.13.0") is False,
#                                reason='NumPy does not include 1.13 features. Update package and add to envvar PYTHONPATH')
#
#using_matplotlib = pytest.mark.skipif(_plugin_import('matplotlib') is False,
#                                reason='Note detecting module matplotlib. Install package if necessary and add to envvar PYTHONPATH')
#
#using_pylibefp = pytest.mark.skipif(_plugin_import('pylibefp') is False,
#                                reason='Not detecting module pylibefp. Install package if necessary and add to envvar PYTHONPATH')
