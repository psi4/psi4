import pytest


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


#def is_numpy_new_enough(version_feature_introduced):
#    if not _plugin_import('numpy'):
#        return False
#    import numpy
#    from pkg_resources import parse_version
#    return parse_version(numpy.version.version) >= parse_version(version_feature_introduced)
#
#
#using_scipy = pytest.mark.skipif(_plugin_import('scipy') is False,
#                                reason='Not detecting module scipy. Install package if necessary and add to envvar PYTHONPATH')

using_psi4 = pytest.mark.skipif(_plugin_import('psi4') is False,
                                 reason='Not detecting module psi4. Install package and add to envvar PYTHONPATH')

#using_psi4_libxc = pytest.mark.skipif(is_psi4_new_enough("1.2a1.dev100") is False,
#                                reason="Psi4 does not include DFT rewrite to use Libxc. Update to development head")
#
#using_psi4_efpmints = pytest.mark.skipif(is_psi4_new_enough("1.2a1.dev507") is False,
#                                reason="Psi4 does not include EFP integrals in mints. Update to development head")
#
#using_psi4_python_integral_deriv = pytest.mark.skipif(is_psi4_new_enough("1000") is False,
#                                reason="Psi4 does not include derivatives of integrals exported to python. Update to development head")

using_psi4_molrec = pytest.mark.skipif(is_psi4_new_enough("1.2a1.dev999") is False,
                                reason="Psi4 does not use the new Molecule parsing. Update to development head")

#using_numpy_113 = pytest.mark.skipif(is_numpy_new_enough("1.13.0") is False,
#                                reason='NumPy does not include 1.13 features. Update package and add to envvar PYTHONPATH')
#
#using_matplotlib = pytest.mark.skipif(_plugin_import('matplotlib') is False,
#                                reason='Note detecting module matplotlib. Install package if necessary and add to envvar PYTHONPATH')

using_pylibefp = pytest.mark.skipif(_plugin_import('pylibefp') is False,
                                reason='Not detecting module pylibefp. Install package if necessary and add to envvar PYTHONPATH')

