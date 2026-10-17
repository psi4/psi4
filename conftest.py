"""Repository-level pytest configuration.

Load command-line plugins here so their options are available before pytest
descends into ``tests/pytests``.  The latter conftest is also installed as the
root conftest of Psi4's packaged test suite.
"""

pytest_plugins = ["addon_report"]
