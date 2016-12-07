"""Dummy setup.py file solely for the purposes of getting an on-the-fly
computed version number into the conda recipe.

"""
from distutils.core import setup


def version_func():
    import subprocess

    command = 'python psi4/versioner.py --formatonly --format={versionlong}'
    process = subprocess.Popen(command.split(), shell=False, stdout=subprocess.PIPE)
    (out, err) = process.communicate()
    return out.strip()

setup(
    version=version_func(),
)
