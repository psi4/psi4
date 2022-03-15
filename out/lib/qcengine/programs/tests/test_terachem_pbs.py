import pytest

from qcengine.programs.terachem_pbs import TeraChemPBSHarness


def test_found_no_tcpb_installed(monkeypatch):
    # Empty sys.path so that no external packages can be found in case tcpb is installed
    with monkeypatch.context() as m:
        m.setattr("sys.path", [])
        harness = TeraChemPBSHarness()
        with pytest.raises(ModuleNotFoundError):
            harness.found(raise_error=True)
        assert harness.found() is False
