import os
import sys
import time

import pytest
from qcelemental.models import AtomicInput

from qcengine import util
from qcengine.exceptions import InputError


def test_model_wrapper():

    with pytest.raises(InputError):
        util.model_wrapper({"bad": "yup"}, AtomicInput)


def test_compute_wrapper_capture():

    oldout = id(sys.stdout)
    with util.compute_wrapper(capture_output=True):

        assert id(sys.stdout) != id(oldout)

    assert id(sys.stdout) == oldout


def test_compute_wrapper_capture_exception():

    oldout = id(sys.stdout)
    with util.compute_wrapper(capture_output=True) as metadata:

        assert id(sys.stdout) != id(oldout)

        raise KeyError("Hello there!")

    assert id(sys.stdout) == oldout

    assert metadata["success"] is False
    assert metadata["error_type"] == "unknown_error"


def test_terminate():

    t = time.time()
    with util.popen(["sleep", "30"]) as proc:

        util.terminate_process(proc["proc"])

    assert (time.time() - t) < 1


def test_tmpdir():

    with util.temporary_directory(child="this") as tmpdir:
        assert str(tmpdir).split(os.path.sep)[-1] == "this"

    with util.temporary_directory() as parentdir:
        with util.temporary_directory(parent=parentdir, child="this") as tmpdir:
            assert str(tmpdir).split(os.path.sep)[-1] == "this"
            assert str(tmpdir).rsplit(os.path.sep, 1)[0] == str(parentdir)

    with util.temporary_directory(suffix="this") as tmpdir:
        assert str(tmpdir).split(os.path.sep)[-1].endswith("this")


def test_disk_files():

    infiles = {"thing1": "hello", "thing2": "world", "other": "everyone"}
    outfiles = {"thing*": None, "other": None}
    with util.temporary_directory(suffix="this") as tmpdir:
        with util.disk_files(infiles=infiles, outfiles=outfiles, cwd=tmpdir):
            pass

    assert outfiles.keys() == {"thing*", "other"}
    assert outfiles["thing*"]["thing1"] == "hello"
    assert outfiles["other"] == "everyone"


def test_popen_tee_output(capsys):
    # Test without passing
    with util.popen(["echo", "hello"]) as proc:
        proc["proc"].wait()
    assert proc["stdout"].strip() == "hello"

    # Test with passing
    with util.popen(["echo", "hello"], pass_output_forward=True) as proc:
        proc["proc"].wait()
    assert proc["stdout"] == "hello\n"
    captured = capsys.readouterr()
    assert captured.out == "hello\n"
