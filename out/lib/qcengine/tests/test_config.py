"""
Tests the QCEngine compute module configuration
"""

import copy

import pydantic
import pytest

import qcengine as qcng
from qcengine.config import NodeDescriptor
from qcengine.util import create_mpi_invocation, environ_context


def test_node_blank():
    node = qcng.config.NodeDescriptor(name="something", hostname_pattern="*")


def test_node_auto():

    desc = {
        "name": "something",
        "hostname_pattern": "*",
        "jobs_per_node": 1,
        "ncores": 4,
        "memory": 10,
        "memory_safety_factor": 0,
    }
    node1 = qcng.config.NodeDescriptor(**desc)
    job1 = qcng.get_config(hostname=node1)
    assert job1.ncores == 4
    assert job1.nnodes == 1
    assert pytest.approx(job1.memory) == 10.0

    desc["jobs_per_node"] = 2
    node2 = qcng.config.NodeDescriptor(**desc)
    job2 = qcng.get_config(hostname=node2)
    assert job2.ncores == 2
    assert pytest.approx(job2.memory) == 5.0


def test_node_environ():

    scratch_name = "myscratch1234"
    with environ_context(env={"QCA_SCRATCH_DIR": scratch_name}):
        description = {"name": "something", "hostname_pattern": "*", "scratch_directory": "$QCA_SCRATCH_DIR"}

        node = qcng.config.NodeDescriptor(**description)
        assert node.scratch_directory == scratch_name


def test_node_skip_environ():
    description = {"name": "something", "hostname_pattern": "*", "scratch_directory": "$RANDOM_ENVIRON"}

    node = qcng.config.NodeDescriptor(**description)
    assert node.scratch_directory is None


@pytest.fixture
def opt_state_basic():
    """
    Capture the options state and temporarily override.
    """

    # Snapshot env
    old_node = copy.deepcopy(qcng.config.NODE_DESCRIPTORS)

    scratch_name = "myscratch1234"
    with environ_context(env={"QCA_SCRATCH_DIR": scratch_name}):

        configs = [
            {
                "name": "dragonstooth",
                "hostname_pattern": "dt*",
                "jobs_per_node": 2,
                "ncores": 12,
                "memory": 120,
                "scratch_directory": "$NOVAR_RANDOM_ABC123",
            },
            {"name": "newriver", "hostname_pattern": "nr*", "jobs_per_node": 2, "ncores": 24, "memory": 240},
            {
                "name": "batchnode",
                "hostname_pattern": "bn*",
                "is_batch_node": True,
                "jobs_per_node": 1,
                "mpiexec_command": "mpirun -n {total_ranks} -N {ranks_per_node}",
                "ncores": 24,
                "memory": 240,
            },
            {
                "name": "default",
                "hostname_pattern": "*",
                "jobs_per_node": 1,
                "memory": 4,
                "memory_safety_factor": 0,
                "ncores": 5,
                "scratch_directory": "$QCA_SCRATCH_DIR",
            },
        ]
        for desc in configs:
            node = qcng.config.NodeDescriptor(**desc)
            qcng.config.NODE_DESCRIPTORS[desc["name"]] = node

        yield

        # Reset env
        qcng.config.NODE_DESCRIPTORS = old_node


def test_node_matching(opt_state_basic):
    node = qcng.config.get_node_descriptor("nomatching")
    assert node.name == "default"

    node = qcng.config.get_node_descriptor("dt149")
    assert node.name == "dragonstooth"

    node = qcng.config.get_node_descriptor("nr149")
    assert node.name == "newriver"


def test_node_env(opt_state_basic):
    node = qcng.config.get_node_descriptor("dt")
    assert node.name == "dragonstooth"
    assert node.scratch_directory is None

    node = qcng.config.get_node_descriptor("nomatching")
    assert node.name == "default"
    assert node.scratch_directory == "myscratch1234"


def test_config_default(opt_state_basic):
    config = qcng.config.get_config(hostname="something")
    assert config.ncores == 5
    assert config.memory == 4

    config = qcng.config.get_config(hostname="dt149")
    assert config.ncores == 6
    assert config.retries == 0
    assert pytest.approx(config.memory, 0.1) == 54


def test_config_local_ncores(opt_state_basic):
    config = qcng.config.get_config(hostname="something", local_options={"ncores": 10, "retries": 3})
    assert config.ncores == 10
    assert config.memory == 4
    assert config.retries == 3


def test_config_local_njobs(opt_state_basic):
    config = qcng.config.get_config(hostname="something", local_options={"jobs_per_node": 5})
    assert config.ncores == 1
    assert pytest.approx(config.memory) == 0.8


def test_config_local_njob_ncore(opt_state_basic):
    config = qcng.config.get_config(hostname="something", local_options={"jobs_per_node": 3, "ncores": 1})
    assert config.ncores == 1
    assert pytest.approx(config.memory, 0.1) == 1.33


def test_config_local_njob_ncore_plus_memory(opt_state_basic):
    config = qcng.config.get_config(hostname="something", local_options={"jobs_per_node": 3, "ncores": 1, "memory": 6})
    assert config.ncores == 1
    assert pytest.approx(config.memory, 0.1) == 6


def test_config_local_nnodes(opt_state_basic):
    # Give a warning that mentions that mpirun is needed if you define a multi-node task
    with pytest.raises(ValueError) as exc:
        qcng.config.get_config(hostname="something", local_options={"nnodes": 10})
    assert "https://qcengine.readthedocs.io/en/stable/environment.html" in str(exc.value)

    # Test with an MPI run command
    local_options = {
        "nnodes": 4,
        "mpiexec_command": "mpirun -n {total_ranks} -N {ranks_per_node} --cpus-per-slot {cores_per_rank}",
    }
    config = qcng.config.get_config(hostname="something", local_options=local_options)
    assert config.use_mpiexec
    assert config.mpiexec_command.startswith("mpirun")
    assert create_mpi_invocation("hello_world", config) == [
        "mpirun",
        "-n",
        "20",
        "-N",
        "5",
        "--cpus-per-slot",
        "1",
        "hello_world",
    ]

    # Change the number of cores per rank
    local_options["cores_per_rank"] = 2
    config = qcng.config.get_config(hostname="something", local_options=local_options)
    assert config.use_mpiexec
    assert config.mpiexec_command.startswith("mpirun")
    assert create_mpi_invocation("hello_world", config) == [
        "mpirun",
        "-n",
        "10",
        "-N",
        "2",
        "--cpus-per-slot",
        "2",
        "hello_world",
    ]


def test_config_validation(opt_state_basic):
    with pytest.raises(pydantic.ValidationError):
        config = qcng.config.get_config(hostname="something", local_options={"bad": 10})


def test_global_repr():
    assert isinstance(qcng.config.global_repr(), str)


def test_batch_node(opt_state_basic):
    # Should always use mpirun
    config = qcng.config.get_config(hostname="bn1")
    assert config.use_mpiexec
    assert config.ncores == 24
    assert config.nnodes == 1

    config = qcng.config.get_config(hostname="bn1", local_options={"nnodes": 2})
    assert config.use_mpiexec
    assert config.ncores == 24
    assert config.nnodes == 2


def test_mpirun_command_errors():
    # Checks for ranks_per_node
    with pytest.raises(ValueError) as exc:
        NodeDescriptor(name="something", hostname_pattern="*", mpiexec_command="mpirun -n {total_ranks}")
    assert "must explicitly state" in str(exc.value)

    with pytest.raises(ValueError) as exc:
        NodeDescriptor(name="something", hostname_pattern="*", mpiexec_command="mpirun -N {ranks_per_node}")
    assert "must contain either" in str(exc.value)
