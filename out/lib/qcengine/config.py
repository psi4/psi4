"""
Creates globals for the qcengine module
"""

import fnmatch
import getpass
import logging
import os
import socket
from typing import Any, Dict, Optional, Union

import pydantic

from .extras import get_information

__all__ = ["get_config", "get_provenance_augments", "global_repr", "NodeDescriptor"]

# Start a globals dictionary with small starting values
_global_values = None
NODE_DESCRIPTORS = {}
LOGGER = logging.getLogger("QCEngine")
LOGGER.setLevel(logging.CRITICAL)


# Generic globals
def get_global(key: Optional[str] = None) -> Union[str, Dict[str, Any]]:
    import cpuinfo
    import psutil

    # TODO (wardlt): Implement a means of getting CPU information from compute nodes on clusters for MPI tasks
    #  The QC code runs on a different node than the node running this Python function, which may have different info

    global _global_values
    if _global_values is None:
        _global_values = {}
        _global_values["hostname"] = socket.gethostname()
        _global_values["memory"] = round(psutil.virtual_memory().available / (1024 ** 3), 3)
        _global_values["username"] = getpass.getuser()

        # Work through VMs and logical cores.
        if hasattr(psutil.Process(), "cpu_affinity"):
            cpu_cnt = len(psutil.Process().cpu_affinity())
        else:
            cpu_cnt = psutil.cpu_count(logical=False)
            if cpu_cnt is None:
                cpu_cnt = psutil.cpu_count(logical=True)

        _global_values["ncores"] = cpu_cnt
        _global_values["nnodes"] = 1

        _global_values["cpuinfo"] = cpuinfo.get_cpu_info()
        try:
            _global_values["cpu_brand"] = _global_values["cpuinfo"]["brand_raw"]
        except KeyError:
            # Remove this if py-cpuinfo is pinned to >=6.0.0
            _global_values["cpu_brand"] = _global_values["cpuinfo"]["brand"]

    if key is None:
        return _global_values.copy()
    else:
        return _global_values[key]


class NodeDescriptor(pydantic.BaseModel):
    """
    Description of an individual node
    """

    # Host data
    hostname_pattern: str
    name: str
    scratch_directory: Optional[str] = None  # What location to use as scratch

    memory: Optional[float] = None
    memory_safety_factor: int = 10  # Percentage of memory as a safety factor

    # Specifications
    ncores: Optional[int] = pydantic.Field(
        None,
        description="""Number of cores accessible to each task on this node
    
    The default value, ``None``, will allow QCEngine to autodetect the number of cores.""",
    )
    jobs_per_node: int = 2
    retries: int = 0

    # Cluster options
    is_batch_node: bool = pydantic.Field(
        False,
        help="""Whether the node running QCEngine is a batch node
    
    Some clusters are configured such that tasks are launched from a special "batch" or "MOM" onto the compute nodes.
    The compute nodes on such clusters often have a different CPU architecture than the batch nodes and 
    often are unable to launch MPI tasks, which has two implications:
        1) QCEngine must make *all* calls to an executable via ``mpirun`` because the executables might not
        be able to run on the batch node. 
        2) QCEngine must run on the batch node to be able to launch tasks on the more than one compute nodes  
    
    ``is_batch_node`` is used when creating the task configuration as a means of determining whether
    ``mpiexec_command`` must always be used even for serial jobs (e.g., getting the version number)
    """,
    )
    mpiexec_command: Optional[str] = pydantic.Field(
        None,
        description="""Invocation for launching node-parallel tasks with MPI
        
        The invocation need not specify the number of nodes, tasks, or cores per node.
        Information about the task configuration will be added to the command by use of
        Python's string formatting. The configuration will be supplied as the following variables:
        
            {nnodes} - Number of nodes
            {ranks_per_node} - Number of MPI ranks per node
            {cores_per_rank} - Number of cores to use for each MPI rank
            {total_ranks} - Total number of MPI ranks
            
        As examples, the ``aprun`` command on Cray systems should be similar to 
        ``aprun -n {total_ranks} -N {ranks_per_node}`` and ``mpirun`` from OpenMPI should
        be similar to ``mpirun -np {total_ranks} -N {ranks_per_node}``.
        
        Programs where each MPI rank can use multiple threads (e.g., QC programs with MPI+OpenMP) can 
        use the {cores_per_rank} option to control the hybrid parallelism. 
        As an example, the Cray ``aprun`` command using this figure could be:
        ``aprun -n {total_ranks} -N {ranks_per_node} -d {cores_per_rank} -j 1``.
        The appropriate number of ranks per node will be determined based on the number of
        cores per node and the number of cores per rank.
        """,
    )

    def __init__(self, **data: Dict[str, Any]):
        data = parse_environment(data)
        super().__init__(**data)

        if self.mpiexec_command is not None:
            # Ensure that the mpiexec_command contains necessary information
            if not ("{total_ranks}" in self.mpiexec_command or "{nnodes}" in self.mpiexec_command):
                raise ValueError("mpiexec_command must contain either {total_ranks} or {nnodes}")
            if "{ranks_per_node}" not in self.mpiexec_command:
                raise ValueError("mpiexec_command must explicitly state the number of ranks per node")

    class Config:
        extra = "forbid"


class TaskConfig(pydantic.BaseModel):
    """Description of the configuration used to launch a task."""

    # Specifications
    ncores: int = pydantic.Field(None, description="Number cores per task on each node")
    nnodes: int = pydantic.Field(None, description="Number of nodes per task")
    memory: float = pydantic.Field(
        None, description="Amount of memory in GiB (2^30 bytes; not GB = 10^9 bytes) per node."
    )
    scratch_directory: Optional[str]  # What location to use as scratch
    retries: int  # Number of retries on random failures
    mpiexec_command: Optional[str]  # Command used to launch MPI tasks, see NodeDescriptor
    use_mpiexec: bool = False  # Whether it is necessary to use MPI to run an executable
    cores_per_rank: int = pydantic.Field(1, description="Number of cores per MPI rank")
    scratch_messy: bool = pydantic.Field(
        False, description="Leave scratch directory and contents on disk after completion."
    )

    class Config:
        extra = "forbid"


def _load_defaults() -> None:
    """
    Pulls the defaults from the QCA folder
    """

    # Find the config
    load_path = None
    test_paths = [os.getcwd(), os.path.join(os.path.expanduser("~"), ".qcarchive")]

    if "DQM_CONFIG_PATH" in os.environ:
        test_paths.insert(0, os.environ["DQM_CONFIG_PATH"])

    for path in test_paths:
        path = os.path.join(path, "qcengine.yaml")
        if os.path.exists(path):
            load_path = path
            break

    if load_path is None:
        LOGGER.info("Could not find 'qcengine.yaml'. Searched the following paths: {}".format(", ".join(test_paths)))
        LOGGER.info("Using default options...")

    else:
        import yaml

        LOGGER.info("Found 'qcengine.yaml' at path: {}".format(load_path))
        with open(load_path, "r") as stream:
            user_config = yaml.load(stream)

        for k, v in user_config.items():
            NODE_DESCRIPTORS[k] = NodeDescriptor(name=k, **v)


def global_repr() -> str:
    """
    A representation of the current global configuration.
    """

    ret = ""
    ret += "Host information:\n"
    ret += "-" * 80 + "\n"

    prov = get_provenance_augments()
    for k in ["username", "hostname", "cpu"]:
        ret += "{:<30} {:<30}\n".format(k, prov[k])

    ret += "\nNode information:\n"
    ret += "-" * 80 + "\n"
    for k, v in get_node_descriptor():
        ret += "  {:<28} {}\n".format(k, v)

        if k in ["scratch_directory", "memory_per_job"]:
            ret += "\n"

    ret += "\nJob information:\n"
    ret += "-" * 80 + "\n"
    for k, v in get_config():
        ret += "  {:<28} {}\n".format(k, v)

    ret += "-" * 80 + "\n"

    return ret


def get_node_descriptor(hostname: Optional[str] = None) -> NodeDescriptor:
    """
    Find the correct NodeDescriptor based off current hostname
    """
    if isinstance(hostname, NodeDescriptor):
        return hostname

    if hostname is None:
        hostname = get_global("hostname")

    # Find a match
    for name, node in NODE_DESCRIPTORS.items():

        if fnmatch.fnmatch(hostname, node.hostname_pattern):
            config = node
            break
    else:
        config = NodeDescriptor(
            name="default", hostname_pattern="*", memory=get_global("memory"), ncores=get_global("ncores")
        )

    return config


def parse_environment(data: Dict[str, Any]) -> Dict[str, Any]:
    """Collects local environment variable values into ``data`` for any keys with RHS starting with ``$``."""
    ret = {}
    for k, var in data.items():
        if isinstance(var, str) and var.startswith("$"):
            var = var.replace("$", "", 1)
            if var in os.environ:
                var = os.environ[var]
            else:
                var = None

        ret[k] = var

    return ret


def get_config(*, hostname: Optional[str] = None, local_options: Dict[str, Any] = None) -> TaskConfig:
    """
    Returns the configuration key for qcengine.
    """

    if local_options is None:
        local_options = {}

    local_options = parse_environment(local_options)
    config = {}

    # Node data
    node = get_node_descriptor(hostname)
    ncores = node.ncores or get_global("ncores")
    config["scratch_directory"] = local_options.pop("scratch_directory", node.scratch_directory)
    config["retries"] = local_options.pop("retries", node.retries)

    # Jobs per node
    jobs_per_node = local_options.pop("jobs_per_node", None) or node.jobs_per_node

    # Handle memory
    memory = local_options.pop("memory", None)
    if memory is None:
        memory = node.memory or get_global("memory")
        memory_coeff = 1 - node.memory_safety_factor / 100
        memory = round(memory * memory_coeff / jobs_per_node, 3)

    config["memory"] = memory

    # Get the number of cores available to each task
    ncores = local_options.pop("ncores", int(ncores / jobs_per_node))
    if ncores < 1:
        raise KeyError("Number of jobs per node exceeds the number of available cores.")

    config["ncores"] = ncores
    config["nnodes"] = local_options.pop("nnodes", 1)

    # Add in the MPI launch command template
    config["mpiexec_command"] = node.mpiexec_command
    config["use_mpiexec"] = node.is_batch_node or config["nnodes"] > 1
    config["cores_per_rank"] = local_options.get("cores_per_rank", 1)

    # Override any settings
    if local_options is not None:
        config.update(local_options)

    # Make sure mpirun command is defined if needed
    if config["use_mpiexec"] and config["mpiexec_command"] is None:
        raise ValueError(
            "You need to define the mpiexec command for this node. "
            "See: https://qcengine.readthedocs.io/en/stable/environment.html"
        )

    return TaskConfig(**config)


def get_provenance_augments() -> Dict[str, str]:
    return {
        "cpu": get_global("cpu_brand"),
        "hostname": get_global("hostname"),
        "username": get_global("username"),
        "qcengine_version": get_information("version"),
    }


def get_logger() -> "Logger":
    return LOGGER


# Pull in the local variables
_load_defaults()
