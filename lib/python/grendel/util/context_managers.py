""" A collection of classes that conform to the Python context manager protocol
"""
import os
from grendel.util.misc import full_path

__all__ = [
    'working_directory',
]

class working_directory(object):
    """ Change to a certain working directory for the execution of a block of code, then change back upon (any kind of) exit.
    """

    def __init__(self, dir, permissions=0, create=False):
        self.directory = full_path(dir)
        if create:
            if not os.path.exists(self.directory):
                os.makedirs(self.directory)
        if not os.path.isdir(self.directory):
            raise OSError("{0} is not a directory.")
        elif not os.access(self.directory, os.R_OK|permissions):
            raise OSError("directory {0} cannot be accessed with sufficient permissions")

    def __enter__(self):
        self.return_directory = full_path(os.curdir)
        os.chdir(self.directory)

    def __exit__(self, exc_type, exc_val, exc_tb):
        os.chdir(self.return_directory)


class acquired_lock(object):
    """ Execute code with a lock acquired before the block of code and released afterwards
    """

    def __init__(self, lock):
        self.lock = lock

    def __enter__(self):
        self.lock.acquire()

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.lock.release()

