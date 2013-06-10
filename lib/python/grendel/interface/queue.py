from __future__ import print_function
from abc import ABCMeta, abstractmethod, abstractproperty

class Queue(object):
    """ Abstract base class for queues of computations to be run.
    """
    __metaclass__ = ABCMeta

    ####################
    # Class Attributes #
    ####################

    quiet = False
    verbose = False

    ##############
    # Attributes #
    ##############

    computations = abstractproperty()
    """ The Computation objects in the queue (regardless of whether they've been run)"""

    ####################
    # Abstract Methods #
    ####################

    @abstractmethod
    def run(self):
        """ Run the queue.
        (Abstract method)
        This runs the whole queue and blocks until the run is complete.  See also `start`
        """
        return NotImplemented

    @abstractmethod
    def start(self):
        """ Start the queue in the background and return to execution.
        Beware that if Python exits, the queue's processes will be orphaned.
        This mode of running may not be supported by all queue types.
        """

    @abstractmethod
    def is_started(self):
        """ True if the queue has been started (i. e. if any of the computations have been started)
        """

    @abstractmethod
    def is_finished(self):
        """ True if all the computations on the queue are finished
        """
        return NotImplemented

    @abstractmethod
    def enqueue(self, *c):
        """ Add computation(s) to the queue.
        """
        return NotImplemented




