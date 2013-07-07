from abc import abstractmethod, ABCMeta
from subprocess import Popen
import threading
from grendel.util.context_managers import acquired_lock

__all__ = [
    "Runner"
]

class Runner(object):
    """ Abstract base class for classes that run a computation (or an equivalent task, such as submit a job to a
    PBS queue).  Not all Computation objects need a runner; for instance XML-parsed computations are constructed in
    an already-run state.
    """
    __metaclass__ = ABCMeta

    ####################
    # Class Attributes #
    ####################

    default_input_filename = "input.dat"
    default_output_filename = "output.dat"

    ##############
    # Attributes #
    ##############

    # Note: Don't change any of these attributes without acquiring the _runner_lock if there is
    #   any chance the runner is running.  (Actually, don't change them unless you know what you
    #   are doing and have a really good reason to do so)
    process = None
    input_file = None
    output_file = None
    stdin = None
    stdout = None
    stderr = None
    environment = None

    working_directory = None
    """ The runner should do pass the working directory to the process to be opened rather than
    changing to that directory. This variable contains that working directory.
    """

    command_sequence = None
    """ The list of strings that gets passed to `subprocess.Popen()`.
    Basically, this is equivalent to the command-line command when each element of the list is separated by a space.
    """

    ######################
    # Private Attributes #
    ######################

    _runner_lock = None

    ##################
    # Initialization #
    ##################

    def __init__(self, input_file=None, output_file=None, **kwargs):
        self.input_file = input_file or self.__class__.default_input_filename
        self.output_file = output_file or self.__class__.default_output_filename
        self.process = None
        self._runner_lock = threading.Lock()
        for kw, val in kwargs.items():
            if hasattr(self.__class__, kw):
                setattr(self, kw, val)
            else:
                raise TypeError("Unknown attribute '{}' given to {} initializer".format(
                    kw,
                    self.__class__.__name__
                ))

    ####################
    # Abstract Methods #
    ####################

    @abstractmethod
    def validate(self):
        """ Determine if we will be able to run.  Raise a ComputationUnavailableError if running will cause a crash
        """
        return NotImplemented

    ###########
    # Methods #
    ###########

    def start(self):
        with acquired_lock(self._runner_lock):
            self.process = Popen(
                self.command_sequence,
                stdin=self.stdin,
                stdout=self.stdout,
                stderr=self.stderr,
                env=self.environment,
                cwd=self.working_directory
            )
        return self.process

    def run(self):
        """ Run in a blocking manner.  Default implementation calls wait() on the process returned by start().
        Overload it if this doesn't make sense.
        """
        self.start().wait()
        return self.process
