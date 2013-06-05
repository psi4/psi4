from grendel.util import SystemInfo
from grendel.interface.runner import Runner
from subprocess import Popen

class LocalQChemRunner(Runner):
    """ Runs Q-Chem on the local machine.
    """

    ##################
    # Initialization #
    ##################

    def __init__(self, **kwargs):
        super(LocalQChemRunner, self).__init__(**kwargs)
        self.command_sequence = [SystemInfo.qchem_executable, self.input_file, self.output_file]


    ###########
    # Methods #
    ###########

    def validate(self):
        # TODO something more sophisticated here?
        return True

