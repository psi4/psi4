from psi4.driver import p4util
from psi4.driver import constants
from psi4.driver.p4util.exceptions import ConvergenceError, ValidationError
from psi4 import core


class smart_solver():
    #TODO: Change to staging procedure
    """Purpose: enable easy extension of "smart" SCF solving capabilities.
    This class is an attribute of a wfn object ... a wfn object is also
            an attribute of this class for easy data transfer.

    Methods:
        __init__: creates objects to hold energy and density history
        smart_iter: one stop method call placed in .iterations() function
                makes decisions about convergence tools

        trailing_conv: is the SCF stuck in trailing convergence?
                if so, turn on SOSCF and some damping.

        init_damp: is the guess case SAD, GWH, or CORE? if so, damp the 
                initial SCF iterations to remove wild energy changes.
                Currently set at 70%, a la ORCA, could fine tune a bit.

    Attributes:
        E_history: simple list of SCF energies in the current set of iterations.
        Drms_history: simple list of Drms values in the current set of iterations.

    Notes:

        The goal with having this class exist at all is to nicely contain
                tools for scf convergence in one spot. Extensibility
                could go something like this:

                    1. add a method to smart_solver class:
                        def great_convergence(self):
                            great convergence tricks

                    2. simply add a call in smart_iter to incorporate it into
                            smart_scf.

                    3. You're done!
    """

    def __init__(self, wfn):
        self.wfn = wfn
        self.E_history = []
        self.Drms_history = []
        self.user_vals = {
            "frac_enabled": self.wfn.frac_enabled,
            "damping_enabled": self.wfn.damping_enabled,
            "soscf_enabled": self.wfn.soscf_enabled,
            "MAXITER": core.get_option('SCF', 'MAXITER')
        }

        #self.stage is where options are directly set by decision methods.
        self.stage = {}

    def smart_iter(self, SCFE, Drms):
        """Called every iteration to update energy history,
        density history, and call decision making methods.
        
        Decision making methods, such as initdamp(), edit a dictionary
        called self.smart_opts_dict. These are 'staged' options, which are
        then validated for conflicts at the end of this method. 
        """

        #Update Energy and Drms history
        self.update_E_history()
        self.update_Drms_history(Drms)
        if not self.initdamp():
            self.trailing_conv()

        #Finally, check for conflicts between set options
        _validate_smart()

    def trailing_conv(self):
        """Detects trailing-convergence and attempts to alleviate the issue
        by engaging SOSCF measures and a conservative amount of damping.

        The criteria for 'trailing convergence' are:
            Greater than 15 SCF cycles have occured,
            the most recent Drms value is less than 1E-4,
            the requested convergence is tighter (lower in magnitude) than
                1E-8
        
        These criteria are somewhat arbitrary but seem reasonable to MMD.

        """

        if (self.wfn.iteration_ > 15 and core.get_options('SCF', 'E_CONVERGENCE') < 1.0E-8
                and self.Drms_history[-1] <= 1E-4):
            self.wfn.soscf_enabled = True
            self.wfn.damping_enabled = True
            self.wfn.damping_percentage = 25
            return True

        else:
            return False

    def initdamp(self):
        """Decides whether to damp initial iterations in SCF.

        The criteria are:
            SCF iteration < 4, 
            and ither the SAD, GWH, or CORE guess was selected.
        """

        #If the iteration is not <= 4, no need to continue.
        if self.wfn.iteration_ > 4:
            return False

        guess_opt = core.get_option('SCF', "GUESS")

        if (guess_opt in {'SAD', 'GWH', 'CORE'}):
            self.stage['DAMPING_PERCENTAGE'] = smart_vals['INIT_DAMP_PERCENT']
            self.stage['DAMPING_ENABLED'] = True
            return True

    def dyn_damp(self, Drms_target=5E-4):
        """
        Offers dynamic damping to hit a target Drms value.
        Default is 5E-4
        """

        if self.wfn.iterations_ < 5:
            return 0
        elif self.Drms_history[-1] > Drms_target:
            self.stage['DAMPING_PERCENTAGE'] = \
                    Drms_target/self.Drms_history[-1]

    def update_E_history(self):
        self.E_history.append(self.wfn.get_energies("Total Energy"))

    def update_Drms_history(self, Drms):
        self.Drms_history.append(Drms)

    def _validate_smart_iter(self):
        """Sanity-check a smart_iter to ensure options make sense. 
        
        If there is a qualitative conflict, such as DIIS and SOSCF both
        requested by the smart_solver, then an error will be printed and the
        original user defined value will be taken. Conflicts like this should 
        not arise if the smart_solver is operating correctly, however it is 
        incldued for defensive coding.

        Raises
        ------

        Returns
        -------
        """
        #Currently no conflicts seem possible, so just return True and push
        #staged values to self.wfn
        for key in self.stage:
            if key.lower == 'damping_percentage':
                self.wfn.damping_percentage = self.stage[key]
            elif key.lower == 'damping_enabled':
                self.wfn.damping_enabled = self.stage[key]
            elif key.lower == 'soscf_enabled':
                self.wfn.soscf_enabled = self.stage[key]

        self.stage = {}
      

#seems silly to have a dict for just one option, but with more capabilities
#I think this makes sense        
gen_values = {'init_damp_percentage': 70}
