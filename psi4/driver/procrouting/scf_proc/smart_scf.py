from psi4.driver import p4util
from psi4.driver import constants
from psi4.driver.p4util.exceptions import ConvergenceError, ValidationError
from psi4 import core
from math import log
        
class smart_solver():
    """Purpose: enable easy extension of "smart" SCF solving capabilities.
    This class is an attribute of a wfn object ... a wfn object is also\
            an attribute of this class for easy data transfer.

    Methods:
        __init__: creates objects to hold energy and density history
        smart_iter: one stop method call placed in .iterations() function\
                makes decisions about convergence tools

        trailing_conv: is the SCF stuck in trailing convergence?\
                if so, turn on SOSCF and some damping.

        init_damp: is the guess case SAD, GWH, or CORE? if so, damp the \
                initial SCF iterations to remove wild energy changes.
                Currently set at 70%, a la ORCA, could fine tune a bit.

        smart_guess: makes a best stab at the guess density in the  case\
                that the user does not specify a guess method.

    Attributes:
        E_history: simple list of SCF energies in the current set of iterations.
        Drms_history: simple list of Drms values in the current set of iterations.
        opt_dict: Dictionary to hold convergence tools enabled at different\
                values of 'smart_level' a la OPTKING's dynamic_level
        smart_level: Level of effort to be expended to get SCF convergence\
                a la OPTKING dynamic_level

    Notes:

        The goal with having this class exist at all is to nicely contain\
                tools for scf convergence in one spot. Extensibility\
                could go something like this:

                    1. add a method to smart_solver class:
                        def great_convergence(self):
                            great convergence tricks

                    2. simply add a call in smart_iter to incorporate it into\
                            smart_scf.

                    3. You're done!
    """

    def __init__(self, wfn):
        self.wfn=wfn
        self._validate_smart()
        self.E_history=[]
        self.Drms_history=[]
        self.smart_level=self.wfn.smart_level #not necessary
        self.opt_dict=smart_opt_dict[self.smart_level] #not necessary
        self.user_vals={
                "frac_enabled":self.wfn.frac_enabled,
                "damping_enabled":self.wfn.damping_enabled,
                "soscf_enabled":self.wfn.soscf_enabled,
                "MAXITER":core.get_option('SCF','MAXITER')
                } 
        pass

    def smart_iter(self,SCFE,Drms):
        """Called every iteration to update energy history,
        density history, and call decision making methods.
        """
        self.update_E_history()
        self.update_Drms_history(Drms)
        if not self.initdamp():
            self.trailing_conv()
        
    def smart_guess(self):
        #!Q: Do we want basis_guess to be default smart behavior?
        do_castup = self.opt_dict["CASTUP"]
        if do_castup:
            castup_basis = self.opt_dict["CASTUP_BASIS"]
            core.set_option('SCF',"BASIS_GUESS",True) 

    def trailing_conv(self):
        """for auto detect trailing convergence and switching on SOSCF
        both 15 and 1E-8 were pulled out of a hat
        damping percentage 25 is included to catch oscillatory cases as well
        """
        if self.wfn.iteration_ > 15 and core.get_options('SCF','E_CONVERGENCE') < 1.0E-8:
            self.wfn.soscf_enabled = True
            self.wfn.damping_enabled = True
            self.wfn.damping_percentage = 25
            return True
        else:
            return False

    def initdamp(self):
        """decides whether to damp initial iterations in SCF and returns 
        True if initial damping occured, False otherwise
        """
        if self.wfn.iteration_ > 4:
            return False
        guess_opt = core.get_option('SCF',"GUESS")
        if (guess_opt in {'SAD','GWH','CORE'}) and self.wfn.iteration_ in range(0,4):
            self.wfn.damping_enabled=self.opt_dict['init_damp']
            self.wfn.damping_percentage=self.opt_dict['init_damp_percentage']
            return True

    def dyn_damp(self,Drms_target=5E-4):
        """
        Offers dynamic damping to hit a on a target Drms value.
        Perhaps 5E-4 is conservative - can be modified with smart_opt_dict
        Default is 5E-4, source: a hat
        """
        if self.wfn.iterations_ < 5:
            return 0 
        elif self.Drms_history[-1] > Drms_target:
            self.wfn.damping_percentage = Drms_target/self.Drms_history[-1]

    def update_E_history(self):
        self.E_history.append(self.wfn.get_energies("Total Energy"))

    def update_Drms_history(self,Drms):
        self.Drms_history.append(Drms)

    def _validate_smart(self):

        """Sanity-check smartSCF options

        Raises
        ------
        ValidationError
            If soscf, damping, DIIS, or other convergence settings don't match
            between user defined settings and smart recommendations, do:
                smart=0 -> no action, smart is desabledchange local.
                smart=1 -> resolve diff by falling to smartscf recommendations, print a notice in output.
                smart>=2 -> resolve diff by falling back to user defined settings
                            for non-smartscf defined values.
        Returns
        -------
        bool
            whether smart was able to resolve differences. Changes local scfiter options if discrepancy. 
        """

        pass

