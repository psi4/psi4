r"""Module to provide mechanism to store and restore option states in driver.

"""
import sys
import PsiMod


class OptionState(object):
    """Class to store the state of a single *option*. If *module* given, the *option*
    value and has_changed value is stored for global, local to *module*, and used by
    *module* scopes; otherwise (used for BASIS keywords), only global scope is stored.
    Class can store, print, and restore option values. ::

        >>> OptionState('SCF_TYPE', 'SCF')

        >>> print(OptionState('DF_BASIS_MP2'))

    """
    def __init__(self, option, module=None):
        self.option = option.upper()
        if module:
            self.module = module.upper()
        else:
            self.module = None

        self.value_global = PsiMod.get_global_option(option)
        self.haschanged_global = PsiMod.has_global_option_changed(option)
        if self.module:
            self.value_local = PsiMod.get_local_option(self.module, option)
            self.haschanged_local = PsiMod.has_local_option_changed(self.module, option)
            self.value_used = PsiMod.get_option(self.module, option)
            self.haschanged_used = PsiMod.has_option_changed(self.module, option)
        else:
            self.value_local = None
            self.haschanged_local = None
            self.value_used = None
            self.haschanged_used = None

    def __str__(self):
        text = ''
        if self.module:
            text += """  ==> %s Option in Module %s <==\n\n""" % (self.option, self.module)
            text += """  Global (has changed?) value: %7s %s\n""" % ('(' + str(self.haschanged_global) + ')', self.value_global)
            text += """  Local (has changed?) value:  %7s %s\n""" % ('(' + str(self.haschanged_local) + ')', self.value_local)
            text += """  Used (has changed?) value:   %7s %s\n""" % ('(' + str(self.haschanged_used) + ')', self.value_used)
        else:
            text += """  ==> %s Option in Global Scope <==\n\n""" % (self.option)
            text += """  Global (has changed?) value: %7s %s\n""" % ('(' + str(self.haschanged_global) + ')', self.value_global)
        text += """\n"""
        return text

    def restore(self):
        PsiMod.set_global_option(self.option, self.value_global)
        if not self.haschanged_global:
            PsiMod.revoke_global_option_changed(self.option)
        if self.module:
            PsiMod.set_local_option(self.module, self.option, self.value_local)
            if not self.haschanged_local:
                PsiMod.revoke_local_option_changed(self.module, self.option)


class OptionsState(object):
    """Class to contain multiple :py:func:`~optproc.OptionsState` objects. ::

           >>> optstash = OptionsState(
               ['SCF', 'DFT_FUNCTIONAL'],
               ['DF_BASIS_SCF'],
               ['SCF', 'SCF_TYPE'],
               ['SCF', 'REFERENCE'])

           >>> print(optstash)

           >>> optstash.restore()

    """
    def __init__(self, *largs):
        self.data = []
        for item in largs:
            if len(item) == 2:
                self.data.append(OptionState(item[1], item[0]))
            elif len(item) == 1:
                self.data.append(OptionState(item[0]))
            else:
                print('ERROR: Each argument to OptionsState should be an array, the first element')
                print('       of which is the module scope and the second element of which is the')
                print('       module name. Bad argument: %s' % (item))
                sys.exit()

    def __str__(self):
        text = ''
        for item in self.data:
            text += str(item)
        return text

    def restore(self):
        for item in self.data:
            item.restore()
