r"""Module to

"""
import sys
import PsiMod

class OptionState(object):
    """Class to

    """
    def __init__(self, option, module):
        self.option = option.upper()
        self.module = module.upper()

        self.value_global = PsiMod.get_global_option(option)
        self.value_local = PsiMod.get_local_option(self.module, option)
        self.value_used = PsiMod.get_option(self.module, option)
        self.haschanged_global = PsiMod.has_global_option_changed(option)
        self.haschanged_local = PsiMod.has_local_option_changed(self.module, option)
        self.haschanged_used = PsiMod.has_option_changed(self.module, option)

    def __str__(self):
        text  = ''
        text += """  ==> %s Option in Module %s <==\n\n""" % (self.option, self.module)
        text += """  Global (has changed?) value: %7s %s\n""" % ('('+str(self.haschanged_global)+')', self.value_global)
        text += """  Local (has changed?) value:  %7s %s\n""" % ('('+str(self.haschanged_local)+')', self.value_local)
        text += """  Used (has changed?) value:   %7s %s\n""" % ('('+str(self.haschanged_used)+')', self.value_used)
        text += """\n"""
        return text

    def restore(self):
        PsiMod.set_global_option(self.option, self.value_global)
        if not self.haschanged_global:
            PsiMod.revoke_global_option_changed(self.option)
        PsiMod.set_local_option(self.module, self.option, self.value_local)
        if not self.haschanged_local:
            PsiMod.revoke_local_option_changed(self.module, self.option)


class OptionsState(object):
    """Class to

    """
    def __init__(self, *largs):
        self.data = []
        for item in largs:
            if len(item) != 2:
                print('ERROR: Each argument to savestate should be an array, the first element of which is the module scope and the second element of which is the module name. Bad argument: %s' % (item))
                sys.exit()
            self.data.append(OptionState(item[1], item[0]))

    def __str__(self):
        text = ''
        for item in self.data:
            text += str(item)
        return text

    def restore(self):
        for item in self.data:
            item.restore()

