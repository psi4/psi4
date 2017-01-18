#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

r"""Module to provide mechanism to store and restore option states in driver.

"""
import sys
from .exceptions import *


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

        self.value_global = core.get_global_option(option)
        self.haschanged_global = core.has_global_option_changed(option)
        if self.module:
            self.value_local = core.get_local_option(self.module, option)
            self.haschanged_local = core.has_local_option_changed(self.module, option)
            self.value_used = core.get_option(self.module, option)
            self.haschanged_used = core.has_option_changed(self.module, option)
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
        core.set_global_option(self.option, self.value_global)
        if not self.haschanged_global:
            core.revoke_global_option_changed(self.option)
        if self.module:
            core.set_local_option(self.module, self.option, self.value_local)
            if not self.haschanged_local:
                core.revoke_local_option_changed(self.module, self.option)


class OptionsState(object):
    """Class to contain multiple :py:func:`~optproc.OptionsState` objects.
    Used in python driver functions to collect several options before altering
    them, then restoring before function return. ::

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
                raise ValidationError('Each argument to OptionsState should be an array, the first element of which is     the module scope and the second element of which is the module name. Bad argument: %s' % (item))

    def __str__(self):
        text = ''
        for item in self.data:
            text += str(item)
        return text

    def restore(self):
        for item in self.data:
            item.restore()
