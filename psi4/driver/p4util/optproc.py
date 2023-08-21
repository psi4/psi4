#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2023 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#
r"""Module to provide mechanism to store and restore option states in driver.

"""

__all__ = [
    "OptionState",
    "OptionsStateCM",
    "OptionsState",
]

import sys
from contextlib import contextmanager
from typing import Iterator, List, Optional

from psi4 import core

from .exceptions import ValidationError


class OptionState():
    """Store the state (value and changed status) of a single `option`.

    Parameters
    ----------
    option
        Name of read_options keyword. All caps.
    module
        Name of read_options module or None if global. All caps.
        If `module` given, the `option` value and has_changed value is stored
        for global, local to `module`, and used by `module` scopes. Otherwise
        (used for BASIS keywords), only global scope is stored.

    Examples
    --------
    >>> OptionState('E_CONVERGENCE', 'SCF')

    >>> print(OptionState('DF_BASIS_MP2'))

    """

    def __init__(self, option: str, module: Optional[str] = None):
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
            text += """  Global (has changed?) value: %7s %s\n""" % ('(' + str(self.haschanged_global) + ')',
                                                                     self.value_global)
            text += """  Local (has changed?) value:  %7s %s\n""" % ('(' + str(self.haschanged_local) + ')',
                                                                     self.value_local)
            text += """  Used (has changed?) value:   %7s %s\n""" % ('(' + str(self.haschanged_used) + ')',
                                                                     self.value_used)
        else:
            text += """  ==> %s Option in Global Scope <==\n\n""" % (self.option)
            text += """  Global (has changed?) value: %7s %s\n""" % ('(' + str(self.haschanged_global) + ')',
                                                                     self.value_global)
        text += """\n"""
        return text

    def restore(self):
        """Restore value and has_changed status to saved condition."""
        core.set_global_option(self.option, self.value_global)
        if not self.haschanged_global:
            core.revoke_global_option_changed(self.option)
        if self.module:
            core.set_local_option(self.module, self.option, self.value_local)
            if not self.haschanged_local:
                core.revoke_local_option_changed(self.module, self.option)


class OptionsState():
    """Store multiple :py:func:`OptionState` objects.
    Use in driver functions to collect several keywords before altering them,
    then restore them before function return.

    Parameters
    ----------
    largs
        Specify which keywords to store value and has_changed state.

    Examples
    --------
    >>> optstash = OptionsState(
            ['DF_BASIS_SCF'],
            ['SCF_TYPE'],
            ['SCF', 'REFERENCE'])

    >>> print(optstash)

    >>> optstash.restore()

    """

    def __init__(self, *largs: List[List[str]]):
        self.data = {}
        for item in largs:
            self.add_option(item)

    def add_option(self, item: List[str]):
        """Store info for another keyword, `item`.

        Parameters
        ----------
        item
            A one-membered list with a global keyword or a two-membered list
            with a module keyword and module.

        """
        if len(item) == 2:
            key = (item[1], item[0])
        elif len(item) == 1:
            key = (item[0], )
        else:
            raise ValidationError(
                'Each argument to OptionsState should be an array, the first element of which is     the module scope and the second element of which is the module name. Bad argument: %s'
                % (item))

        if key in self.data:
            raise ValidationError(
                'Malformed options state, duplicate key adds of "{}". This should not happen, please raise a issue on github.com/psi4/psi4'.format(key))
        else:
            self.data[key] = OptionState(*key)

    def __str__(self):
        text = ''
        for key, item in self.data.items():
            text += str(item)
        return text

    def restore(self):
        """Restore value and has_changed status of each keyword to saved condition."""
        for key, item in self.data.items():
            item.restore()


@contextmanager
def OptionsStateCM(osd) -> Iterator[None]:
    """Return a context manager that will collect the state (value and changed
    status) of a list of keywords `osd` that can initialize
    :py:class:`OptionsState` on entry to the with-statement and restore the
    collected state when exiting the with-statement.

    """
    oso = OptionsState(osd)
    yield
    oso.restore()
