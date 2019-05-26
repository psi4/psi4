#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2019 The Psi4 Developers.
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

import collections

import numpy as np


class QCAspect(collections.namedtuple('QCAspect', 'lbl units data comment doi glossary')):
    """Facilitates the storage of quantum chemical results by labeling them with basic metadata."""

    def __new__(cls, lbl, units, data, comment='', doi=None, glossary=''):
        return super(QCAspect, cls).__new__(cls, lbl, units, data, comment, doi, glossary)

    def __str__(self, label=''):
        width = 40
        text = []
        text.append('-' * width)
        text.append('{:^{width}}'.format('QCAspect ' + self.lbl, width=width))
        if label:
            text.append('{:^{width}}'.format(label))
        text.append('-' * width)
        text.append('Data:     {}'.format(self.data))
        text.append('Units:    [{}]'.format(self.units))
        text.append('doi:      {}'.format(self.doi))
        text.append('Comment:  {}'.format(self.comment))
        text.append('Glossary: {}'.format(self.glossary))
        text.append('-' * width)
        return ('\n'.join(text))

    def to_dict(self):
        dicary = dict(self._asdict())  # dict, not OrderedDict
        for d in ['doi', 'comment', 'glossary']:
            dicary.pop(d)
        if isinstance(self.data, (np.ndarray, np.number)):
            if self.data.dtype == np.complex:
                dicary['data'] = [dicary['data'].real.tolist(), dicary['data'].imag.tolist()]
            else:
                dicary['data'] = dicary['data'].tolist()
        elif isinstance(self.data, (complex, np.complex)):
            dicary['data'] = [self.data.real, self.data.imag]

        return dicary


def print_variables(qcvars=None):
    """Form a printable representation of qcvariables.

    Parameters
    ----------
    qcvars : dict of QCAspect, optional
        Group of QCAspect objects to print. If `None`, will use `qcdb.pe.active_qcvars`.

    Returns
    -------
    str
        Printable string representation of label, data, and unit in QCAspect-s.

    """
    text = []
    text.append('\n  Variable Map:')
    text.append('  ----------------------------------------------------------------------------')

    if qcvars is None:
        qcvars = pe.active_qcvars

    if len(qcvars) == 0:
        text.append('  (none)')
        return '\n'.join(text)

    largest_key = max(len(k) for k in qcvars) + 2  # for quotation marks
    for k, qca in sorted(qcvars.items()):
        if k != qca.lbl:
            raise ValidationError('Huh? {} != {}'.format(k, qca.label))

        if isinstance(qca.data, np.ndarray):
            data = np.array_str(qca.data, max_line_width=120, precision=8, suppress_small=True)
            data = '\n'.join('        ' + ln for ln in data.splitlines())
            text.append("""  {:{keywidth}} => {:{width}} [{}]""".format(
                '"' + k + '"', '', qca.units, keywidth=largest_key, width=20))
            text.append(data)
        else:
            text.append("""  {:{keywidth}} => {:{width}.{prec}f} [{}]""".format(
                '"' + k + '"', qca.data, qca.units, keywidth=largest_key, width=20, prec=12))

    text.append('')
    return '\n'.join(text)
