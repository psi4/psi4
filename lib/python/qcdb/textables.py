import itertools
try:
    from collections import OrderedDict
except ImportError:
    from oldpymodules import OrderedDict
from qcdb.modelchems import Method, BasisSet, Error, methods, bases, errors

mc_archive = {'mtd': methods, 'bas': bases, 'err': errors}


# define helper functions for formatting table cells
def val(kw):
    return r"""%s""" % (kw['matelem'])


def graphics(kw):
    mc = '-'.join([kw[bit] for bit in ['mtd', 'opt', 'bas']])
    return r"""\includegraphics[width=6.67cm,height=3.5mm]{%s%s.pdf}""" % (kw['plotpath'], mc)


def lmtdbas(kw):
    return """%-20s""" % (methods[kw['mtd']].latex + '/' + bases[kw['bas']].latex)


def label(kw):
    return """ %-20s""" % (mc_archive[kw['target']][kw[kw['target']]].latex)


def table_generic(dbse, serrors,
    mtd, bas, columnplan, rowplan=['bas', 'mtd'],
    opt=['CP'], err=['mae'], sset=['default'],
    landscape=False, standalone=True, subjoin=True,
    plotpath='', theme='', filename=None):
    """
    Arrays *mtd* and *bas* contain the keys to the qcdb.Method and
    qcdb.BasisSet objects that span all those that the table may
    encompass. If method and basis are to be scanned over, the arrays
    should be in the desired order.

    """
    def table_header(kw, abbr, head1, head0, head2):  # TODO caption, label
        """Form table header"""
        text.append('')
        text.append(r"""\begingroup""")
        text.append(r"""\squeezetable""")
        text.append(r"""\begin{%s}[h!tp]""" % ('sidewaystable' if landscape else 'table'))
        text.append(r"""\renewcommand{\baselinestretch}{1}""")
    #    text += r"""\\caption{$errorhash{$error} of the interaction energy for databases and their subsets with basis set $basishash{$basis}.""" + '\n'
    #    text += r"""\\label{tbl:qc-merge4dbse-".$error."-".$basis."-r".$round."}}""" + '\n'
        text.append(r"""\caption{""")
        text.append(r"""\label{tbl:qcdb-%s-%s}}""" % (theme, '-'.join([kw[bit] for bit in tag])))
        text.append(r"""\begin{ruledtabular}""")
        text.append(r"""\begin{tabular}{%s}""" % (abbr))
        text.append(head1)
        text.append(head0)
        text.append(head2)
        text.append(hline)

    def table_footer():
        """Form table footer"""
        text.append(r"""\end{tabular}""")
        text.append(r"""\end{ruledtabular}""")
        text.append(r"""\footnotetext[1]{Errors with respect to Gold Standard (see Sec. II D for plot details). Guide lines are at 0, 0.3, and 1.0 kcal/mol overbound ($-$) and underbound ($+$).}""")
        text.append(r"""\end{%s}""" % ('sidewaystable' if landscape else 'table'))
        text.append(r"""\endgroup""")
        text.append(r"""\clearpage""")
        text.append('')

    def matelem(dict_row, dict_col):
        """Return merge of index dictionaries *dict_row* and *dict_col* (precedence) with error string from serrors appended at key 'matelem'."""
        kw = dict(dict_row, **dict_col)
        kw['matelem'] = serrors['-'.join([kw[bit] for bit in ['mtd', 'opt', 'bas']])][kw['sset']][kw['dbse']][kw['err']]
        return kw

    # form LaTeX reference tag
    tag = []
    for key in ['dbse', 'sset', 'mtd', 'opt', 'bas', 'err']:
        if len(locals()[key]) == 1 or (key == rowplan[0] and not subjoin):
            tag.append(key)
    tag = set(tag)
    for col in columnplan:
        tag -= set(col[4].keys())

    # form column headers
    start = 1
    stop = 1
    head0 = ''
    for index in range(2, len(columnplan)):
        if columnplan[index][1] == columnplan[index - 1][1]:
            stop = index
        else:
            head0 += r"""\cline{%d-%d}""" % (start + 1, stop + 1)
            start = index
            stop = index
        if index + 1 == len(columnplan):
            head0 += r"""\cline{%d-%d}""" % (start + 1, stop + 1)

    abbr = ''.join([col[0] for col in columnplan])
    h1 = [(k, len(list(g))) for k, g in itertools.groupby([col[1] for col in columnplan])]
    head1 = ' & '.join([r"""\multicolumn{%d}{c}{\textbf{%s}}""" % (repeat, label) for (label, repeat) in h1]) + r""" \\"""
    h2 = [(k, len(list(g))) for k, g in itertools.groupby([col[2] for col in columnplan])]
    head2 = ' & '.join([r"""\multicolumn{%d}{c}{\textbf{%s}}""" % (repeat, label) for (label, repeat) in h2]) + r""" \\"""

    # form table body
    text = []
    nH = len(rowplan)
    hline = r"""\hline"""
    kw = {'plotpath': plotpath, 'sset': sset[0], 'dbse': dbse[0], 'err': err[0],
          'mtd': mtd[0], 'opt': opt[0], 'bas': bas[0]}

    if standalone:
        text.append(r"""""")
        text.append(r"""\documentclass[aip,jcp,preprint,superscriptaddress,floatfix]{revtex4-1}""")
        text.append(r"""\usepackage{bm}""")
        text.append(r"""\usepackage{dcolumn}""")
        text.append(r"""\usepackage{rotating}""")
        text.append(r"""\begin{document}""")

    if nH == 1:
        subjoin = True

    if subjoin:
        table_header(kw, abbr, head1, head0, head2)
        if text[-1] != hline:
            text.append(hline)

    for hier0 in locals()[rowplan[0]]:
        kw[rowplan[0]] = hier0
        kw['target'] = rowplan[0]
        if nH > 1:

            if not subjoin:
                table_header(kw, abbr, head1, head0, head2)
            if text[-1] != hline:
                text.append(hline)
            text.append(r"""\textbf{%s} \\""" % (mc_archive[rowplan[0]][hier0].latex))

            for hier1 in locals()[rowplan[1]]:
                kw[rowplan[1]] = hier1
                kw['target'] = rowplan[1]
                if nH > 2:
                    text.append(r"""\enspace\textbf{%s} \\""" % (mc_archive[rowplan[1]][hier1].latex))

                    for hier2 in locals()[rowplan[2]]:
                        kw[rowplan[2]] = hier2
                        kw['target'] = rowplan[2]

                        text.append(r"""\enspace\enspace""" + ' & '.join([col[3](matelem(kw, col[4])) for col in columnplan]) + r""" \\""")
                else:
                    text.append(r"""\enspace""" + ' & '.join([col[3](matelem(kw, col[4])) for col in columnplan]) + r""" \\""")
            if not subjoin:
                table_footer()
        else:
            text.append(' & '.join([col[3](matelem(kw, col[4])) for col in columnplan]) + r""" \\""")

    if subjoin:
        table_footer()

    if standalone:
        text.append(r"""\end{document}""")

    text = '\n'.join(text)
    print text


if __name__ == "__main__":

    serrors = {'MP2-CP-adtz': {'hb': {'S22': OrderedDict([('mae', '    0.38'), ('mape', '     2.5')]), 'HBC1': OrderedDict([('mae', '    0.28'), ('mape', '     2.5')]), 'NBC1': None, 'HSG': OrderedDict([('mae', '    0.28'), ('mape', '     1.9')]), 'DB4': OrderedDict([('mae', '    0.31'), ('mape', '     2.3')])}, 'default': {'S22': OrderedDict([('mae', '    0.89'), ('mape', '    19.2')]), 'HBC1': OrderedDict([('mae', '    0.28'), ('mape', '     2.5')]), 'NBC1': OrderedDict([('mae', '    1.10'), ('mape', '   139.0')]), 'HSG': OrderedDict([('mae', '    0.18'), ('mape', '     9.6')]), 'DB4': OrderedDict([('mae', '    0.61'), ('mape', '    42.6')])}, 'mxdd': {'S22': OrderedDict([('mae', '    1.13'), ('mape', '    27.1')]), 'HBC1': None, 'NBC1': OrderedDict([('mae', '    1.10'), ('mape', '   139.0')]), 'HSG': OrderedDict([('mae', '    0.16'), ('mape', '    10.9')]), 'DB4': OrderedDict([('mae', '    0.80'), ('mape', '    59.0')])}}, 'CCSD-CP-adz': {'hb': {'S22': OrderedDict([('mae', '    2.42'), ('mape', '    18.0')]), 'HBC1': OrderedDict([('mae', '    2.09'), ('mape', '    15.3')]), 'NBC1': None, 'HSG': OrderedDict([('mae', '    2.29'), ('mape', '    14.8')]), 'DB4': OrderedDict([('mae', '    2.27'), ('mape', '    16.1')])}, 'default': {'S22': OrderedDict([('mae', '    1.82'), ('mape', '    32.3')]), 'HBC1': OrderedDict([('mae', '    2.09'), ('mape', '    15.3')]), 'NBC1': OrderedDict([('mae', '    1.15'), ('mape', '   147.4')]), 'HSG': OrderedDict([('mae', '    1.14'), ('mape', '    85.9')]), 'DB4': OrderedDict([('mae', '    1.55'), ('mape', '    70.2')])}, 'mxdd': {'S22': OrderedDict([('mae', '    1.54'), ('mape', '    39.0')]), 'HBC1': None, 'NBC1': OrderedDict([('mae', '    1.15'), ('mape', '   147.4')]), 'HSG': OrderedDict([('mae', '    0.95'), ('mape', '    97.7')]), 'DB4': OrderedDict([('mae', '    1.21'), ('mape', '    94.7')])}}, 'MP2-CP-atz': {'hb': {'S22': OrderedDict([('mae', '    0.70'), ('mape', '     5.0')]), 'HBC1': OrderedDict([('mae', '    0.52'), ('mape', '     4.1')]), 'NBC1': None, 'HSG': OrderedDict([('mae', '    0.61'), ('mape', '     4.0')]), 'DB4': OrderedDict([('mae', '    0.61'), ('mape', '     4.4')])}, 'default': {'S22': OrderedDict([('mae', '    0.86'), ('mape', '    17.0')]), 'HBC1': OrderedDict([('mae', '    0.52'), ('mape', '     4.1')]), 'NBC1': OrderedDict([('mae', '    0.98'), ('mape', '   126.7')]), 'HSG': OrderedDict([('mae', '    0.24'), ('mape', '    13.5')]), 'DB4': OrderedDict([('mae', '    0.65'), ('mape', '    40.3')])}, 'mxdd': {'S22': OrderedDict([('mae', '    0.94'), ('mape', '    22.6')]), 'HBC1': None, 'NBC1': OrderedDict([('mae', '    0.98'), ('mape', '   126.7')]), 'HSG': OrderedDict([('mae', '    0.18'), ('mape', '    15.1')]), 'DB4': OrderedDict([('mae', '    0.70'), ('mape', '    54.8')])}}, 'CCSD-CP-atz': {'hb': {'S22': OrderedDict([('mae', '    1.41'), ('mape', '    10.5')]), 'HBC1': OrderedDict([('mae', '    1.11'), ('mape', '     8.4')]), 'NBC1': None, 'HSG': OrderedDict([('mae', '    1.39'), ('mape', '     9.0')]), 'DB4': OrderedDict([('mae', '    1.30'), ('mape', '     9.3')])}, 'default': {'S22': OrderedDict([('mae', '    1.27'), ('mape', '    23.5')]), 'HBC1': OrderedDict([('mae', '    1.11'), ('mape', '     8.4')]), 'NBC1': OrderedDict([('mae', '    0.97'), ('mape', '   123.2')]), 'HSG': OrderedDict([('mae', '    0.77'), ('mape', '    57.9')]), 'DB4': OrderedDict([('mae', '    1.03'), ('mape', '    53.3')])}, 'mxdd': {'S22': OrderedDict([('mae', '    1.20'), ('mape', '    29.6')]), 'HBC1': None, 'NBC1': OrderedDict([('mae', '    0.97'), ('mape', '   123.2')]), 'HSG': OrderedDict([('mae', '    0.67'), ('mape', '    66.1')]), 'DB4': OrderedDict([('mae', '    0.95'), ('mape', '    73.0')])}}, 'CCSD-CP-adtz': {'hb': {'S22': OrderedDict([('mae', '    1.03'), ('mape', '     7.6')]), 'HBC1': OrderedDict([('mae', '    0.76'), ('mape', '     5.9')]), 'NBC1': None, 'HSG': OrderedDict([('mae', '    1.02'), ('mape', '     6.6')]), 'DB4': OrderedDict([('mae', '    0.94'), ('mape', '     6.7')])}, 'default': {'S22': OrderedDict([('mae', '    1.05'), ('mape', '    20.0')]), 'HBC1': OrderedDict([('mae', '    0.76'), ('mape', '     5.9')]), 'NBC1': OrderedDict([('mae', '    0.90'), ('mape', '   113.9')]), 'HSG': OrderedDict([('mae', '    0.62'), ('mape', '    46.7')]), 'DB4': OrderedDict([('mae', '    0.83'), ('mape', '    46.6')])}, 'mxdd': {'S22': OrderedDict([('mae', '    1.06'), ('mape', '    25.7')]), 'HBC1': None, 'NBC1': OrderedDict([('mae', '    0.90'), ('mape', '   113.9')]), 'HSG': OrderedDict([('mae', '    0.55'), ('mape', '    53.3')]), 'DB4': OrderedDict([('mae', '    0.84'), ('mape', '    64.3')])}}, 'MP2-CP-adz': {'hb': {'S22': OrderedDict([('mae', '    1.64'), ('mape', '    12.2')]), 'HBC1': OrderedDict([('mae', '    1.41'), ('mape', '    10.4')]), 'NBC1': None, 'HSG': OrderedDict([('mae', '    1.42'), ('mape', '     9.1')]), 'DB4': OrderedDict([('mae', '    1.49'), ('mape', '    10.6')])}, 'default': {'S22': OrderedDict([('mae', '    0.97'), ('mape', '    16.2')]), 'HBC1': OrderedDict([('mae', '    1.41'), ('mape', '    10.4')]), 'NBC1': OrderedDict([('mae', '    0.70'), ('mape', '    95.9')]), 'HSG': OrderedDict([('mae', '    0.46'), ('mape', '    30.5')]), 'DB4': OrderedDict([('mae', '    0.88'), ('mape', '    38.3')])}, 'mxdd': {'S22': OrderedDict([('mae', '    0.65'), ('mape', '    18.0')]), 'HBC1': None, 'NBC1': OrderedDict([('mae', '    0.70'), ('mape', '    95.9')]), 'HSG': OrderedDict([('mae', '    0.30'), ('mape', '    34.1')]), 'DB4': OrderedDict([('mae', '    0.55'), ('mape', '    49.3')])}}}
    #print serrors['CCSD-CP-adtz']['hb']['HSG']['mape']

    columnplan = [
       ['l', r"""Method \& Basis Set""", '', label, {}],
       ['d', r'S22', 'HB', val, {'sset': 'hb', 'dbse': 'S22'}],
       ['d', r'S22', 'MX/DD', val, {'sset': 'mxdd', 'dbse': 'S22'}],
       ['d', r'S22', 'TT', val, {'sset': 'default', 'dbse': 'S22'}],
       ['d', r'Overall', 'HB', val, {'sset': 'hb', 'dbse': 'DB4'}],
       ['d', r'Overall', 'MX/DD', val, {'sset': 'mxdd', 'dbse': 'DB4'}],
       ['d', r'Overall', 'TT', val, {'sset': 'default', 'dbse': 'DB4'}]]

    table_generic(columnplan=columnplan, dbse=['DB4'], serrors=serrors, mtd=['MP2', 'CCSD'], bas=['adz', 'atz'], theme='test', subjoin=False)
