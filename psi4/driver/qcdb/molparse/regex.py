import re

NUCLEUS = r"""(?:
   (?P<gh1>@)|(?P<gh2>Gh\())?                # optional ghost: @stuff or Gh(stuff) ...
        (                                    # mandatory element: AEuser or Zuser
            (?P<label1>
                (?P<A>\d+)?                  # optional mass number, A
                (?P<E>[A-Z]{1,3})            # mandatory atomic symbol, E
                (?P<user1>(_\w+)|(\d+))?) |  # optional user label: Enumber or E_alphanum
            (?P<label2>
                (?P<Z>\d{1,3})               # mandatory atomic number, Z
                (?P<user2>(_\w+))?)          # optional user label: Z_alphanum
        )
        (?:@(?P<mass>\d+\.\d+))?             # optional mass value [u]
   (?(gh2)\)                                 # ... ghost
          )"""

NUMBER = r"""(
    (?:[-+]?\d*\.\d+(?:[DdEe][-+]?\d+)?) |   # .num with optional sign, exponent, wholenum
    (?:[-+]?\d+\.\d*(?:[DdEe][-+]?\d+)?) |   # num. with optional sign, exponent, decimals
    (?:[-+]?\d+(?:[DdEe][-+]?\d+)?)          # num with optional sign, exponent
         )"""

SEP = r"""[\t ,]+"""
ENDL = r"""[\t ,]*$"""

CHGMULT = r"""(?P<chg>""" + NUMBER + r')' + SEP + r"""(?P<mult>\d+)"""
CARTXYZ = r'(?P<x>' + NUMBER + r')' + SEP + r'(?P<y>' + NUMBER + r')' + SEP + r'(?P<z>' + NUMBER + r')'
