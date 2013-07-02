"""
Static module containing various information about elements
"""

from grendel.chemistry.element import Element
from grendel.util.units import AtomicMassUnit

MASS_UNIT = AtomicMassUnit


Elements = dict()

H = Element("H", 1,
    [
        ( 1.00782503207, 0.999885 ),
        ( 2.0141017778, 0.000115 ),
        ( 3.0160492777, None ),
    ]
)
Elements['H'] = H

He = Element("He", 2,
    [
        (4.00260325415, 0.99999866),
        (3.0160293191, 1.34e-5),
    ]

)
Elements['He'] = He

B = Element("B", 5,
    [
        ( 10.0129370, 19.9 ),
        ( 11.0093054, 80.1 )
    ]
)
Elements['B'] = B

C = Element("C", 6,
    [
        ( 11.0114336, None ),
        ( 12.000000, 0.9893 ),
        ( 13.0033548378, 0.0107 ),
        ( 14.003241989, None ),
    ]
)
Elements['C'] = C

N = Element("N", 7,
    [
        ( 14.0030740048, 0.99636 ),
        ( 15.0001088982, 0.00364 )
    ]
)
Elements['N'] = N

O = Element("O", 8,
    [
        ( 15.99491461956, 0.99757 ),
        ( 16.99913170, 3.8e-4 ),
        ( 17.9991610, 2.05e-3 ),
    ]
)
Elements['O'] = O

F = Element("F", 9,
    [
        ( 18.99840322, 1.0 )
    ]
)
Elements['F'] = F

Ne = Element("Ne", 10,
    [
        ( 19.9924401754, 0.9048 ),
        ( 20.99384668, 0.0027 ),
        ( 21.991385114, 0.0925 )
    ]
)
Elements['Ne'] = Ne

Ge = Element("Ge", 32,
    [
        ( 69.9242474, 0.2038 ),
        ( 71.9220758, 0.2731 ),
        ( 72.9234589, 0.0776 ),
        ( 73.9211778, 0.3672 ),
        ( 75.9214026, 0.0783 )
    ]

)
Elements['Ge'] = Ge

#TODO Write a script to parse the rest of this data into the format above from some reputable source online.  Perhaps this page: http://www.nndc.bnl.gov/nudat2/reCenter.jsp?
# Here's a useful command:
# curl -d 'nuc=&spnuc=zan&z=Ne&a=&n=&zmin=&zmax=&amin=&amax=&nmin=&nmax=&evenz=any&evena=any&evenn=any&dmed=disabled&dmn=ANY&eled=disabled&elmin=0&elmax=40&jled=disabled&jlv=&plv=ANY&tled=disabled&tlmin=0&utlow=FS&tlmax=1E10&utupp=GY&ord=zalt&out=wp&unc=nds&sub=Search' -e 'http://www.nndc.bnl.gov/nudat2/indx_sigma.jsp' -O http://www.nndc.bnl.gov/nudat2/sigma_searchi.jsp

