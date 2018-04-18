import numpy as np

from .util import filter_comments

def load_hessian(shess, dtype):

    # list o'lines w/o comments or blanks
    shess = filter_comments(shess)
    lhess = list(filter(None, map(str.strip, shess.splitlines())))

    if dtype in ['fcmfinal', 'cfour']:
        nat = int(lhess[0].split()[0])
        ndof = 3 * nat
        datastr = '\n'.join(lhess[1:])
        nhess = np.fromstring(datastr, sep=' ')
        nhess = nhess.reshape(ndof, ndof)
    else:
        raise ValidationError('Unknown dtype: {}'.format(dtype))

    return nhess


def to_string(hess, handle, dtype='psi4'):

    nat = hess.shape[0] // 3
    assert hess.shape == (3 * nat, 3 * nat)

    if dtype in ['fcmfinal', 'cfour', 'psi4', 'intder']:
        second_number = (6 * nat) if dtype == 'intder' else (3 * nat)
        header = '{:5}{:5}'.format(nat, second_number)

        np.savetxt(handle, hess.reshape((-1, 3)), fmt='%20.10f', delimiter='', newline='\n', header=header, comments='')

        # Bounty! a Psi4 mug or similar gear to anyone who trace the `6 * nat` above to a pre-PSI/CCQC source.
        #   See discussion starting https://github.com/psi4/psi4/pull/953#issuecomment-381447849
