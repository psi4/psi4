import os
import uuid
import shutil
import socket
import pathlib
import subprocess


def dftd3_subprocess(dftd3rec):  # dftd3rec@i -> dftd3rec@io
    """Minimal localized DFTD3 call and harvest from and into `dftd3rec`.

    Required Input Fields
    ---------------------
    command : list
        Command and arguments to execute.
    dftd3par : str
        Parameter file contents defining job run.
    dftd3_geometry : str
        XYZ of real atoms file contents defining job geometry.

    Optional Input Fields
    ---------------------
    scratch_location : str, optional
        Override the default scratch location.
        Note that this is PARENT dir.
    scratch_messy : bool, optional
        If present and `True`, scratch is left behind.

    Output Fields
    -------------
    stdout : str
        Main output file that gets written to stdout.

    Optional Output Fields
    ----------------------
    output_dftd3_gradient : str, optional
        If `-grad` present in `command`, contents of file with real atom gradient.

    """
    try:
        dftd3rec['command']
        dftd3rec['dftd3par']
        dftd3rec['dftd3_geometry']
    except KeyError as err:
        raise KeyError('Required field ({}) missing from ({})'.format(str(err), list(dftd3rec.keys()))) from err

    current_directory = os.getcwd()

    # move ~/.dftd3par.<hostname> out of the way so it won't interfere
    # TODO remove str() below as soon as py35 dropped
    defaultfile = os.path.join(str(pathlib.Path.home()), '.dftd3par.' + socket.gethostname())
    defmoved = False
    if os.path.isfile(defaultfile):
        os.rename(defaultfile, defaultfile + '_hide')
        defmoved = True

    # find environment by merging PSIPATH and PATH environment variables
    # * filter out None values as subprocess will fault on them
    lenv = {
        'HOME': os.environ.get('HOME'),
        'PATH': os.pathsep.join([os.path.abspath(x) for x in os.environ.get('PSIPATH', '').split(os.pathsep) if x != '']) + \
                os.pathsep + os.environ.get('PATH'),
        'LD_LIBRARY_PATH': os.environ.get('LD_LIBRARY_PATH')
        }
    lenv = {k: v for k, v in lenv.items() if v is not None}

    # set up unique scratch directory and move in
    if 'scratch_location' in dftd3rec:
        basedir = dftd3rec['scratch_location']
    else:
        basedir = os.environ['HOME'] + os.sep
    dftd3_tmpdir = basedir + 'dftd3_' + str(uuid.uuid4())[:8]
    if not os.path.exists(dftd3_tmpdir):
        os.mkdir(dftd3_tmpdir)
    os.chdir(dftd3_tmpdir)

    # write governing inputs
    # * old patched DFTD3 took parameters from file 'dftd3_parameters'
    with open('.dftd3par.local', 'w') as handle:
        handle.write(dftd3rec['dftd3par'])
    with open('dftd3_geometry.xyz', 'w') as handle:
        handle.write(dftd3rec['dftd3_geometry'])

    # call `dftd3` program
    try:
        spcall = subprocess.Popen(dftd3rec['command'], stdout=subprocess.PIPE, env=lenv)
    except OSError as err:
        raise OSError('Command (`{}`) failed with PATH ({})'.format(' '.join(dftd3rec['command']), lenv['PATH'])) from err

    # recover output data
    out, err = spcall.communicate()
    dftd3rec['stdout'] = out.decode('utf-8')

    try:
        with open('dftd3_gradient', 'r') as handle:
            dftd3rec['output_dftd3_gradient'] = handle.read()
    except (OSError, FileNotFoundError) as err:
        pass

    # clean up files and remove scratch directory
    if 'scratch_messy' not in dftd3rec or dftd3rec['scratch_messy'] is False:
        os.chdir('..')
        try:
            shutil.rmtree(dftd3_tmpdir)
        except OSError as err:
            raise OSError('Unable to remove dftd3 temporary directory: {}'.format(dftd3_tmpdir)) from err

    if defmoved is True:
        os.rename(defaultfile + '_hide', defaultfile)

    os.chdir(current_directory)

    return dftd3rec
