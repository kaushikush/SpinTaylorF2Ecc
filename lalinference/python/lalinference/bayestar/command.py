#
# Copyright (C) 2013-2017  Leo Singer
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
"""
Functions that support the command line interface.
"""

from __future__ import print_function
import argparse
from distutils.dir_util import mkpath
from distutils.errors import DistutilsFileError
import glob
import inspect
import itertools
import logging
import os
import sys
import tempfile
import matplotlib
from matplotlib import cm
from ..plot import cmap
from ..util import sqlite


# Set no-op Matplotlib backend to defer importing anything that requires a GUI
# until we have determined that it is necessary based on the command line
# arguments.
if 'matplotlib.pyplot' in sys.modules:
    from matplotlib import pyplot as plt
    plt.switch_backend('Template')
else:
    matplotlib.use('Template', warn=False, force=True)


# FIXME: Remove this after all Matplotlib monkeypatches are obsolete.
import matplotlib
import distutils.version
mpl_version = distutils.version.LooseVersion(matplotlib.__version__)


def get_version():
    from .. import InferenceVCSInfo as vcs_info
    return vcs_info.name + ' ' + vcs_info.version


class GlobAction(argparse._StoreAction):
    """Generate a list of filenames from a list of filenames and globs."""

    def __call__(self, parser, namespace, values, *args, **kwargs):
        values = list(
            itertools.chain.from_iterable(glob.iglob(s) for s in values))
        if values:
            super(GlobAction, self).__call__(
                parser, namespace, values, *args, **kwargs)
        nvalues = getattr(namespace, self.dest)
        nvalues = 0 if nvalues is None else len(nvalues)
        if self.nargs == argparse.OPTIONAL:
            if nvalues > 1:
                msg = 'expected at most one file'
            else:
                msg = None
        elif self.nargs == argparse.ONE_OR_MORE:
            if nvalues < 1:
                msg = 'expected at least one file'
            else:
                msg = None
        elif self.nargs == argparse.ZERO_OR_MORE:
            msg = None
        elif int(self.nargs) != nvalues:
            msg = 'expected exactly %s file' % self.nargs
            if self.nargs != 1:
                msg += 's'
        else:
            msg = None
        if msg is not None:
            msg += ', but found '
            msg += '{} file'.format(nvalues)
            if nvalues != 1:
                msg += 's'
            raise argparse.ArgumentError(self, msg)


waveform_parser = argparse.ArgumentParser(add_help=False)
group = waveform_parser.add_argument_group(
    'waveform options', 'Options that affect template waveform generation')
# FIXME: The O1 uberbank high-mass template, SEOBNRv2_ROM_DoubleSpin, does
# not support frequencies less than 30 Hz.
group.add_argument(
    '--f-low', type=float, metavar='Hz', default=30,
    help='Low frequency cutoff [default: %(default)s]')
group.add_argument(
    '--f-high-truncate', type=float, default=0.95,
    help='Truncate waveform at this fraction of the maximum frequency of the '
    'PSD [default: %(default)s]')
group.add_argument(
    '--waveform', default='o2-uberbank',
    help='Template waveform approximant (e.g., TaylorF2threePointFivePN) '
    '[default: O2 uberbank mass-dependent waveform]')
del group


prior_parser = argparse.ArgumentParser(add_help=False)
group = prior_parser.add_argument_group(
    'prior options', 'Options that affect the BAYESTAR likelihood')
group.add_argument(
    '--min-distance', type=float, metavar='Mpc',
    help='Minimum distance of prior in megaparsecs '
    '[default: infer from effective distance]')
group.add_argument(
    '--max-distance', type=float, metavar='Mpc',
    help='Maximum distance of prior in megaparsecs '
    '[default: infer from effective distance]')
group.add_argument(
    '--prior-distance-power', type=int, metavar='-1|2', default=2,
    help='Distance prior '
    '[-1 for uniform in log, 2 for uniform in volume, default: %(default)s]')
group.add_argument(
    '--cosmology', default=False, action='store_true',
    help='Use cosmological comoving volume prior [default: %(default)s]')
group.add_argument(
    '--disable-snr-series', dest='enable_snr_series', action='store_false',
    help='Disable input of SNR time series (WARNING: UNREVIEWED!) '
    '[default: enabled]')
del group


skymap_parser = argparse.ArgumentParser(add_help=False)
group = skymap_parser.add_argument_group(
    'sky map output options', 'Options that affect sky map output')
group.add_argument(
    '--nside', '-n', type=int, default=-1,
    help='HEALPix resolution [default: auto]')
group.add_argument(
    '--chain-dump', default=False, action='store_true',
    help='For MCMC methods, dump the sample chain to disk [default: no]')
del group


class MatplotlibFigureType(argparse.FileType):

    def __init__(self):
        super(MatplotlibFigureType, self).__init__('wb')

    @staticmethod
    def __show():
        from matplotlib import pyplot as plt
        return plt.show()

    def __save(self):
        from matplotlib import pyplot as plt
        _, ext = os.path.splitext(self.string)
        ext = ext.lower()
        program, _ = os.path.splitext(os.path.basename(sys.argv[0]))
        cmdline = ' '.join([program] + sys.argv[1:])
        metadata = {'Title': cmdline}
        if ext == '.png':
            metadata['Software'] = get_version()
        elif ext in {'.pdf', '.ps', '.eps'}:
            metadata['Creator'] = get_version()
        return plt.savefig(self.string, metadata=metadata)

    def __call__(self, string):
        from matplotlib import pyplot as plt
        if string == '-':
            plt.switch_backend(matplotlib.rcParamsOrig['backend'])
            return self.__show
        else:
            with super(MatplotlibFigureType, self).__call__(string):
                pass
            plt.switch_backend('agg')
            self.string = string
            return self.__save


class HelpChoicesAction(argparse.Action):

    def __init__(self,
                 option_strings,
                 choices=(),
                 dest=argparse.SUPPRESS,
                 default=argparse.SUPPRESS):
        name = option_strings[0].replace('--help-', '')
        super(HelpChoicesAction, self).__init__(
            option_strings=option_strings,
            dest=dest,
            default=default,
            nargs=0,
            help='show support values for --' + name + ' and exit')
        self._name = name
        self._choices = choices

    def __call__(self, parser, namespace, values, option_string=None):
        print('Supported values for --' + self._name + ':')
        for choice in self._choices:
            print(choice)
        parser.exit()


def type_with_sideeffect(type):
    def decorator(sideeffect):
        def func(value):
            ret = type(value)
            sideeffect(ret)
            return ret
        return func
    return decorator


@type_with_sideeffect(str)
def colormap(value):
    from matplotlib import rcParams
    rcParams['image.cmap'] = value


@type_with_sideeffect(float)
def figwidth(value):
    from matplotlib import rcParams
    rcParams['figure.figsize'][0] = float(value)


@type_with_sideeffect(float)
def figheight(value):
    from matplotlib import rcParams
    rcParams['figure.figsize'][1] = float(value)


@type_with_sideeffect(int)
def dpi(value):
    from matplotlib import rcParams
    rcParams['figure.dpi'] = rcParams['savefig.dpi'] = float(value)


@type_with_sideeffect(int)
def transparent(value):
    from matplotlib import rcParams
    rcParams['savefig.transparent'] = bool(value)


figure_parser = argparse.ArgumentParser(add_help=False)
colormap_choices = sorted(cm.cmap_d.keys())
group = figure_parser.add_argument_group(
    'figure options', 'Options that affect figure output format')
group.add_argument(
    '-o', '--output', metavar='FILE.{pdf,png}',
    default='-', type=MatplotlibFigureType(),
    help='name of output file [default: plot to screen]')
group.add_argument(
    '--colormap', default='cylon', choices=colormap_choices,
    type=colormap, metavar='CMAP',
    help='name of matplotlib colormap [default: %(default)s]')
group.add_argument(
    '--help-colormap', action=HelpChoicesAction, choices=colormap_choices)
group.add_argument(
    '--figure-width', metavar='INCHES', type=figwidth, default='8',
    help='width of figure in inches [default: %(default)s]')
group.add_argument(
    '--figure-height', metavar='INCHES', type=figheight, default='6',
    help='height of figure in inches [default: %(default)s]')
group.add_argument(
    '--dpi', metavar='PIXELS', type=dpi, default=300,
    help='resolution of figure in dots per inch [default: %(default)s]')
# FIXME: the savefig.transparent rcparam was added in Matplotlib 1.4,
# but we have to support Matplotlib 1.2 for Scientific Linux 7.
if mpl_version >= '1.4':
    group.add_argument(
        '--transparent', const='1', default='0', nargs='?', type=transparent,
        help='Save image with transparent background [default: false]')
del colormap_choices
del group


# Defer loading SWIG bindings until version string is needed.
class VersionAction(argparse._VersionAction):
    def __call__(self, parser, namespace, values, option_string=None):
        self.version = get_version()
        super(VersionAction, self).__call__(
            parser, namespace, values, option_string)


@type_with_sideeffect(str)
def loglevel_type(value):
    try:
        value = int(value)
    except ValueError:
        value = value.upper()
    logging.basicConfig(level=value)


class LogLevelAction(argparse._StoreAction):

    def __init__(
            self, option_strings, dest, nargs=None, const=None, default=None,
            type=None, choices=None, required=False, help=None, metavar=None):
        # FIXME: this broke because of internal changes in the Python standard
        # library logging module between Python 2.7 and 3.6. We should not rely
        # on these undocumented module variables in the first place.
        try:
            logging._levelNames
        except AttributeError:
            metavar = '|'.join(logging._levelToName.values())
        else:
            metavar = '|'.join(
                _ for _ in logging._levelNames.keys() if isinstance(_, str))
        type = loglevel_type
        super(LogLevelAction, self).__init__(
            option_strings, dest, nargs=nargs, const=const, default=default,
            type=type, choices=choices, required=required, help=help,
            metavar=metavar)


class ArgumentParser(argparse.ArgumentParser):
    """
    An ArgumentParser subclass with some sensible defaults.

    - Any ``.py`` suffix is stripped from the program name, because the
      program is probably being invoked from the stub shell script.

    - The description is taken from the docstring of the file in which the
      ArgumentParser is created.

    - If the description is taken from the docstring, then whitespace in
      the description is preserved.

    - A ``--version`` option is added that prints the version of LALInference.
    """
    def __init__(self,
                 prog=None,
                 usage=None,
                 description=None,
                 epilog=None,
                 parents=[],
                 formatter_class=None,
                 prefix_chars='-',
                 fromfile_prefix_chars=None,
                 argument_default=None,
                 conflict_handler='error',
                 add_help=True):
        if prog is None:
            prog = os.path.basename(sys.argv[0]).replace('.py', '')
        if description is None:
            parent_frame = inspect.currentframe().f_back
            description = parent_frame.f_locals.get('__doc__', None)
            if formatter_class is None:
                formatter_class = argparse.RawDescriptionHelpFormatter
        if formatter_class is None:
            formatter_class = argparse.HelpFormatter
        super(ArgumentParser, self).__init__(
                 prog=prog,
                 usage=usage,
                 description=description,
                 epilog=epilog,
                 parents=parents,
                 formatter_class=argparse.RawDescriptionHelpFormatter,
                 prefix_chars=prefix_chars,
                 fromfile_prefix_chars=fromfile_prefix_chars,
                 argument_default=argument_default,
                 conflict_handler=conflict_handler,
                 add_help=add_help)
        self.register('action', 'glob', GlobAction)
        self.register('action', 'loglevel', LogLevelAction)
        self.register('action', 'version', VersionAction)
        self.add_argument('--version', action='version')
        self.add_argument('-l', '--loglevel', action='loglevel', default='INFO')


class DirType(object):
    """Factory for directory arguments."""

    def __init__(self, create=False):
        self._create = create

    def __call__(self, string):
        if self._create:
            try:
                mkpath(string)
            except DistutilsFileError as e:
                raise argparse.ArgumentTypeError(e.message)
        else:
            try:
                os.listdir(string)
            except OSError as e:
                raise argparse.ArgumentTypeError(e)
        return string


class SQLiteType(argparse.FileType):
    """Open an SQLite database, or fail if it does not exist.
    FIXME: use SQLite URI when we drop support for Python < 3.4.
    See: https://docs.python.org/3.4/whatsnew/3.4.html#sqlite3

    Here is an example of trying to open a file that does not exist for
    reading (mode='r'). It should raise an exception:

    >>> import tempfile
    >>> filetype = SQLiteType('r')
    >>> filename = tempfile.mktemp()
    >>> # Note, simply check or a FileNotFound error in Python 3.
    >>> filetype(filename)
    Traceback (most recent call last):
      ...
    argparse.ArgumentTypeError: ...

    If the file already exists, then it's fine:
    >>> import sqlite3
    >>> filetype = SQLiteType('r')
    >>> with tempfile.NamedTemporaryFile() as f:
    ...     with sqlite3.connect(f.name) as db:
    ...         _ = db.execute('create table foo (bar char)')
    ...     filetype(f.name)
    <sqlite3.Connection object at ...>

    Here is an example of opening a file for writing (mode='w'), which should
    overwrite the file if it exists. Even if the file was not an SQLite
    database beforehand, this should work:

    >>> filetype = SQLiteType('w')
    >>> with tempfile.NamedTemporaryFile(mode='w') as f:
    ...     print('This is definitely not an SQLite file.', file=f)
    ...     f.flush()
    ...     with filetype(f.name) as db:
    ...         db.execute('create table foo (bar char)')
    <sqlite3.Cursor object at ...>

    Here is an example of opening a file for appending (mode='a'), which should
    NOT overwrite the file if it exists. If the file was not an SQLite database
    beforehand, this should raise an exception.

    >>> import pytest
    >>> filetype = SQLiteType('a')
    >>> with tempfile.NamedTemporaryFile(mode='w') as f:
    ...     print('This is definitely not an SQLite file.', file=f)
    ...     f.flush()
    ...     with filetype(f.name) as db:
    ...         db.execute('create table foo (bar char)')
    Traceback (most recent call last):
      ...
    sqlite3.DatabaseError: ...

    And if the database did exist beforehand, then opening for appending
    (mode='a') should not clobber existing tables.

    >>> filetype = SQLiteType('a')
    >>> with tempfile.NamedTemporaryFile() as f:
    ...     with sqlite3.connect(f.name) as db:
    ...         _ = db.execute('create table foo (bar char)')
    ...     with filetype(f.name) as db:
    ...         db.execute('select count(*) from foo').fetchone()
    (0,)
    """

    def __init__(self, mode):
        if mode not in 'arw':
            raise ValueError('Unknown file mode: {}'.format(mode))
        self.mode = mode

    def __call__(self, string):
        try:
            return sqlite.open(string, self.mode)
        except OSError as e:
            raise argparse.ArgumentTypeError(e)


def _sanitize_arg_value_for_xmldoc(value):
    if hasattr(value, 'read'):
        return value.name
    elif isinstance(value, tuple):
        return tuple(_sanitize_arg_value_for_xmldoc(v) for v in value)
    elif isinstance(value, list):
        return [_sanitize_arg_value_for_xmldoc(v) for v in value]
    else:
        return value


def register_to_xmldoc(xmldoc, parser, opts, **kwargs):
    from glue.ligolw.utils import process
    params = {key: _sanitize_arg_value_for_xmldoc(value)
              for key, value in opts.__dict__.items()}
    return process.register_to_xmldoc(xmldoc, parser.prog, params, **kwargs)


start_msg = '\
Waiting for input on stdin. Type control-D followed by a newline to terminate.'
stop_msg = 'Reached end of file. Exiting.'


def iterlines(file, start_message=start_msg, stop_message=stop_msg):
    """Iterate over non-emtpy lines in a file."""
    is_tty = os.isatty(file.fileno())

    if is_tty:
        print(start_message, file=sys.stderr)

    while True:
        # Read a line.
        line = file.readline()

        if not line:
            # If we reached EOF, then exit.
            break

        # Strip off the trailing newline and any whitespace.
        line = line.strip()

        # Emit the line if it is not empty.
        if line:
            yield line

    if is_tty:
        print(stop_message, file=sys.stderr)


from lalinference.bayestar.deprecation import warn
warn('ligo.skymap.tool')
