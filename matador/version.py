from subprocess import check_output
from matador.utils.print_utils import print_warning
from os.path import realpath, dirname
from os import getcwd, chdir
from pkg_resources import require

__version__ = require('matador')[0].version
