from subprocess import check_output
from os.path import realpath, dirname
from os import getcwd, chdir

cwd = getcwd()
chdir(dirname(realpath(__file__)))
__version__ = check_output(["git", "describe", "--tags"]).strip()
__version__ += '-' + check_output(["git", "rev-parse", "--abbrev-ref", "HEAD"])
chdir(cwd)
