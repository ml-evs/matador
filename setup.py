from setuptools import setup, find_packages
from subprocess import check_output
from glob import glob
try:
    __version__ = check_output(['git', 'describe', '--tags']).decode('utf-8').strip()
    __version__ += '+' + (check_output(['git', 'rev-parse', '--abbrev-ref', 'HEAD'])
                          .decode('utf-8').strip())
except:
    __version__ = 'xxx'

setup(name='matador',
      version=__version__,
      description='MATerial and Atomic Database Of Refined structures.',
      long_description=open('README.md').read(),
      url='https://github.com/ml-evs/matador',
      author='Matthew Evans',
      author_email='me388@cam.ac.uk',
      license='MIT',
      packages=find_packages(),
      install_requires=[
          'numpy>=1.10',
          'scipy>=0.18',
          'pymongo>=3.4',
          'periodictable>=1.4',
          'psutil',
          'seekpath==1.8.1',
          'networkx>=2.0',
          'spglib'],
      scripts=glob('bin/*') + glob('scripts/*'),
      package_data={'matador.scrapers': ['words', 'nouns'],
                    'matador.scrapers.tests': ['data/*'],
                    'matador.tests': ['data/*']},
      test_suite='matador.tests',
      extras_require={
          'pdffit': ['diffpy.Structure',
                     'diffpy.srfit'],
          'cifs': ['ase'],
          'scrape_oqmd': ['mysqlclient'],
          'docs': ['sphinx'],
          'viz': ['ase', 'nglview'],
          'plotting': ['matplotlib>=2.0',
                       'python-ternary==1.0.3',
                       'seaborn'],
          'stats': ['ascii_graph>=1.2'],
          'progressbars': ['progressbar2'],
          'encapsulation': ['pyairss']
      },
      zip_safe=False)
