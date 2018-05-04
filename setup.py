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
      python_requires='>=3.5',
      install_requires=[
          'numpy>=1.10',
          'scipy>=0.18',
          'periodictable>=1.4',
          'psutil',
          'seekpath==1.8.1',
          'spglib',
      ],
      scripts=glob('bin/*') + glob('scripts/*'),
      test_suite='matador.tests',
      include_package_data=True,
      extras_require={
          'networks': ['networkx>=2.0'],
          'cifs': ['ase'],
          'database': ['pymongo>=3.4', 'ascii_graph>=1.2'],
          'docs': ['sphinx', 'sphinx-argparse', 'sphinx-rtd-theme', 'sphinxcontrib'],
          'viz': ['ase', 'nglview'],
          'plotting': ['matplotlib>=2.0', 'python-ternary==1.0.3', 'seaborn']
      },
      entry_points={'console_scripts': ['matador = matador.cli.matador_cli:main',
                                        'dispersion = matador.cli.dispersion:main',
                                        'run3 = matador.cli.run3:main']},
      zip_safe=False)
