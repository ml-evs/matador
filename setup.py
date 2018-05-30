# coding: utf-8
# Distributed under the terms of the MIT License.


import subprocess as sp
from setuptools import setup, find_packages
from glob import glob

try:
    __version__ = sp.check_output(['git', 'describe', '--tags']).decode('utf-8').strip()
    __version__ += '+' + (sp.check_output(['git', 'rev-parse', '--abbrev-ref', 'HEAD'])
                          .decode('utf-8').strip())
except sp.CalledProcessError:
    __version__ = '0.9a1'

with open('requirements/requirements.txt', 'r') as f:
    requirements = [line.strip() for line in f.readlines()]

extra_requirements = dict()
for subreq in ['docs', 'test', 'plotting', 'viz', 'db', 'network']:
    with open('requirements/{}_requirements.txt'.format(subreq), 'r') as f:
        extra_requirements[subreq] = [line.strip() for line in f.readlines()]

setup(name='matador',
      version=__version__,
      description='MATerial and Atomic Database Of Refined structures.',
      long_description=open('README.rst').read(),
      url='https://bitbucket.org/ml-evs/matador',
      author='Matthew Evans',
      author_email='me388@cam.ac.uk',
      license='MIT',
      packages=find_packages(),
      python_requires='>=3.5',
      install_requires=requirements,
      scripts=glob('bin/*') + glob('scripts/*'),
      test_suite='matador.tests',
      include_package_data=True,
      extras_require=extra_requirements,
      entry_points={'console_scripts': ['matador = matador.cli.cli:main',
                                        'dispersion = matador.cli.dispersion:main',
                                        'run3 = matador.cli.run3:main']},
      zip_safe=False)
