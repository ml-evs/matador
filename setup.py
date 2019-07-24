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

with open('requirements/pip_requirements.txt', 'r') as f:
    requirements += [line.strip() for line in f.readlines()]

extra_requirements = dict(all=[])
for subreq in ['docs', 'test', 'plotting', 'viz', 'crystal', 'pip']:
    with open('requirements/{}_requirements.txt'.format(subreq), 'r') as f:
        extra_requirements[subreq] = [line.strip() for line in f.readlines()]
        extra_requirements['all'] += extra_requirements[subreq]


setup(name='matador',
      version=__version__,
      summary='MATerial and Atomic Database Of Refined structures.',
      description_file='README.rst',
      url='https://bitbucket.org/ml-evs/matador',
      author='Matthew Evans',
      author_email='matthew@ml-evs.science',
      maintainer='Matthew Evans',
      maintainer_email='matthew@ml-evs.science',
      license='MIT',
      packages=find_packages(),
      python_requires='>=3.6',
      install_requires=requirements,
      scripts=glob('bin/*') + glob('scripts/*'),
      test_suite='matador.tests',
      include_package_data=True,
      extras_require=extra_requirements,
      entry_points={'console_scripts': ['matador = matador.cli.cli:main',
                                        'dispersion = matador.cli.dispersion:main',
                                        'run3 = matador.cli.run3:main']},
      zip_safe=False)
