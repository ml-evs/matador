# coding: utf-8
# Distributed under the terms of the MIT License.

import os
from setuptools import setup, find_packages
from glob import glob

from matador import __version__

with open('requirements/requirements.txt', 'r') as f:
    requirements = [line.strip() for line in f.readlines()]

with open('requirements/pip_requirements.txt', 'r') as f:
    requirements += [line.strip() for line in f.readlines()]

extra_requirements = dict(all=[])
req_files = glob('requirements/*.txt')
for _file in req_files:
    if _file not in ['requirements/requirements.txt', 'requirements/pip_requirements.txt']:
        with open(_file, 'r') as f:
            subreq = _file.split('/')[-1].split('_')[0]
            extra_requirements[subreq] = [line.strip() for line in f.readlines()]
            extra_requirements['all'] += extra_requirements[subreq]

with open("README.rst", "r") as f:
    long_description = f.read()


setup(name='matador-db',
      version=__version__,
      description='MATerial and Atomic Databases Of Refined structures.',
      long_description=long_description,
      url='https://github.com/ml-evs/matador',
      author='Matthew Evans',
      author_email='matador@ml-evs.science',
      maintainer='Matthew Evans',
      maintainer_email='matador@ml-evs.science',
      license='MIT',
      packages=find_packages(),
      python_requires='>=3.6',
      install_requires=requirements,
      scripts=[script for script in glob('scripts/*') if os.path.isfile(script)],
      test_suite='matador.tests',
      include_package_data=True,
      extras_require=extra_requirements,
      classifiers=[
          "Development Status :: 4 - Beta",
          "Intended Audience :: Science/Research",
          "License :: OSI Approved :: MIT License",
          "Programming Language :: Python :: 3",
          "Programming Language :: Python :: 3.6",
          "Programming Language :: Python :: 3.7",
          "Topic :: Scientific/Engineering",
          "Topic :: Scientific/Engineering :: Chemistry",
          "Topic :: Scientific/Engineering :: Physics"
      ],
      entry_points={'console_scripts': ['matador = matador.cli.cli:main',
                                        'dispersion = matador.cli.dispersion:main',
                                        'run3 = matador.cli.run3:main']},
      zip_safe=False)
