from setuptools import setup, find_packages
from subprocess import check_output

__version__ = check_output(["git", "describe", "--tags"]).decode("utf-8").strip()
__version__ += '-' + check_output(["git", "rev-parse", "--abbrev-ref", "HEAD"]).decode("utf-8").strip()

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
          'pymongo>=3.3.0',
          'periodictable>=1.4',
          'spglib>=1.9',
          'matplotlib>=1.5',
          'bokeh>=0.12',
          'python-ternary'],
      scripts=['bin/matador', 'bin/run3', 'bin/oddjob'],
      dependency_links=["https://github.com/marcharper/python-ternary/tarball/master"],
      package_data={'matador.scrapers': ['words', 'nouns'],
                    'matador.scrapers.tests': ['data/*'],
                    'matador.tests': ['data/*']},
      test_suite="matador.tests",
      extras_require={
          'pdffit': ["diffpy.Structure", "diffpy.srfit"],
          'scrape_oqmd': ["mysqlclient"],
          'docs': ["sphinx"],
      },
      zip_safe=False)
