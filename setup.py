import os
from setuptools import setup, find_packages

here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the relevant file
with open(os.path.join(here, 'README.rst')) as f:
      long_description = f.read()

def get_version(string):
      """ Parse the version number variable __version__ from a script. """
      import re
      version_re = r"^__version__ = ['\"]([^'\"]*)['\"]"
      version_str = re.search(version_re, string, re.M).group(1)
      return version_str


setup(name='svim',
      version=get_version(open('src/svim/svim').read()),
      description='A structural variant caller for long reads.',
      long_description=long_description,
      url='https://github.com/eldariont/svim',
      author='David Heller',
      author_email='heller_d@molgen.mpg.de',
      license='GPLv3',
      classifiers=[
      'Development Status :: 5 - Production/Stable',
      'Environment :: Console',
      'Intended Audience :: Science/Research',
      'Topic :: Scientific/Engineering :: Bio-Informatics',
      'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
      'Programming Language :: Python :: 3.6'
      ],
      keywords='svim SV PacBio structural variation caller',
      packages = find_packages("src"),
      package_dir = {"": "src"},
      data_files = [("", ["LICENSE"])],
      zip_safe=False,
      install_requires=['pysam', 'numpy', 'scipy', 'matplotlib', 'edlib', 'pyspoa'],
      scripts=['src/svim/svim'])
