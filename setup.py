#from distutils.core import setup
from setuptools import setup
from distutils.extension import Extension
import sys

# fix problems with pythons terrible import system
import os
file_dir = os.path.dirname(os.path.realpath(__file__))

SRC_DIR = 'protocol'

import protocol
version = protocol.__version__
AUTHOR = 'Collin Tokheim'
EMAIL = 'fake@gmail.com'
URL = 'https://github.com/KarchinLab/protocol'
DESCRIPTION = 'Evaluation protocol for computational methods predicting cancer driver genes'
PACKAGES = [SRC_DIR,]
setup(name='protocol',
        version=version,
        description=DESCRIPTION,
        author=AUTHOR,
        author_email=EMAIL,
        url=URL,
        packages=PACKAGES,
        install_requires=['numpy', 'scipy', 'pandas', 'PyYAML',
                          'matplotlib', 'seaborn'],
        entry_points={
            'console_scripts':[
                'driver_protocol = protocol.protocol:cli_main',
                'standard_plots = protocol.standard_plots:cli_main',
            ]
        },
        #long_description=open('README.md').read(),
        classifiers=['Topic :: Scientific/Engineering :: Bio-Informatics',
                    'Environment :: Console',
                    'Intended Audience :: Developers',
                    'Intended Audience :: Science/Research'],
        use_2to3=True,
)
