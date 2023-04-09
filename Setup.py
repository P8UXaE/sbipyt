
# setup.py

'''
The setup script is the centre of all activity in building, distributing, and installing modules
using the Distutils. The main purpose of the setup script is to describe your module distribution 
to the Distutils, so that the various commands that operate on your modules do the right thing

you can set up your required packages in the environment created for the work, then it would not affect 
the versions of the packages that you are working in your base.
'''

from setuptools import setup, find_packages
from setuptools.command.install import install
import subprocess
import sys
import os

class InstallCommand(install):
    """Custom install command to create virtual environment and install dependencies."""
    def run(self):
        # Create virtual environment
        subprocess.run([sys.executable, '-m', 'venv', 'myenv'])

        # Activate virtual environment and install pandas
        subprocess.run('source myenv/bin/activate', shell=True)
        #subprocess.run(['myenv/bin/python', '-m', 'pip', 'install', 'pandas'])
        
        # Run regular install command
        install.run(self)
        
setup(name='Protein_binding_site_prediction',      
      version='1.0',      
      description='The required packages needed for the biter program',      
      author='Pau Pujol Vives, Junhua Ye',      
      author_email='paupujolvives@gmail.com, junhuay00@gmail.com',      
      url='https://github.com/P8UXaE/sbipyt',      
      # packages=['proteinclass'], 
      packages=find_packages(),
      install_requires=[
        'numpy==1.23.5',
        'pandas>=1.2.0',
        'matplotlib>=3.3.0',
        'tqdm',
        'scipy',
        'torch_geometric',
        'torch==1.13.1',
        'Biopython',
        'scikit-learn',
        'typing_extensions',
        'numba==0.56.4',
    ],
        cmdclass={
        'install': InstallCommand,
    },
)

