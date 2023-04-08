
# setup.py


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
        
setup(name='Protein binding site prediction',      
      version='1.0',      
      description='Description of my project',      
      author='Pau Pujol Vives, Junhua Ye',      
      author_email='paupujolvives@gmail.com, junhuay00@gmail.com',      
      url='http://www.url_of_my_project.org',      
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


'''
# Create the virtual environment
subprocess.run(['python3.9', '-m', 'venv', 'myenv'])

subprocess.run('source myenv/bin/activate && pip install pandas', shell=True)

# Install packages inside the virtual environment
subprocess.run(['pip', 'install', 'pandas'])
'''
