# setup.py

from setuptools import setup, find_packages

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
    ]
)

