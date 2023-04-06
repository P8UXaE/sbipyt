# setup.py

from setuptools import setup 
setup(name='Protein binding site prediction',      
      version='1.0',      
      description='Description of my project',      
      author='Pau Pujol Vives, Junhua Ye',      
      author_email='junhuay00@gmail.com',      
      url='http://www.url_of_my_project.org',      
      # packages=['proteinclass'], 
      install_requires=[
        'numpy>=1.19.0',
        'pandas>=1.2.0',
        'matplotlib>=3.3.0',
        'tqdm',
        'itertools',
        'scipy',
        'torch_geometric',
        'torch>=2.0.0',
        'Biopython',
        'scikit-learn',
        
    ]
      
)

