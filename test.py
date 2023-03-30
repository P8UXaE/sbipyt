import proteinclass
from tqdm import tqdm
import os
import sys


mol = proteinclass.readMod2('scPDB/1iki_1/protein.mol2')
print(mol.geometric_binding())