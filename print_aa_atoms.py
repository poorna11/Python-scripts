# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 20:54:44 2014

@author: Poorna
"""

from Bio.PDB.PDBParser import PDBParser

p= PDBParser(PERMISSIVE=1)
structure_id = '2AW7'
structure = p.get_structure('2AW7', 'E://Leontis//Python scripts//CIF//2AW7.pdb')
model = structure[0]
chain=model['I']
residue= chain[10]
print residue
for atom in residue:
    print atom

    
    



