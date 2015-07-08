# -*- coding: utf-8 -*-
"""
Created on Wed Nov 26 12:44:30 2014 @author: Poorna
"""

"""Detect and plot RNA base- amino acid interactions."""
from fr3d.cif.reader import Cif
import numpy as np
#from fr3d.definitions import aa_fg
import csv

def get_structure(filename):
    with open(filename, 'rb') as raw:
        structure = Cif(raw).structure()
        """All RNA bases are placed in the standard orientation. All Hydrogen
 atoms are inferred. Rotation matrix is calculated for each base."""
        #structure.infer_hydrogens()
        return structure

def atom_dist(aa_residue1,aa_residue2):
    """Calculates atom to atom distance of part "aa_part" of neighboring amino acids 
    of type "aa" from each atom of base"""
    min_aa = 4
    n = 0
    for aa1_atom in aa_residue1.atoms():
        for aa2_atom in aa_residue2.atoms():
            distance = np.subtract(aa1_atom.coordinates(), aa2_atom.coordinates())
            distance_aa = np.linalg.norm(distance)
            if distance_aa <= min_aa:
                n = n + 1
    if n >= 2:
        return True
             

def find_neighbors(aa1, aa1_part, aa2, aa2_part, dist_cent_cutoff):
    """Finds all amino acids of type "aa" for which center of "aa_part" is within
    specified distance of center of bases of type "base" and returns superposed bases"""
    
    count = 0
    aas = list(aa2)
    aa2List_len = None
    new_aa2List_len = None
    list_aa1_aa2 = []
    
    for aa1_residue in aa1:
        aa1_center = aa1_residue.centers[aa1_part]
        if aa1_center is None:
                continue
        
        
        for aa2_residue in aas:
            aa2_center = aa2_residue.centers[aa2_part]
            if aa2_center is None:
                continue
            
            dist_vector = np.subtract(aa2_center,aa1_center)
            dist_scalar = np.linalg.norm(dist_vector)
            if dist_scalar <= dist_cent_cutoff and \
            atom_dist(aa1_residue, aa2_residue):
                count = count + 1
                print aa1_residue, aa2_residue
                
                tup1= (aa1_residue.unit_id(),aa2_residue.unit_id())
                list_aa1_aa2.append(tup1)
                
    if aa2List_len == new_aa2List_len:
        print "No neighbors detected"
              
    return count, list_aa1_aa2
    
def test_stacking(aa1_center, aa2_center):
    """Detects stacking interaction between amino acids and RNA bases"""
    aa2_x = aa2_center[0]
    aa2_y = aa2_center[1]
    aa1_x = aa1_center[0]
    aa1_y = aa1_center[1]
    a = aa1_x - 3
    b = aa1_x + 3
    c = aa1_y - 3
    d = aa1_y + 3
    return a <= aa2_x <= b and c <= aa2_y <= d
            
def csv_output(result_list):
    with open('E:\\FRET\\protein-protein_%s.csv' % PDB, 'wb') as csvfile:
        fieldnames = ['Protein1 Chain ID', 'AA1 residue','AA1 residue number','Protein2 Chain ID', 'AA2 residue','AA2 residue number']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for aa1_residue, aa2_residue in result_list:
                    aa1_component = str(aa1_residue).split("|")
                    aa2_component = str(aa2_residue).split("|")
                    writer.writerow({'Protein1 Chain ID': aa1_component[2], 'AA1 residue':aa1_component[3],'AA1 residue number': aa1_component[4],'Protein2 Chain ID':aa2_component[2],'AA2 residue': aa2_component[3],'AA2 residue number': aa2_component[4]})    

"""Inputs a list of PDBs of interest to generate super-imposed plots"""   
PDB_List = ['3LVL']
aa1_seq_list = ['ALA','VAL','ILE','LEU','ARG','LYS','HIS','ASP','GLU','ASN','GLN','THR','SER','TYR','TRP','PHE','PRO','CYS','MET']
#aa_list = ['ALA','VAL','ILE','LEU','ARG','LYS','HIS','ASP','GLU','ASN','GLN','THR','SER','TYR','TRP','PHE','PRO','CYS','MET']
aa2_seq_list = ['ALA','VAL','ILE','LEU','ARG','LYS','HIS','ASP','GLU','ASN','GLN','THR','SER','TYR','TRP','PHE','PRO','CYS','MET']

"""Inputs base, amino acid, aa_part of interest and cut-off distance for subsequent functions"""   
if __name__=="__main__":
    
    for PDB in PDB_List:
        structure = get_structure('E:\\FRET\\IscU figures\\%s.cif' % PDB)
        result_aa = []
        for aa1_seq in aa1_seq_list:
            for aa2_seq in aa2_seq_list:
                aa1_part = 'aa_fg'
                aa2_part = 'aa_fg'
        
                chain1 = 'A'
                chain2 = 'B'
            
                aa1 = structure.residues(chain = chain1, sequence = aa1_seq, symmetry = '1_555')
                aa2 = structure.residues(chain = chain2, sequence = aa2_seq, symmetry = '1_555')
                count, list_aa1_aa2 = find_neighbors(aa1, aa1_part, aa2, aa2_part, 10)
                result_aa.extend(list_aa1_aa2)
        
        csv_output(result_aa)