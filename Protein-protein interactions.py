# -*- coding: utf-8 -*-
"""
Created on Wed Nov 26 12:44:30 2014 @author: Poorna
"""

"""Detect and plot RNA base- amino acid interactions."""
from fr3d.cif.reader import Cif
import numpy as np
#from fr3d.definitions import aa_fg
from fr3d.definitions import aa_backconnect
from fr3d.definitions import ChainNames
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
    
    #count = 0
    aa1s = list (aa1)
    aa2s = list(aa2)
    list_aa1_aa2 = []
    aa2List_len = None
    new_aa2List_len = None
    
    
    for aa1_residue in aa1s:
        aa1_center = aa1_residue.centers[aa1_part]
        if aa1_center is None:
                continue
        #The following statemnt checks if the aa sequence is a key in the dictionary aa_backconnect
        if aa1_residue.sequence not in aa_backconnect:
            continue
                
        for aa2_residue in aa2s:
            if aa2_residue.sequence not in aa_backconnect:
                 continue
            
            if aa1_residue.chain == aa2_residue.chain:
                continue
        
            aa2_center = aa2_residue.centers[aa2_part]
            
            if aa2_center is None:
                continue
            
            #count_interactions = open('E:\\Leontis\\Python scripts\\%s.csv' % PDB, 'wb')
            
            dist_vector = np.subtract(aa2_center,aa1_center)
            dist_scalar = np.linalg.norm(dist_vector)
            if dist_scalar <= dist_cent_cutoff and \
            atom_dist(aa1_residue, aa2_residue):
            
                #print aa1_residue, aa2_residue
                
                tup1= (aa1_residue, aa2_residue)
                sorted_tup1 = tuple(sorted(tup1))
                
                if sorted_tup1 not in list_aa1_aa2:
                    list_aa1_aa2.append(sorted_tup1)
                                               
                new_aa2List_len = len(list_aa1_aa2)
                
    if aa2List_len == new_aa2List_len:
        print "No neighbors detected"
              
    return list_aa1_aa2


def chainwise_interactions(count_list, PDB):
    chain_count = {}
    
    for aa1, aa2 in count_list:
        chain1 = ChainNames[PDB][aa1.chain]
        chain2 = ChainNames[PDB][aa2.chain]
        key = chain1, chain2
        sorted_key = tuple(sorted(key))
        #sorted is an in-built function that sorts the strings of the key in an order 
        #and produces a list of strings. Tuple converts that list to a tuple

        if sorted_key not in chain_count:
            #count = 1
            chain_count[sorted_key]= 1
            #tup = (chain1, chain2, count)
            #chainwise_output.append(tup)
        else:
            chain_count[sorted_key] += 1
            # a= a+1 is same as a +=1
            
    #print chain_count
    return chain_count

def chainwise_interactions_csv(result_dict):
    with open('E:\\Leontis\\Python scripts\\Outputs\\Chainwise-%s.csv' % PDB, 'wb') as csvfile:
        fieldnames = ['Chain1 ID', 'Chain2 ID', 'No. of interactions']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for (aa1_chain, aa2_chain), count in result_dict.items():
            writer.writerow({'Chain1 ID': aa1_chain, 'Chain2 ID':aa2_chain,'No. of interactions':count})    
                    
def csv_output(result_list):
    with open('E:\\Leontis\\Python scripts\\Outputs\\%s.csv' % PDB, 'wb') as csvfile:
        fieldnames = ['Protein1 Chain ID', 'AA1 residue','AA1 residue number','Protein2 Chain ID', 'AA2 residue','AA2 residue number']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for aa1_residue, aa2_residue in result_list:
                    aa1_component = str(aa1_residue.unit_id()).split("|")
                    aa2_component = str(aa2_residue.unit_id()).split("|")
                                        
                    writer.writerow({'Protein1 Chain ID': ChainNames[PDB][aa1_component[2]], 'AA1 residue':aa1_component[3],'AA1 residue number': aa1_component[4],'Protein2 Chain ID':ChainNames[PDB][aa2_component[2]],'AA2 residue': aa2_component[3],'AA2 residue number': aa2_component[4]})    

"""Inputs a list of PDBs of interest to generate super-imposed plots"""   
PDB_List = ['3I8G']
aa_seq_list = ['ALA','VAL','ARG']
#aa_seq_list = ['MET', 'VAL']
#aa2_seq_list = ['ALA','VAL','ILE','LEU','ARG','LYS','HIS','ASP','GLU','ASN','GLN','THR','SER','TYR','TRP','PHE','PRO','CYS','MET']

"""Inputs base, amino acid, aa_part of interest and cut-off distance for subsequent functions"""   
if __name__=="__main__":
    
    for PDB in PDB_List:
        structure = get_structure('E:\\Leontis\\Python scripts\\CIF\\%s.cif' % PDB)
        result_aa = []
        count_chain = []
        aa1_part = 'aa_fg'
        aa2_part = 'aa_fg'
        
        #chain1 = 'A'
        #chain2 = 'B'
            
        aa1 = structure.residues(sequence = aa_seq_list, symmetry = '1_555')
        aa2 = structure.residues(sequence = aa_seq_list, symmetry = '1_555')
                
        list_aa1_aa2 = find_neighbors(aa1, aa1_part, aa2, aa2_part, 10)
                
        result_aa.extend(list_aa1_aa2)
                       
        #csv_output(result_aa)
        #chainwise_interactions_csv(chainwise_interactions(result_aa, PDB))
        print chainwise_interactions(result_aa, PDB)