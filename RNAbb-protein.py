# -*- coding: utf-8 -*-
"""
Created on Wed Nov 26 12:44:30 2014 @author: Poorna
Name: RNA-protein detection
"""

"""Detect and plot RNA base- amino acid interactions."""
from fr3d.cif.reader import Cif
from fr3d.definitions import aa_fg
from fr3d.definitions import nt_sugar
from fr3d.definitions import nt_phosphate
from fr3d.definitions import tilt_cutoff
from fr3d.definitions import planar_atoms
import numpy as np
import csv


def get_structure(filename):
    with open(filename, 'rb') as raw:
        structure = Cif(raw).structure()
        """All RNA bases are placed in the standard orientation. All Hydrogen
 atoms are inferred. Rotation matrix is calculated for each base."""
        structure.infer_hydrogens()
        return structure


def atom_dist_basepart(base_residue, aa_residue, base_atoms):
    """Calculates atom to atom distance of part "aa_part" of neighboring amino acids 
    of type "aa" from each atom of base. Only returns a pair of aa/nt if two 
    or more atoms are within the cutoff distance"""
    min_distance = 4
    n = 0
    for base_atom in base_residue.atoms(name=base_atoms):
        for aa_atom in aa_residue.atoms(name=aa_fg[aa_residue.sequence]):
            distance = np.subtract(base_atom.coordinates(), aa_atom.coordinates())
            distance = np.linalg.norm(distance)
            if distance <= min_distance:
                n = n+1                
    if n>=1:
        return True
                     
def find_neighbors(bases, amino_acids, aa_part, dist_cent_cutoff):
    """Finds all amino acids of type "aa" for which center of "aa_part" is within
    specified distance of center of bases of type "base" and returns superposed bases"""
    #count_total = 0
    count_pair = 0
    list_aa_coord = [] 
    list_base_coord = [] 
    aas = list(amino_acids)
    aaList_len = None
    new_aaList_len = None
    list_base_aa = []
    cationic_aa = set (["TYR", "HIS", "ARG", "ASN", "GLN", "LYS"])

    for base_residue in bases:
        
        base_seq = base_residue.sequence
        if base_part == 'sugar':
            base_atoms = nt_sugar[base_seq]
        elif base_part == 'phosphate':
            base_atoms = nt_phosphate[base_seq]
                    
        try:
            base_center = base_residue.centers[tuple(base_atoms)]
                    
            if not base_center.any():
                continue
        except:
            print "Incomplete residue", base_residue.unit_id()
            continue
        
        aaList_len = len(list_aa_coord)
        new_aaList_len = 0
        for aa_residue in aas:
            aa_center = aa_residue.centers[aa_part]
            if not aa_center.any():
                continue
            #print base_center, aa_center, aa_residue.unit_id()            
            dist_vector = np.subtract(base_center, aa_center)
            dist_scalar = np.linalg.norm(dist_vector)
            
            if dist_scalar <= dist_cent_cutoff and \
            atom_dist_basepart(base_residue, aa_residue, base_atoms): 
                count_pair = count_pair + 1
                
                rotation_matrix = base_residue.rotation_matrix
                
                base_coordinates = {}
                for base_atom in base_residue.atoms():
                    base_key = base_atom.name
                    base_coordinates[base_key]= translate_rotate(base_atom, base_center, rotation_matrix)
                    
                aa_coordinates = {}                           
                for atom in aa_residue.atoms():
                    key = atom.name
                    aa_coordinates[key]= translate_rotate(atom, base_center, rotation_matrix)
                
                if base_part == 'sugar':                                                        
                    interaction = ribose_stacking(base_residue, aa_residue,
                aa_coordinates, base_coordinates)
                    if interaction is not None:
                        base_aa = annotate(base_residue, aa_residue, interaction)
                        list_base_aa.append(base_aa)
                    
                elif base_part == 'phosphate' and aa_residue.sequence in cationic_aa:
                    interaction = "Ionic"  
                    if interaction is not None:
                        base_aa = annotate(base_residue, aa_residue, interaction)
                        list_base_aa.append(base_aa)
                  
                    for base_atom in base_residue.atoms():
                        list_base_coord.append(base_coordinates)
                    for aa_atom in aa_residue.atoms():
                        list_aa_coord.append(aa_coordinates)

        new_aaList_len = len(list_base_aa)
                           
        new_aaList_len = len(list_aa_coord)
        #list_base_residue.append(base_residue)
    try:
        if aaList_len == new_aaList_len:
            
            print 'No neighbors detected with %s' % aa_residue.sequence
    except:
       print "done"
         
    return list_base_aa, list_aa_coord, list_base_coord 
    

def annotate(base_residue, aa_residue, interaction):
    base_aa = (base_residue, aa_residue, interaction)
    return base_aa
    
def ribose_stacking(base_residue, aa_residue, aa_coordinates, base_coordinates):
    base_x = 0
    base_y = 0
    base_z = 0
    n = 0 

    squared_xy_dist_list = []
    z_dist =[]
    
    """Defines different sets of amino acids"""
    stacked_aa = set (["TRP", "TYR", "PHE", "HIS", "ARG", "ASN", "GLN"])
    non_planar_aa = set (["LEU", "ILE", "PRO", "VAL", "MET"])       
    
    for base_atom in base_residue.atoms(name=nt_sugar[base_residue.sequence]):
        base_key = base_atom.name
        #print base_key, base_coordinates[base_key]
        n = n+1        
        base_x= base_x + base_coordinates[base_key][0]
        base_y= base_y + base_coordinates[base_key][0]
        base_z= base_z + base_coordinates[base_key][0]
    
    base_center_x = base_x/n
    base_center_y = base_y/n
    base_center_z = base_z/n
    
    for aa_atom in aa_residue.atoms(name=aa_fg[aa_residue.sequence]):
        key = aa_atom.name
        aa_x= np.subtract(aa_coordinates[key][0], base_center_x)
        aa_y= np.subtract(aa_coordinates[key][1], base_center_y)
        aa_z = np.subtract(aa_coordinates[key][2], base_center_z)
        squared_xy_dist = (aa_x**2) + (aa_y**2)
        squared_xy_dist_list.append(squared_xy_dist)
        z_dist.append(aa_z)
        
    mean_z = np.mean(aa_z)
    mean_z = abs(mean_z)
    #print "Before", base_residue.unit_id(), aa_residue.unit_id(), min(squared_xy_dist_list), mean_z
    if min(squared_xy_dist_list) <= 2  and mean_z <= 5:
        #print "After", base_residue.unit_id(), aa_residue.unit_id(), min(squared_xy_dist_list), mean_z
        if aa_residue.sequence in stacked_aa:
            print "stacking?", base_residue.unit_id(), aa_residue.unit_id(), min(squared_xy_dist_list), mean_z
            return stacking_angle(base_residue, aa_residue)
        elif aa_residue.sequence in non_planar_aa:
            return stacking_tilt(aa_residue, aa_coordinates)
         
def calculate_angle (base_residue, aa_residue):
    vec1 = vector_calculation(base_residue)
    vec2 = vector_calculation(aa_residue)
                
    angle = angle_between_planes(vec1, vec2)
    return angle

def stacking_angle (base_residue, aa_residue):
    vec1 = vector_calculation(base_residue)
    vec2 = vector_calculation(aa_residue)
         
    angle = angle_between_planes(vec1, vec2)
    print "angle", base_residue.unit_id(), aa_residue.unit_id(),angle
    if angle <=0.65 or 2.25 <= angle <= 3.15:
        return "stacked"
    else:
        return None

def stacking_tilt(aa_residue, aa_coordinates):
    baa_dist_list = []     
        
    for aa_atom in aa_residue.atoms(name=aa_fg[aa_residue.sequence]):
        key = aa_atom.name
        aa_z = aa_coordinates[key][2]
        baa_dist_list.append(aa_z)        
    max_baa = max(baa_dist_list)
    min_baa = min(baa_dist_list)
    diff = max_baa - min_baa
    #print aa_residue.unit_id(), diff
    if diff <= tilt_cutoff[aa_residue.sequence]:
        return "stacked"
    else:
        return None
    
def vector_calculation(residue):
    key = residue.sequence
    P1 = residue.centers[planar_atoms[key][0]]
    P2 = residue.centers[planar_atoms[key][1]]
    P3 = residue.centers[planar_atoms[key][2]]
    #print P1, P2, P3
    vector = np.cross((P2 - P1),(P3-P1))
    return vector

def angle_between_planes (vec1, vec2):
    cosang = np.dot(vec1, vec2)
    sinang = np.linalg.norm(np.cross(vec1, vec2))
    angle = np.arctan2(sinang, cosang)
    return angle

def translate_rotate(atom, reference, rotation_matrix):
     atom_coord = atom.coordinates()
     dist_translate = np.subtract(atom_coord, reference)
     dist_aa_matrix = np.matrix(dist_translate)
     #transposed_rotation = rotation_matrix.transpose()
     rotated_atom = dist_aa_matrix * rotation_matrix
     coord_array = np.array(rotated_atom)
     a = coord_array.flatten()
     coord = a.tolist()    
     return coord
                
def text_output(result_list):
    with open('E:\\Leontis\\Python scripts\\Outputs\\proteinRNA_%s.txt' % PDB, 'wb') as target:
        for result in result_list:
            target.write(str(result))
            target.write("\r\n")
            target.close

def csv_output(result_list):
    with open('E:\\Leontis\\Python scripts\\Outputs\\aa-fg_nt-sugar_%s.csv' % PDB, 'wb') as csvfile:
        fieldnames = ['RNA Chain ID', 'RNA residue','RNA residue number','Protein Chain ID', 'AA residue','AA residue number', 'Interaction']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        for base_residue, aa_residue, interaction in result_list:
            base = base_residue.unit_id()
            aa = aa_residue.unit_id()
            #print base, aa, interaction               
            base_component = str(base).split("|")
            aa_component = str(aa).split("|")
            writer.writerow({'RNA Chain ID': base_component[2], 'RNA residue':base_component[3],'RNA residue number': base_component[4],'Protein Chain ID':aa_component[2],'AA residue': aa_component[3],'AA residue number': aa_component[4], 'Interaction': interaction})
        

"""Inputs a list of PDBs of interest to generate super-imposed plots"""   
PDB_List = ['2AW7']
base_seq_list = ['A','U','C','G']
#base_seq_list = ['A']
aa_list = ['ALA','VAL','ILE','LEU','ARG','LYS','HIS','ASP','GLU','ASN','GLN','THR','SER','TYR','TRP','PHE','PRO','CYS','MET']
#aa_list = ['TYR', 'TRP', 'PHE', 'ARG']

"""Inputs base, amino acid, aa_part of interest and cut-off distance for subsequent functions"""   
if __name__=="__main__":
    for PDB in PDB_List:
        structure = get_structure('E:\\Leontis\\Python scripts\\CIF\\%s.cif' % PDB)
        result_nt_aa = []
        
        aa_part = 'aa_fg'
        base_part = 'sugar'
                 
        bases = structure.residues(sequence= base_seq_list)
        amino_acids = structure.residues(sequence=aa_list)
                
        list_base_aa, list_aa_coord, list_base_coord = find_neighbors(bases, amino_acids, aa_part, 10)
               
        #making the list of resultant RNA-aa pairs
        result_nt_aa.extend(list_base_aa)
        
        #writing out output files                
        csv_output(result_nt_aa)
