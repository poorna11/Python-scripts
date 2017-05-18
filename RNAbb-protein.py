# -*- coding: utf-8 -*-
"""
Created on Wed Nov 26 12:44:30 2014 @author: Poorna
Name: RNA-protein detection
"""

"""Detect and plot RNA base- amino acid interactions."""
from fr3d.cif.reader import Cif
from fr3d.definitions import RNAconnections
from fr3d.definitions import Ribophos_connect
from fr3d.definitions import aa_connections
from fr3d.definitions import aa_backconnect
from fr3d.definitions import aa_fg
from fr3d.definitions import nt_sugar
from fr3d.definitions import nt_phosphate
from fr3d.definitions import tilt_cutoff
from fr3d.definitions import planar_atoms
import numpy as np
import csv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



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
                #print "rotation matrix", base_residue, rotation_matrix
                
                base_coordinates = {}
                for base_atom in base_residue.atoms():
                    base_key = base_atom.name
                    base_coordinates[base_key]= translate_rotate(base_atom, base_center, rotation_matrix)
                    
                #for base_atoms in base_residue.atoms(name=nt_sugar[base_residue.sequence]):
                    
                aa_coordinates = {}                           
                for atom in aa_residue.atoms():
                    key = atom.name
                    aa_coordinates[key]= translate_rotate(atom, base_center, rotation_matrix)
                
                if base_part == 'sugar':                                                        
                    interaction = ribose_stacking(base_residue, aa_residue,
                aa_coordinates, base_coordinates)
                    
                elif base_part == 'phosphate' and atom_dist_basepart(base_residue, aa_residue, base_atoms):
                    interaction = "Ionic"                
                
                base_aa = annotate(base_residue, aa_residue, interaction)
                   
                if base_aa is not None:                      
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
    #return list_aa_coord, list_base_coord, count, list_base_aa

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
    stacked_aa = set (["TRP", "TYR", "PHE", "HIS", "ARG", "ASN", "GLN", "LEU", "ILE", "PRO"])
    
    for base_atom in base_residue.atoms(name=nt_sugar[base_residue.sequence]):
        base_key = base_atom.name
        #print base_key, base_coordinates[base_key]
        n = n+1        
        base_x= base_x + base_coordinates[base_key][0]
        base_y= base_y + base_coordinates[base_key][0]
        base_z= base_z + base_coordinates[base_key][0]
    #print "number of atoms", n
    base_center_x = base_x/n
    base_center_y = base_y/n
    base_center_z = base_z/n
    
    #print "ribose center",base_coord
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
    print base_residue.unit_id(), aa_residue.unit_id(), min(squared_xy_dist_list), mean_z
    if min(squared_xy_dist_list) <= 5:
        print base_residue.unit_id(), aa_residue.unit_id(), min(squared_xy_dist_list), mean_z
        if aa_residue.sequence in stacked_aa and mean_z <= 3.5:
            #print "stacking?", base_residue.unit_id(), aa_residue.unit_id(), min(squared_xy_dist_list), mean_z
            return stacking_angle(base_residue, aa_residue, min(squared_xy_dist_list))
        else:
            return stacking_tilt(aa_residue, aa_coordinates)
         
def calculate_angle (base_residue, aa_residue):
    vec1 = vector_calculation(base_residue)
    vec2 = vector_calculation(aa_residue)
                
    angle = angle_between_planes(vec1, vec2)
    return angle

def stacking_angle (base_residue, aa_residue, min_dist):
    vec1 = vector_calculation(base_residue)
    vec2 = vector_calculation(aa_residue)
    
    stacked_aa = set (["TRP", "TYR", "PHE", "HIS", "ARG", "ASN", "GLN"])       
    angle = angle_between_planes(vec1, vec2)
    #print base_residue.unit_id(), aa_residue.unit_id(), min_dist, angle
    if aa_residue.sequence in stacked_aa:
        if angle <=0.5 or 2.6 <= angle <= 3.15:
            return "stacked"    

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
    non_planar_aa = set (["LEU", "ILE", "PRO", "VAL", "MET"])       
    if aa_residue.sequence in non_planar_aa and diff <= tilt_cutoff[aa_residue.sequence]:
        return "stacked"
    
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
     #print rotated_atom
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
    with open('E:\\Leontis\\Python scripts\\Outputs\\aa-fg_nt-backbone_%s.csv' % PDB, 'wb') as csvfile:
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
        
        """for base_residue, aa_residue,interaction in result_list:
                    base_component = str(base_residue).split("|")
                    aa_component = str(aa_residue).split("|")
                    writer.writerow({'RNA Chain ID': base_component[2], 'RNA residue':base_component[3],'RNA residue number': base_component[4],'Protein Chain ID':ChainNames[PDB][aa_component[2]],'AA residue': aa_component[3],'AA residue number': aa_component[4], 'Interaction': interaction})"""
                    
               
   
def draw_base(base_seq, ax):
    """Connects atoms to draw neighboring bases and amino acids for 3D plots"""
     #creates lists of rotated base coordinates
    for basecoord_list in list_base_coord:
        new_base_x = []
        new_base_y = []
        new_base_z = [] 
        
        back_base_x = []
        back_base_y = []
        back_base_z = []
        

        try:
            for atomname in RNAconnections[base_seq]:
                coord_base = []
                coord_base= basecoord_list[atomname]
                new_base_x.append(coord_base[0])
                new_base_y.append(coord_base[1])
                new_base_z.append(coord_base[2])
            base_lines= ax.plot(new_base_x, new_base_y, new_base_z, label= 'Base')
            #ax.scatter(basecenter[0], basecenter[1], basecenter[2], zdir='y', color='b', marker='o')
            #ax.scatter(x = 0, y= 0, z= 0, color='b', marker='o')
            plt.setp(base_lines, 'color', 'b', 'linewidth', 1.0)
    
            for atomname in Ribophos_connect[base_seq]:
                back_base=[]           
                back_base= basecoord_list[atomname]
                back_base_x.append(back_base[0])
                back_base_y.append(back_base[1])
                back_base_z.append(back_base[2])
            base_lines= ax.plot(back_base_x, back_base_y, back_base_z, label= 'Base')
            plt.setp(base_lines, 'color', 'g', 'linewidth', 1.0)
            #ax.text(9, 1, 1, base_residue)
        except:
            print "Missing residues"
            continue

def draw_aa(aa, ax):
    #Connects atoms to draw neighboring bases and amino acids for 3D plots
    for aacoord_list in list_aa_coord:
        new_aa_x=[]
        new_aa_y=[]
        new_aa_z=[]
        
        back_aa_x=[]
        back_aa_y=[]
        back_aa_z=[]
            
        try:
            for atomname in aa_connections[aa]:
                coord_aa=[]           
                coord_aa= aacoord_list[atomname]
                new_aa_x.append(coord_aa[0])
                new_aa_y.append(coord_aa[1])
                new_aa_z.append(coord_aa[2])
            aa_lines= ax.plot(new_aa_x, new_aa_y, new_aa_z, label= 'Amino acid')
            plt.setp(aa_lines, 'color', 'r', 'linewidth', 1.0)
        
            for atomname in aa_backconnect[aa]:
                back_aa=[]           
                back_aa= aacoord_list[atomname]
                back_aa_x.append(back_aa[0])
                back_aa_y.append(back_aa[1])
                back_aa_z.append(back_aa[2])
            aa_lines= ax.plot(back_aa_x, back_aa_y, back_aa_z, label= 'Amino acid')
            plt.setp(aa_lines, 'color', 'y', 'linewidth', 1.0)
        except:
            print "Missing residues"
            continue
        
def draw_aa_cent(aa, aa_part, ax):
    #Connects atoms to draw neighboring bases and amino acids for 3D plots
    for aacoord_list in list_aa_coord:
        new_aa_x=[]
        new_aa_y=[]
        new_aa_z=[]
        
        aa_center_x = 0
        aa_center_y = 0
        aa_center_z = 0
        n = 0
        
        if aa_part == 'aa_fg':
            connections = aa_connections
        elif aa_part == 'aa_backbone':
            connections = aa_backconnect
        try:
            for atomname in connections[aa]:
                coord_aa=[]           
                coord_aa= aacoord_list[atomname]
                new_aa_x.append(coord_aa[0])
                new_aa_y.append(coord_aa[1])
                new_aa_z.append(coord_aa[2])
                
                aa_center_x = aa_center_x + coord_aa[0]
                aa_center_y = aa_center_y + coord_aa[1]
                aa_center_z = aa_center_z + coord_aa[2]
                n = n + 1
            ax.scatter(aa_center_x/n, aa_center_y/n, aa_center_z/n, c= 'r', marker = 'o')
        except:
            print "Missing residues"
            continue
                
"""Inputs a list of PDBs of interest to generate super-imposed plots"""   
PDB_List = ['2AW7']
base_seq_list = ['A','U','C','G']
#base_seq_list = ['A']
aa_list = ['ALA','VAL','ILE','LEU','ARG','LYS','HIS','ASP','GLU','ASN','GLN','THR','SER','TYR','TRP','PHE','PRO','CYS','MET']
#aa_list = ['TYR', 'TRP', 'PHE', 'ARG']

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')

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
               
        """ 3D plots of base-aa interactions
        for base, aa, interaction in list_base_aa:
            base_seq = base.sequence
            aa= aa.sequence
                   
            draw_base(base_seq, ax)
            draw_aa(aa, ax)
            #draw_aa_cent(aa, aa_part, ax)
        
            ax.set_xlabel('X Axis')
            ax.set_ylabel('Y Axis')
            ax.set_zlabel('Z Axis')
            ax.set_xlim3d(10, -15)
            ax.set_ylim3d(10, -15)
            ax.set_zlim3d(10, -15)
            plt.title('%s with ' % base_seq +'%s' % aa + ' %s' % aa_part)
            plt.show()
                      """
        #making the list of resultant RNA-aa pairs
        result_nt_aa.extend(list_base_aa)
        
        #writing out output files                
        #text_output(result_nt_aa)
            
        csv_output(result_nt_aa)
