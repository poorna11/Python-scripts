# -*- coding: utf-8 -*-
"""
Created on Wed Nov 26 12:44:30 2014 @author: Poorna
"""

"""Detect and plot RNA base- amino acid interactions."""
from fr3d.cif.reader import Cif
from fr3d.definitions import RNAconnections
from fr3d.definitions import RNAbaseheavyatoms
from fr3d.definitions import Ribophos_connect
from fr3d.definitions import aa_connections
from fr3d.definitions import aa_backconnect
from fr3d.definitions import aa_fg
from fr3d.definitions import nt_backbone
from fr3d.definitions import tilt_cutoff
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

def atom_dist(base_residue,aa_residue):
    """Calculates atom to atom distance of part "aa_part" of neighboring amino acids 
    of type "aa" from each atom of base"""
    min_baa = 4
    for base_atom in base_residue.atoms():
        base_coord = base_atom.coordinates()        
        for atom in aa_residue.atoms():
            aa_atom = atom.coordinates()
            dist_baa = np.subtract(aa_atom,base_coord)
            baa_scalar = np.linalg.norm(dist_baa)
            #print baa_scalar
            if baa_scalar <= min_baa:
                return True

def atom_dist_basepart(base_residue, aa_residue, atom_names):
    """Calculates atom to atom distance of part "aa_part" of neighboring amino acids 
    of type "aa" from each atom of base"""
    min_distance = 3.5
    for atom in base_residue.atoms(name=atom_names):
        for aa_atom in aa_residue.atoms(name=aa_fg[aa_residue.sequence]):
            # aa_atom = atom.coordinates()
            distance = np.subtract(atom.coordinates(), aa_atom.coordinates())
            distance = np.linalg.norm(distance)
            if distance <= min_distance:
                return distance
        
                     
def find_neighbors(PDB, bases, base_atoms, amino_acids, aa, aa_part, dist_cent_cutoff):
    """Finds all amino acids of type "aa" for which center of "aa_part" is within
    specified distance of center of bases of type "base" and returns superposed bases"""
    #count_total = 0
    count = 0
    #count_stack = 0
    list_aa_coord = [] 
    list_base_coord = [] 
    aas = list(amino_acids)
    aaList_len = None
    new_aaList_len = None
    list_base_aa = []
    
    target = open('E:\\Leontis\\Python scripts\\RNAprotein-count_%s.txt' % PDB, 'a')
    for base_residue in bases:
        try:
            base_center = base_residue.centers[tuple(base_atoms)]
        
            if base_center is None:
                continue
        except:
            print "Incomplete residue", base_residue.unit_id()
            continue
        
        aaList_len = len(list_aa_coord)
        new_aaList_len = 0
        for aa_residue in aas:
            aa_center = aa_residue.centers[aa_part]
            if aa_center is None:
                continue
            
            dist_vector = np.subtract(base_center, aa_center)
            dist_scalar = np.linalg.norm(dist_vector)
            base_seq = base_residue.sequence
            if dist_scalar <= dist_cent_cutoff and \
            atom_dist_basepart(base_residue, aa_residue, base_atoms):
                count = count + 1
                
                rotation_matrix = base_residue.rotation_matrix
                
                base_coordinates = {}
                for base_atom in base_residue.atoms():
                    base_key = base_atom.name
                    base_coordinates[base_key]= translate_rotate(base_atom, base_center, rotation_matrix)
                    # base_coordinates is a list of the Atoms
                
                aa_coordinates = {}                           
                for atom in aa_residue.atoms():
                    key = atom.name
                    aa_coordinates[key]= translate_rotate(atom, base_center, rotation_matrix)
                    #print key, translate_rotate(atom, base_center, base_residue)
                                
                edge = detect_edge(base_residue, aa_residue, aa_coordinates, base_atoms)
                if edge == "WC" or edge == "Sugar" or edge == "Hoogsteen":                
                    #print base_residue, aa_residue, edge
                    
                    dist = atom_dist_basepart(base_residue, aa_residue, base_atoms)
                    tup1= (base_residue.unit_id(),aa_residue.unit_id(), dist, edge)
                    list_base_aa.append(tup1)
                        
                    for base_atom in base_residue.atoms():
                        list_base_coord.append(base_coordinates)
                    for aa_atom in aa_residue.atoms():
                        list_aa_coord.append(aa_coordinates)
                   
                    
                new_aaList_len = len(list_aa_coord)
        #list_base_residue.append(base_residue)
    if aaList_len == new_aaList_len:
        print "No neighbors detected in %s" % PDB +' with %s' % aa                 
   
    print "%d neighbors" % count + ' detected in %s' % PDB + ' with %s' % aa 
    
    result = str(count) + ' between %s' % base_part + ' of %s' % base_seq + ' and %s' % aa_part + ' of %s' % aa 
    target.write(str(result))
    target.write("\n")
    target.close
    
    return list_aa_coord, list_base_coord, count, list_base_aa

def test_stacking(base_residue, aa_residue, aa_coordinates):
    
    """Detects stacking interaction between amino acids and RNA bases"""
    #creates a circle around the rotated base center, which is now at (0,0,0)
    squared_xy_dist_list = []
    for aa_atom in aa_residue.atoms(name=aa_fg[aa_residue.sequence]):
        key = aa_atom.name
        aa_x= aa_coordinates[key][0]
        aa_y= aa_coordinates[key][1]
        squared_xy_dist = (aa_x**2) + (aa_y**2)
        squared_xy_dist_list.append(squared_xy_dist)
    #print aa_residue, min(squared_xy_dist_list)
    if min(squared_xy_dist_list) <= 3:
        
        baa_dist_list = []     
        
        for aa_atom in aa_residue.atoms(name=aa_fg[aa_residue.sequence]):
            key = aa_atom.name
            aa_z = aa_coordinates[key][2]
            baa_dist_list.append(aa_z)        
        max_baa = max(baa_dist_list)
        min_baa = min(baa_dist_list)
        #print 'max distance: %s' % max_baa + ' min distance: %s' % min_baa
        diff = max_baa - min_baa
        
        return diff <= tilt_cutoff[aa_residue.sequence]
    else:
        return False
    

def translate_rotate(atom, reference, rotation_matrix):
     atom_coord = atom.coordinates()
     dist_translate = np.subtract(atom_coord, reference)
     dist_aa_matrix = np.matrix(dist_translate)
     #rotation_matrix = base_residue.rotation_matrix
     #transposed_rotation = rotation_matrix.transpose()
     rotated_atom = dist_aa_matrix * rotation_matrix
     coord_array = np.array(rotated_atom)
     a = coord_array.flatten()
     coord = a.tolist()    
     return coord
                
def detect_edge(base_residue, aa_residue, aa_coordinates, atom_names):
    #Detects edge of nucletide interacting with part "aa_part" of neighboring amino acids of type "aa"
    squared_xy_dist_list = []
    for aa_atom in aa_residue.atoms(name=aa_fg[aa_residue.sequence]):
        key = aa_atom.name
        aa_x= aa_coordinates[key][0]
        aa_y= aa_coordinates[key][1]
        squared_xy_dist = (aa_x**2) + (aa_y**2)
        squared_xy_dist_list.append(squared_xy_dist)
        
    #print aa_residue, min(squared_xy_dist_list)
    if min(squared_xy_dist_list) > 3:
        n = 0
        aa_z = 0
        for aa_atom in aa_residue.atoms(name=aa_fg[aa_residue.sequence]):
            key = aa_atom.name
            aa_x= aa_x + aa_coordinates[key][0]
            aa_y= aa_y + aa_coordinates[key][1]
            aa_z= aa_z + aa_coordinates[key][2]
            n = n+1
        aa_center_x = aa_x/n
        aa_center_y = aa_y/n    
        aa_center_z = aa_z/n
        print base_residue, aa_residue, (aa_center_x, aa_center_y, aa_center_z)
        if 2.9 <= aa_center_x <= 8.0 and -3.9 <= aa_center_y <= 9.2 and -5.0 <= aa_center_z <= 3.3:
            edge = "WC"
            return edge
        elif 0.5 <= aa_center_x <= 5.9 and -7.0 <= aa_center_y <= -3.1 and -3.4 <= aa_center_z <= 2.0:
            edge = "Sugar"
            return edge
        elif -8.9 <= aa_center_x <= 2.4 and -1.5 <= aa_center_y <= 8.6 and -4.0 <= aa_center_z <= 3.9:
            edge = "Hoogsteen"
            return edge
            
def text_output(result_list):
    with open('E:\\Leontis\\Python scripts\\proteinRNA_%s.txt' % PDB, 'wb') as target:
        for result in result_list:
            target.write(str(result))
            target.write("\r\n")
            target.close
        
def csv_output(result_list):
    with open('E:\\Leontis\\Python scripts\\proteinRNA_%s.csv' % PDB, 'wb') as csvfile:
        fieldnames = ['RNA Chain ID', 'RNA residue','RNA residue number','Protein Chain ID', 'AA residue','AA residue number','aa_nt distance','edge']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for base_residue, aa_residue, dist, edge in result_list:
                    base_component = str(base_residue).split("|")
                    aa_component = str(aa_residue).split("|")
                    writer.writerow({'RNA Chain ID': base_component[2], 'RNA residue':base_component[3],'RNA residue number': base_component[4],'Protein Chain ID':aa_component[2],'AA residue': aa_component[3],'AA residue number': aa_component[4],'aa_nt distance':dist, 'edge':edge})
               
   
def draw_base(base_seq, ax):
    """Connects atoms to draw neighboring bases and amino acids for 3D plots"""
     #creates lists of rotated base coordinates
    for basecoord_list in list_base:
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
    for aacoord_list in list_aa:
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
    for aacoord_list in list_aa:
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
PDB_List = ['1S72']
base_seq_list = ['G']
#aa_list = ['ALA','VAL','ILE','LEU','ARG','LYS','HIS','ASP','GLU','ASN','GLN','THR','SER','TYR','TRP','PHE','PRO','CYS','MET']
aa_list = ['LYS']

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

"""Inputs base, amino acid, aa_part of interest and cut-off distance for subsequent functions"""   
if __name__=="__main__":
    for PDB in PDB_List:
        structure = get_structure('E:\\Leontis\\Python scripts\\%s.cif' % PDB)
        result_nt_aa = []
        
        for base_seq in base_seq_list:
            for aa in aa_list:
                aa_part = 'aa_fg'
                base_part = 'base'
                
                residue_atoms = []
                if base_part == 'base':
                    residue_atoms = RNAbaseheavyatoms[base_seq]
                elif base_part == 'nt_backbone':
                    residue_atoms = nt_backbone[base_seq]
                    
                 
                bases = structure.residues(sequence= base_seq)
                amino_acids = structure.residues(sequence=aa)
                
                list_aa, list_base, count, list_base_aa = find_neighbors(PDB, bases, residue_atoms, amino_acids, aa, aa_part, 7)
                
                # 3D plots of base-aa interactions
                draw_base(base_seq, ax)
                draw_aa(aa, ax)
                #draw_aa_cent(aa, aa_part, ax)
                
                ax.set_xlabel('X Axis')
                ax.set_ylabel('Y Axis')
                ax.set_zlabel('Z Axis')
                ax.set_xlim3d(10, -15)
                ax.set_ylim3d(10, -15)
                ax.set_zlim3d(10, -15)
                #plt.title('%s with ' % base_seq +'%s' % aa + ' %s' % aa_part)
                plt.show()
                
                #making the list of resultant RNA-aa pairs
                result_nt_aa.extend(list_base_aa)
        
        #writing out output files                
        #text_output(result_nt_aa)
        csv_output(result_nt_aa)