# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 11:40:12 2016

@author: Poorna
"""

import csv

SASA = {}
Phosphate = []
Sugar= []
     
Phosphate = ["P", "OP1", "OP2"]
Sugar = ["C1'", "C2'", "C3'", "O4'", "C5'", "O5'", "O3'"]
Backbone = Sugar + Phosphate

with open('E:\\Leontis\\Python scripts\\G524.csv', 'rb') as csvfile:
     filereader = csv.reader(csvfile )
     
     for row in filereader:
         key = row[2]
         if key not in SASA:
             SASA[key] = 0
        
         
         atomname = row[1]
         atomarea = float(row[4])
        
         if atomname not in set(Backbone):
            SASA[key] += atomarea
                  
with open('SASA-G24.csv', 'wb') as csvoutput:
    writer = csv.writer(csvoutput)
    for key, value in SASA.items():
        writer.writerow([key, value])
     
    
     
         