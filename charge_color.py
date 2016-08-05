# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 13:38:52 2016

@author: Poorna
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Nov 21 11:11:23 2014

@author: Poorna
"""
from cmd import *
import csv
import numpy as np

charge_magnitude=[]


f = open("E:/FRET/Charge Conservation/resultcharge_edited.csv", "rb")
for i,row in enumerate(csv.reader(f)):
    if i==0:
        continue
    resi = (row[0], float(row[2]), float(row[3]))   
    
    charge_magnitude.append(resi)
    
f.close()

        
for i, magnitude, sign in charge_magnitude:
    
        
    #print "magnitude", magnitude
    #print "sign", sign
    if 0 <= sign <= 1:
        r = (0.7-sign)
        g = (0.7-sign)
        b = 0.7
        
        #color_value_i = abs((magnitude + 0.3)*np.array([r,g,b]))
        color_value_i = (magnitude + 0.3)*np.array([r,g,b])
        color_value_i = color_value_i.tolist()
        
        for index, item in enumerate(color_value_i):
            if item < 0:
                color_value_i[index] = 0
                      
                        
    elif -1 <= sign < 0:
        r = 0.7
        g = (0.7+sign)
        b = (0.7+sign)
                
        color_value_i = ((magnitude + 0.3)*np.array([r,g,b]))
        color_value_i = color_value_i.tolist()
        
        for index, item in enumerate(color_value_i):
            if item < 0:
                color_value_i[index] = 0
    
    color_name = str("color%s"% i)
    #print color_name, color_value_i
    set_color (color_name, color_value_i)
    resnum = str(i)
    color(color_name, "IscU and resi "+resnum)
           
