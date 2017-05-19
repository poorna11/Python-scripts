# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 16:18:10 2016

@author: Poorna
"""
from Bio import AlignIO
from FormalCharge import charge_magnitude
from FormalCharge import charge_sign
import numpy as np
import csv
import matplotlib.pyplot as plt


alignment = AlignIO.read(open("E:\\FRET\\Charge Conservation\\IscU.fasta"), "fasta")
"""alignment = AlignIO.read(sys.argv[1], "fasta")
ref_id = sys.argv[2]"""

#ref_id = input('Please enter the ID of the reference sequence:')
reference_key = None
sequence = {}
charge_code = []
net_charge_mag = []
net_charge_sign = []

for record in alignment:
    key = record.id
    
    if reference_key is None:
        reference_key = key
    sequence[key] = record.seq

for key, value in sequence.iteritems():
    value = value.upper()
    seq = list(value)
    #print seq
    
    list_charge_magnitude = []
    list_charge_sign = []
    list_unique_id = []
    list_ref_percent = []
                    
    for aa_code in seq:
                
        if aa_code in charge_magnitude:
            m = charge_magnitude[aa_code]
            s = charge_sign[aa_code]
                        
            list_charge_magnitude.append(m)
            list_charge_sign.append(s)
            
        else:
            print "Invalid amino acid sequence"
            continue
        
    
    #charge_sequence[key] = (list_charge_magnitude, list_charge_sign)
   
    net_charge_mag.append(list_charge_magnitude)
    net_charge_sign.append(list_charge_sign)
    
    
    if key == reference_key:
        ref_seq = seq
        #print "Reference:", ref_seq
        
 
result1 = np.nanmean(net_charge_mag, axis = 0)
result2 = np.nanmean(net_charge_sign, axis = 0)

result = zip(ref_seq, result1,result2)
#2print "Result:", result

out = open('resultdraft.csv', 'wb')
 
writer = csv.writer(out, delimiter=",", quotechar='"', quoting=csv.QUOTE_ALL)
for item in result:
    writer.writerow(item)
    
out.close()

#Getting user input for assigning correct residue numbers

res_start = input('Please enter the starting position in the reference sequence: ')
plot_start = res_start

plot_magnitude = []
plot_sign = []

output = open('resultdraft.csv', 'rb')

edited_output= open('resultcharge_final.csv', 'wb')
fieldnames = ['Residue number', 'Residue', 'Charge Magnitude', 'Charge Sign']
writer = csv.DictWriter(edited_output, fieldnames=fieldnames)
writer.writeheader()

for row in csv.reader(output):
    seq = row[0]
    
    if seq != "-":
        writer.writerow({'Residue number': res_start,'Residue': row[0], 
        'Charge Magnitude': row[1], 'Charge Sign': row[2]})
        
        plot_magnitude.append(float(row[1]))
        plot_sign.append(float(row[2]))
        res_start = res_start + 1


#Plotting the charge magnitude and sign


fig = plt.figure()
ax = fig.add_subplot(111)

#Setting plot ranges
plot_end = len(plot_sign)+plot_start
x = np.arange(plot_start, plot_end)
position = x.tolist()

# set axes range
plt.xlim(plot_start-2, plot_end+2)
plt.ylim(-1.2, 1.2)

#Coloring the output plot
shade_x = np.arange(plot_start-5,plot_end+5)
ax.fill_between(shade_x, 0.2, -0.2, facecolor='khaki', alpha=0.5)

ax.scatter(position, plot_sign, s=50, c="black", marker = "o", alpha=1.0, label = 'Charge Sign')
ax.plot(position, plot_sign, c= "black")

ax.scatter(position, plot_magnitude, s=90, c="green", marker = "*", alpha=1.0, label = 'Charge Magnitude')
ax.plot(position, plot_magnitude, c= "green")

#Drawing the threshold lines for positive and negative charge conservation

ax.axhline(y=0, linewidth=2, color='black')

ax.axhline(y=-0.75, linewidth=2, color='red')
ax.axhline(y=0.75, linewidth=2, color='blue')


plt.legend(loc='lower right');
plt.show()

output.close()
edited_output.close()

  
