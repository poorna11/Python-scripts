# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 16:18:10 2016

@author: Poorna
"""
from FormalCharge import charge_magnitude
from FormalCharge import charge_sign
import numpy as np
import csv
import matplotlib.pyplot as plt

sequence = {}
charge_sequence = {}
charge_code = []
net_charge_mag = []
net_charge_sign = []

with open('E:\\FRET\\Charge Conservation\\IscU_alignments_CDD.txt') as f:
    
    lineslist = f.readlines()
for line in lineslist:
    
    line=line.upper()    
    string = line.split()
    #print string
    key = string[1]
    #print key
    if key not in sequence:
        sequence[key]= string[3]
    else:
        sequence[key]+=string[3]

with open('E:\\FRET\\Charge Conservation\\IscU_alignments_CDD.txt') as f:
    f.seek(0)
    first_line = f.readline()
    ref = first_line.split()
    ref_key = ref[1]
    
#print "Sequence:", sequence
#print len(sequence)
for key, value in sequence.iteritems():
    seq = list(value)
    if key == ref_key:
        ref_seq = seq
        #print ref_seq
    #print key, len(seq)
    list_charge_magnitude = []
    list_charge_sign = []
    
    for aa_code in seq:
        if aa_code in charge_magnitude:
            m = charge_magnitude[aa_code]
            s = charge_sign[aa_code]
            
            list_charge_magnitude.append(m)
            list_charge_sign.append(s)
        else:
            print "Invalid amino acid sequence"
            continue
    charge_sequence[key] = (list_charge_magnitude, list_charge_sign)
    #print charge_sequence[key]
    #len(list_charge_magnitude)
    
    net_charge_mag.append(list_charge_magnitude)
    net_charge_sign.append(list_charge_sign)
    
#print net_charge


result1 = np.nanmean(net_charge_mag, axis = 0)
result2 = np.nanmean(net_charge_sign, axis = 0)
result = zip(ref_seq, result1,result2)
#print "Result:", result
out = open('resultcharge.csv', 'wb')
writer = csv.writer(out, delimiter=",", quotechar='"', quoting=csv.QUOTE_ALL)
for item in result:
    writer.writerow(item)
    
out.close()

plot_magnitude = []
plot_sign = []
output = open('resultcharge.csv', 'rb')
for row in csv.reader(output):
    seq = row[0]
    if seq != "-":
        plot_magnitude.append(float(row[1]))
        plot_sign.append(float(row[2]))

"""Plotting the charge magnitude and sign"""

# set axes range
fig = plt.figure()
ax = fig.add_subplot(111)
plt.xlim(-0.5, len(plot_sign)+5)
plt.ylim(-1.2, 1.2)

x = np.arange(1, len(plot_sign)+1)
position = x.tolist()

ax.scatter(position, plot_sign, s=50, c="black", marker = "o", alpha=1.0, label = 'Charge Sign')
ax.plot(position, plot_sign, c= "black")
ax.scatter(position, plot_magnitude, s=90, c="green", marker = "*", alpha=1.0, label = 'Charge Magnitude')
ax.plot(position, plot_magnitude, c= "green")

ax.axhline(y=0)
ax.axhline(linewidth=2, color='black')

plt.legend(loc='lower right');
plt.show()

output.close()

    
