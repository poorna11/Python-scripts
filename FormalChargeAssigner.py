# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 16:18:10 2016

@author: Poorna
"""
from FormalCharge import charge_magnitude
from FormalCharge import charge_sign
from FormalCharge import unique_id
import numpy as np
import csv
import matplotlib.pyplot as plt

sequence = {}
charge_sequence = {}
charge_code = []
net_charge_mag = []
net_charge_sign = []
net_unique_id = []

with open('E:\\FRET\\Charge Conservation\\IscU_alignments_CDD.txt') as f:
    
    lineslist = f.readlines()
    reference_key = None
    
for line in lineslist:
    line=line.upper()    
    string = line.split()
    

    if len(string)<3:
        print string
        string = string[0], '###', string[1]
        print string
       
    elif len(string) == 5:
        string = string[0]+string[1], string[2], string[3], string[4]
        
    key = string[0]
    
    if reference_key is None:
        reference_key = key
        #print reference_key
    #print key
    if key not in sequence:
        sequence[key]= string[2]
    else:
        sequence[key]+=string[2]

for key, value in sequence.iteritems():
    #n = 1
    
    seq = list(value)
    #print seq
 
    if key == reference_key:
        ref_seq = seq
        
        
    list_charge_magnitude = []
    list_charge_sign = []
    list_unique_id = []
                    
    for aa_code in seq:
                
        if aa_code in charge_magnitude:
            m = charge_magnitude[aa_code]
            s = charge_sign[aa_code]
            u = unique_id[aa_code]
            
            list_charge_magnitude.append(m)
            list_charge_sign.append(s)
            list_unique_id.append(u)
        else:
            print "Invalid amino acid sequence"
            continue
    """print "Reference:",ref_seq
    net = zip(seq, ref_seq)
    
    for aa_code, ref_code in net:
        print (aa_code, ref_code)
        if aa_code == ref_code:
            n = n+1
    
    print n, len(net)
        
    #charge_sequence[key] = (list_charge_magnitude, list_charge_sign)
"""    
    net_charge_mag.append(list_charge_magnitude)
    net_charge_sign.append(list_charge_sign)
    net_unique_id.append(list_unique_id)
    


result1 = np.nanmean(net_charge_mag, axis = 0)
result2 = np.nanmean(net_charge_sign, axis = 0)

consensus_unique_id = np.nanmean(net_unique_id, axis = 0)
consensus_unique_id = (consensus_unique_id).flatten()

#print consensus_unique_id
consensus_seq = []

for id in consensus_unique_id:
    int_id = round(id,0)
    for key, value in unique_id.iteritems():
        if value==int_id:
            consensus_seq.append(key)

result = zip(ref_seq, result1,result2, consensus_seq)
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
fieldnames = ['Residue number', 'Residue', 'Charge Magnitude', 'Charge Sign', 
'Consensus seq']
writer = csv.DictWriter(edited_output, fieldnames=fieldnames)
writer.writeheader()

for row in csv.reader(output):
    seq = row[0]
    
    if seq != "-":
        writer.writerow({'Residue number': res_start,'Residue': row[0], 
        'Charge Magnitude': row[1], 'Charge Sign': row[2], 'Consensus seq': row[3]})
        
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
