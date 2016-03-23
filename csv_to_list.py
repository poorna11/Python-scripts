# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 14:38:06 2015

@author: Poorna
"""

import csv
your_list = []
with open('E:\\Leontis\\Python scripts\\CIF\\3J9M.csv', 'rb') as f:
    reader = csv.reader(f)
    for entry in reader:
        your_list.extend(entry)

print your_list