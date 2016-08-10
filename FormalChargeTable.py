# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 16:36:40 2016

@author: Poorna
"""
import numpy 
charge_magnitude = {}
charge_sign = {}
unique_id = {}
nan = numpy.NAN
charge_magnitude = {'G': 0, 'A':0, 'L':0, 'I':0, 'V':0, 'F':0, 'M':0, 'P':0, 'W':0,
 'S':0.5, 'T':0.5, 'Y':0.5, 'C':0.5, 'N':0.5, 'Q':0.5, 'D':1, 'E':1, 'K':1, 'R':1, 'H':0.5
 , '-': nan}
 
charge_sign = {'G': 0, 'A':0, 'L':0, 'I':0, 'V':0, 'F':0, 'M':0, 'P':0, 'W':0,
 'S':0, 'T':0, 'Y':0, 'C':0, 'N':0, 'Q':0, 'D':-1, 'E':-1, 'K':1, 'R':1, 'H':0
 , '-': nan}

unique_id = {'G': 1, 'A':2, 'L':3, 'I':4, 'V':5, 'F':6, 'M':7, 'P':8, 'W':9,
 'S':10, 'T':11, 'Y':12, 'C':13, 'N':14, 'Q':15, 'D':16, 'E':17, 'K':18, 'R':19,
 'H':20 , '-': nan}
