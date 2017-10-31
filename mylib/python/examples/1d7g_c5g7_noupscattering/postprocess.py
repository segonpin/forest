#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 12:35:58 2017

@author: Antoni Vidal
"""
import sys
sys.path.insert(0, '../..')

from c5g7 import oneDc5g7_noupscattering
from prob.parse_utils import get_calutation_time, get_sweeps, get_iterations
from prob.parse_utils import latexRow


files = []
files.append(oneDc5g7_noupscattering(1.26, 'PI'))
files.append(oneDc5g7_noupscattering(1.26, 'kr', 1, 3))
files.append(oneDc5g7_noupscattering(1.26, 'kr', 1, 5))
files.append(oneDc5g7_noupscattering(1.26, 'kr', 1, 10))

files.append(oneDc5g7_noupscattering(1.50, 'PI'))
files.append(oneDc5g7_noupscattering(1.50, 'kr', 1, 3))
files.append(oneDc5g7_noupscattering(1.50, 'kr', 1, 5))
files.append(oneDc5g7_noupscattering(1.50, 'kr', 1, 10))

files.append(oneDc5g7_noupscattering(2.00, 'PI'))
files.append(oneDc5g7_noupscattering(2.00, 'kr', 1, 3))
files.append(oneDc5g7_noupscattering(2.00, 'kr', 1, 5))
files.append(oneDc5g7_noupscattering(2.00, 'kr', 1, 10))     

#########################################

dr = [0.8953, 0.9458, 0.9713]
time = []
it = []
sw =[]
for i, f in enumerate(files):
    log_file = f + '.log'
    sw.append(get_sweeps(log_file))
    it.append(get_iterations(log_file))
    time.append(get_calutation_time(log_file))
    
    
    
print latexRow(['\delta','Method', 'M', 'N', 'Time (s)'])
for i in range(len(files)):
    if (i% 4) == 0 :
        method = 'Power Iteration'
    else:
        method = 'Krylov-Schur' 
    print latexRow([round(dr[i/(len(files)/len(dr))], 3),method, it[i], sw[i], time[i]])
    


#print latexRow(['PI', get_keff(files[0] + ".out.xml"), get_calutation_time(file_pi), 
#          get_iterations(file_pi), get_sweeps(file_pi)/3])
#print latexRow(['Arnoldi 1', get_keff(files[1] + ".out.xml"), get_calutation_time(file_kr1), 
#          get_iterations(file_kr1), get_sweeps(file_kr1)/3])
#print latexRow(['Arnoldi 4', get_keff(files[2] + ".out.xml"), get_calutation_time(file_kr4), 
#          get_iterations(file_kr4), get_sweeps(file_kr4)/3])
#print ''


#file_pi = files[3] + '.log'
#file_kr1 = files[4] + '.log'
#file_kr4 = files[5] + '.log'
#
#print 'FISSION DENSITY'
#print latexRow(['PI', get_keff(files[3] + ".out.xml"), get_calutation_time(file_pi), 
#          get_iterations(file_pi), get_sweeps(file_pi)/3])
#print latexRow(['Arnoldi 1', get_keff(files[4] + ".out.xml"), get_calutation_time(file_kr1), 
#          get_iterations(file_kr1), get_sweeps(file_kr1)/3])
#print latexRow(['Arnoldi 4', get_keff(files[5] + ".out.xml"), get_calutation_time(file_kr4), 
#          get_iterations(file_kr4), get_sweeps(file_kr4)/3])

