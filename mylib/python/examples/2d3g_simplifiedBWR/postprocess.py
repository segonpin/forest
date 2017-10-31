#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 12:35:58 2017

@author: Antoni Vidal
"""

import xml.etree.ElementTree as ET

def parseFileInline(filename, begin=''):
    """ Parse a File and return what is below a begin title.
    It must be defined a the begin title and the end word or/and the maximum
    number of lines to read.
    """
    f = open(filename)
    line = ' '
    out = []
    while (line != ''):
        line = f.readline()
        # Begin Found
        if line[0:len(begin)] == begin:
            return line[len(begin):-1].split()[0].strip('.')
    f.close()       
    return out

def get_calutation_time(filename) :
    """ Get the calculation time from an log file"""
    
    title= "FOREST:TIME  ::Total time calculated:"
    time = parseFileInline(filename, begin=title)
    return round(float(time),1)

def get_sweeps(filename) :
    """ Get the total number of sweeps from a log file"""
    title= "FOREST::  Total number of sweeps:"
    sweeps = parseFileInline(filename, begin=title)
    return int(sweeps)

def get_iterations(filename) :
    """ Get the total number of iterations from a log file"""
    title= "FOREST::  Total number of iterations:"
    iterations = parseFileInline(filename, begin=title)
    return int(iterations)


def get_keff(filename):
    """ Get the keff from a out.xml file"""
    tree = ET.parse(filename)
    keff = tree.getroot()[0].text
    return round(float(keff), 6)

def latexRow(array):
    """Print an array as a latex table row.
    >>> latexRow(['Table', 2.0, 3.0])
    'Table & 2.0 & 3.0 \\\\\\\\'
    """
    out = ''
    for i in xrange(len(array)):
        if (i == len(array)-1):
            out += str(array[i]) + " \\\\"
        else:
            out += str(array[i]) + " & "
    return out


files = ["2d3g_simplifiedBWR_PI",  "2d3g_simplifiedBWR_Kr1", "2d3g_simplifiedBWR_Kr4",
         "2d3g_simplifiedBWR_PI_FD", "2d3g_simplifiedBWR_Kr1_FD", "2d3g_simplifiedBWR_Kr4_FD"]

## RUNNING ALL EXAMPLES
import os

#os.system("python " + files[0] + ".py")
#os.system("python " + files[1] + ".py")
#os.system("python " + files[2] + ".py")
#os.system("python " + files[3] + ".py")
#os.system("python " + files[4] + ".py")
#os.system("python " + files[5] + ".py")
####################################



file_pi = files[0] + '.log'
file_kr1 = files[1] + '.log'
file_kr4 = files[2] + '.log'

print 'STANDARD'
print latexRow(['PI', get_keff(files[0] + ".out.xml"), get_calutation_time(file_pi), 
          get_iterations(file_pi), get_sweeps(file_pi)/3])
print latexRow(['Arnoldi 1', get_keff(files[1] + ".out.xml"), get_calutation_time(file_kr1), 
          get_iterations(file_kr1), get_sweeps(file_kr1)/3])
print latexRow(['Arnoldi 4', get_keff(files[2] + ".out.xml"), get_calutation_time(file_kr4), 
          get_iterations(file_kr4), get_sweeps(file_kr4)/3])
print ''


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

