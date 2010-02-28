#!/usr/bin/env python

##################################################
# Nicolas Poilvert, September 2009               #
#                                                #
# Computes the mean quantum conductance of a set #
# of conductance files in the current working di #
# rectory.                                       #
# The format of a conductance file is the follo- #
# wing :                                         #
# <comment> on the first line                    #
# <energy>    <quantum conductance>              #
# see : http://www.wannier.org/user_guide.html   #
# for more details about the format              #
##################################################

import os
import sys
import numpy

# the current working directory should contain the
# conductance files
directory = os.listdir('./')
for file in directory:
  if 'compute_mean_qc' in file:
    directory.remove(file)

# total number of quantum conductance files
number_of_measurements = len(directory)

# some informations about the quantum conductance files
starting_energy = []
ending_energy   = []
energy_step     = []
number_of_lines = []

# intermediate list containing all the conductance files
# that are not to be considered for computing the mean
to_remove       = []

for file in directory:
  try:
    file_unit   = open(file,'r')
    lines       = file_unit.readlines()
    second_line = lines[1].split()
    starting_energy.append( float(second_line[0]) )
    third_line  = lines[2].split()
    energy_step.append( float(third_line[0]) - float(second_line[0]) )
    last_line   = lines[-1].split()
    ending_energy.append( float(last_line[0]) )
    number_of_lines.append(len(lines))
  except:
    print ""
    print " Problem encountered when extracting informations from file : '%s' " %(file,)
    print " Disregarding this file for computing the mean Quantum Conductance "
    print ""
    to_remove.append(file)
  finally:
    file_unit.close()

for file in to_remove:
  directory.remove(file)

number_of_measurements = len(directory)

# making sure that all the conductance files have the exact same
# energy range in order to compute a meaningful average
for i in xrange(len(starting_energy)):
  if (abs(starting_energy[i]-starting_energy[0]) >= 0.000001):
    print ""
    print " Starting energies are not consistent among files "
    print ""
    sys.exit(0)
for i in xrange(len(ending_energy)):
  if (abs(ending_energy[i]-ending_energy[0]) >= 0.000001):
    print ""
    print " Ending energies are not consistent among files "
    print ""
    sys.exit(0)
for i in xrange(len(energy_step)):
  if (abs(energy_step[i]-energy_step[0]) >= 0.000001):
    print ""
    print " Energy steps are not consistent among files "
    print ""
    sys.exit(0)
for i in xrange(len(number_of_lines)):
  if (number_of_lines[i]-number_of_lines[0] != 0):
    print ""
    print " Inconsistent number of lines among files "
    print ""
    sys.exit(0)

import linecache
from math import sqrt

# actually computing the mean conductance
output_file  = open('Mean-QC.dat','w')
output_file.write(' #  Energy      average QC     average QC - sigma     average QC + sigma \n')
num_of_lines = number_of_lines[0]
qc_values    = numpy.zeros(number_of_measurements)

for i in range(1,num_of_lines):
  iterator = 0
  for file in directory:
    qc_values[iterator] = float(linecache.getline(file,i+1).split()[1])
    iterator += 1
  mean_qc  = qc_values.mean()
  variance = sqrt(qc_values.var())/sqrt(number_of_measurements-1)
  output_file.write(' %10.6f   %10.6f        %10.6f             %10.6f \n' 
                      %(starting_energy[0]+(i-1)*energy_step[0],
                        mean_qc,
                        mean_qc-variance,
                        mean_qc+variance))
  linecache.clearcache()

output_file.close()

