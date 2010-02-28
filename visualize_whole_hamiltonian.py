#!/usr/bin/env python

#########################################################
# Nicolas Poilvert, January 2010                        #
#                                                       #
# This program will take a set of Wannier90 Hamiltonian #
# matrices (H_L, H_LC, H_C, H_CR and possibly H_R)  and #
# reconstruct the entire matrix.                        #
# As an output, one ends up with a Wannier90-like total #
# Hamiltonian matrix that is visualized with matplotlib #
#                                                       #
# Note: one needs to be in the directory containing the #
# Wannier90  Hamiltonian  matrices for  the program  to #
# work.                                                 #
#########################################################

import os
import sys
import numpy
import optparse
import matplotlib.pyplot as plt

def main():

  # parsing the command-line arguments
  parser = optparse.OptionParser()

  parser.add_option('-f','--file',
                    dest="seedname",
                    help="seedname of the Wannier90 Hamiltonian matrix files")
  parser.add_option('-w','--wallval',
                    dest="wall_value",
                    type="float",
                    default=15.0,
                    help="value of the pixels that separate the matrices")

  options, remainder = parser.parse_args()

  # if the program is called without options it prints the help menu
  if len(sys.argv[1:])==0:
    parser.print_help()
    sys.exit(0)

  # calling the function that will process the data
  visualize_hamiltonian(options.seedname,options.wall_value)

  return

def visualize_hamiltonian(seedname,wall_value):

  #################################################################
  # 1) trying to determine which Hamiltonian matrices are present #
  #################################################################
  #
  cwd_content = os.listdir('./')
  file_names  = [seedname+'_htL.dat',seedname+'_htLC.dat',seedname+'_htC.dat',
                 seedname+'_htCR.dat']
  for file in file_names:
    if file not in cwd_content:
      print "                                            "
      print " %s is missing. The program cannot continue " %(file)
      print " without this matrix.                       "
      print "                                            "
      sys.exit(0)
  if seedname+'_htR.dat' not in cwd_content:
    print "                                            "
    print " %s is missing, the program assumes that it " %(seedname+'_htR.dat')
    print " can be replaced by %s.                     " %(seedname+'_htL.dat')
    print "                                            "

  ###############################################################################
  # 2) now we will determine the needed matrix dimensions and load the matrices #
  ###############################################################################
  #
  ######
  # HL #
  ######
  datafile = open(seedname+'_htL.dat','r')
  lines    = datafile.readlines()
  ll       = int( lines[1].split()[0] )         # size of the HL00 and HL01 square matrices
  param    = ( len( lines ) - 3 )/2             # parameter needed to extract HL00 and HL01
  HL00     = numpy.zeros((ll,ll),dtype='float') # extracting HL00
  iterator = 0
  for i in xrange(2,2+param):
    line = lines[i].split()
    for item in line:
      col = iterator/ll
      row = iterator%ll
      HL00[row,col] = float(item)
      iterator += 1
  HL01     = numpy.zeros((ll,ll),dtype='float') # extracting HL01
  iterator = 0
  for i in xrange(3+param,3+2*param):
    line = lines[i].split()
    for item in line:
      col = iterator/ll
      row = iterator%ll
      HL01[row,col] = float(item)
      iterator += 1
  datafile.close()
  #######
  # HLC #
  #######
  datafile = open(seedname+'_htLC.dat','r')
  lines    = datafile.readlines()
  lc       = int( lines[1].split()[1] )         # number of columns in HLC
  HLC      = numpy.zeros((ll,lc),dtype='float') # extracting HLC
  iterator = 0
  for i in xrange(2,len(lines)):
    line = lines[i].split()
    for item in line:
      col = iterator/ll
      row = iterator%ll
      HLC[row,col] = float(item)
      iterator += 1
  datafile.close()
  ######
  # HC #
  ######
  datafile = open(seedname+'_htC.dat','r')
  lines    = datafile.readlines()
  cc       = int( lines[1].split()[0] )         # size of HC
  HC       = numpy.zeros((cc,cc),dtype='float') # extracting HC
  iterator = 0
  for i in xrange(2,len(lines)):
    line = lines[i].split()
    for item in line:
      col = iterator/cc
      row = iterator%cc
      HC[row,col] = float(item)
      iterator += 1
  datafile.close()
  #######
  # HCR #
  #######
  datafile = open(seedname+'_htCR.dat','r')
  lines    = datafile.readlines()
  cr       = int( lines[1].split()[0] )         # number of rows in HLC
  rr       = int( lines[1].split()[1] )         # number of columns in HLC
  HCR      = numpy.zeros((cr,rr),dtype='float') # extracting HCR
  iterator = 0
  for i in xrange(2,len(lines)):
    line = lines[i].split()
    for item in line:
      col = iterator/cr
      row = iterator%cr
      HCR[row,col] = float(item)
      iterator += 1
  datafile.close()
  ######
  # HR #
  ######
  if seedname+'_htR.dat' in cwd_content:
    datafile = open(seedname+'_htR.dat','r')
    lines    = datafile.readlines()
    rr       = int( lines[1].split()[0] )         # size of the HR00 and HR01 square matrices
    param    = ( len( lines ) - 3 )/2             # parameter needed to extract HR00 and HR01
    HR00     = numpy.zeros((rr,rr),dtype='float') # extracting HR00
    iterator = 0
    for i in xrange(2,2+param):
      line = lines[i].split()
      for item in line:
        col = iterator/rr
        row = iterator%rr
        HR00[row,col] = float(item)
        iterator += 1
    HR01     = numpy.zeros((rr,rr),dtype='float') # extracting HR01
    iterator = 0
    for i in xrange(3+param,3+2*param):
      line = lines[i].split()
      for item in line:
        col = iterator/rr
        row = iterator%rr
        HR01[row,col] = float(item)
        iterator += 1
    datafile.close()
  else:
    rr = ll
    HR00 = numpy.zeros((rr,rr),dtype='float')
    HR00 = HL00
    HR01 = numpy.zeros((rr,rr),dtype='float')
    HR01 = HL01

  #####################################################
  # 3) we need to build the global Hamiltonian matrix #
  #####################################################
  #
  size = 2*ll + 2*rr + cc + 4 # total size of the global Hamiltonian matrix
  Htot = numpy.zeros((size,size),dtype='float')
  # inserting HL blocks
  Htot[0:ll,0:ll] = HL00
  Htot[0:ll,ll+1:2*ll+1] = HL01
  Htot[ll+1:2*ll+1,0:ll] = numpy.transpose(HL01)
  Htot[ll+1:2*ll+1,ll+1:2*ll+1] = HL00
  # inserting HLC blocks
  Htot[ll+1:2*ll+1,2*ll+2:2*ll+2+lc] = HLC
  Htot[2*ll+2:2*ll+2+lc,ll+1:2*ll+1] = numpy.transpose(HLC)
  # inserting HC block
  Htot[2*ll+2:2*ll+2+cc,2*ll+2:2*ll+2+cc] = HC
  # inserting HCR blocks
  Htot[2*ll+2+cc-cr:2*ll+2+cc,2*ll+2+cc+1:2*ll+2+cc+rr+1] = HCR
  Htot[2*ll+2+cc+1:2*ll+2+cc+rr+1,2*ll+2+cc-cr:2*ll+2+cc] = numpy.transpose(HCR)
  # inserting HR blocks
  Htot[2*ll+2+cc+1:2*ll+2+cc+rr+1,2*ll+2+cc+1:2*ll+2+cc+rr+1] = HR00
  Htot[2*ll+2+cc+1:2*ll+2+cc+rr+1,2*ll+2+cc+rr+2:size] = HR01
  Htot[2*ll+2+cc+rr+2:size,2*ll+2+cc+1:2*ll+2+cc+rr+1] = numpy.transpose(HR01)
  Htot[2*ll+cc+rr+4:size,2*ll+cc+rr+4:size] = HR00
  # inserting delimiters
  Htot[ll,0:2*ll+2] = wall_value
  Htot[0:2*ll+2,ll] = wall_value
  Htot[2*ll+1,0:2*ll+cc+3] = wall_value
  Htot[0:2*ll+cc+3,2*ll+1] = wall_value
  Htot[2*ll+cc+3+rr,2*ll+cc+2:size] = wall_value
  Htot[2*ll+cc+2:size,2*ll+cc+3+rr] = wall_value
  Htot[2*ll+cc+2,2*ll+1:size] = wall_value
  Htot[2*ll+1:size,2*ll+cc+2] = wall_value

  #####################################################
  # 4) we now visualize the global Hamiltonian matrix #
  #####################################################
  #
  plt.matshow(Htot,cmap=plt.cm.gray)
  plt.title(' Hamiltonian matrix, size = (%6i,%6i) ' %(size-4,size-4))
  plt.show()

  return

if __name__=="__main__":
  main()

