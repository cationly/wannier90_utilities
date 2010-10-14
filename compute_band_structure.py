#!/usr/bin/env python

"""
 Author
 ======
        Nicolas Poilvert

 Date
 ====
        October 2010

 This function reads a Wannier90 output  file that stores
 a  lead  matrix (files ending  in "_htB.dat", "_htL.dat"
 and "_htR.dat")  and  computes the band structure of the
 lead. For this, sub-blocks  of  the H00 and H01 matrices
 are  extracted  and  used  to  diagonalize a k-dependent
 matrix whose eigenvalues  are exactly  the band energies
 at  this  k  point.  This  way  of  computing  the  band 
 structure  is similar  to  the  "cut"  mode  implemented
 in Wannier90 but here  no other knowledge  than the lead
 matrix is required.

 :Inputs:  a string corresponding to the full path to the
           file containing the lead matrix

           an integer giving the number of unit cells per
           principal layer

           an integer giving the total number of points
           in the k space grid for sampling the band
           energies

 :Outputs: a file containing the band energies in XY for-
           mat

           if matplotlib is present on the user's comput-
           er, the plot is directly shown on screen.

"""

import os
import sys
import time
try:
  import numpy as np
except:
  print ""
  print " Error in compute_band_structure.py :"
  print " Needed package 'Numpy' could not be found."
  print " Exiting the program..."
  print ""
  sys.exit(1)
import optparse

#---------------------------------------------------#
# Main program which parses the user's command line #
#---------------------------------------------------#

def main():
  #
  # command-line options
  parser = optparse.OptionParser()
  parser.add_option('-m', '--matrix',     
                    dest="lead_matrix",
                    help="full path to the file containing the lead matrix")
  parser.add_option('-c', '--cells',
                    dest="cell_number",
                    type="int",
                    help="number of lead unit cell(s) in one principal layer")
  parser.add_option('-n', '--pointnumber',
                    dest="point_number",
                    type="int",
                    default=314,
                    help="number of points in the interval [-pi/a,pi/a] for band structure")
  options, remainder = parser.parse_args()
  #
  # if the program is called without options it prints the help menu
  if len(sys.argv[1:])==0:
    parser.print_help()
    sys.exit(0)
  #
  # loading the H00 and H01 matrices from the given file
  (H00,H01) = read_lead(options.lead_matrix)
  #
  # calling the routine that computes the band structure
  t_start = time.time()
  (grid,eigs) = compute_band_structure(H00,H01,options.cell_number,options.point_number)
  t_stop  = time.time()
  #
  # dumping the eigenvalues in an XY formatted file
  print_to_file(grid,eigs)
  #
  # plotting the band structure if matplotlib is installed
  try:
    import matplotlib.pyplot as plt
    plot_bands(grid,eigs)
  except:
    print ""
    print " Warning in compute_band_structure.py :"
    print " Module 'matplotlib' could not be found on this"
    print " computer. Dropping the bands plot."
    print ""
  #
  # timing information
  print ""
  elapsed_time = t_stop - t_start
  if elapsed_time < 1.0:
    print " Time to compute band structure : %6.2f ms" %(elapsed_time*1000)
  else:
    print " Time to compute band structure : %6.2f s " %(elapsed_time)
  print ""

  return

#------------------------------------------------------#
# Here is the routine that computes the band structure #
# from the knowledge of H00 and H01                    #
#------------------------------------------------------#

def compute_band_structure(H00,H01,cell_num,pt_num):
  """
  This routine computes the band structure from the
  knowledge of H00, H01 and the number of lead unit
  cells inside ONE principal layer

  :Inputs:  a 2D numpy array representing H00
            a 2D numpy array representing H01
            an integer giving the number of cells per PL
            an integer giving the number of points in the
            k-space grid

  :Outputs: a 1D numpy array corresponding to the k space grid
            a 2D numpy array corresponding to the band energies
            for each grid value (the rows of that array
            represent the bands)
  """
  #
  # making sure that the number of unit cells is compatible
  # with the size of H00 and H01
  if H00.shape!=H01.shape:
    print ""
    print " Error in compute_band_structure :"
    print " H00 and H01 have different shapes"
    print " H00 is a (%i,%i) matrix" %(H00.shape[0],H00.shape[1])
    print " H01 is a (%i,%i) matrix" %(H01.shape[0],H01.shape[1])
    print " Those matrices should have similar dimensions"
    print " Exiting the program..."
    print ""
    sys.exit(1)
  wf_num = H00.shape[0]/cell_num
  if wf_num*cell_num!=H00.shape[0]:
    print ""
    print " Error in compute_band_structure :"
    print " The dimension of H00 is not divisible by"
    print " the number of cells provided."
    print " Exiting the program..."
    print ""
    sys.exit(1)
  #
  # setting up the grid in k space (we only go from
  # 0 to pi/a because of time reversal symmetry)
  from math import pi, cos, sin
  grid = np.linspace(0.0,pi,pt_num,endpoint=True)
  #
  # extracting the fundamental blocks corresponding to the
  # "on-site" , "nearest neighbours", etc... (between unit
  # cells), see matrices H_0i below
  block = []
  for i in xrange(cell_num):
    block.append( H00[:wf_num,i*wf_num:(i+1)*wf_num].astype(complex) )
  block.append( H01[:wf_num,:wf_num].astype(complex) )
  #
  # computing the band structure from the following formula :
  #
  # H_k = H_00 + (H_01 * exp(i.k.R_1) + H_01^+ * exp(-i.k.R_1))
  #            + (H_02 * exp(i.k.R_2) + H_02^+ * exp(-i.k.R_2))
  #            + (H_03 * exp(i.k.R_3) + H_03^+ * exp(-i.k.R_3))
  #            + ...
  # where :
  #        - H_0i is the i-th (wf_num*wf_num) block of matrix
  #          elements between cell 0 and cell i (to the right)
  #        - H_0i^+ means taking the hermitian conjugate of
  #          matrix H_0i
  #        - R_i is nothing else than i*ka where k is the k
  #          vector and a the length of a UNIT CELL of lead
  #          (not the PL length)
  eigs = np.zeros((block[0].shape[0],grid.shape[0]),dtype="float")
  for k in xrange(grid.shape[0]):
    k_hamiltonian = np.zeros((block[0].shape[0],block[0].shape[1]),dtype="complex")
    k_hamiltonian += block[0]
    for j in xrange(1,len(block)):
      k_hamiltonian += block[j]*complex(cos(j*grid[k]),sin(j*grid[k])) + \
          np.conjugate(block[j].transpose())*complex(cos(j*grid[k]),-sin(j*grid[k]))
    k_eigenvalues = np.zeros(k_hamiltonian.shape[0],dtype="float")
    k_eigenvalues = np.real( np.linalg.eigvalsh( k_hamiltonian ) )
    eigs[:,k] = k_eigenvalues

  return (grid,eigs)

#-------------------------------------------#
# Utility routines used in the main program #
#-------------------------------------------#

def read_lead(path):
  """
  This function extracts the H00 and H01 matrices
  from a `*_htL.dat`, `*_htR.dat` or `*_htB.dat` Wannier90
  formatted file.

  The H00 and H01 matrices correspond respectively
  to the **on-site** and **off-diagonal** matrix sub-blocks
  of a lead principal layer.

  :Input(s):  a string corresponding to the name (full
              path) of the file containing the principal
              layer matrices

  :Output(s): a tuple with two rank-2 numpy arrays in the
              following order (H00,H01)
  """
  #
  # tests whether the path is actually existing and
  # if there is a file to read
  if not os.path.exists(path):
    print ""
    print " Error in read_lead :"
    print " It seems that '%s' does not exist" %(path)
    print " Exiting the program..."
    print ""
    sys.exit(1)
  elif not os.path.isfile(path):
    print ""
    print " Error in read_lead :"
    print " It seems that '%s' is not a file" %(path)
    print " Exiting the program..."
    print ""
    sys.exit(1)
  temp_file = open(path,'r')
  lines = temp_file.readlines()
  temp_file.close()

  starting_line_for_H00 = 1
  starting_line_for_H01 = int( (len(lines)+1)/2 )
  #
  # extracting the size of H00
  try:
    H00_size = int(lines[starting_line_for_H00].split()[0])
  except:
    print ""
    print " Error in read_lead :"
    print " could not extract the size of H00 on line 1"
    print " in file %s" %(path)
    print " Exiting the program..."
    print ""
    sys.exit(1)
  #
  # extracting the size of H01
  try:
    H01_size = int(lines[starting_line_for_H01].split()[0])
  except:
    print ""
    print " Error in read_lead :"
    print " could not extract the size of H01 on line %i" %(starting_line_for_H01)
    print " in file %s" %(path)
    print " Exiting the program..."
    print ""
    sys.exit(1)

  H00 = np.zeros((H00_size,H00_size),dtype='float')
  H01 = np.zeros((H01_size,H01_size),dtype='float')
  #
  # extracting the H00 matrix
  iterator = 0
  line_number = starting_line_for_H00+1
  try:
    for i in xrange(starting_line_for_H00+1,starting_line_for_H01):
      line = lines[i].split()
      for item in line:
        col = iterator/H00_size
        row = iterator%H00_size
        H00[row,col] = float(item)
        iterator += 1
      line_number += 1
  except:
    print ""
    print " Error in read_lead :"
    print " could not extract matrix H00 on line %i" %(line_number)
    print " in file %s" %(path)
    print " Exiting the program..."
    print ""
    sys.exit(1)
  #
  # extracting the H01 matrix
  iterator = 0
  line_number = starting_line_for_H01+1
  try:
    for i in xrange(starting_line_for_H01+1,len(lines)):
      line = lines[i].split()
      for item in line:
        col = iterator/H01_size
        row = iterator%H01_size
        H01[row,col] = float(item)
        iterator += 1
      line_number += 1
  except:
    print ""
    print " Error in read_lead :"
    print " could not extract matrix H01 on line %i" %(line_number)
    print " in file %s" %(path)
    print " Exiting the program..."
    print ""
    sys.exit(1)

  return (H00,H01)

def print_to_file(grid,eigs):
  """
  This function takes a grid of k points between -pi/a and pi/a and
  the band energies on that grid and dumps those in a file in XY
  format

  :Inputs:  a 1D numpy array of floats (the grid)
            a 2D numpy array whose rows contain the band
            energies given at the grid points

  :Output:  a file containing the grid as X axis and the band energies
            as Y values. The file is in XY format
  """
  #
  # checking for dimension match between grid and eigs
  if grid.shape[0]!=eigs.shape[1]:
    print ""
    print " Error in compute_band_structure.py :"
    print " There was a mismatch between the size of the grid"
    print " of points in k space and the size of the band energy"
    print " arrays."
    print " Cannot print band energies to file"
    print ""
    sys.exit(1)
  #
  # if no dimension problems, the band energies are printed to file
  outfile = open('./band_energies.dat','w')
  for i in xrange(grid.shape[0]):
    outfile.write(' %10.7f   ' %(grid[i]))
    for j in xrange(eigs.shape[0]):
      outfile.write('%10.7f   ' %(eigs[j,i]))
    outfile.write('\n')
  outfile.close()

  return

def plot_bands(grid,eigs):
  """
  This function uses matplotlib to plot the band energies on screen

  :Inputs:  a 1D numpy array corresponding to the 'x' values 'grid'
            a 2D numpy array corresponding to the bands 'eigs'

  :Output:  a matplotlib graph on screen
  """
  import matplotlib.pyplot as plt
  #
  # building the plot
  plt.figure()
  plt.xticks(np.array([grid.min(),grid.max()]), 
            ('0', r'$\frac{\pi}{a}$'),
            fontsize=20)
  plt.ylabel(r'Band Energies (eV)', fontsize=18)
  for i in xrange(eigs.shape[0]):
    plt.plot(grid,eigs[i,:],color=(float(i)/eigs.shape[0],0.5,0.5),label='band %i' %(i))
  plt.legend()
  plt.show()

  return

if __name__=="__main__":
  main()

