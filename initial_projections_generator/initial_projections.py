#!/usr/bin/env python

#################################################
# Nicolas Poilvert, April 2009                  #
#                                               #
# This program computes the initial projections #
# necessary for a Wannier function calculation. #
# It  reads  the  atomic  positions from a file #
# written in xyz format, and the periodicity of #
# the  system  from  the user. Then it tries to #
# find an initial guess for the projections ba- #
# sed on the nature  of the chemical species in #
# the system.                                   #
#                                               #
# INPUT :                                       #
#   - a file containing the atomic positions in #
#     xyz format                                #
#   - the periodicity of the system (1D,2D,3D)  #
#   - the length  of  the primitive lattice vec #
#     in the periodic directions                #
#                                               #
# OUTPUT :                                      #
#   - a Wannier90 master input file             #
#################################################

import os
import re
import sys
import os.path
import optparse

import chemistry

from chemistry import ProjectionError

def main():
  #
  # parsing  command line arguments
  #
  parser = optparse.OptionParser()

  parser.add_option('-f', '--file',     
    dest="atomic_positions_file",
    help="name of the file containing the atomic positions in xyz format and in Angstroms")

  parser.add_option('-d', '--debug',    
    action="store_true",
    dest="debug_flag",
    default=False,
    help="flag to switch on and off debugging mode")

  options, remainder = parser.parse_args()
  # if the program is called without options it prints the help menu
  if len(sys.argv[1:])==0:
    parser.print_help()
    sys.exit(0)

  #####
  # 1 #
  #####
  ###################################
  # extracting the atomic positions #
  # from the file                   #
  ###################################
  #
  file  = open(options.atomic_positions_file,'r')

  lines = file.readlines()

  float_regex        = r'[+-]?[0-9]*\.[0-9]*[eE]?[+-]?[0-9]*'
  atomic_coordinates = r'(\w+)\s+('+float_regex+r')\s+('+float_regex+r')\s+('+float_regex+r')'

  atoms              = {}
  iterator           = 0
  for line in lines:
    #
    match = re.search(atomic_coordinates, line)
    #
    if match:
      #
      atoms[iterator] = [match.group(1),float(match.group(2)), \
                                        float(match.group(3)), \
                                        float(match.group(4))]
      iterator       += 1

  file.close()

  # printing the atomic positions to the debug file
  #
  if options.debug_flag:
    #
    debug("\n")
    debug(" atomic positions : \n")
    debug("\n")
    dummy = ""
    for key in atoms.keys():
      #
      dummy = dummy + "%s %9.6f %9.6f %9.6f \n" %(atoms[key][0], \
                                                  atoms[key][1], \
                                                  atoms[key][2], \
                                                  atoms[key][3])
    debug(dummy)

  # checking whether all the chemical species of the system are supported
  # by the "chemistry.py" module

  for key in atoms.keys():
    #
    if atoms[key][0] not in chemistry.supported_elements:
      #
      print "                                             "
      print " At least one chemical species of the system "
      print " is not supported by the Chemistry.py module "
      print " Exiting ...                                 "
      print "                                             "
      print " Supported species are : %s " %(chemistry.supported_elements,)
      print "                                             "
      sys.exit(0)

  # asking for the seedname
  #
  seedname = raw_input("\n What is the seedname for this system ?\n | \n seedname = ")

  #####
  # 2 #
  #####
  #####################################################
  # asking the user for the periodicity of the system #
  #####################################################
  #
  dimensionality = int(raw_input("\n Is the system 1D, 2D or 3D periodic ?\n 1=1D\n 2=2D\n 3=3D\n | \n dimensionality = "))

  if dimensionality not in [1,2,3]:
    #
    print "                                  "
    print " dimensionality must be 1, 2 or 3 "
    print " You entered : %i                 " %(dimensionality,)
    print " Exiting...                       "
    print "                                  "
    sys.exit(0)

  if dimensionality==1:
    #
    one_dim_dir = str(raw_input("\n Is the system periodic in x, y or z ?\n | \n periodic direction = "))

    if one_dim_dir=='x':
      a_dir = 1
    elif one_dim_dir=='y':
      a_dir = 2
    elif one_dim_dir=='z':
      a_dir = 3
    else:
      print "                                   "
      print " direction '%s' was not recognised " %(one_dim_dir,)
      print "                                   "
      sys.exit(2)

    a_length = float(raw_input("\n What is the length a of the periodic direction in Angstrom(s) ?\n | \n a = "))

  elif dimensionality==2:
    #
    two_dim_dir_x = str(raw_input("\n What is the first  periodic direction (x, y or z) ?\n | \n first  periodic direction = "))
    two_dim_dir_y = str(raw_input("\n What is the second periodic direction (x, y or z) ?\n | \n second periodic direction = "))

    if two_dim_dir_x=='x':
      a_dir = 1
    elif two_dim_dir_x=='y':
      a_dir = 2
    elif two_dim_dir_x=='z':
      a_dir = 3
    else:
      print "                                   "
      print " direction '%s' was not recognised " %(two_dim_dir_x,)
      print "                                   "
      sys.exit(2)

    if two_dim_dir_y=='x':
      b_dir = 1
    elif two_dim_dir_y=='y':
      b_dir = 2
    elif two_dim_dir_y=='z':
      b_dir = 3
    else:
      print "                                   "
      print " direction '%s' was not recognised " %(two_dim_dir_y,)
      print "                                   "
      sys.exit(2)

    a_length = float(raw_input("\n What is the length a of the first  periodic direction in Angstrom(s) ?\n | \n a = "))
    b_length = float(raw_input("\n What is the length b of the second periodic direction in Angstrom(s) ?\n | \n b = "))

  else:
    #
    three_dim_dir_x = str(raw_input("\n What is the first  periodic direction (x, y or z) ?\n | \n first  periodic direction = "))
    three_dim_dir_y = str(raw_input("\n What is the second periodic direction (x, y or z) ?\n | \n second periodic direction = "))
    three_dim_dir_z = str(raw_input("\n What is the third  periodic direction (x, y or z) ?\n | \n third  periodic direction = "))

    if three_dim_dir_x=='x':
      a_dir = 1
    elif three_dim_dir_x=='y':
      a_dir = 2
    elif three_dim_dir_x=='z':
      a_dir = 3
    else:
      print "                                   "
      print " direction '%s' was not recognised " %(three_dim_dir_x,)
      print "                                   "
      sys.exit(2)

    if three_dim_dir_y=='x':
      b_dir = 1
    elif three_dim_dir_y=='y':
      b_dir = 2
    elif three_dim_dir_y=='z':
      b_dir = 3
    else:
      print "                                   "
      print " direction '%s' was not recognised " %(three_dim_dir_y,)
      print "                                   "
      sys.exit(2)

    if three_dim_dir_z=='x':
      c_dir = 1
    elif three_dim_dir_z=='y':
      c_dir = 2
    elif three_dim_dir_z=='z':
      c_dir = 3
    else:
      print "                                   "
      print " direction '%s' was not recognised " %(three_dim_dir_z,)
      print "                                   "
      sys.exit(2)

    a_length = float(raw_input("\n What is the length a of the first  periodic direction in Angstrom(s) ?\n | \n a = "))
    b_length = float(raw_input("\n What is the length b of the second periodic direction in Angstrom(s) ?\n | \n b = "))
    c_length = float(raw_input("\n What is the length c of the third  periodic direction in Angstrom(s) ?\n | \n c = "))

  # printing the atomic positions to the debug file
  #
  if options.debug_flag:
    #
    debug("\n")
    debug(" periodicity of the system : \n")
    debug("\n")
    if dimensionality==1:
      #
      debug("one_dim_dir = '%s'    \n" %(one_dim_dir,))
      debug("a_dir       =  %i     \n" %(a_dir,))
      debug("a_length    = %11.5f A\n" %(a_length,))
    elif dimensionality==2:
      #
      debug("two_dim_dir_x = '%s'    \n" %(two_dim_dir_x,))
      debug("two_dim_dir_y = '%s'    \n" %(two_dim_dir_y,))
      debug("a_dir         =  %i     \n" %(a_dir,))
      debug("b_dir         =  %i     \n" %(b_dir,))
      debug("a_length      = %11.5f A\n" %(a_length,))
      debug("b_length      = %11.5f A\n" %(b_length,))
    else:
      #
      debug("three_dim_dir_x = '%s'    \n" %(three_dim_dir_x,))
      debug("three_dim_dir_y = '%s'    \n" %(three_dim_dir_y,))
      debug("three_dim_dir_z = '%s'    \n" %(three_dim_dir_z,))
      debug("a_dir           =  %i     \n" %(a_dir,))
      debug("b_dir           =  %i     \n" %(b_dir,))
      debug("c_dir           =  %i     \n" %(c_dir,))
      debug("a_length        = %11.5f A\n" %(a_length,))
      debug("b_length        = %11.5f A\n" %(b_length,))
      debug("c_length        = %11.5f A\n" %(c_length,))

  #####
  # 3 #
  #####
  #########################################
  # creating a dictionnary containing the #
  # periodic images of the atoms          #
  #########################################
  #
  periodic_atoms = {}

  periodic_atoms_query = {}

  reverse_periodic_atoms_query = {}

  if dimensionality==1:
    #
    iterator = 0
    #
    for key in atoms.keys():
      #
      for alpha in xrange(-1,2):
        #
        if not (alpha==0):
          #
          dummy        = []
          # nature and coordinates of the original atom
          dummy.append(atoms[key][0])
          dummy.append(atoms[key][1])
          dummy.append(atoms[key][2])
          dummy.append(atoms[key][3])
          # computing the coordinate in the periodic direction
          dummy[a_dir] = dummy[a_dir] + alpha * a_length
          periodic_atoms[iterator] = dummy
          # periodic_atoms_query
          periodic_atoms_query[iterator] = "%i,%i" %(key,int(alpha*a_dir))
          # reverse_periodic_atoms_query
          reverse_periodic_atoms_query["%i,%i" %(key,int(alpha*a_dir))] = iterator
          iterator    += 1

  elif dimensionality==2:
    #
    iterator = 0
    #
    for key in atoms.keys():
      #
      for alpha in xrange(-1,2):
        #
        for beta in xrange(-1,2):
          #
          if not (alpha==0 and beta==0):
            #
            dummy        = []
            dummy.append(atoms[key][0])
            dummy.append(atoms[key][1])
            dummy.append(atoms[key][2])
            dummy.append(atoms[key][3])
            dummy[a_dir] = dummy[a_dir] + alpha * a_length
            dummy[b_dir] = dummy[b_dir] + beta  * b_length
            periodic_atoms[iterator] = dummy
            # periodic_atoms_query
            periodic_atoms_query[iterator] = "%i,%i,%i" %(key,int(alpha*a_dir),int(beta*b_dir))
            # reverse_periodic_atoms_query
            reverse_periodic_atoms_query["%i,%i,%i" %(key,int(alpha*a_dir),int(beta*b_dir))] = iterator
            iterator    += 1

  else:
    #
    iterator = 0
    #
    for key in atoms.keys():
      #
      for alpha in xrange(-1,2):
        #
        for beta in xrange(-1,2):
          #
          for gamma in xrange(-1,2):
            #
            if not (alpha==0 and beta==0 and gamma==0):
              #
              dummy        = []
              dummy.append(atoms[key][0])
              dummy.append(atoms[key][1])
              dummy.append(atoms[key][2])
              dummy.append(atoms[key][3])
              dummy[a_dir] = dummy[a_dir] + alpha * a_length
              dummy[b_dir] = dummy[b_dir] + beta  * b_length
              dummy[c_dir] = dummy[c_dir] + gamma * c_length
              periodic_atoms[iterator] = dummy
              # periodic_atoms_query
              periodic_atoms_query[iterator] = "%i,%i,%i,%i" %(key,int(alpha*a_dir),\
                                                                   int(beta*b_dir), \
                                                                   int(gamma*c_dir))
              # reverse_periodic_atoms_query
              reverse_periodic_atoms_query["%i,%i,%i,%i" %(key,int(alpha*a_dir),\
                                                               int(beta*b_dir), \
                                                               int(gamma*c_dir))] = iterator
              iterator    += 1

  # printing periodic atomic positions in debug_file if asked
  #
  if options.debug_flag:
    #
    debug("\n")
    debug(" periodic atomic positions : \n")
    debug("\n")
    dummy = ""
    for key in periodic_atoms.keys():
      #
      dummy = dummy + "%s %9.6f %9.6f %9.6f \n" %(periodic_atoms[key][0], \
                                                  periodic_atoms[key][1], \
                                                  periodic_atoms[key][2], \
                                                  periodic_atoms[key][3])
    debug(dummy)

  #####
  # 4 #
  #####
  ####################################################
  # now we have atoms and periodic_atoms. We need to #
  # find the nearest neighbours for each atom in the #
  # supercell                                        #
  ####################################################
  #
  from math import sqrt

  neighbours = {}

  for key in atoms.keys():
    #
    local_neighbours = {}
    #
    specie1 = atoms[key][0]
    x1 = atoms[key][1]
    y1 = atoms[key][2]
    z1 = atoms[key][3]
    #
    # loop over the other atoms in the supercell
    #
    supercell_neighbours = []
    #
    for key_bis in atoms.keys():
      #
      # we exclude the case where key_bis is key!
      #
      if key!=key_bis:
        #
        specie2 = atoms[key_bis][0]
        x2 = atoms[key_bis][1]
        y2 = atoms[key_bis][2]
        z2 = atoms[key_bis][3]
        #
        distance = sqrt( (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2 )
        #
        if distance <= chemistry.critical_length(specie1,specie2):
          #
          supercell_neighbours.append(key_bis)
    #
    local_neighbours['supercell_neighbours'] = supercell_neighbours
    #
    # loop over the periodic atoms
    #
    periodic_neighbours = []
    #
    for key_bis in periodic_atoms.keys():
      #
      specie2 = periodic_atoms[key_bis][0]
      x2 = periodic_atoms[key_bis][1]
      y2 = periodic_atoms[key_bis][2]
      z2 = periodic_atoms[key_bis][3]
      #
      distance = sqrt( (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2 )
      #
      if distance <= chemistry.critical_length(specie1,specie2):
        #
        periodic_neighbours.append(key_bis)
    #
    local_neighbours['periodic_neighbours'] = periodic_neighbours
    #
    neighbours[key] = local_neighbours

  # printing list of neighbours in debug_file if asked
  #
  if options.debug_flag:
    #
    debug("\n")
    debug(" atoms' neighbours : \n")
    debug("\n")
    dummy = ""
    for key in neighbours.keys():
      #
      dummy = dummy + " atom %5i : \n" %(key,)
      dummy = dummy + " | \n"
      dummy = dummy + " => supercell neighbours : \n"
      dummy = dummy + "      %s \n" %(neighbours[key]["supercell_neighbours"],)
      dummy = dummy + " | \n"
      dummy = dummy + " => periodic  neighbours : \n"
      dummy = dummy + "      %s \n" %(neighbours[key]["periodic_neighbours"],)
      dummy = dummy + "\n"
    debug(dummy)

  #####
  # 5 #
  #####
  #############################################################
  # At  this  point we have  a dictionnary for the atoms, the #
  # periodic images and the neighbours of each supercell atom #
  # Now we need to build lists of forbidden links             #
  #############################################################
  #
  forbidden_supercell_links = []

  for key in neighbours.keys():
    #
    neighbour_list = neighbours[key]['supercell_neighbours']
    number         = len(neighbour_list)
    #
    if number > 0:
      #
      for i in xrange(number):
        #
        if [key,neighbour_list[i]] not in forbidden_supercell_links:
          #
          forbidden_supercell_links.append([neighbour_list[i],key])

  forbidden_periodic_links = []

  for key in neighbours.keys():
    #
    neighbour_list = neighbours[key]['periodic_neighbours']
    number         = len(neighbour_list)
    #
    if number > 0:
      #
      for i in xrange(number):
        #
        if [key,neighbour_list[i]] not in forbidden_periodic_links:
          #
          reverse_pair = get_reverse_pair(key,neighbour_list[i], \
                                          periodic_atoms_query,  \
                                          reverse_periodic_atoms_query)
          forbidden_periodic_links.append(reverse_pair)

  # printing list of forbidden links in debug_file if asked
  #
  if options.debug_flag:
    #
    debug("\n")
    debug(" supercell forbidden links : \n")
    debug("\n")
    dummy = "%s \n" %(forbidden_supercell_links,)
    dummy = dummy + "\n"
    debug(dummy)
    debug(" periodic  forbidden links : \n")
    debug("\n")
    dummy = "%s \n" %(forbidden_periodic_links,)
    dummy = dummy + "\n"
    debug(dummy)

  #####
  # 6 #
  #####
  ###########################################################
  # it is now time to build the Wannier initial projections #
  ###########################################################
  #
  try:
    #
    projections = chemistry.find_projections(atoms,periodic_atoms,neighbours, \
                                                   forbidden_supercell_links, \
                                                   forbidden_periodic_links)

  except ProjectionError, Error:
    #
    print "                                    "
    print " Error when computing projections : "
    print "                                    "
    print Error.message
    sys.exit(2)

  number_of_projections = len(projections)

  #
  # just to check the projections, I create a file with
  # both the atomic positions and the projections in an
  # xyz format
  #
  check_file = open("check_projections.xyz","w")
  #
  check_file.write("%6i \n" %(len(atoms)+number_of_projections,))
  check_file.write("%5i atoms, %5i projections \n" %(len(atoms),number_of_projections))
  for key in atoms.keys():
    #
    check_file.write("%s %12.7f %12.7f %12.7f \n" %(atoms[key][0], \
                                                    atoms[key][1], \
                                                    atoms[key][2], \
                                                    atoms[key][3]))
  for i in xrange(number_of_projections):
    #
    float_regex = r'\s+[+-]?[0-9]*\.[0-9]*[eE]?[+-]?[0-9]*|\s*[+-]?[0-9]*\.[0-9]*[eE]?[+-]?[0-9]*'
    pattern     = r'c=('+float_regex+r'),('+float_regex+r'),('+float_regex+r'):'
    #
    match = re.search(pattern, projections[i])
    #
    if match:
      #
      check_file.write("X %12.7f %12.7f %12.7f \n" %(float(match.group(1)), \
                                                     float(match.group(2)), \
                                                     float(match.group(3))))
  check_file.close()

  #####
  # 7 #
  #####
  ############################################
  # creating the Wannier90 master input file #
  ############################################
  #
  wannier_master_input_string = """\
num_bands       = 'at least %i'
num_wann        = %i
num_iter        = 5000
conv_tol        = 1.0d-10
conv_window     = 20
num_dump_cycles = 200

dis_win_max     = 'highest band energy in nscf calculation'
dis_froz_max    = 'Fermi energy of the nscf calculation'
dis_num_iter    = 5000
dis_conv_tol    = 1.0d-10
dis_conv_window = 20
dis_mix_ratio   = 0.5
guiding_centres = .true.

translate_home_cell = .true.

num_print_cycles = 10

#Begin Kpoint_Path
#G 0.00  0.00  0.00    A 1.00  0.00  0.00
#End Kpoint_Path

#restart = plot
#bands_plot = .false.

iprint = 4

hr_plot = .true.

## SYSTEM PARAMETERS ##

begin unit_cell_cart
ang
 x.xx  0.00  0.00
 0.00  y.yy  0.00
 0.00  0.00  z.zz
end unit_cell_cart

""" %(number_of_projections,   \
      number_of_projections)

  atomic_positions_string = "begin atoms_cart \nAng\n"

  # I prefer not to put the atomic positions automatically
  # since the Wannier code is very picky on the fact  that
  # those positions should be EXACTLY equal to the ones in
  # the original scf and nscf files.
  #
  #for key in atoms.keys():
  #  #
  #  atomic_positions_string += "%s %12.7f %12.7f %12.7f \n" %(atoms[key][0], \
  #                                                            atoms[key][1], \
  #                                                            atoms[key][2], \
  #                                                            atoms[key][3])

  atomic_positions_string += "end atoms_cart \n \n"

  wannier_master_input_string += atomic_positions_string

  projections_string = "begin projections \nAng\n"

  for i in xrange(number_of_projections):
    #
    projections_string += projections[i]

  projections_string += "end projections \n \n"

  wannier_master_input_string += projections_string

  k_points_string = """\
## KPOINTS ##

gamma_only = .true.
mp_grid    =  1 1 1

begin kpoints
0.0000  0.0000  0.0000
end kpoints
"""

  wannier_master_input_string += k_points_string

  wannier_file = open("./%s.win" %(seedname,),"w")

  wannier_file.write(wannier_master_input_string)

  wannier_file.close()

  # explaining the user what needs to be changed in the seedname.win file
  #
  print "                                                     "
  print " Remember to change the following in the %s.win file " %(seedname,)
  print " |                                                   "
  print " 1) number of bands                                  "
  print " 2) the upper values for the frozen and outer windows"
  print " 3) the supercell parameters                         "
  print " 4) copy the atomic positions from nscf file         "
  print " 5) everything concerning k points                   "
  print "                                                     "

  return

def get_reverse_pair(supercell_atom,periodic_atom,periodic_atoms_query,reverse_periodic_atoms_query):
  """This function takes a supercell_atom number and a
periodic atom number and returns a tuple that cor
responds to the reverse pair of atoms.

e.g: atom 4 of the supercell
     and atom 32 of the periodic images

     atom 32 originates  from atom 7  of the supercell 
     displaced in directions +1 and -3.
     So the reverse pair will correspond to 7 and  the
     periodic atom of 4 displaced in directions -1 and
     +3 that is for example atom 27.
     So all in all (4,32) => (7,27)"""

  get_string = periodic_atoms_query[periodic_atom]

  split_string = get_string.split(",")

  new_string = "%i" %(supercell_atom,)

  for i in xrange(1,len(split_string)):
    #
    new_string = new_string + "," + str(-1*int(split_string[i]))

  new_periodic_atom = int(reverse_periodic_atoms_query[new_string])

  new_supercell_atom = int(split_string[0])

  return [new_supercell_atom,new_periodic_atom]

def debug(string_to_file):

  # creating a debug file if it does not
  # exist already
  #
  if not os.path.exists("./debug_file"):
    #
    debug_file = open("debug_file","w")
  else:
    #
    debug_file = open("debug_file","a")

  debug_file.write(string_to_file)

  debug_file.close()

  return

if __name__ == '__main__':
  main()

