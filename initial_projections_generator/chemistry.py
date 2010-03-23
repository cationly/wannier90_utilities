#!/usr/bin/env python

from math import sqrt
from math import cos
from math import sin
import numpy

##################################################
# This module will contains all the informations #
# and functions related somehow to some chemical #
# knowledge about the species in the system.     #
##################################################

supported_elements = ['C','H','N', \
                      'O','F']

###
###### function : "critical_length"
###

def critical_length(specie1,specie2):
  #
  # gives  a cutoff distance  in A beyond which
  # the  two species  are considered non-bonded
  # to each other. This cutoff distance depends
  # on the nature of the 2 species.
  #
  # INPUTS : 2 strings that caracterizes the
  #          two chemical species to look at

  # possible bonds of carbon
  #
  if (specie1=='C' and specie2=='C'):
    #
    cutoff_distance = 1.72
  
  elif (specie1=='C' and specie2=='H') \
    or (specie1=='H' and specie2=='C'):
    #
    cutoff_distance = 1.30

  elif (specie1=='C' and specie2=='N') \
    or (specie1=='N' and specie2=='C'):
    #
    cutoff_distance = 2.20

  elif (specie1=='C' and specie2=='O') \
    or (specie1=='O' and specie2=='C'):
    #
    cutoff_distance = 2.20

  elif (specie1=='C' and specie2=='F') \
    or (specie1=='F' and specie2=='C'):
    #
    cutoff_distance = 1.85

  # possible bonds of hydrogen
  #
  elif (specie1=='H' and specie2=='H'):
    #
    cutoff_distance = 0.94

  elif (specie1=='H' and specie2=='N') \
    or (specie1=='N' and specie2=='H'):
    #
    cutoff_distance = 1.21

  elif (specie1=='H' and specie2=='O') \
    or (specie1=='O' and specie2=='H'):
    #
    cutoff_distance = 1.16

  elif (specie1=='H' and specie2=='F') \
    or (specie1=='F' and specie2=='H'):
    #
    cutoff_distance = 1.12

  # possible bonds of nitrogen
  #
  elif (specie1=='N' and specie2=='N'):
    #
    cutoff_distance = 1.55

  elif (specie1=='N' and specie2=='O') \
    or (specie1=='O' and specie2=='N'):
    #
    cutoff_distance = 1.67

  elif (specie1=='N' and specie2=='F') \
    or (specie1=='F' and specie2=='N'):
    #
    cutoff_distance = 1.64

  # possible bonds of oxygen
  #
  elif (specie1=='O' and specie2=='O'):
    #
    cutoff_distance = 1.69

  elif (specie1=='O' and specie2=='F') \
    or (specie1=='F' and specie2=='O'):
    #
    cutoff_distance = 1.62

  # possible bonds of fluorine
  #
  elif (specie1=='F' and specie2=='F'):
    #
    cutoff_distance = 1.62

  return cutoff_distance

###
###### function : "electronegativity"
###

def electronegativity(element):
  #
  # this little function returns the Pauling
  # electronegativity of the element
  #
  # from Wikipedia on electronegativity
  #
  # INPUT :  a string that caracterizes the
  #          chemical species
  
  dictionnary = {'C' : 2.55, \
                 'H' : 2.20, \
                 'N' : 3.04, \
                 'O' : 3.44, \
                 'F' : 3.98}

  if element in dictionnary.keys():
    #
    value = dictionnary[element]

  else:
    #
    value = 1.0

  return value

###
###### function : "find_projections"
###

def find_projections(atoms,
                     periodic_atoms,
                     neighbours,
                     forbidden_supercell_links,
                     forbidden_periodic_links):

    #
    # This function takes the list of nearest neighbours
    # for each atom of a supercell and finds the non-red
    # ondant wannier initial projections
    #
    
    projections = []

    for atom in atoms.keys():
        #
        # extracting the nature of the atom
        #
        nature = atoms[atom][0]
        #
        # computing the total number of nearest neighbours
        #
        number_of_neighbours = len(neighbours[atom]["supercell_neighbours"]) + \
                               len(neighbours[atom]["periodic_neighbours"])
        #
        # extracting the atomic coordinate of the atom
        #
        x0 = atoms[atom][1]
        y0 = atoms[atom][2]
        z0 = atoms[atom][3]
        #
        # extracting the nature and coodinates of each nearest neighbour
        # along with their "supercell" or "periodic" status and their 
        # index as a supercell or periodic atom
        #
        atom_neighbours = []
        #
        if len(neighbours[atom]["supercell_neighbours"]) > 0:
            #
            for i in xrange(len(neighbours[atom]["supercell_neighbours"])):
                #
                current_neighbour_number = neighbours[atom]["supercell_neighbours"][i]
                current_neighbour = atoms[current_neighbour_number]+["supercell"]+[current_neighbour_number]
                atom_neighbours.append(current_neighbour)
        #
        if len(neighbours[atom]["periodic_neighbours"]) > 0:
            #
            for i in xrange(len(neighbours[atom]["periodic_neighbours"])):
                #
                current_neighbour_number = neighbours[atom]["periodic_neighbours"][i]
                current_neighbour = periodic_atoms[current_neighbour_number]+["periodic"]+[current_neighbour_number]
                atom_neighbours.append(current_neighbour)
        #
        ##########################
        # projections for Carbon #
        ##########################
        #
        if nature=='C':
            #
            if number_of_neighbours > 4 or number_of_neighbours < 2:
                #
                # Error message corresponding to a weird number
                # of nearest neighbours
                #
                msg = " A carbon atom has been found with %i neighbours ! \n" %(number_of_neighbours,) + \
                      " Expecting 2, 3 or 4 neighbours for Carbon...      \n" +                          \
                      " Coordinates : %9.6f %9.6f %9.6f                   \n" %(atoms[atom][1],          \
                                                                                atoms[atom][2],          \
                                                                                atoms[atom][3])
                raise ProjectionError, msg
            else:
                #        
                if number_of_neighbours==4:
                    #
                    #    4  s-type  projections  on  mid-bonds
                    # corrected by elements' electronegativities
                    #
                    for j in xrange(len(atom_neighbours)):
                        #
                        if atom_neighbours[j][-2]=="supercell":
                            couple = [atom,atom_neighbours[j][-1]]
                            possible_bond = True
                            if couple in forbidden_supercell_links:
                                possible_bond = False
                        elif atom_neighbours[j][-2]=="periodic":
                            couple = [atom,atom_neighbours[j][-1]]
                            possible_bond = True
                            if couple in forbidden_periodic_links:
                                possible_bond = False
                        #
                        nature1 = atom_neighbours[j][0]
                        x1 = atom_neighbours[j][1]
                        y1 = atom_neighbours[j][2]
                        z1 = atom_neighbours[j][3]
                        #
                        if possible_bond:
                            #
                            local_projection = single_bond(nature,x0,y0,z0,nature1,x1,y1,z1)
                            projections.append(local_projection)
                #
                elif number_of_neighbours==3:
                    #
                    #  3 s-type projections on mid-bonds corrected
                    # by electronegativity + p_z orthogonal to 3 neighbour
                    # plane
                    #
                    for j in xrange(len(atom_neighbours)):
                        #
                        if atom_neighbours[j][-2]=="supercell":
                            couple = [atom,atom_neighbours[j][-1]]
                            possible_bond = True
                            if couple in forbidden_supercell_links:
                                possible_bond = False
                        elif atom_neighbours[j][-2]=="periodic":
                            couple = [atom,atom_neighbours[j][-1]]
                            possible_bond = True
                            if couple in forbidden_periodic_links:
                                possible_bond = False
                        #
                        nature1 = atom_neighbours[j][0]
                        x1 = atom_neighbours[j][1]
                        y1 = atom_neighbours[j][2]
                        z1 = atom_neighbours[j][3]
                        #
                        if possible_bond:
                            #
                            local_projection = single_bond(nature,x0,y0,z0,nature1,x1,y1,z1)
                            projections.append(local_projection)
                    #
                    # p_z orbital
                    #
                    direction = get_plane_normal(atom_neighbours)
                    #
                    projections.append("c=%11.7f,%11.7f,%11.7f:pz:z=%3.2f,%3.2f,%3.2f \n" %(x0,y0,z0,     \
                                                                                         direction[0], \
                                                                                         direction[1], \
                                                                                         direction[2]))
                #
                elif number_of_neighbours==2:
                    #
                    #  2 s-type projections on mid-bonds corrected by 
                    # electronegativity + p_x + p_y in a plane orthogonal
                    # to the line joining the 2 neighbours and the carbon
                    # atom
                    #
                    for j in xrange(len(atom_neighbours)):
                        #
                        if atom_neighbours[j][-2]=="supercell":
                            couple = [atom,atom_neighbours[j][-1]]
                            possible_bond = True
                            if couple in forbidden_supercell_links:
                                possible_bond = False
                        elif atom_neighbours[j][-2]=="periodic":
                            couple = [atom,atom_neighbours[j][-1]]
                            possible_bond = True
                            if couple in forbidden_periodic_links:
                                possible_bond = False
                        #
                        nature1 = atom_neighbours[j][0]
                        x1 = atom_neighbours[j][1]
                        y1 = atom_neighbours[j][2]
                        z1 = atom_neighbours[j][3]
                        #
                        if possible_bond:
                            #
                            local_projection = single_bond(nature,x0,y0,z0,nature1,x1,y1,z1)
                            projections.append(local_projection)
                    #
                    # p_x and p_y orbitals
                    #
                    directions = get_line_normals(atom_neighbours)
                    #
                    projections.append("c=%11.7f,%11.7f,%11.7f:pz:z=%3.2f,%3.2f,%3.2f \n" %(x0,y0,z0,     \
                                                                                         directions[0], \
                                                                                         directions[1], \
                                                                                         directions[2]))
                    projections.append("c=%11.7f,%11.7f,%11.7f:pz:z=%3.2f,%3.2f,%3.2f \n" %(x0,y0,z0,     \
                                                                                         directions[3], \
                                                                                         directions[4], \
                                                                                         directions[5]))
        #
        ############################
        # projections for Hydrogen #
        ############################
        #
        if nature=='H':
            #
            if number_of_neighbours > 1 or number_of_neighbours < 1:
                #
                # Error message corresponding to a weird number
                # of nearest neighbours
                #
                msg = " An Hydrogen atom has been found with %i neighbours ! \n" %(number_of_neighbours,) + \
                      " Expecting only 1 neighbour for Hydrogen...           \n" +                          \
                      " Coordinates : %9.6f %9.6f %9.6f                      \n" %(atoms[atom][1],          \
                                                                                   atoms[atom][2],          \
                                                                                   atoms[atom][3])
                raise ProjectionError, msg
            else:
                #
                if number_of_neighbours==1:
                    #
                    #  1 s-type projection on mid-bond corrected by 
                    # electronegativity
                    #
                    for j in xrange(len(atom_neighbours)):
                        #
                        if atom_neighbours[j][-2]=="supercell":
                            couple = [atom,atom_neighbours[j][-1]]
                            possible_bond = True
                            if couple in forbidden_supercell_links:
                                possible_bond = False
                        elif atom_neighbours[j][-2]=="periodic":
                            couple = [atom,atom_neighbours[j][-1]]
                            possible_bond = True
                            if couple in forbidden_periodic_links:
                                possible_bond = False
                        #
                        nature1 = atom_neighbours[j][0]
                        x1 = atom_neighbours[j][1]
                        y1 = atom_neighbours[j][2]
                        z1 = atom_neighbours[j][3]
                        #
                        if possible_bond:
                            #
                            local_projection = single_bond(nature,x0,y0,z0,nature1,x1,y1,z1)
                            projections.append(local_projection)
        #
        ############################
        # projections for Fluorine #
        ############################
        #
        if nature=='F':
            #
            if number_of_neighbours > 1 or number_of_neighbours < 1:
                #
                # Error message corresponding to a weird number
                # of nearest neighbours
                #
                msg = " A  Fluorine atom has been found with %i neighbours ! \n" %(number_of_neighbours,) + \
                      " Expecting only 1 neighbour for Fluorine...           \n" +                          \
                      " Coordinates : %9.6f %9.6f %9.6f                      \n" %(atoms[atom][1],          \
                                                                                   atoms[atom][2],          \
                                                                                   atoms[atom][3])
                raise ProjectionError, msg
            else:
                #
                if number_of_neighbours==1:
                    #
                    #  1 s-type projection on mid-bond corrected by 
                    # electronegativity
                    #
                    for j in xrange(len(atom_neighbours)):
                        #
                        if atom_neighbours[j][-2]=="supercell":
                            couple = [atom,atom_neighbours[j][-1]]
                            possible_bond = True
                            if couple in forbidden_supercell_links:
                                possible_bond = False
                        elif atom_neighbours[j][-2]=="periodic":
                            couple = [atom,atom_neighbours[j][-1]]
                            possible_bond = True
                            if couple in forbidden_periodic_links:
                                possible_bond = False
                        #
                        nature1 = atom_neighbours[j][0]
                        x1 = atom_neighbours[j][1]
                        y1 = atom_neighbours[j][2]
                        z1 = atom_neighbours[j][3]
                        #
                        if possible_bond:
                            #
                            local_projection = single_bond(nature,x0,y0,z0,nature1,x1,y1,z1)
                            projections.append(local_projection)
                    #
                    # now the 3 lone pair of fluorine such that the 4
                    # projections creates a regular tetrahedron
                    #
                    # building a first normalized vector in the
                    # direction : Fluorine to first nearest neighbour
                    #
                    x1 = atom_neighbours[0][1] - x0
                    y1 = atom_neighbours[0][2] - y0
                    z1 = atom_neighbours[0][3] - z0
                    length = sqrt( x1**2 + y1**2 + z1**2 )
                    #
                    n1 = numpy.array([x1/length,y1/length,z1/length])
                    #
                    # now that we have a first direction, n1, we
                    # can use the "get_line_normals" function to
                    # compute two linearly ind. orthogonal dir. to
                    # n1
                    #
                    two_atoms = [['F',x0,y0,z0],[atom_neighbours[0][0], \
                                                 atom_neighbours[0][1], \
                                                 atom_neighbours[0][2], \
                                                 atom_neighbours[0][3]]]
                    directions = get_line_normals(two_atoms)
                    n2 = numpy.array([0.0,0.0,0.0])
                    n3 = numpy.array([0.0,0.0,0.0])
                    n2[0] = float(directions[0])
                    n2[1] = float(directions[1])
                    n2[2] = float(directions[2])
                    n3[0] = float(directions[3])
                    n3[1] = float(directions[4])
                    n3[2] = float(directions[5])
                    #
                    # in the basis n1, n2, n3 we find the atomic
                    # positions of the 3 lone pairs of Fluorine
                    #
                    angle1 = 1.9094 # 109.4 degres from a regular tetrahedron
                    angle2 = 2.0944 # 2pi/3 rotation angle about n1
                    #
                    first_pair  = length/2.1*(cos(angle1)*n1 + sin(angle1)*n2)
                    second_pair = length/2.1*(cos(angle1)*n1 + sin(angle1)*cos(angle2)*n2 + sin(angle1)*sin(angle2)*n3)
                    third_pair  = length/2.1*(cos(angle1)*n1 + sin(angle1)*cos(2*angle2)*n2 + sin(angle1)*sin(2*angle2)*n3)
                    #
                    projections.append("c=%11.7f,%11.7f,%11.7f:s \n" %(first_pair[0]+x0,first_pair[1]+y0,first_pair[2]+z0))
                    projections.append("c=%11.7f,%11.7f,%11.7f:s \n" %(second_pair[0]+x0,second_pair[1]+y0,second_pair[2]+z0))
                    projections.append("c=%11.7f,%11.7f,%11.7f:s \n" %(third_pair[0]+x0,third_pair[1]+y0,third_pair[2]+z0))
        #
        ##########################
        # projections for Oxygen #
        ##########################
        #
        if nature=='O':
            #
            if number_of_neighbours > 2 or number_of_neighbours < 1:
                #
                # Error message corresponding to a weird number
                # of nearest neighbours
                #
                msg = " An Oxygen atom has been found with %i neighbours ! \n" %(number_of_neighbours,) + \
                      " Expecting only 1 or 2 neighbours for Oxygen...     \n" +                          \
                      " Coordinates : %9.6f %9.6f %9.6f                    \n" %(atoms[atom][1],          \
                                                                                 atoms[atom][2],          \
                                                                                 atoms[atom][3])
                raise ProjectionError, msg
            else:
                #
                if number_of_neighbours==1:
                    #
                    #  1 s-type projection on mid-bond corrected by 
                    # electronegativity
                    #
                    for j in xrange(len(atom_neighbours)):
                        #
                        if atom_neighbours[j][-2]=="supercell":
                            couple = [atom,atom_neighbours[j][-1]]
                            possible_bond = True
                            if couple in forbidden_supercell_links:
                                possible_bond = False
                        elif atom_neighbours[j][-2]=="periodic":
                            couple = [atom,atom_neighbours[j][-1]]
                            possible_bond = True
                            if couple in forbidden_periodic_links:
                                possible_bond = False
                        #
                        nature1 = atom_neighbours[j][0]
                        x1 = atom_neighbours[j][1]
                        y1 = atom_neighbours[j][2]
                        z1 = atom_neighbours[j][3]
                        #
                        if possible_bond:
                            #
                            local_projection = single_bond(nature,x0,y0,z0,nature1,x1,y1,z1)
                            projections.append(local_projection)
                    #
                    # now we need a p_z orthogonal to the plane of the
                    # oxygen, the nearest neighbour to oxygen and a 
                    # neighbour to the nearest neighbour! Plus 2 lone
                    # pairs in plane
                    #
                    three_atoms = []
                    # the oxygen atom itself
                    three_atoms.append(['O',x0,y0,z0])
                    # the only nearest neighbour to it
                    three_atoms.append([atom_neighbours[0][0],atom_neighbours[0][1], \
                                        atom_neighbours[0][2],atom_neighbours[0][3]])
                    # the index of the nearest neighbour
                    neighbour_number = atom_neighbours[0][-1]
                    # now we find a nearest neighbour to the nearest neighbour
                    # which is NOT the oxygen atom!
                    if atom_neighbours[0][-2]=="supercell":
                        #
                        if len(neighbours[neighbour_number]["supercell_neighbours"]) + \
                           len(neighbours[neighbour_number]["periodic_neighbours"]) > 1:
                            #
                            if len(neighbours[neighbour_number]["periodic_neighbours"]) > 0:
                                #
                                second_neighbour = neighbours[neighbour_number]["periodic_neighbours"][0]
                                three_atoms.append([periodic_atoms[second_neighbour][0],periodic_atoms[second_neighbour][1], \
                                                    periodic_atoms[second_neighbour][2],periodic_atoms[second_neighbour][3]])
                            else:
                                #
                                iterator = 0
                                #
                                while neighbours[neighbour_number]["supercell_neighbours"][iterator]==atom:
                                    #
                                    iterator += 1
                                #
                                second_neighbour = neighbours[neighbour_number]["supercell_neighbours"][iterator]
                                three_atoms.append([atoms[second_neighbour][0],atoms[second_neighbour][1], \
                                                    atoms[second_neighbour][2],atoms[second_neighbour][3]])
                    #
                    normal = [0.0,0.0,0.0]
                    if len(three_atoms)==3:
                        #
                        direction = get_plane_normal(three_atoms)
                        normal = [direction[0],direction[1],direction[2]]
                    #
                    else:
                        #
                        directions = get_line_normals(three_atoms)
                        normal = [directions[0],directions[1],directions[2]]
                    #
                    projections.append("c=%11.7f,%11.7f,%11.7f:pz:z=%3.2f,%3.2f,%3.2f \n" %(x0,y0,z0,  \
                                                                                         normal[0], \
                                                                                         normal[1], \
                                                                                         normal[2]))
                    #
                    # now the two remaining lone pairs
                    #
                    x1 = atom_neighbours[0][1]
                    y1 = atom_neighbours[0][2]
                    z1 = atom_neighbours[0][3]
                    #
                    u0 = x0 - x1
                    v0 = y0 - y1
                    w0 = z0 - z1
                    #
                    length = sqrt( u0**2 + v0**2 + w0**2 )
                    #
                    u0 = u0 / length
                    v0 = v0 / length
                    w0 = w0 / length
                    #
                    theta = 109.4*3.1415926535/180
                    #
                    u1 = v0*normal[2]-w0*normal[1]
                    v1 = w0*normal[0]-u0*normal[2]
                    w1 = u0*normal[1]-v0*normal[0]
                    #
                    lone_pair_1 = [0.6*(cos(theta/2.0)*u0+sin(theta/2.0)*u1), \
                                   0.6*(cos(theta/2.0)*v0+sin(theta/2.0)*v1), \
                                   0.6*(cos(theta/2.0)*w0+sin(theta/2.0)*w1)]
                    lone_pair_2 = [0.6*(cos(theta/2.0)*u0-sin(theta/2.0)*u1), \
                                   0.6*(cos(theta/2.0)*v0-sin(theta/2.0)*v1), \
                                   0.6*(cos(theta/2.0)*w0-sin(theta/2.0)*w1)]
                    projections.append("c=%11.7f,%11.7f,%11.7f:s \n" %(x0+lone_pair_1[0],y0+lone_pair_1[1],z0+lone_pair_1[2]))
                    projections.append("c=%11.7f,%11.7f,%11.7f:s \n" %(x0+lone_pair_2[0],y0+lone_pair_2[1],z0+lone_pair_2[2]))
                #
                elif number_of_neighbours==2:
                    #
                    #  2 s-type projections on mid-bonds corrected by 
                    # electronegativity
                    #
                    for j in xrange(len(atom_neighbours)):
                        #
                        if atom_neighbours[j][-2]=="supercell":
                            couple = [atom,atom_neighbours[j][-1]]
                            possible_bond = True
                            if couple in forbidden_supercell_links:
                                possible_bond = False
                        elif atom_neighbours[j][-2]=="periodic":
                            couple = [atom,atom_neighbours[j][-1]]
                            possible_bond = True
                            if couple in forbidden_periodic_links:
                                possible_bond = False
                        #
                        nature1 = atom_neighbours[j][0]
                        x1 = atom_neighbours[j][1]
                        y1 = atom_neighbours[j][2]
                        z1 = atom_neighbours[j][3]
                        #
                        if possible_bond:
                            #
                            local_projection = single_bond(nature,x0,y0,z0,nature1,x1,y1,z1)
                            projections.append(local_projection)
                    #
                    # now 2 lone pairs on both sides of a bisecting plane
                    #
                    atom1_x = atom_neighbours[0][1]
                    atom1_y = atom_neighbours[0][2]
                    atom1_z = atom_neighbours[0][3]
                    atom2_x = atom_neighbours[1][1]
                    atom2_y = atom_neighbours[1][2]
                    atom2_z = atom_neighbours[1][3]
                    length  = sqrt( (atom2_x-atom1_x)**2 + (atom2_y-atom1_y)**2 + (atom2_z-atom1_z)**2 )
                    m_x = (atom2_x-atom1_x) / length
                    m_y = (atom2_y-atom1_y) / length
                    m_z = (atom2_z-atom1_z) / length
                    
                    direction = get_plane_normal([['O',x0,y0,z0], \
                                                  ['X',atom1_x,atom1_y,atom1_z], \
                                                  ['X',atom2_x,atom2_y,atom2_z]])

                    n_x = direction[0]
                    n_y = direction[1]
                    n_z = direction[2]
                    nn = numpy.array([n_x,n_y,n_z])
                    p_x = m_y*n_z-m_z*n_y
                    p_y = m_z*n_x-m_x*n_z
                    p_z = m_x*n_y-m_y*n_x
                    pp = numpy.array([p_x,p_y,p_z])
                    angle = 109.4*3.1415926535/(2.0*180)
                    # length of 0.6 Angstrom from the oxygen
                    lone_pair_1 = numpy.array([x0,y0,z0])+ 0.6*(sin(angle)*nn-cos(angle)*pp)
                    lone_pair_2 = numpy.array([x0,y0,z0])+ 0.6*(-sin(angle)*nn-cos(angle)*pp)
                    projections.append("c=%11.7f,%11.7f,%11.7f:s \n" %(lone_pair_1[0],lone_pair_1[1],lone_pair_1[2]))
                    projections.append("c=%11.7f,%11.7f,%11.7f:s \n" %(lone_pair_2[0],lone_pair_2[1],lone_pair_2[2]))
        #
        ############################
        # projections for Nitrogen #
        ############################
        #
        if nature=='N':
            #
            if number_of_neighbours > 3 or number_of_neighbours < 1:
                #
                # Error message corresponding to a weird number
                # of nearest neighbours
                #
                msg = " A Nitrogen atom has been found with %i neighbours !  \n" %(number_of_neighbours,) + \
                      " Expecting only 1, 2 or 3 neighbours for Nitrogen...  \n" +                          \
                      " Coordinates : %9.6f %9.6f %9.6f                      \n" %(atoms[atom][1],          \
                                                                                   atoms[atom][2],          \
                                                                                   atoms[atom][3])
                raise ProjectionError, msg
            else:
                #
                if number_of_neighbours==3:
                    #
                    #  3 s-type projections on mid-bonds corrected by 
                    # electronegativity
                    #
                    for j in xrange(len(atom_neighbours)):
                        #
                        if atom_neighbours[j][-2]=="supercell":
                            couple = [atom,atom_neighbours[j][-1]]
                            possible_bond = True
                            if couple in forbidden_supercell_links:
                                possible_bond = False
                        elif atom_neighbours[j][-2]=="periodic":
                            couple = [atom,atom_neighbours[j][-1]]
                            possible_bond = True
                            if couple in forbidden_periodic_links:
                                possible_bond = False
                        #
                        nature1 = atom_neighbours[j][0]
                        x1 = atom_neighbours[j][1]
                        y1 = atom_neighbours[j][2]
                        z1 = atom_neighbours[j][3]
                        #
                        if possible_bond:
                            #
                            local_projection = single_bond(nature,x0,y0,z0,nature1,x1,y1,z1)
                            projections.append(local_projection)
                    #
                    # now a lone pair "atop" the plane defined by the three neighbours
                    #
                    direction = get_plane_normal(atom_neighbours)
                    #
                    atom1_x = x0-atom_neighbours[0][1]
                    atom1_y = y0-atom_neighbours[0][2]
                    atom1_z = z0-atom_neighbours[0][3]
                    #
                    alpha = atom1_x*direction[0]+atom1_y*direction[1]+atom1_z*direction[2]
                    #
                    if alpha > 0:
                        #
                        lone_pair = [x0+0.65*direction[0],y0+0.65*direction[1],z0+0.65*direction[2]]
                    else:
                        #
                        lone_pair = [x0-0.65*direction[0],y0-0.65*direction[1],z0-0.65*direction[2]]
                    #
                    projections.append("c=%11.7f,%11.7f,%11.7f:s \n" %(lone_pair[0], \
                                                                    lone_pair[1], \
                                                                    lone_pair[2]))
                elif number_of_neighbours==1:
                    #
                    #  1 s-type projection on mid-bond corrected by 
                    # electronegativity
                    #
                    for j in xrange(len(atom_neighbours)):
                        #
                        if atom_neighbours[j][-2]=="supercell":
                            couple = [atom,atom_neighbours[j][-1]]
                            possible_bond = True
                            if couple in forbidden_supercell_links:
                                possible_bond = False
                        elif atom_neighbours[j][-2]=="periodic":
                            couple = [atom,atom_neighbours[j][-1]]
                            possible_bond = True
                            if couple in forbidden_periodic_links:
                                possible_bond = False
                        #
                        nature1 = atom_neighbours[j][0]
                        x1 = atom_neighbours[j][1]
                        y1 = atom_neighbours[j][2]
                        z1 = atom_neighbours[j][3]
                        #
                        if possible_bond:
                            #
                            local_projection = single_bond(nature,x0,y0,z0,nature1,x1,y1,z1)
                            projections.append(local_projection)
                    #
                    # 2 p orbitals orthogonal to the bond between nitrogen
                    # and its nearest neighbour plus one s opposite to the 
                    # bond
                    #
                    n_x = atom_neighbours[0][1]-x0
                    n_y = atom_neighbours[0][2]-y0
                    n_z = atom_neighbours[0][3]-z0
                    #
                    length = sqrt( n_x**2 + n_y**2 + n_z**2 )
                    #
                    n_x = n_x / length
                    n_y = n_y / length
                    n_z = n_z / length
                    #
                    projections.append("c=%11.7f,%11.7f,%11.7f:s \n" %(x0-0.65*n_x,y0-0.65*n_y,z0-0.65*n_z))
                    #
                    two_atoms = [['N',x0,y0,z0],['X',atom_neighbours[0][1],atom_neighbours[0][2],atom_neighbours[0][3]]]
                    directions = get_line_normals(two_atoms)
                    #
                    projections.append("c=%11.7f,%11.7f,%11.7f:pz:z=%3.2f,%3.2f,%3.2f \n" %(x0,y0,z0,      \
                                                                                         directions[0], \
                                                                                         directions[1], \
                                                                                         directions[2]))
                    projections.append("c=%11.7f,%11.7f,%11.7f:pz:z=%3.2f,%3.2f,%3.2f \n" %(x0,y0,z0,      \
                                                                                         directions[3], \
                                                                                         directions[4], \
                                                                                         directions[5]))
                #
                elif number_of_neighbours==2:
                    #
                    #  2 s-type projections on mid-bonds corrected by 
                    # electronegativity
                    #
                    for j in xrange(len(atom_neighbours)):
                        #
                        if atom_neighbours[j][-2]=="supercell":
                            couple = [atom,atom_neighbours[j][-1]]
                            possible_bond = True
                            if couple in forbidden_supercell_links:
                                possible_bond = False
                        elif atom_neighbours[j][-2]=="periodic":
                            couple = [atom,atom_neighbours[j][-1]]
                            possible_bond = True
                            if couple in forbidden_periodic_links:
                                possible_bond = False
                        #
                        nature1 = atom_neighbours[j][0]
                        x1 = atom_neighbours[j][1]
                        y1 = atom_neighbours[j][2]
                        z1 = atom_neighbours[j][3]
                        #
                        if possible_bond:
                            #
                            local_projection = single_bond(nature,x0,y0,z0,nature1,x1,y1,z1)
                            projections.append(local_projection)
                    #
                    # now we need a pz and a lone pair
                    #
                    if atom_neighbours[0][0]=='H':
                        #
                        good = 1
                    elif atom_neighbours[1][0]=='H':
                        #
                        good = 0
                    else:
                        #
                        l1 = sqrt((atom_neighbours[0][1]-x0)**2+(atom_neighbours[0][2]-y0)**2+ \
                                  (atom_neighbours[0][3]-z0)**2)
                        l2 = sqrt((atom_neighbours[1][1]-x0)**2+(atom_neighbours[1][2]-y0)**2+ \
                                  (atom_neighbours[1][3]-z0)**2)
                        if l1 > l2:
                            #
                            good = 1
                        else:
                            #
                            good = 0
                    #
                    three_atoms = []
                    # the nitrogen atom itself
                    three_atoms.append(['N',x0,y0,z0])
                    # the only nearest neighbour that is closest
                    three_atoms.append([atom_neighbours[good][0],atom_neighbours[good][1], \
                                        atom_neighbours[good][2],atom_neighbours[good][3]])
                    # the index of the nearest neighbour
                    neighbour_number = atom_neighbours[good][-1]
                    # now we find a nearest neighbour to the nearest neighbour
                    # which is NOT the nitrogen atom!
                    if atom_neighbours[good][-2]=="supercell":
                        #
                        if len(neighbours[neighbour_number]["supercell_neighbours"]) + \
                           len(neighbours[neighbour_number]["periodic_neighbours"]) > 1:
                            #
                            if len(neighbours[neighbour_number]["periodic_neighbours"]) > 0:
                                #
                                second_neighbour = neighbours[neighbour_number]["periodic_neighbours"][0]
                                three_atoms.append([periodic_atoms[second_neighbour][0],periodic_atoms[second_neighbour][1], \
                                                    periodic_atoms[second_neighbour][2],periodic_atoms[second_neighbour][3]])
                            else:
                                #
                                iterator = 0
                                #
                                while neighbours[neighbour_number]["supercell_neighbours"][iterator]==atom:
                                    #
                                    iterator += 1
                                #
                                second_neighbour = neighbours[neighbour_number]["supercell_neighbours"][iterator]
                                three_atoms.append([atoms[second_neighbour][0],atoms[second_neighbour][1], \
                                                    atoms[second_neighbour][2],atoms[second_neighbour][3]])
                    #
                    normal = [0.0,0.0,0.0]
                    if len(three_atoms)==3:
                        #
                        direction = get_plane_normal(three_atoms)
                        normal = [direction[0],direction[1],direction[2]]
                    #
                    else:
                        #
                        directions = get_line_normals(three_atoms)
                        normal = [directions[0],directions[1],directions[2]]
                    #
                    projections.append("c=%11.7f,%11.7f,%11.7f:pz:z=%3.2f,%3.2f,%3.2f \n" %(x0,y0,z0,  \
                                                                                         normal[0], \
                                                                                         normal[1], \
                                                                                         normal[2]))                    
                    #
                    # the lone pair
                    #
                    atom1_x = atom_neighbours[0][1]
                    atom1_y = atom_neighbours[0][2]
                    atom1_z = atom_neighbours[0][3]
                    atom2_x = atom_neighbours[1][1]
                    atom2_y = atom_neighbours[1][2]
                    atom2_z = atom_neighbours[1][3]
                    length  = sqrt( (atom2_x-atom1_x)**2 + (atom2_y-atom1_y)**2 + (atom2_z-atom1_z)**2 )
                    m_x = (atom2_x-atom1_x) / length
                    m_y = (atom2_y-atom1_y) / length
                    m_z = (atom2_z-atom1_z) / length
                    
                    direction = get_plane_normal([['N',x0,y0,z0], \
                                                  ['X',atom1_x,atom1_y,atom1_z], \
                                                  ['X',atom2_x,atom2_y,atom2_z]])

                    n_x = direction[0]
                    n_y = direction[1]
                    n_z = direction[2]
                    nn = numpy.array([n_x,n_y,n_z])
                    p_x = m_y*n_z-m_z*n_y
                    p_y = m_z*n_x-m_x*n_z
                    p_z = m_x*n_y-m_y*n_x
                    pp = numpy.array([p_x,p_y,p_z])
                    projections.append("c=%11.7f,%11.7f,%11.7f:s \n" %(x0-0.65*pp[0],y0-0.65*pp[1],z0-0.65*pp[2]))

    return projections

def get_line_normals(atom_neighbours):
    #
    # get two linearly independent orthogonal directions to a line
    #
    #
    # vector linking atom 0 to atom 1
    #
    x0 = atom_neighbours[1][1] - atom_neighbours[0][1]
    y0 = atom_neighbours[1][2] - atom_neighbours[0][2]
    z0 = atom_neighbours[1][3] - atom_neighbours[0][3]
    #
    length = sqrt( x0**2 + y0**2 + z0**2 )
    #
    x0 = x0 / length
    y0 = y0 / length
    z0 = z0 / length
    #
    abs_vector = [abs(x0),abs(y0),abs(z0)]
    #
    if abs_vector.index(max(abs_vector))==0:
        #
        first_normal = [-(y0+z0)/x0,1.0,1.0]
    elif abs_vector.index(max(abs_vector))==1:
        #
        first_normal = [1.0,-(x0+z0)/y0,1.0]
    elif abs_vector.index(max(abs_vector))==2:
        #
        first_normal = [1.0,1.0,-(x0+y0)/z0]
    #
    # normalize first_normal
    #
    length = sqrt( first_normal[0]**2 + first_normal[1]**2 + first_normal[2]**2 )
    #
    first_normal[0] = first_normal[0] / length
    first_normal[1] = first_normal[1] / length
    first_normal[2] = first_normal[2] / length
    #
    # second normal vector
    #
    second_normal = [y0*first_normal[2]-z0*first_normal[1], \
                     z0*first_normal[0]-x0*first_normal[2], \
                     x0*first_normal[1]-y0*first_normal[0]]
    #
    directions = first_normal + second_normal
    #
    return directions

def get_plane_normal(atom_neighbours):
    #
    # get the normal to the plane defined by
    # the 3 atoms in atom_neighbours
    #
    # vector linking atom 0 to atom 1
    #
    x0 = atom_neighbours[1][1] - atom_neighbours[0][1]
    y0 = atom_neighbours[1][2] - atom_neighbours[0][2]
    z0 = atom_neighbours[1][3] - atom_neighbours[0][3]
    #
    # vector linking atom 0 to atom 2
    #
    x1 = atom_neighbours[2][1] - atom_neighbours[0][1]
    y1 = atom_neighbours[2][2] - atom_neighbours[0][2]
    z1 = atom_neighbours[2][3] - atom_neighbours[0][3]
    #
    # computing the vectorial product of
    # the above two vectors
    #
    x2 = y0*z1 - z0*y1
    y2 = z0*x1 - x0*z1
    z2 = x0*y1 - y0*x1
    #
    # normalizing the vector
    #
    length = sqrt( x2**2 + y2**2 + z2**2 )
    #
    x2 = x2 / length
    y2 = y2 / length
    z2 = z2 / length
    #
    return [x2,y2,z2]

def single_bond(nature,x0,y0,z0,nature1,x1,y1,z1):
    #
    # returns an s-type projection between 2 atoms
    # at a distance from the atoms that takes into
    # account their electronegativity
    #
    alpha = electronegativity(nature) / ( electronegativity(nature) + electronegativity(nature1) )
    #
    x = alpha*x0 + (1.0-alpha)*x1
    y = alpha*y0 + (1.0-alpha)*y1
    z = alpha*z0 + (1.0-alpha)*z1
    #
    projection = "c=%11.7f,%11.7f,%11.7f:s \n" %(x,y,z)
    #
    return projection

class ProjectionError(Exception):
    """Exception Raised when a Projection problem is encountered"""
    def __init__(self, message):
        self.message = message

