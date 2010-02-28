This set of programs has been designed to enable an automatic
computation of the Wannier initial guesses. It  is taking 1D,
2D  or  3D  periodicity into account fully and should produce
high-quality projections.

For this to happen, the code  is restricted to  a subclass of
elements  for  which  the chemical environment  can be easily
obtained (mostly elements entering the composition of organic
compounds).

The code takes  a file  as input. This file should contain the
atomic positions in a CARTESIAN set of axes with the coordina-
tes expressed in ANGSTROM.

Then the code will  ask a few questions to the user concerning
the periodicity of the system.

 - namely : is the system 1D, 2D or 3D periodic
 - then, once the dimensionality has been entered
   the code asks for the axes of periodicity. Those
   should be either 'x', 'y' or 'z'
 - following this the code finally asks for the
   length of the period in each direction (for example
   5 Angstroms in 'x')

!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!
!!!! When the code extracts  the atomic positions !!!!
!!!! it also makes  some checks  on the nature of !!!!
!!!! the elements to see  if there are some  non- !!!!
!!!! supported elements. If there is one then the !!!!
!!!! code STOPS                                   !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

HOW TO EXTEND THE CAPABILITIES OF THE CODE ?
============================================

Basically the only thing a user might want, would be to extend
the number of supported elements.
To  do just  that follow the steps below (the  user  will ONLY
have to focus on the module "chemistry.py")

1) add the symbol of the element that you want
   to the "supported_elements" list at the beginning
   of the "chemistry.py" module.

2) then you will need to change the "critical_length"
   function by adding the values for the cutoff distances
   between this new element and ALL the previous ones.
   Keep the structure of the function intact.

3) add the electronegativity of the element to the
   "electronegativity" function

4) Last but not least, you will need to implement all
   the possible chemical environments your new atom
   might encounter in a structure. Go in the
   "find_projections" function and look at the implemented
   projections for the other supported elements. Then
   implement the projections for your element.

