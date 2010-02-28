Introduction :
==============

This folder contains some Python scripts for making life easier
when using the Wannier90 code (http://www.wannier.org/).

Content of the directory :
==========================

- compute_mean_qc.py, this program computes the average (i.e. the
mean) quantum conductance from a set of quantum conductance files
given by Wannier90. The  energy  interval needs  to be indentical
for all the files.

- visualize_whole_hamiltonian.py, which allows  one  to read some
Hamiltonian matrices produced by the Wannier90 and combines  them
into one  big Hamiltonian matrix  in order to "visualize" it. The
matrices are given by the names :
 -> <seedname>_htL.dat
 -> <seedname>_htLC.dat
 -> <seedname>_htC.dat
 -> <seedname>_htCR.dat
 -> <seedname>_htR.dat
  where "seedname" is the root of the *.win file (the * part).

- initial_projections_generator, which is  a directory containing
some Python routines  used  to generate  a Wannier90 master input
file with some pre-computed quasi-optimal initial projections.
This utility should be extremely useful since it basically allows
one  to alleviate  the burden  of finding initial projections for
the wannierisation.
However, only a handful of elements are supported by the program.

Dependencies :
==============

- compute_mean_qc.py, only needs a working version of Numpy.

- visualize_whole_hamiltonian.py, needs Numpy and Matplotlib.

- initial_projections_generator, only needs Numpy.

