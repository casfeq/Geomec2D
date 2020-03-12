# Geomec2D

This program implements a Finite Volume Method for discretization and solution of two-dimensional consolidation problems as part of my master's thesis entitled "Analysis of Numerical Schemes in Collocated and Staggered Grids for Problems of Poroelasticity". It solves the problems presented and solved by Terzaghi [3] and Mandel [2]. The governing equations are discretized within the FVM and the resulting linear system of equations is solved with a LU Factorization found in PETSc suite [1].

The GUI uses Zenity package (which is by default installed in Ubuntu). PETsc is also required. To run this program, first edit the "CMakeLists.txt" file with the location of PETsc installation. Then call the "geomec.sh" file.

Written by FERREIRA, C. A. S.

PETSc:
https://www.mcs.anl.gov/petsc/

Zenity:
https://help.gnome.org/users/zenity/

[1] BALAY et al. PETSc User Manual. Technical Report, Argonne National Laboratory, 2017.
[2] MANDEL, J. Consolidation Des Sols (Étude Mathématique). Géotechnique, v. 3, n. 7, pp. 287-
299, 1953.
[3] TERZAGHI, K. Erdbaumechanik auf Bodenphysikalischer Grundlage. Franz Deuticke, Leipzig,
1925.
