# GeomecFVLib

This code implements a Finite Volume Method for discretization and solution of two-dimensional consolidation problems as part of a master's thesis entitled "Analysis Numerical Schemes in Collocated and Staggered Grids for Poroelasticity Problems". It solves the problems presented and solved by Terzaghi [3] and Mandel [2]. The governing equations are discretized within the FVM and the resulting linear system of equations is solved with a LU Factorization found in PETSc [1].

Written by FERREIRA, C. A. S.

Florianópolis, 2019.

[1] BALAY et al. PETSc User Manual. Technical Report, Argonne National Laboratory, 2017.
[2] MANDEL, J. Consolidation Des Sols (Étude Mathématique). Géotechnique, v. 3, n. 7, pp. 287-
299, 1953.
[3] TERZAGHI, K. Erdbaumechanik auf Bodenphysikalischer Grundlage. Franz Deuticke, Leipzig,
1925.
