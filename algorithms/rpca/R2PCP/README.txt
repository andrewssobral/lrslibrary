      Riemannian Robust Principle Component Pursuit (R2PCP)
       http://www.uni-graz.at/imawww/ifb/r2pcp/index.html


________________________________________________________________
OVERVIEW

Robust principle component pursuit (RPCP) aims at decomposing a given data 
matrix additively into a low-rank matrix and a sparse matrix. R2PCP
solves such a matrix decomposition problem by formulating a least-
squares problem subject to rank- and cardinality constraints. The 
resulting constrained optimization problem is tackled by an (inexact)
alternating minimization scheme on the respective matrix manifolds.
In particular, the low-rank matrix subproblem is resolved by a
Riemannian optimization step. In addition, based on particular 
applications, a heuristic trimming procedure is incorporated to adjust 
the rank and cardinality dynamically.

The R2PCP package was implemented in MATLAB. 


_________________________________________________________________
REFERENCE

Michael Hintermueller and Tao Wu, "Robust Principal Component Pursuit 
via Inexact Alternating Minimization on Matrix Manifolds", Journal 
of Mathematical Imaging and Vision, to appear. 
DOI:10.1007/s10851-014-0527-y


_________________________________________________________________
DIRECTORY STRUCTURE

demo_synthetic_data.m
	Demonstration of R2PCP on synthetic data.
demo_airport_video.m
	Demonstration of R2PCP for background-foreground separation 
	on the "airport" video.
demo_lobby_video.m
	Demonstration of R2PCP for background-foreground separation 
	on the "lobby" video.
README.txt
        This file.
Riem_RPCP (dir)
        Source codes of R2PCP.
video_lobby (dir)
	Raw data for the "airport" video.
video_lobby (dir)
	Raw data for the "lobby" video.
utilities (dir)
	Utility functions used in demonstrations.


_________________________________________________________________
ACKNOWLEDGEMENT

The R2PCP package was developed by Michael Hintermueller from 
Department of Mathematics at the Humboldt-University of Berlin and 
Tao Wu from Institute for Mathematics and Scientific Computing at 
the University of Graz. This work is supported by the Austria 
Science Fund (FWF) under START-program Y305 "Interfaces and Free 
Boundaries" and the SFB F32 "Mathematical Optimization and 
Applications in Biomedical Sciences".


_________________________________________________________________
DISCLAIMER

The R2PCP package (including code modifications) may only be used 
for NON-COMMERCIAL RESEARCH purposes. For inquiries concerning a 
different use, please contact Prof. Michael Hintermueller from the 
Humboldt-University of Berlin (hint"at"math.hu-berlin.de).

Your comments are welcome. Please keep track of bugs or missing/
confusing instructions and report them to

Michael Hintermueller                 <hint"at"math.hu-berlin.de>
Tao Wu                                <tao.wu"at"uni-graz.at>

The algorithms contained in the R2PCP package were implemented by
Tao Wu.
