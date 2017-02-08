.. _credits:

============================
Credits and Acknowledgements
============================

CVX was designed by Michael Grant and Stephen Boyd, with input from Yinyu Ye; and was
implemented by Michael Grant [GBY06]_. It incorporates ideas from earlier works
by Löfberg [Löf04]_, Dahl and [DV04]_, Wu and Boyd [WB00]_,
and many others. The modeling language follows the spirit of `AMPL
<http://www.ampl.com>`_ or `GAMS <http://www.gams.com>`_; unlike these
packages, however, CVX was designed from the beginning to fully exploit
convexity. The specific method for implementing CVX in Matlab draws
heavily from `YALMIP <http://users.isy.liu.se/johanl/yalmip>`_. 

We wish to thank the following people for their contributions:
Toh Kim Chuan, Laurent El Ghaoui, Arpita Ghosh,
Siddharth Joshi, Johan Löberg, Almir Mutapcic, Michael Overton and his
students, Art Owen, Rahul Panicker, Imre Polik, Joëlle Skaf, Lieven
Vandenberghe, Argyris Zymnis. We are also grateful to the many students
in several universities who have (perhaps unwittingly) served as beta
testers by using CVX in their classwork. We thank Igal Sason for
catching many typos in an earlier version of this document, and
generally helping us to improve its clarity. 

We would like to thank
`Gurobi Optimization <http://gurobi.com>`_ and `MOSEK ApS <http://mosek.com>`_
for their generous assistance as we developed the interfaces
to their commercial products.

.. raw:: html

	<h2>References</h2>

.. [AG00] F. Alizadeh and D. Goldfarb.
	Second-order cone programming.
	*Mathematical Programming, Series B*, 95:3-51, 2001.
	http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.23.5133

.. [BKVH05] S. Boyd, S. J. Kim, L. Vandenberghe, and A. Hassibi,.
	A tutorial on geometric programming.
	*Optimization and Engineering*, 8(1):67-127, 2007.
	http://stanford.edu/~boyd/papers/gp_tutorial.html

.. [BV04] S. Boyd and L. Vandenberghe.
 	*Convex Optimization*.
 	Cambridge University Press, 2004.
	http://stanford.edu/~boyd/cvxbook.html
	
.. [Cru02] C. Crusius.
	*A Parser/Solver for Convex Optimization Problems*. 
	Ph.D. thesis, Information Systems Laboratory, 
	Department of Electrical Engineering, Stanford University, 2002.

.. [DV04] J. Dahl and L. Vandenberghe,
	CVXOPT: A Python package for convex optimization (version 1.1.5).		
	http://abel.ee.ucla.edu/cvxopt/

.. [GBY06] M. Grant and S. Boyd and Y. Ye.
	Disciplined convex programming.
	In *Global Optimization: from Theory to Implementation*, 
	Nonconvex Optimization and Its Applications,
	L. Liberti and N. Maculan, *eds.*, Springer, 2006.
	http://stanford.edu/~boyd/disc_cvx_prog.html
		
.. [Gra04] M. Grant.
	*Disciplined Convex Programming*.
	Ph.D. thesis, Information Systems Laboratory, 
	Department of Electrical Engineering, Stanford University, 2004.
	http://stanford.edu/~boyd/disc_cvx_prog.html
	
.. [Löf04] J. Löfberg.
 	YALMIP: a toolbox for modeling and optimization in MATLAB.
	*Proceedings of the 2004 International Symposium on Computer Aided Control Systems Design*,
	IEEE Press, September 2004, pp. 284-289.
	http://users.isy.liu.se/johanl/yalmip/
	
.. [Owen06] A. Owen.
	A robust hybrid of lasso and ridge regression.
	Technical report, Department of Statistics, Stanford University, October 2006.
	http://www-stat.stanford.edu/~owen/reports/hhu.pdf
	
.. [Stu99] J.F. Sturm, 
	Using SeDuMi 1.02, a MATLAB toolbox for optimization over symmetric cones.
	*Optimization Methods and Software*, 11-12:625-633, 1999.
	Special issue on Interior Point Methods (CD supplement with software).
	http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.49.6954

.. [TTT03] R.H. Tütüncü, K.C. Toh, and M.J. Todd.
	Solving semidefinite-quadratic-linear programs using SDPT3.
	*Mathematical Programming, Series B*, 95:189-217, 2003.
	http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.127.4807

.. [WB00] S.P. Wu and S. Boyd.
	SDPSOL: A parser/solver for semidefinite programs with matrix structure.
	In *Recent Advances in LMI Methods for Control*, 
	L. El Ghaoui and S.I. Niculescu, *eds.*, SIAM, pp. 79-91, 2000.
	http://www.stanford.edu/~boyd/sdpsol.html

