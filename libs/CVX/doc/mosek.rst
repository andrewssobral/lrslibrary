.. _mosek:

====================
Using MOSEK with CVX
====================

About MOSEK
-----------

`MOSEK ApS <http://www.mosek.com>`_ is widely considered the leader in commercial software
for nonlinear convex optimization. The company is led by CEO 
`Erling Andersen <http://www.linkedin.com/in/edandersen>`_, and its  board is chaired by 
Stanford Professor `Yinyu Ye <http://www.stanford.edu/~yyye/>`_. Both are internationally
recognized for their contributions  to the field of convex optimization, and remain active 
in research and publication. With its existing support for integer variables and the 
addition of semidefinite  programming capability in version 7.0, the MOSEK solver can 
address a wider variety of CVX models than any other solver.

Using CVX with MOSEK requires a CVX Professional license. Please visit 
`Licensing`_ to learn more about licensing
options. Academic users can obtain a free CVX Professional license by 
submitting an `Academic License Request`_.

.. note::

	If you intend to use CVX with *both* MOSEK and Gurobi, please follow the 
	directions on the page :ref:`gurobi`.

Using the bundled version of MOSEK
----------------------------------

.. note:: 

    The bundled version of MOSEK can only be used within CVX. If you wish to use MOSEK
    outside of CVX as well, you will need a standalone MOSEK installation.

The simplest way to use MOSEK with CVX is by installing the appropriate 
CVX Professional Solver Bundle, available `here <http://cvxr.com/cvx/download>`_, 
along with a MOSEK-enabled CVX Professional license. Please see :ref:`install`
and :ref:`licinstall`
for general installation instructions. Once the CVX Professional license has been
properly installed, MOSEK will be enabled.

Using CVX with separate MOSEK installation
------------------------------------------

If you wish to use CVX with a separate installation of MOSEK 6.0 or 7.0,
follow these steps after you have successfully installed MOSEK:

* Make sure that MATLAB can locate your current installation of MOSEK. 
  If you have already been using the ``mosekopt`` command within MATLAB, there
  is no further configuration needed. Otherwise, you will need to modify your MATLAB
  search path so it can find your MOSEK installation. For information, please see 
  the relevant page in your MOSEK documentation:
  
  - `MOSEK 7: Installation <http://docs.mosek.com/7.0/toolbox/Installation.html>`_
  - `MOSEK 6: Insatllation <http://docs.mosek.com/6.0/toolbox/node006.html>`_

* If you have not done so yet, download and install CVX and a CVX Professional license
  according to the instructions in :ref:`install` and :ref:`licinstall`.
  
* If CVX and your license had already been installed, simply re-run ``cvx_setup`` so
  that CVX can rebuild its solver list and include MOSEK.
  
If successful, the output of ``cvx_setup`` should show that MOSEK is among the list
of available solvers. If it fails to find MOSEK, it will offer diagnostic information
that you can use to correct the problem. If those remedies fail, feel free to contact
`CVX Support`_.

Selecting MOSEK as your default solver
--------------------------------------

Even if MOSEK is successfully added to your solver list, it will not automatically
be selected as your default solver. To change this, type the following two commands
on the MATLAB command line:

::

    cvx_solver mosek
    cvx_save_prefs
    
The first command changes the active solver to MOSEK, but only for the current session.
The second line saves that change to CVX's preference file, so that MOSEK will be 
selected as the active solver every time you start MATLAB.

Obtaining support for CVX and MOSEK
------------------------------------

If you encounter problems using CVX and MOSEK, please contact 
`CVX Support`_ first instead of MOSEK ApS.
If we can reproduce your problem, we will determine whether or not it is an
issue that is unique to CVX or needs to be forwarded to MOSEK ApS for further
analysis.

.. _CVX Support: http://support.cvxr.com/
.. _CVX Sales: mailto:sales@cvxr.com   
.. _Licensing: http://cvxr.com/cvx/licensing
.. _Academic License Request: http://cvxr.com/cvx/academic
