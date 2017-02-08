.. _gurobi:

=====================
Using Gurobi with CVX
=====================

About Gurobi
------------

`Gurobi Optimization <http://www.gurobi.com>`_ was founded in 2008 by some of the most
experienced and respected members of the optimization community. The Gurobi solver quickly became an industry performance leader in linear, quadratic, and mixed-integer programming. Gurobi is a fantastic solver for use with CVX, particularly with the new integer and binary variable capability added in CVX 2.0.

Using CVX with Gurobi requires both a CVX Professional license and a
Gurobi license. Please visit `Licensing`_ for information about
CVX licensing, and Gurobi's `Licensing Overview <http://www.gurobi.com/products/licensing-and-pricing/licensing-overview>`_
page for information about Gurobi licenses.

Academic users can obtain both licenses at no charge. An academic
CVX Professional license is obtained by submitting an `Academic License Request`_. 
For instructions on obtaining an academic license for Gurobi, please see
Gurobi's `Academic Licenses <http://www.gurobi.com/products/licensing-and-pricing/academic-licensing>`_ page.

.. _gurobilic:

Using the bundled version of Gurobi
-----------------------------------

.. note:: 

    The bundled version of Gurobi can only be used within CVX. If you wish to use Gurobi
    outside of CVX as well, you will need a standalone Gurobi installation.

If you wish to use CVX with the bundled version of Gurobi, you will need three things:

* A CVX Professional Solver Bundle, available `here <http://cvxr.com/cvx/download>`_.
* A Gurobi license code, which is composed of 32 hexidecimal digits in the format
  ``xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx``. If you purchase a commercial CVX+Gurobi
  package, you will receive this code in an email from CVX Research. If you are an
  academic user, you will receive it directly from Gurobi.
* A CVX Professional license, saved to a convenient location on your local disk.
  
Installation proceeds as follows:

* First, install CVX in the standard manner according to the directions found in
  :ref:`install`. Do not attempt to install either license at this stage.
* Next, retrieve your Gurobi license key by running the command ``cvx_grbgetkey`` *{code}*,
  where *{code}* is the 32-digit Gurobi key. The command will look something like this::

    cvx_grbgetkey xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
    
  *Important note for academic users:* this step must be run from a computer
  connected to your university network (a VPN is usually sufficient). Please
  consult `this page <http://www.gurobi.com/documentation/5.5/quick-start-guide/node9>`_
  of the Gurobi documentation for details.
* Finally, install your CVX Professional license according to the directions
  found in :ref:`licinstall`. 
      
If you complete these steps successfully, ``cvx_setup`` will add Gurobi to its solver
list. If you have an academic or dual-solver CVX Professional license, the MOSEK solver
will be added to the solver list as well. If for some reason installation fails, the
output of ``cvx_setup`` will provide diagnostic information that you can use to rectify
the problem. If you are still unable to complete the installation, feel free to contact
`CVX Support`_.

.. _gurobistandalone:

Using CVX with a standalone Gurobi installation
-----------------------------------------------

If you wish to use CVX with a standalone installation of Gurobi, you will
need the following four things:

* A Gurobi installation package, or a pre-existing Gurobi installation. 
  CVX works with Gurobi 5.0 or later, but use
  of the latest version is always recommended.
* A Gurobi license code or key file, if you are installing Gurobi for the first time.
* A standard CVX package, available `here <http://cvxr.com/cvx/download>`_.
  You do not need the Professional Solver bundle.
* A CVX Professional license, saved to a convenient location on your local disk.

Installation proceeds as follows:

* *Install Gurobi*. See `Downloading and Installation <http://www.gurobi.com/documentation/5.5/quick-start-guide/node1>`_
  from the Gurobi Quick Start Guide.
* *Install the Gurobi license.* See `How to Obtain and Install a Gurobi License <http://www.gurobi.com/documentation/5.5/quick-start-guide/node5>`_.
* *Connect your Gurobi installation to MATLAB.* See `Setting up Gurobi for MATLAB 
  <http://www.gurobi.com/documentation/5.5/quick-start-guide/node120>`_.
* *Install CVX and/or the CVX Professional license*. See :ref:`install` and :ref:`licinstall`.
  Even if these have already been installed, you *must* at least re-run ``cvx_setup``
  so that CVX can locate Gurobi and add it to your solver list.
  
If you complete these steps successfully, ``cvx_setup`` will show that Gurobi has been
recognized and added to the solver list. If for some reason installation fails, the
output of ``cvx_setup`` will provide diagnostic information that you can use to rectify
the problem. If you are still unable to complete the installation, feel free to contact
`CVX Support`_.

Selecting Gurobi as your default solver
---------------------------------------

Even if Gurobi is successfully added to your solver list, it will not automatically
be selected as your default solver. To change this, type the following two commands
on the MATLAB command line:

::

    cvx_solver gurobi
    cvx_save_prefs
    
The first command changes the active solver to Gurobi, but only for the current session.
The second line saves that change to CVX's preference file, so that Gurobi will be 
selected as the active solver every time you start MATLAB.

Obtaining support for CVX and Gurobi
------------------------------------

If you encounter problems using CVX and Gurobi, please contact 
`CVX Support`_ first instead of Gurobi Optimization.
If we can reproduce your problem, we will determine whether or not it is an
issue that is unique to CVX or needs to be forwarded to Gurobi for further
analysis.

.. _CVX Support: http://support.cvxr.com/
.. _CVX Sales: mailto:sales@cvxr.com   
.. _Licensing: http://cvxr.com/cvx/licensing
.. _Academic License Request: http://cvxr.com/cvx/academic
