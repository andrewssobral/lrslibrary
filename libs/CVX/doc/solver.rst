.. _solvers:

=======
Solvers
=======

.. _supported-solvers:

Supported solvers
-----------------

This version of CVX supports four solvers, each with different capabilities:

============================================================= ==== ==== ====== ===== ====  ========= ======== ========
 Solver name                                                   LP   QP   SOCP   SDP   GP    Integer   MATLAB   Octave 
============================================================= ==== ==== ====== ===== ====  ========= ======== ========
`SeDuMi <http://sedumi.ie.lehigh.edu>`_                        Y    Y    Y       Y    E     N         Y        *soon*     
`SDPT3 <http://www.math.nus.edu.sg/~mattohkc/sdpt3.html>`_     Y    Y    Y       Y    E     N         Y        *soon*     
`Gurobi <http://gurobi.com>`_                                  Y    Y    Y       N    N     Y         P        *soon*     
`MOSEK <http://mosek.com>`_                                    Y    Y    Y       Y*   E     Y         P        *soon*     
`GLPK <http://www.gnu.org/software/glpk/>`_                    Y    N    N       N    N     Y         N        *soon*     
============================================================= ==== ==== ====== ===== ====  ========= ======== ========

(key: Y = Yes, N = No, E = Experimental, P = CVX Professional license required, * = Mosek 7 or later is required.)

Each solver has different capabilities and different levels of performance. For instance,
SeDuMi [Stu99]_, SDPT3 [TTT03]_, and MOSEK 7 support all of the continuous (non-integer) models 
that CVX itself supports, while Gurobi is more limited, in that it does not support semidefinite
constraints; and GLPK is limited even further. On the other hand, Gurobi, GLPK, and
MOSEK support integer consraints, while SeDuMi and SDPT3 do not.

SeDuMi and SDPT3 are included with the standard CVX distribution, so you do not need
to download an additional solver to start using CVX. We have also entered into contractual
arrangements with the developers of Gurobi and MOSEK that allow us to ship their binaries
with CVX as well, but using those solvers requires a CVX Professional license. Due to
license differences, we are *not* able to supply GLPK with CVX. However,
	
If you are having difficulty with one solver, *please try another*. No one solver performs
better than the others on *every* model CVX can generate---including commercial solvers.
That said, if you encounter a problem that one solver can handle well and another 
cannot, please send us a bug report (see :ref:`support`) and we will forward the
results to the solver's authors.

We have created special sections in this user guide for using Gurobi and MOSEK with CVX:

* Gurobi: :ref:`gurobi`
* Mosek:  :ref:`mosek`

Support for GLPK should be considered experimental, and has been provided primarly to support
upcoming Octave capability (that is *not* ready yet.)

.. _solver-selection:

Selecting a solver
------------------

The default solver is currently SDPT3. We have found that SeDuMi is faster for most
problems, but unfortunately not as reliable. None of the solvers are perfect, however,
and you may find for your application that another solver is preferred.

To see which solver is currently selected, simply type

::

    cvx_solver

To change the current solver, simply follow the ``cvx_solver`` with the name of your
chosen solver. For example, to select SeDuMi, type

::

    cvx_solver sedumi

The ``cvx_solver`` command is case insensitive, so ``cvx_solver SeDuMi`` 
will work just fine as well.

If you issue this command inside a model---that is, between ``cvx_begin`` and
``cvx_end`` it will change the solver *only* for that model; the next model will
use the previous choice. If, only the other hand, you issue a ``cvx_solver`` command
*outside* of a model, it will change the solver used for the remainder of your Matlab
session (or until you change it again).

If you would like to change the default solver *permanently*---that is, so that it remains
the default even if you quit and re-start Matlab---then make sure it is set properly, 
and then issue the command

::

	cvx_save_prefs
	
This command saves not only your solver choice, but also your settings for ``cvx_expert``,
``cvx_power_warning``, and ``cvx_precision`` as well.	

.. _solver-output:   
    
Controlling screen output
-------------------------

Once you gain confidence in using CVX and start incorporating it
into your larger algorithms and programs, you are likely going to want
to silence the messages it delivers to the screen. To do so, simply add
the ``quiet`` keyword to the ``cvx_begin`` command; that is,

::

    cvx_begin quiet
        ...
    cvx_end

Previous versions of CVX utilized a separate ``cvx_quiet`` command
and that command is still available in this version as well, if you
prefer it. Entering ``cvx_quiet true`` suppresses screen output from the
solver, while entering ``cvx_quiet false`` restores the screen output.
If you enter these commands within a model---that is, between
``cvx_begin`` and ``cvx_end``---it will affect only that model. If you
enter it *outside* of a model, it will affect all subsequent models.
Entering cvx_quiet with no arguments returns the current setting.

.. _interpreting:

Interpreting the results
------------------------

After a complete CVX specification has been entered and the
cvx_end command issued, the solver is called to generate a numerical
result. It proceeds to replace the variables in your model with the
computed numerical values, and creates the variable cvx_optval
containing the value of the objective function. It also summarizes the
result of its efforts in the form of a string named ``cvx_status``. The
possible values of ``cvx_status`` are as follows:

``Solved``
    A complementary (primal and dual) solution has been found. The
    primal and dual variables are replaced with their computed values,
    and the the optimal value of the problem is placed in cvx_optval
    (which, by convention, is :math:`0` for feasibility problems).

``Unbounded``
    The solver has determined that the problem is unbounded. The value
    of ``cvx_optval`` is set to ``-Inf`` for minimizations, and ``+Inf``
    for maximizations. (Feasibility problems, by construction, never
    produce an ``Unbounded`` status.) The values of any dual variables
    are replaced with ``NaN``, as the dual problem is in fact
    infeasible.

    For unbounded problems, CVX stores an *unbounded direction* into
    the problem variables. This is is a *direction* along which the
    feasible set is unbounded, and the optimal value approaches
    :math:`\pm\infty`. It is important to understand that this value is
    very likely *not* a feasible point. If a feasible point is required,
    the problem should be re-solved as a feasibility problem by omitting
    the objective. Mathematically speaking, given an unbounded direction
    :math:`v` and a feasible point :math:`x`, :math:`x+tv` is feasible
    for all :math:`t\geq0`, and the objective tends to :math:`-\infty`
    (for minimizations; :math:`+\infty` for maximizations) as 
    :math:`t\rightarrow+\infty` itself.

``Infeasible``
    The problem has been proven to be infeasible through the discovery
    of an unbounded direction. The values of the variables are filled
    with ``NaN``, and the value of ``cvx_optval`` is set to ``+Inf``
    for minimizations and feasibility problems, and ``-Inf`` for
    maximizations.

    Associated with a provably infeasible problem is an *unbounded dual
    direction*. Appropriate components of this direction are stored in
    the dual variables. Similarly to the ``Unbounded`` case, it is
    important to understand that the unbounded dual direction is very
    likely not a feasible dual point.

``Inaccurate/Solved``, ``Inaccurate/Unbounded``, ``Inaccurate/Infeasible``
    These three status values indicate that the solver was unable to
    make a determination to within the default numerical tolerance.
    However, it determined that the results obtained satisfied a
    "relaxed" tolerance leve and therefore may still be suitable for
    further use. If this occurs, you should test the validity of the
    computed solution before using it in further calculations. See
    :ref:`solver-precision` for a more advanced
    discussion of solver tolerances and how to make adjustments.
    
``Suboptimal``
    This status is possible only for *mixed-integer* problems. It is
    returned when the branching algorithm has discovered at least one
    feasible integer solution, but it was unable to continue the search
    process to global optimality. This will occur if the solver is 
    required to terminate due to a time limit or a forced interruption
    (for example, if the user types `Ctrl-C`.)     

``Failed``
    The solver failed to make sufficient progress towards a solution,
    even to within the "relaxed" tolerance setting. The values of
    cvx_optval and primal and dual variables are filled with
    ``NaN``. This result can occur because of numerical problems
    within SeDuMi, often because the problem is particularly "nasty" in
    some way (*e.g.*, a non-zero duality gap).

``Overdetermined``
    The presolver has determined that the problem has more equality
    constraints than variables, which means that the coefficient matrix
    of the equality constraints is singular. In practice, such problems
    are often, but not always, infeasible. Unfortunately, solvers
    typically cannot handle such problems, so a precise conclusion
    cannot be reached. The situations that most commonly produce an
    Overdetermined result are discussed in :ref:`overdetermined`.
   
.. _solver-precision:

Controlling precision
----------------------

.. note::

	We consider the modification of solver precision to be an advanced feature, to be
	used sparingly, if at all---and only after you have become 
	comfortable building models in CVX.

Numerical methods for convex optimization are not exact; they compute
their results to within a predefined numerical precision or tolerance.
Upon solution of your model, the tolerance level the solver has achieved
is returned in the ``cvx_slvtol`` variable. Attempts to interpret this
tolerance level in any absolute sense are not recommended. For one
thing, each solver computes it differently. For another, it depends
heavily on the considerable transformations that CVX applies to your
model before delivering it to the solver. So while you may find its
value interesting we strongly discourage dependence upon it within your
applications.

The tolerance levels that CVX selects by default have been inherited
from some of the underlying solvers being used, with minor modifications.
CVX actually considers *three* different tolerance levels
:math:`\epsilon_{\text{solver}}\leq\epsilon_{\text{standard}}\leq\epsilon_{\text{reduced}}`
when solving a model:

-  The *solver tolerance* :math:`\epsilon_{\text{solver}}` is the level
   requested of the solver. The solver will stop as soon as it achieves
   this level, or until no further progress is possible.
-  The *standard tolerance* :math:`\epsilon_{\text{standard}}` is the
   level at which CVX considers the model solved to full precision.
-  The *reduced tolerance* :math:`\epsilon_{\text{reduced}}` is the
   level at which CVX considers the model "inaccurately" sovled,
   returning a status with the ``Inaccurate/`` prefix. If this tolerance
   cannot be achieved, CVX returns a status of ``Failed``, and the
   values of the variables should not be considered reliable.

(See :ref:`interpreting` for more information about the
status messages.) Typically,
:math:`\epsilon_{\text{solver}}=\epsilon_{\text{standard}}`, but setting
:math:`\epsilon_{\text{standard}}<\epsilon_{\text{solver}}` has a useful
interpretation: it allows the solver to search for more accurate
solutions without causing an ``Inaccurate/`` or ``Failed`` condition if
it cannot do so. The default values of
:math:`[\epsilon_{\text{solver}},\epsilon_{\text{standard}},\epsilon_{\text{reduced}}]`
are set to :math:`[ \epsilon^{1/2}, \epsilon^{1/2}, \epsilon^{1/4} ]`,
where :math:`\epsilon=2.22\times10^{-16}` is the machine precision. This
should be quite sufficient for most applications.

If you wish to modify the tolerances, you may do so using the
``cvx_precision`` command. There are three ways to invoke this command.
Called with no arguments, it will print the current tolerance levels
to the screen; or if called as a function, it will return those levels
in a 3-element row vector.

Calling ``cvx_precision`` with a string argument allows you to select
from a set of predefined precision modes:

-  ``cvx_precision low``:
   :math:`[ \epsilon^{3/8}, \epsilon^{1/4}, \epsilon^{1/4} ]`
-  ``cvx_precision medium``:
   :math:`[ \epsilon^{1/2}, \epsilon^{3/8}, \epsilon^{1/4} ]`
-  ``cvx_precision default``:
   :math:`[ \epsilon^{1/2}, \epsilon^{1/2}, \epsilon^{1/4} ]`
-  ``cvx_precision high``:
   :math:`[ \epsilon^{3/4}, \epsilon^{3/4}, \epsilon^{3/8} ]`
-  ``cvx_precision best``: :math:`[ 0, \epsilon^{1/2}, \epsilon^{1/4} ]`

In function mode, these calls look like ``cvx_precision('low')``, etc.
Note that the ``best`` precision settings sets the solver target to
zero, which means that the solver continues as long as it can make
progress. It will often be slower than ``default``, but it is just as
reliable, and sometimes produces more accurate solutions.

Finally, the ``cvx_precision`` command can be called with a scalar, a
length-2 vector, or a length-3 vector. If you pass it a scalar, it will
set the solver and standard tolerances to that value, and it will
compute a default reduced precision value for you. Roughly speaking,
that reduced precision will be the square root of the standard
precision, with some bounds imposed to make sure that it stays
reasonable. If you supply two values, the smaller will be used for the
solver and standard tolerances, and the larger for the reduced
tolerance. If you supply three values, their values will be sorted, and
each tolerance will be set separately.

The ``cvx_precision`` command can be used either *within* a CVX
model or *outside* of it; and its behavior differs in each case. If you
call it from within a model, *e.g.*,

::

    cvx_begin
        cvx_precision high
        ...
    cvx_end

then the setting you choose will apply only until ``cvx_end`` is
reached. If you call it outside a model, *e.g.*,

::

    cvx_precision high
    cvx_begin
        ...
    cvx_end

then the setting you choose will apply *globally*; that is, to any
subsequent models that are created and solved. The local approach should
be preferred in an application where multiple models are constructed and
solved at different levels of precision.

If you call ``cvx_precision`` in function mode, either with a string or
a numeric value, it will return as its output the *previous* precision
vector---the same result you would obtain if you called it with no
arguments. This may seem confusing at first, but this is done so that
you can save the previous value in a variable, and restore it at the end
of your calculations; e.g.,

::

    cvxp = cvx_precision( 'high' );
    cvx_begin
        ...
    cvx_end
    cvx_precision( cvxp );

This is considered good coding etiquette in a larger application where
multiple CVX models at multiple precision levels may be employed. Of
course, a simpler but equally courteous approach is to call
``cvx_precision`` within the CVX model, as described above, so that
its effect lasts only for that model.

.. _solver-settings:

Advanced solver settings
------------------------

.. warning::

	This is an **advanced topic** for users who have a deep understanding of the 
	underlying solver they are using, or who have received specific advice from 
	the solver's developer for improving performance. Improper use of the
	``cvx_solver_settings`` command can cause unpredictable results.

Solvers can be tuned and adjusted in a variety of ways. Solver vendors attempt to select
default settings that will provide good performance across a broad range of
problems. But no solver, and no choice of settings, will perform well for every
possible model. On occasion, it may be worthwhile to give a particular special instructions
to improve its performance for a specific application. Unfortunately, such settings differ
from solver to solver, so there is no way for CVX to provide this ability in a verifiable,
reliable, global fashion.

Nevertheless, using the new ``cvx_solver_settings`` command, you can customize a solver's
settings when a specific model demands it. We cannot emphasize enough that this is an
*expert* feature to be employed by experienced modelers only. Indeed, if you are an
expert, you understand that these warnings are essential:

- CVX does not check the correctness of the settings you supply. If the solver rejects the
  settings, CVX will fail until you change or remove those settings.
- There is no guarantee that altering the settings will improve performance in any
  way; indeed, it can make the performance worse.
- CVX Research provides *no* documentation on the specific settings available for each
  solver; you will have to consult the solver's own documentation for this.
- The settings set here *override* any default values CVX may have chosen for each solver.
  Thus in certain cases, using this feature this may actually confuse CVX and cause it to
  misinterpret the results. For this reason, we cannot support all possible 
  combinations of custom settings.  
- Unless you have turned off solver output completely, CVX will warn you if any custom 
  settings are in effect every time you solve model.
  
With this warning out of the way, let us introduce ``cvx_solver_settings``. Typing

::

	cvx_solver_settings
	
at the command prompt provides a listing of the custom settings that have been provided
for the active solver. Custom settings are *specific to each solver*. Typing

::

	cvx_solver_settings -all
	
will provide a full list of the custom settings provided for *all* solvers.

To create a new custom setting for the current solver, use this syntax:

::

	cvx_solver_settings( '{name}', {value} )
	
``{name}`` must be a valid MATLAB variable/field name. ``{value}`` can be *any* valid Matlab
object; CVX does not check its value in any way.

To clear all custom settings for the active solver, type

::

	cvx_solver_settings -clear

To clear just a single setting, type

::

	cvx_solver_settings -clear {<name>}
	
To clear all settings for all solvers, type

::

	cvx_solver_settings -clearall
	
The settings created by the ``cvx_solver_settings`` command enjoy the same scope as
``cvx_solver``, ``cvx_precision``, and so forth. For instance, if you use this command
*within* a model---between ``cvx_begin`` and ``cvx_end``---the changes will apply only
to that particular model. If you issue the command *outside* of a particular model, the
change will persist through the end of the MATLAB session (or until you change it again).
Finally, if you use the ``cvx_save_prefs`` command, any custom settings you have added
will be saved and restored the next time you start Matlab.



	
