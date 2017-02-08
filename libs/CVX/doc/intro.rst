.. _introduction:

============
Introduction
============

What is CVX?
------------

.. index:: CVX, LP, DCP, MIDCP, QP, SOCP, SDP, GP, Matlab

.. index::
	see: Linear program; LP
	see: Disciplined convex program; DCP
	see: Mixed-integer disciplined convex program; MIDCP
	see: Mixed-integer problems; MIDCP
	see: Quadratic program; QP
	see: Second-order cone program; SOCP
	see: Semidefinite program; SDP
	see: Geometric program; GP

CVX is a modeling system for constructing and solving
:ref:`disciplined convex programs <what-is-dcp>` (DCPs).
CVX supports a number of standard problem types, including linear and
quadratic programs (LPs/QPs), second-order cone programs (SOCPs), and
semidefinite programs (SDPs). CVX can also solve much more complex convex
optimization problems, including many involving nondifferentiable
functions, such as :math:`\ell_1` norms. You can use CVX to conveniently
formulate and solve constrained norm minimization, entropy maximization,
determinant maximization, and many other convex programs. As of version
2.0, CVX also solves :ref:`mixed integer disciplined convex programs <what-is-midcp>` (MIDCPs)
as well, with an appropriate integer-capable solver.

To use CVX effectively, you need to know at least a bit about convex
optimization. For background on convex optimization, see the book
`Convex Optimization <http://www.stanford.edu/~boyd/cvxbook>`_ [BV04]_ or the
`Stanford course EE364A <http://www.stanford.edu/class/ee364a>`_.

CVX is implemented in `Matlab <http://mathworks.com>`_, effectively
turning Matlab into an optimization modeling language. Model
specifications are constructed using common Matlab operations and
functions, and standard Matlab code can be freely mixed with these
specifications. This combination makes it simple to perform the
calculations needed to form optimization problems, or to process the
results obtained from their solution. For example, it is easy to compute
an optimal trade-off curve by forming and solving a family of
optimization problems by varying the constraints. As another example,
CVX can be used as a component of a larger system that uses convex
optimization, such as a branch and bound method, or an engineering
design framework.

.. index:: SDP mode, GP mode

CVX provides special modes to simplify the construction of problems
from two specific problem classes. In 
:ref:`semidefinite programming (SDP) mode <sdp-mode>`, 
CVX applies a matrix interpretation to the
inequality operator, so that linear matrix inequalities (LMIs) and
SDPs may be expressed in a more natural form. In 
:ref:`geometric programming (GP) mode <gp-mode>`, 
CVX accepts all of the special functions and
combination rules of geometric programming, including monomials,
posynomials, and generalized posynomials, and transforms such problems
into convex form so that they can be solved efficiently. For background
on geometric programming, see this 
`tutorial paper <http://www.stanford.edu/~boyd/papers/gp_tutorial.html>`_ [BKVH05]_.

.. index:: SeDuMi, SDPT3, Gurobi, MOSEK, Solvers

.. index::
	single: Solvers; SeDuMi
	single: Solvers; SDPT3
	single: Solvers; Gurobi
	single: Solvers; MOSEK
	single: Solvers; commercial

Previous versions of CVX supported two free SQLP solvers,
`SeDuMi <http://sedumi.ie.lehigh.edu>`_ [Stu99]_ and
`SDPT3 <http://www.math.nus.edu.sg/~mattohkc/sdpt3.html>`_ [TTT03]_. These 
solvers are included with the CVX distribution. Starting with version 2.0,
CVX supports two *commercial* solvers as well,
`Gurobi <http://gurobi.com>`_ and `MOSEK <http://mosek.com>`_. For 
more information, see :ref:`solvers`.

The ability
to use CVX with commercial solvers is a new capability that we have decided
to include under a new CVX Professional license model. Academic users will
be able to utilize these features at no charge, but commercial users will require
a paid CVX Professional license. For more details, see :ref:`licensing`.

What's new?
~~~~~~~~~~~

If you browse the source code and documentation, 
you will find indications of support for Octave with CVX. However:

.. note:: 

  Unfortunately, for average end users (this means you!), Octave
  will *not* work. The *currently released* versions of Octave,
  including versions 3.8.0 and earlier, do not
  support CVX. Please do not waste your time by trying!

We are working hard with the Octave team on final updates
to bring CVX to Octave, and we anticipate version 3.8.1 or 3.9.0 will
be ready. We add this here to warn you *not* to interpret the mentions
of Octave in the code as a hidden code to try it yourself!

.. index:: DCP

.. _what-is-dcp:

What is disciplined convex programming?
---------------------------------------

*Disciplined convex programming* is a methodology for constructing
convex optimization problems proposed by 
Michael Grant, Stephen Boyd, and Yinyu Ye [GBY06]_, [Gra04]_.
It is meant to support the formulation and construction of optimization 
problems that the user intends *from the outset* to be convex.

.. index:: DCP ruleset

Disciplined convex programming
imposes a set of conventions or rules, which we call :ref:`the DCP ruleset <dcp>`.
Problems which adhere to the ruleset can be rapidly and automatically
verified as convex and converted to solvable form. Problems that violate
the ruleset are rejected---even when the problem is convex. That is not
to say that such problems cannot be solved using DCP; they just need to
be rewritten in a way that conforms to the DCP ruleset.

A detailed description of the DCP ruleset is given in :ref:`dcp`.
It is extremely important for anyone who intends to
actively use CVX to understand it. The ruleset is simple to learn, and
is drawn from basic principles of convex analysis. In return for
accepting the restrictions imposed by the ruleset, we obtain
considerable benefits, such as automatic conversion of problems to
solvable form, and full support for nondifferentiable functions. In
practice, we have found that disciplined convex programs closely
resemble their natural mathematical forms.

.. index:: MIDCP

.. _what-is-midcp:

Mixed integer problems
~~~~~~~~~~~~~~~~~~~~~~

With version 2.0, CVX now supports *mixed integer* disciplined convex programs (MIDCPs).
A MIDCP is a model that obeys the same convexity rules as standard DCPs, except
that one or more of its variables is constrained to take on integral values. In other
words, if the integer constraints are removed, the result is a standard DCP.

Unlike a true DCP, a mixed integer problem is *not* convex. Finding the global optimum
requires the combination of a traditional convex optimization algorithm with an exhaustive
search such as a branch-and-bound algorithm. Some CVX solvers do not include this second
piece and therefore do not support MIDCPs; see :ref:`solvers` for more information.
What is more, even the best solvers cannot
guarantee that every moderately-sized MIDCP can be solved in a reasonable amount of time.

Mixed integer disciplined convex programming represents new territory for the 
CVX modeling framework---and for the supporting solvers as well. While solvers
for mixed integer linear and quadratic programs (MILP/MIQP) are reasonably mature,
support for more general convex nonlinearities is a relatively new
development. We anticipate that MIDCP support will improve over time.

What CVX is *not*
------------------

CVX is *not* meant to be a tool for checking if your problem is convex.
You need to know a bit about convex optimization to effectively use CVX;
otherwise you are the proverbial monkey at the typewriter, hoping to
(accidentally) type in a valid disciplined convex program. If you are
not certain that your problem is convex *before* you enter it into CVX,
you are using the tool improperly, and your efforts will likely fail.

CVX is *not* meant for very large problems, so if your problem is very
large (for example, a large image processing or machine learning problem), CVX is unlikely
to work well (or at all). For such problems you will likely need to
directly call a solver, or to develop your own methods, to get the
efficiency you need.

For such problems CVX can play an important role, however. Before
starting to develop a specialized large-scale method, you can use CVX to
solve scaled-down or simplified versions of the problem, to rapidly
experiment with exactly what problem you want to solve. For image
reconstruction, for example, you might use CVX to experiment with
different problem formulations on :math:`50 \times 50` pixel images.

CVX *will* solve many medium and large scale problems, provided they
have exploitable structure (such as sparsity), and you avoid ``for``
loops, which can be slow in Matlab, and functions like ``log`` and ``exp`` that
require successive approximation. If you encounter difficulties in
solving large problem instances, consider posting your model to the
`CVX Forum <http://ask.cvxr.com>`_; the CVX community
may be able to suggest an equivalent formulation that CVX 
can process more efficiently.

.. index::
	single: License
	single: License; free
	single: License; commercial
	single: License; academic
	single: Academic licensing
	single: CVX Professional
	
.. _licensing:

Licensing
---------

CVX is free for use in both academic and commercial settings when paired with
a free solver---including the versions of SeDuMi and SDPT3 that are included with the 
package.

With version 2.0, we have added the ability to connect CVX to *commercial* solvers as well. 
This new functionality is released under a *CVX Professional* product tier
which we intend to license to commercial users for a fee, and offer to academic users
at no charge. The licensing structure is as follows:

* *All users* are free to use the standard features of CVX *at no charge*.
  This includes the ability to construct and solve any of the models 
  supported by the free solvers SeDuMi and SDPT3.
* *Commercial users* who wish to solve CVX models using Gurobi or MOSEK will
  need to purchase a CVX Professional license. Please send an email to
  `CVX Research <mailto:sales@cvxr.com>`_ for inquiries.
  for an availability schedule and pricing details.
* *Academic users* may utilize the CVX Professional capability *at no charge*.
  To obtain an academic license, please visit the
  `Academic licenses <http://cvxr.com/cvx/academic>`_ page on the 
  CVX Research web site.
  
The bulk of CVX remains open source under a slightly modified version of the GPL Version
2 license. A small number of files that support the CVX Professional functionality remain
closed source. If those files are removed, the modified package remains *fully functional*
using the free solvers, SeDuMi and SDPT3. Users
may freely modify, augment, and redistribute this free version of CVX, as long as all
modifications are themselves released under the same license. This includes adding support
for new solvers released under a free software license such as the GPL. 
For more details, please see the full :ref:`Licensing <licensing2>` section.

