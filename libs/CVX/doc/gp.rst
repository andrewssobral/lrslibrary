.. _gp-mode:

==========================
Geometric programming mode
==========================

Geometric programs (GPs) are special mathematical programs that can be
converted to convex form using a change of variables. The convex form of
GPs can be expressed as DCPs, but CVX also provides a special mode
that allows a GP to be specified in its native form. CVX will
automatically perform the necessary conversion, compute a numerical
solution, and translate the results back to the original problem.

To utilize GP mode, you must begin your CVX specification with the
command ``cvx_begin gp`` or ``cvx_begin GP`` instead of simply
``cvx_begin``. For example, the following code, found in the example
library at :file:`gp/max_volume_box.m`, determines the maximum volume box
subject to various area and ratio constraints:

::

    cvx_begin gp
        variables w h d
        maximize( w * h * d )
        subject to
            2*(h*w+h*d) <= Awall;
            w*d <= Afloor;
            alpha <= h/w >= beta;
            gamma <= d/w <= delta;
    cvx_end

As the example illustrates, CVX supports the construction of
monomials and posynomials using addition, multiplication, division (when
appropriate), and powers. In addition, CVX supports the construction
of *generalized geometric programs* (GGPs), by permitting the use of
*generalized posynomials* wherever posynomials are permitted in standard
GP. More information about generalized geometric programs is provided in
this
`tutorial <http://www.stanford.edu/~boyd/papers/gp_tutorial.html>`_.

The solvers used in this version of CVX do not support geometric
programming natively. Instead, they are solved using the successive
approximation technique described in :ref:`successive`.
This means that solving GPs can be slow, but for small and medium sized problems, the method
works well.

In the remainder of this section, we will describe specific rules that
apply when constructing models in GP mode.

Top-level rules
---------------

CVX supports three types of geometric programs:

-  A *minimization problem*, consisting of a generalized posynomial
   objective and zero or more constraints.
-  A *maximization problem*, consisting of a *monomial* objective and
   zero or more constraints.
-  A *feasibility problem*, consisting of one or more constraints.

The asymmetry between minimizations and maximizations---specifically,
that only monomial objectives are allowed in the latter---is an
unavoidable artifact of the geometry of GPs and GGPs.

Constraints
-----------

Three types of constraints may be specified in geometric programs:

-  An *equality constraint*, constructed using ``==``, where both sides
   are monomials.
-  A *less-than inequality constraint* ``<=`` where the left side is a
   generalized posynomial and the right side is a monomial.
-  A *greater-than inequality constraint* ``>=`` where the left side is
   a monomial and the right side is a generalized posynomial.

As with DCPs, non-equality constraints are not permitted; and while
strict inequalities ``<``, ``>`` are supported, they are treated as
non-strict inequalities and should therefore be avoided.

Expressions
-----------

The basic building blocks of generalized geometric programming are
monomials, posynomials, and generalized posynomials. A valid monomial is

-  a declared variable;
-  the product of two or more monomials;
-  the ratio of two monomials;
-  a monomial raised to a real power; or
-  a call to one of the following functions with monomial arguments:
   ``prod``, ``cumprod``, ``geo_mean``, ``sqrt``.

A valid posynomial expression is

-  a valid monomial;
-  the sum of two or more posynomials;
-  the product of two or more posynomials;
-  the ratio of a posynomial and a monomial;
-  a posynomial raised to a positive integral power; or
-  a call to one of the following functions with posynomial arguments:
   ``sum``, ``cumsum``, ``mean``, ``prod``, ``cumprod``.

A valid generalized posynomial expression is

-  a valid posynomial;
-  the sum of two or more generalized posynomials;
-  the product of two or more generalized posynomials;
-  the ratio of a generalized posynomial and a monomial;
-  a generalized posynomial raised to a positive real power; or
-  a call to one of the following functions with arguments that are
   generalized posynomials: ``sum``, ``cumsum``, ``mean``, ``prod``,
   ``cumprod``, ``geo_mean``, ``sqrt``, ``norm``, ``sum_largest``,
   ``norm_largest``.

It is entirely possible to create and manipulate arrays of monomials,
posynomials, and/or generalized posynomials in CVX, in which case
these rules extend in an obvious manner. For example, the product of two
monomial matrices produces a matrix whose entries are polynomials (or monomials
in special cases).