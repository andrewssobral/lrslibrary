.. index:: Examples

.. _quickstart:

=============
A quick start
=============

Once you have installed CVX (see :ref:`install`), you can start using it by
entering a CVX *specification* into a Matlab script or function, or
directly from the command prompt. To delineate CVX specifications
from surrounding Matlab code, they are preceded with the statement
``cvx_begin`` and followed with the statement ``cvx_end``. A
specification can include any ordinary Matlab statements, as well as
special CVX-specific commands for declaring primal and dual
optimization variables and specifying constraints and objective
functions.

Within a CVX specification, optimization variables have no numerical
value; instead, they are special Matlab objects. This enables Matlab to
distinguish between ordinary commands and CVX objective functions
and constraints. As CVX reads a problem specification, it builds an
internal representation of the optimization problem. If it encounters a
violation of the rules of disciplined convex programming (such as an
invalid use of a composition rule or an invalid constraint), an error
message is generated. When Matlab reaches the ``cvx_end`` command, it
completes the conversion of the CVX specification to a canonical
form, and calls the underlying core solver to solve it.

If the optimization is successful, the optimization variables declared
in the CVX specification are converted from objects to ordinary
Matlab numerical values that can be used in any further Matlab
calculations. In addition, CVX also assigns a few other related
Matlab variables. One, for example, gives the status of the problem (i.e.,
whether an optimal solution was found, or the problem was determined to
be infeasible or unbounded). Another gives the optimal value of the
problem. Dual variables can also be assigned.

This processing flow will become clearer as we introduce a number of
simple examples. We invite the reader to actually follow along with
these examples in Matlab, by running the ``quickstart`` script found in
the ``examples`` subdirectory of the CVX distribution. For example,
if you are on Windows, and you have installed the CVX distribution
in the directory ``D:\Matlab\cvx``, then you would type

::

    cd D:\Matlab\cvx\examples
    quickstart

at the Matlab command prompt. The script will automatically print key
excerpts of its code, and pause periodically so you can examine its
output. (Pressing "Enter" or "Return" resumes progress.)

.. index:: Least squares

Least squares
-------------

We first consider the most basic convex optimization problem,
least-squares (also known as linear regression). In a least-squares problem, we seek
:math:`x \in \mathbf{R}^n` that minimizes :math:`\|Ax-b\|_2`, where
:math:`A\in \mathbf{R}^{m \times n}` is skinny and full rank (i.e.,
:math:`m\geq n` and :math:`\operatorname*{\textbf{Rank}}(A)=n`). Let us
create the data for a small test problem in Matlab:

::

    m = 16; n = 8;
    A = randn(m,n);
    b = randn(m,1);

Then the least-squares solution :math:`x=(A^TA)^{-1}A^Tb` is
easily computed using the backslash operator:

::

    x_ls = A \ b;

Using CVX, the same problem can be solved as follows:

::

    cvx_begin
        variable x(n)
        minimize( norm(A*x-b) )
    cvx_end

(The indentation is used for purely stylistic reasons and is optional.)
Let us examine this specification line by line:

-  ``cvx_begin`` creates a placeholder for the new CVX
   specification, and prepares Matlab to accept variable declarations,
   constraints, an objective function, and so forth.
-  ``variable x(n)`` declares ``x`` to be an optimization variable of
   dimension :math:`n`. CVX requires that all problem variables be
   declared before they are used in the objective function or
   constraints.
-  ``minimize( norm(A*x-b) )`` specifies the objective function to be
   minimized.
-  ``cvx_end`` signals the end of the CVX specification, and causes
   the problem to be solved.

Clearly there is no reason to use
CVX to solve a simple least-squares problem. But this example serves
as sort of a "Hello world!" program in CVX; i.e., the simplest code
segment that actually does something useful.

When Matlab reaches the ``cvx_end`` command, the least-squares problem
is solved, and the Matlab variable ``x`` is overwritten with the
solution of the least-squares problem, i.e., :math:`(A^TA)^{-1}A^Tb`. Now
``x`` is an ordinary length-:math:`n` numerical vector, identical to
what would be obtained in the traditional approach, at least to within
the accuracy of the solver. In addition, several additional Matlab
variables are created; for instance,

-  ``cvx_optval`` contains the value of the objective function;
-  ``cvx_status`` contains a string describing the status of the
   calculation (see :ref:`interpreting`).

All of these quantities---``x``, ``cvx_optval``, and ``cvx_status``,
*etc.*---may now be freely used in other Matlab statements, just like
any other numeric or string values. [1]_

There is not much room for error in specifying a simple least-squares
problem, but if you make one, you will get an error or warning message.
For example, if you replace the objective function with

::

    maximize( norm(A*x-b) );

which asks for the norm to be maximized, you will get an error message
stating that a convex function cannot be maximized (at least in
disciplined convex programming):

::

    ??? Error using ==> maximize
    Disciplined convex programming error:
    Objective function in a maximization must be concave.

.. index::
	single: Least squares; bound-constrained
	single: Bound-constrained least squares
	    
Bound-constrained least squares
-------------------------------

Suppose we wish to add some simple upper and lower bounds to the
least-squares problem above: *i.e*.,

.. math::

   \begin{array}{ll}
       \mbox{minimize} & \|Ax-b\|_2\\
       \mbox{subject to} & l \preceq x \preceq u
   \end{array}

where :math:`l` and :math:`u` are given data vectors with the same
dimension as :math:`x`. The vector inequality
:math:`u \preceq v` means componentwise, i.e., :math:`u_i \leq v_i` for
all :math:`i`. We can no longer use the simple backslash notation to
solve this problem, but it can be transformed into a quadratic program
(QP) which can be solved without difficulty with a standard QP solver. [2]_

Let us provide some numeric values for ``l`` and ``u``:

::

    bnds = randn(n,2);
    l = min( bnds, [], 2 );
    u = max( bnds, [], 2 );

If you have the `Matlab Optimization
Toolbox <http://www.mathworks.com/products/optimization>`_, you can use ``quadprog``
to solve the problem as follows:

::

    x_qp = quadprog( 2*A'*A, -2*A'*b, [], [], [], [], l, u );

This actually minimizes the square of the norm, which is the same as
minimizing the norm itself. In contrast, the CVX specification is
given by

::

    cvx_begin
        variable x(n)
        minimize( norm(A*x-b) )
        subject to
            l <= x <= u
    cvx_end

Two new lines of CVX code have been added to the CVX specification:

-  The ``subject to`` statement does nothing---CVX provides this
   statement simply to make specifications more readable. As with
   indentation, it is optional.
-  The line ``l <= x <= u`` represents the :math:`2n` inequality
   constraints.

As before, when the ``cvx_end`` command is reached, the problem is
solved, and the numerical solution is assigned to the variable ``x``.
Incidentally, CVX will *not* transform this problem into a QP by
squaring the objective; instead, it will transform it into an SOCP. The
result is the same, and the transformation is done automatically.

In this example, as in our first, the CVX specification is longer
than the Matlab alternative. On the other hand, it is easier to read the
CVX version and relate it to the original problem. In contrast, the
``quadprog`` version requires us to know in advance the transformation
to QP form, including the calculations such as ``2*A'*A`` and
``-2*A'*b``. For all but the simplest cases, a CVX specification is
simpler, more readable, and more compact than equivalent Matlab code to
solve the same problem.

Other norms and functions
-------------------------

.. index:: linprog (MATLAB function)

Now let us consider some alternatives to the least-squares problem. Norm
minimization problems involving the :math:`\ell_\infty` or
:math:`\ell_1` norms can be reformulated as LPs, and solved using a
linear programming solver such as ``linprog`` in the Matlab Optimization
Toolbox; see, *e.g.*, Section 6.1 of `Convex
Optimization <http://www.stanford.edu/~boyd/cvxbook>`_. However,
because these norms are part of CVX's base library of functions,
CVX can handle these problems directly.

For example, to find the value of :math:`x` that minimizes the Chebyshev
norm :math:`\|Ax-b\|_\infty`, we can employ the ``linprog`` command from
the Matlab Optimization Toolbox:

::

    f    = [ zeros(n,1); 1          ];
    Ane  = [ +A,         -ones(m,1)  ; ...
             -A,         -ones(m,1) ];
    bne  = [ +b;         -b         ];
    xt   = linprog(f,Ane,bne);
    x_cheb = xt(1:n,:);

With CVX, the same problem is specified as follows:

::

    cvx_begin
        variable x(n)
        minimize( norm(A*x-b,Inf) )
    cvx_end

The code based on ``linprog``, and the CVX specification above will
both solve the Chebyshev norm minimization problem, i.e., each will
produce an :math:`x` that minimizes :math:`\|Ax-b\|_\infty`. Chebyshev
norm minimization problems can have multiple optimal points, however, so
the particular :math:`x`'s produced by the two methods can be different.
The two points, however, must have the same value of
:math:`\|Ax-b\|_\infty`.

Similarly, to minimize the :math:`\ell_1` norm :math:`\|\cdot\|_1`, we
can use ``linprog`` as follows:

::

    f    = [ zeros(n,1); ones(m,1);  ones(m,1)  ];
    Aeq  = [ A,          -eye(m),    +eye(m)    ];
    lb   = [ -Inf(n,1);  zeros(m,1); zeros(m,1) ];
    xzz  = linprog(f,[],[],Aeq,b,lb,[]);
    x_l1 = xzz(1:n,:) - xzz(n+1:end,:);

The CVX version is, not surprisingly,

::

    cvx_begin
        variable x(n)
        minimize( norm(A*x-b,1) )
    cvx_end

CVX automatically transforms both of these problems into LPs, not
unlike those generated manually for ``linprog``.

The advantage that automatic transformation provides is magnified if we
consider functions (and their resulting transformations) that are less
well-known than the :math:`\ell_\infty` and :math:`\ell_1` norms. For
example, consider the norm

.. math:: 

	\| Ax-b\|_{\mathrm{lgst},k} = |Ax-b|_{[1]}+ \cdots + |Ax-b|_{[k]},

where :math:`|Ax-b|_{[i]}` denotes the :math:`i`\ th largest element of
the absolute values of the entries of :math:`Ax-b`. This is indeed a
norm, albeit a fairly esoteric one. (When :math:`k=1`, it reduces to the
:math:`\ell_\infty` norm; when :math:`k=m`, the dimension of
:math:`Ax-b`, it reduces to the :math:`\ell_1` norm.) The problem of
minimizing :math:`\| Ax-b\|_{\mathrm{lgst},k}` over :math:`x` can be
cast as an LP, but the transformation is by no means obvious so we will
omit it here. But this norm is provided in the base CVX library, and
has the name ``norm_largest``, so to specify and solve the problem using
CVX is easy:

::

    k = 5;
    cvx_begin
        variable x(n);
        minimize( norm_largest(A*x-b,k) );
    cvx_end

Unlike the :math:`\ell_1`, :math:`\ell_2`, or :math:`\ell_\infty` norms,
this norm is not part of the standard Matlab distribution. Once you have
installed CVX, though, the norm is available as an ordinary Matlab
function outside a CVX specification. For example, once the code
above is processed, ``x`` is a numerical vector, so we can type

::

    cvx_optval
    norm_largest(A*x-b,k)

The first line displays the optimal value as determined by CVX; the
second recomputes the same value from the optimal vector ``x`` as
determined by CVX.

The list of supported nonlinear functions in CVX goes well beyond
``norm`` and ``norm_largest``. For example, consider the Huber penalty
minimization problem

.. math::

   \begin{array}{ll}
       \text{minimize} & \sum_{i=1}^m \phi( (Ax-b)_i )
   \end{array},

with variable :math:`x \in \mathbf{R}^n`, where :math:`\phi` is the
Huber penalty function

.. math:: 

	\phi(z) = \begin{cases} |z|^2 & |z|\leq 1 \\ 2|z|-1 & |z|\geq 1\end{cases}.

The Huber penalty function is convex, and has been provided in the
CVX function library. So solving the Huber penalty minimization
problem in CVX is simple:

::

    cvx_begin
        variable x(n);
        minimize( sum(huber(A*x-b)) );
    cvx_end

CVX automatically transforms this problem into an SOCP, which the
core solver then solves. (The CVX user, however, does not need to
know how the transformation is carried out.)

Other constraints
-----------------

We hope that, by now, it is not surprising that adding the simple
bounds :math:`l\preceq x\preceq u` to the problems above
is as simple as inserting the line ``l <= x <= u``
before the ``cvx_end`` statement in each CVX specification. In fact,
CVX supports more complex constraints as well. For example, let us
define new matrices ``C`` and ``d`` in Matlab as follows,

::

    p = 4;
    C = randn(p,n);
    d = randn(p,1);

Now let us add an equality constraint and a nonlinear inequality
constraint to the original least-squares problem:

::

    cvx_begin
        variable x(n);
        minimize( norm(A*x-b) );
        subject to
            C*x == d;
            norm(x,Inf) <= 1;
    cvx_end

Both of the added constraints conform to the DCP rules, and so are
accepted by CVX. After the ``cvx_end`` command, CVX converts
this problem to an SOCP, and solves it.

Expressions using comparison operators (``==``, ``>=``, *etc.*) behave
quite differently when they involve CVX optimization variables, or
expressions constructed from CVX optimization variables, than when
they involve simple numeric values. For example, because ``x`` is a
declared variable, the expression ``C*x==d`` causes a constraint to be
included in the CVX specification, and returns no value at all. On
the other hand, outside of a CVX specification, if ``x`` has an
appropriate numeric value---for example immediately after the
``cvx_end`` command---that same expression would return a vector of
``1``\ s and ``0``\ s, corresponding to the truth or falsity of each
equality. [3]_ Likewise, within a CVX specification, the statement
``norm(x,Inf)<=1`` adds a nonlinear constraint to the specification;
outside of it, it returns a ``1`` or a ``0`` depending on the numeric
value of ``x`` (specifically, whether its :math:`\ell_\infty`-norm is
less than or equal to, or more than, :math:`1`).

Because CVX is designed to support convex optimization, it must be
able to verify that problems are convex. To that end, CVX adopts
certain rules that govern how constraint and objective
expressions are constructed. For example, CVX requires that the
left- and right-hand sides of an equality constraint be affine. So a
constraint such as

::

    norm(x,Inf) == 1;

results in the following error:

::

    ??? Error using ==> cvx.eq
    Disciplined convex programming error:
    Both sides of an equality constraint must be affine.

Inequality constraints of the form :math:`f(x) \leq g(x)` or
:math:`g(x) \geq f(x)` are accepted only if :math:`f` can be verified as
convex and :math:`g` verified as concave. So a constraint such as

::

    norm(x,Inf) >= 1;

results in the following error:

::

    ??? Error using ==> cvx.ge
    Disciplined convex programming error:
    The left-hand side of a ">=" inequality must be concave.

The specifics of the construction rules are discussed in more detail in
:ref:`dcp`. These rules are relatively intuitive if
you know the basics of convex analysis and convex optimization.

An optimal trade-off curve
--------------------------

For our final example in this section, let us show how traditional
Matlab code and CVX specifications can be mixed to form and solve
multiple optimization problems. The following code solves the problem of
minimizing :math:`\|Ax-b\|_2 +\gamma \|x\|_1`, for a logarithmically
spaced vector of (positive) values of :math:`\gamma`. This gives us
points on the optimal trade-off curve between :math:`\|Ax-b\|_2` and
:math:`\|x\|_1`. An example of this curve is given in the figure below.

::

    gamma = logspace( -2, 2, 20 );
    l2norm = zeros(size(gamma));
    l1norm = zeros(size(gamma));
    fprintf( 1, '   gamma       norm(x,1)    norm(A*x-b)\n' );
    fprintf( 1, '---------------------------------------\n' );
    for k = 1:length(gamma),
        fprintf( 1, '%8.4e', gamma(k) );
        cvx_begin
            variable x(n);
            minimize( norm(A*x-b)+gamma(k)*norm(x,1) );
        cvx_end
        l1norm(k) = norm(x,1);
        l2norm(k) = norm(A*x-b);
        fprintf( 1, '   %8.4e   %8.4e\n', l1norm(k), l2norm(k) );
    end
    plot( l1norm, l2norm );
    xlabel( 'norm(x,1)' );
    ylabel( 'norm(A*x-b)' );
    grid on

.. figure:: tradeoff.pdf

   An example trade-off curve from the ``quickstart.m`` demo.

The ``minimize`` statement above illustrates one of the construction
rules to be discussed in :ref:`dcp`. A basic
principle of convex analysis is that a convex function can be multiplied
by a nonnegative scalar, or added to another convex function, and the
result is then convex. CVX recognizes such combinations and allows
them to be used anywhere a simple convex function can be---such as an
objective function to be minimized, or on the appropriate side of an
inequality constraint. So in our example, the expression

::

    norm(A*x-b)+gamma(k)*norm(x,1)

is recognized as convex by CVX, as long as ``gamma(k)`` is positive
or zero. If ``gamma(k)`` were negative, then this expression becomes the
sum of a convex term and a concave term, which causes CVX to
generate the following error:

::

    ??? Error using ==> cvx.plus
    Disciplined convex programming error:
    Addition of convex and concave terms is forbidden.

.. [1]
   If you type ``who`` or ``whos`` at the command prompt, you may see
   other, unfamiliar variables as well. Any variable that begins with
   the prefix ``cvx_`` is reserved for internal use by ``CVX`` itself,
   and should not be changed.
   
.. [2]
   There are also a number of solvers specifically designed to solve bound-constrained
   least-squares problems, such as `BCLS by Michael Friedlander <http://www.cs.ubc.ca/~mpf/bcls/>`_.  
   
.. [3]
   In fact, immediately after the ``cvx_end`` command above, you would
   likely find that most if not all of the values returned would be
   ``0``. This is because, as is the case with many numerical
   algorithms, solutions are determined only to within some nonzero
   numeric tolerance. So the equality constraints will be satisfied
   closely, but often not exactly.

