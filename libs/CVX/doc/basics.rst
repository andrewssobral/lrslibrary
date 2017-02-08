==========
The Basics
==========

.. _begin-end:

``cvx_begin`` and ``cvx_end``
-----------------------------

All CVX models must be preceded by the command ``cvx_begin`` and
terminated with the command ``cvx_end``. All variable declarations,
objective functions, and constraints should fall in between.
The cvx_begin command may include one more more modifiers:

``cvx_begin quiet``
	Prevents the model from producing any screen output while it is being solved.
	
``cvx_begin sdp``
	Invokes :ref:`semidefinite programming mode <sdp-mode>`.
	
``cvx_begin gp``
	Invokes :ref:`geometric programming mode <gp-mode>`.
	
These modifiers may be combined when appropriate; for instance, ``cvx_begin sdp quiet`` 
invokes SDP mode and silences the solver output.

.. _variables:

Variables
---------

All variables must be declared using the
``variable`` command (or ``variables`` command; see below)
before they can be used in constraints or an objective function.
A ``variable`` command includes the name of the variable, an optional
dimension list, and one or more keywords that provide additional
information about the content or structure of the variable.

Variables can be real or complex scalars, vectors, matrices, or
:math:`n`-dimensional arrays. For instance,

::

    variable X
    variable Y(20,10)
    variable Z(5,5,5)

declares a total of 326 (scalar) variables: a scalar ``X``, a 20x10 matrix ``Y`` (containing 200 scalar variables),
and a 5x5x5 array ``Z`` (containing 125 scalar variables).

Variable declarations can also include one or more *keywords* to denote various
structures or conditions on the variable. For instance, to declare a complex
variable, use the ``complex`` keyword:

::

	variable w(50) complex
	
Nonnegative variables and symmetric/Hermitian positive semidefinite (PSD) matrices can 
be specified with the ``nonnegative`` and ``semidefinite`` keywords, respectively:

::

    variable x(10) nonnegative
    variable Z(5,5) semidefinite
    variable Q(5,5) complex semidefinite
    
In this example, ``x`` is a nonnegative vector, and ``Z`` is a real symmetric PSD matrix
and ``Q``is a complex Hermitian PSD matrix. As we will see below, ``hermitian semidefinite``
would be an equivalent choice for this third case.

For MIDCPs, the ``integer`` and ``binary`` keywords are used to declare integer
and binary variables, respectively:

::

    variable p(10) integer
    variable q binary
    
A variety of keywords are available to help construct variables with
*matrix structure* such as symmetry or bandedness.
For example, the code segment

::

    variable Y(50,50) symmetric
    variable Z(100,100) hermitian toeplitz

declares Y to be a real :math:`50 \times 50` symmetric matrix
variable, and Z a :math:`100 \times 100` Hermitian
Toeplitz matrix variable. (Note that the ``hermitian`` keyword also specifies
that the matrix is complex.) The currently supported structure
keywords are:

::

	banded(lb,ub)      diagonal           hankel             hermitian
	skew_symmetric     symmetric          toeplitz           tridiagonal
	lower_bidiagonal   lower_hessenberg   lower_triangular   
	upper_bidiagonal   upper_hankel       upper_hessenberg   upper_triangular
	
The underscores can actually be omitted; so, for example, ``lower triangular``
is acceptable as well. These keywords are self-explanatory with a couple of exceptions:

``banded(lb,ub)``
   the matrix is banded with a lower bandwidth lb
   and an upper bandwidth ub. If both lb and ub are zero,
   then a diagonal matrix results. ub can be omitted, in which case
   it is set equal to lb. For example, ``banded(1,1)`` (or
   ``banded(1)``) is a tridiagonal matrix.

``upper_hankel``
   The matrix is Hankel (i.e., constant along
   antidiagonals), and zero below the central antidiagonal, i.e., for
   :math:`i+j>n+1`.

When multiple keywords are supplied, the resulting matrix structure is
determined by intersection. For example, ``symmetric tridiagonal`` is a
valid combination. That said, CVX does reject combinations such as
``symmetric lower_triangular`` when a more reasonable alternative
exists---``diagonal``, in this case. Furthermore, if the keywords fully
conflict, such that \emph{no} non-zero matrix that satisfies all keywords, 
an error will result.

Matrix-specific keywords can be applied to
:math:`n`-dimensional arrays as well: each 2-dimensional "slice" of the
array is given the stated structure. So for instance, the declaration

::

	variable R(10,10,8) hermitian semidefinite
	
constructs 8 :math:`10\times 10` complex Hermitian PSD matrices, stored in the 
2-D slices of ``R``.

As flexible as the ``variable`` statement may be, it can only be used to declare 
a single variable, which can be inconvenient if you have a lot of variables to
declare. For this reason, the ``variables`` statement is provided which
allows you to declare multiple variables; i.e.,

:: 

	variables x1 x2 x3 y1(10) y2(10,10,10);

The one limitation of the ``variables`` command is that it cannot
declare complex, integer, or structured variables.
These must be declared one at a time, using the singular ``variable``
command.

.. _objectives:

Objective functions
-------------------

Declaring an objective function requires the use of the ``minimize`` or
``maximize`` function, as appropriate. (For the benefit of our users whose
English favors it, the synonyms ``minimise`` and ``maximise`` are
provided as well.) The objective function in a call to minimize must
be convex; the objective function in a call to maximize must be
concave; for instance:

::

	minimize( norm( x, 1 ) )
	maximize( geo_mean( x ) )

At most one objective function may be declared in a
CVX specification, and it must have a scalar value.

If no objective function is specified, the problem is interpreted as a
*feasibility problem*, which is the same as performing a minimization with
the objective function set to zero. In this case, cvx_optval is
either 0, if a feasible point is found, or +Inf, if the
constraints are not feasible.

.. _constraints:

Constraints
-----------

The following constraint types are supported in CVX:

-  Equality ``==`` constraints, where both the left- and right-hand
   sides are affine expressions.
-  Less-than ``<=`` inequality constraints, where the left-hand
   expression is convex, and the right-hand expression is concave.
-  Greater-than ``>=`` constraints, where the left-hand expression is
   concave, and the right-hand expression is convex.
   
The non-equality operator ``~=`` may *never* be 
used in a constraint; in any case, such constraints are rarely convex. 
The latest version of CVX now allows you to chain inequalities
together; *e.g.*, ``l <= x <= u``.  (Previous versions did not allow chained 
inequalities.)

Note the important distinction between the single equals ``=``, which is an
assignment, and the double equals ``==``, which denotes equality; for more 
on this distinction, see :ref:`assignment` below.

Strict inequalities ``<`` and ``>`` are accepted as well, but they are interpreted 
identically to their nonstrict counterparts. We strongly discourage their use, and
a future version of CVX may remove them altogether. For the reasoning behind
this, please see the fuller discussion in :ref:`strict`.

Inequality and equality constraints are applied in an elementwise fashion,
matching the behavior of MATLAB itself. For instance, if ``A`` and ``B`` are
:math:`m \times n` arrays, then ``A<=B`` is
interpreted as :math:`mn` (scalar) inequalities ``A(i,j)<=B(i,j)``. When one
side or the other is a scalar, that value is replicated; for instance, 
``A>0`` is interpreted as ``A(i,j)>=0``.  

The elementwise treatment of inequalities is altered in
:ref:`semidefinite programming mode <sdp-mode>`; see that section for more details.

CVX also supports a *set membership* constraint; see :ref:`sets` below.

.. _functions:

Functions
---------

The base CVX function library includes a variety of convex,
concave, and affine functions which accept CVX variables or
expressions as arguments. Many are common Matlab functions such as
sum, trace, diag, sqrt, max, and min,
re-implemented as needed to support CVX; others are new functions
not found in Matlab. A complete list of the functions in the base
library can be found in :ref:`funcref`. It is
also possible to add your own new functions; see
:ref:`newfunc`.

An example of a function in the base library is the quadratic-over-linear
function ``quad_over_lin``:

.. math::

	f:\mathbf{R}^{n}\times\mathbf{R} \rightarrow \mathbf{R}, \quad
	f(x,y) = \begin{cases} x^T x / y & y > 0 \\ +\infty & y \leq 0 \end{cases}
	
(The function also accepts complex :math:`x`, but we'll consider
real :math:`x` to keep things simple.) The quadratic-over-linear
function is convex in :math:`x` and :math:`y`, and so can be used as an
objective, in an appropriate constraint, or in a more complicated
expression. We can, for example, minimize the quadratic-over-linear
function of :math:`(Ax-b,c^Tx+d)` using

::

    minimize( quad_over_lin( A * x - b, c' * x + d ) );

inside a CVX specification, assuming x is a vector optimization
variable, A is a matrix, b and c are vectors, and d is a
scalar. CVX recognizes this objective expression as a convex
function, since it is the composition of a convex function (the
quadratic-over-linear function) with an affine function.

You can also use the function ``quad_over_lin`` *outside* a CVX
specification. In this case, it just computes its (numerical) value,
given (numerical) arguments. If ``c'*x+d`` is positive, then the
result is numerically equivalent tp

::

	( ( A * x - b )' * ( A * x - b ) ) / ( c' * x + d )

However, the ``quad_over_lin`` function also
performs a domain check, so it returns ``Inf`` if ``c'*x+d`` is zero or negative.

.. _sets:

Set membership
--------------

CVX supports the definition and use of convex sets. The base
library includes the cone of positive semidefinite :math:`n \times n`
matrices, the second-order or Lorentz cone, 
and various norm balls. A complete list of sets supplied in the base library
is given in :ref:`sets-ref`.

Unfortunately, the Matlab language does not have a set membership
operator, such as ``x in S``, to denote :math:`x \in S`. So in CVX,
we use a slightly different syntax to require that an expression is in a
set. To represent a set we use a *function* that returns an unnamed
variable that is required to be in the set. Consider, for example,
:math:`\mathbf{S}^n_+`, the cone of symmetric positive semidefinite
:math:`n \times n` matrices. In CVX, we represent this by the
function semidefinite(n), which returns an unnamed new variable,
that is constrained to be positive semidefinite. To require that the
matrix expression X be symmetric positive semidefinite, we use the
syntax 

::

	X == semidefinite(n)

The literal meaning of this is that
``X`` is constrained to be equal to some unnamed variable, which is
required to be an :math:`n \times n` symmetric positive semidefinite
matrix. This is, of course, equivalent to saying that X must itself be
symmetric positive semidefinite.

As an example, consider the constraint that a (matrix) variable X
is a correlation matrix, i.e., it is symmetric, has unit diagonal
elements, and is positive semidefinite. In CVX we can declare such a
variable and impose these constraints using

::

    variable X(n,n) symmetric;
    X == semidefinite(n);
    diag(X) == 1;

The second line here imposes the constraint that X be positive
semidefinite. (You can read "==\ " here as "is" or "is in", so the
second line can be read as X is positive semidefinite'.) The
lefthand side of the third line is a vector containing the diagonal
elements of X, whose elements we require to be equal to one.

If this use of equality constraints to represent set membership remains
confusing or simply aesthetically displeasing, we have created a
"pseudo-operator" ``<In>`` that you can use in its place. So, for
example, the semidefinite constraint above can be replaced by

::

    X <In> semidefinite(n);

This is exactly equivalent to using the equality constraint operator,
but if you find it more pleasing, feel free to use it. Implementing this
operator required some Matlab trickery, so don't expect to be able to
use it outside of CVX models.

Sets can be combined in affine expressions, and we can constrain an
affine expression to be in a convex set. For example, we can impose a
constraint of the form

::

    A*X*A'-X <In> B*semidefinite(n)*B';

where X is an :math:`n \times n` symmetric variable matrix, and
A and B are :math:`n \times n` constant matrices. This
constraint requires that :math:`AXA^T-X=BYB^T`, for some
:math:`Y \in \mathbf{S}^n_+`.

CVX also supports sets whose elements are ordered lists of
quantities. As an example, consider the second-order or Lorentz cone,

.. math:: 

	\mathbf{Q}^m = \left\{\, (x,y) \in \mathbf{R}^m\times\mathbf{R}\,~|~\, \| x \|_2 \leq y \,\right\} = \operatorname{\textbf{epi}}\|\cdot\|_2,

where :math:`\operatorname{\textbf{epi}}` denotes the epigraph of a
function. An element of :math:`\mathbf{Q}^m` is an ordered list, with
two elements: the first is an :math:`m`-vector, and the second is a
scalar. We can use this cone to express the simple least-squares problem
from the section `Least squares <#least-squares>`_ (in a fairly
complicated way) as follows:

.. math::

   \begin{array}{ll}
       \text{minimize}   & y \\
       \text{subject to} & ( A x - b, y ) \in \mathbf{Q}^m.
   \end{array}

CVX uses Matlab's cell array facility to mimic this notation:

::

    cvx_begin
        variables x(n) y;
        minimize( y );
        subject to
            { A*x-b, y } <In> lorentz(m);
    cvx_end

The function call ``lorentz(m)`` returns an unnamed variable (i.e., a pair
consisting of a vector and a scalar variable), constrained to lie in the
Lorentz cone of length ``m``. So the constraint in this specification
requires that the pair ``{ A*x-b, y }`` lies in the appropriately-sized
Lorentz cone.

.. _dual-variables:

Dual variables
--------------

When a disciplined convex program is solved, the associated *dual
problem* is also solved. (In this context, the original problem is
called the *primal problem*.) The optimal dual variables, each of which
is associated with a constraint in the original problem, give valuable
information about the original problem, such as the sensitivities with
respect to perturbing the constraints (*c.f.* `Convex
Optimization <http://www.stanford.edu/~boyd/cvxbook>`_, chapter 5). To
get access to the optimal dual variables in CVX, you simply declare
them, and associate them with the constraints. Consider, for example,
the LP

.. math::

   \begin{array}{llcll}
       \mbox{minimize} & c^Tx \\
       \mbox{subject to} & Ax \preceq b,
   \end{array}

with variable :math:`x\in\mathbf{R}^n`, and :math:`m` inequality
constraints.  To associate 
the dual variable :math:`y` with the inequality
constraint :math:`Ax\preceq b` in this LP, we use the following
syntax:

::

    n = size(A,2);
    cvx_begin
        variable x(n);
        dual variable y;
        minimize( c' * x );
        subject to
            y : A * x <= b;
    cvx_end

The line

::

    dual variable y

tells CVX that y will represent the dual variable, and the line

::

    y : A * x <= b;

associates it with the inequality constraint. Notice how the colon ``:``
operator is being used in a different manner than in standard Matlab,
where it is used to construct numeric sequences like ``1:10``. This new
behavior is in effect only when a dual variable is present, so there
should be no confusion or conflict. No dimensions are given for ``y``;
they are automatically determined from the constraint with which it is
associated. For example, if :math:`m=20`, typing ``y`` at the Matlab
command prompt immediately before cvx_end yields

::

    y =
        cvx dual variable (20x1 vector)

It is not necessary to place the dual variable on the left side of the
constraint; for example, the line above can also be written in this way:

::

    A * x <= b : y;

In addition, dual variables for inequality constraints will always be
nonnegative, which means that the sense of the inequality can be
reversed without changing the dual variable's value; i.e.,

::

    b >= A * x : y;

yields an identical result. For *equality* constraints, on the other
hand, swapping the left- and right- hand sides of an equality constraint
will *negate* the optimal value of the dual variable.

After the ``cvx_end`` statement is processed, and assuming the
optimization was successful, CVX assigns numerical values to ``x``
and ``y``---the optimal primal and dual variable values, respectively.
Optimal primal and dual variables for this LP must satisfy the
*complementary slackness conditions*

.. math:: 

	y_i ( b - A x )_i = 0, \quad i=1,\dots,m.

You can check this in Matlab with the line

::

    y .* (b-A*x)

which prints out the products of the entries of ``y`` and ``b-A*x``,
which should be nearly zero. This line must be executed *after* the
``cvx_end`` command (which assigns numerical values to ``x`` and ``y``);
it will generate an error if it is executed inside the CVX
specification, where ``y`` and ``b-A*x`` are still just abstract
expressions.

If the optimization is *not* successful, because either the problem is
infeasible or unbounded, then ``x`` and ``y`` will have different
values. In the unbounded case, ``x`` will contain an *unbounded
direction*; *i.e.*, a point :math:`x` satisfying

.. math:: 

	c^T x = -1, \quad A x \preceq 0,

and ``y`` will be filled with ``NaN`` values, reflecting the fact that
the dual problem is infeasible. In the infeasible case, x is filled
with ``NaN`` values, while y contains an *unbounded dual direction*;
*i.e.*, a point :math:`y` satisfying

.. math:: 

	b^T y = -1, \quad A^T y = 0, \quad y \succeq 0

Of course, the precise interpretation of primal and dual points and/or
directions depends on the structure of the problem. See references such
as `Convex Optimization <http://www.stanford.edu/~boyd/cvxbook>`_ for
more on the interpretation of dual information.

CVX also supports the declaration of *indexed* dual variables. These
prove useful when the *number* of constraints in a model (and,
therefore, the number of dual variables) depends upon the parameters
themselves. For more information on indexed dual variables, see
:ref:`indexed-dual`.

.. _assignment:

Assignment and expression holders
---------------------------------

Anyone with experience with C or Matlab understands the difference between the
single-equal *assignment* operator ``=`` and the double-equal *equality* operator ``==``.
This distinction is vitally important in CVX as well, and CVX takes steps to ensure
that assignments are not used improperly. For instance, consider the following code snippet:

::

	variable X(n,n) symmetric;
	X = semidefinite(n);

At first glance, the statement ``X = semidefinite(n);`` may look like it
constrains ``X`` to be positive semidefinite. But since the assignment operator is
used, ``X`` is actually *overwritten* by the anonymous semidefinite variable instead.
Fortunately, CVX forbids declared variables from being overwritten in this way; when
``cvx_end`` is reached, this model would issue the following error:

::

    ??? Error using ==> cvx_end
    The following cvx variable(s) have been overwritten:
       X
    This is often an indication that an equality constraint was
    written with one equals '=' instead of two '=='. The model
    must be rewritten before cvx can proceed.

We hope that this check will prevent at least some typographical errors
from having frustrating consequences in your models.

Despite this warning, assignments can be genuinely useful, so we encourage their
use with appropriate care. For instance, consider the following excerpt:

::

    variables x y
    z = 2 * x - y;
    square( z ) <= 3;
    quad_over_lin( x, z ) <= 1;

The construction ``z = 2 * x - y`` is *not* an equality constraint; it
is an assignment. It is storing an intermediate calculation
``2 * x - y``, which is an affine expression, which is then used later
in two different constraints. We call ``z`` an *expression holder* to
differentiate it from a formally declared CVX variable. 

Often it will be useful to accumulate an array of expressions into a
single Matlab variable. Unfortunately, a somewhat technical detail of
the Matlab object model can cause problems in such cases. Consider this
construction:

::

    variable u(9);
    x(1) = 1;
    for k = 1 : 9,
        x(k+1) = sqrt( x(k) + u(k) );
    end

This seems reasonable enough: ``x`` should be a vector whose first value
is ``1``, and whose subsequent values are concave CVX expressions.
But if you try this in a CVX model, Matlab will give you a rather
cryptic error:

::

    ??? The following error occurred converting from cvx to double:
    Error using ==> double
    Conversion to double from cvx is not possible.

The reason this occurs is that the Matlab variable ``x`` is initialized
as a numeric array when the assignment ``x(1)=1`` is made; and Matlab
will not permit CVX objects to be subsequently inserted into numeric
arrays.

The solution is to explicitly *declare* ``x`` to be an expression holder
before assigning values to it. We have provided keywords expression
and expressions for just this purpose, for declaring a single or
multiple expression holders for future assignment. Once an expression
holder has been declared, you may freely insert both numeric and CVX
expressions into it. For example, the previous example can be corrected
as follows:

::

    variable u(9);
    expression x(10);
    x(1) = 1;
    for k = 1 : 9,
        x(k+1) = sqrt( x(k) + u(k) );
    end

CVX will accept this construction without error. You can then use
the concave expressions ``x(1)``, ..., ``x(10)`` in any appropriate ways;
for example, you could maximize ``x(10)``.

The differences between a variable object and an expression
object are quite significant. A variable object holds an
optimization variable, and cannot be overwritten or assigned in the
CVX specification. (After solving the problem, however, CVX will
overwrite optimization variables with optimal values.) An expression
object, on the other hand, is initialized to zero, and should be thought
of as a temporary place to store CVX expressions; it can be assigned
to, freely re-assigned, and overwritten in a CVX specification.

Of course, as our first example shows, it is not always *necessary* to
declare an expression holder before it is created or used. But doing so
provides an extra measure of clarity to models, so we strongly recommend
it.
