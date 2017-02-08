.. _dcp:

===============
The DCP ruleset
===============

CVX enforces the conventions dictated by the disciplined convex
programming ruleset, or *DCP ruleset* for short. CVX will issue an
error message whenever it encounters a violation of any of the rules, so
it is important to understand them before beginning to build models. The
rules are drawn from basic principles of convex analysis, and are easy
to learn, once you've had an exposure to convex analysis and convex
optimization.

The DCP ruleset is a set of sufficient, but not necessary, conditions
for convexity. So it is possible to construct expressions that violate
the ruleset but are in fact convex. As an example consider the entropy
function, :math:`-\sum_{i=1}^n x_i \log x_i`, defined for :math:`x>0`,
which is concave. If it is expressed as

::

    - sum( x .* log( x ) )

CVX will reject it, because its concavity does not follow from any
of the composition rules. (Specifically, it violates the no-product rule
described in :ref:`expressions`.) Problems involving entropy,
however, can be solved, by explicitly using the entropy function,

::

    sum( entr( x ) )

which is in the base CVX library, and thus recognized as concave by
CVX. If a convex (or concave) function is not recognized as convex
or concave by CVX, it can be added as a new atom; see
:ref:`newfunc`.

As another example consider the function
:math:`\sqrt{x^2+1}=\|[x~1]\|_2`, which is convex. If it is written as

::

    norm( [ x 1 ] )

(assuming ``x`` is a scalar variable or affine expression) it will be
recognized by CVX as a convex expression, and therefore can be used
in (appropriate) constraints and objectives. But if it is written as

::

    sqrt( x^2 + 1 )

CVX will reject it, since convexity of this function does not follow
from the CVX ruleset.

A taxonomy of curvature
-----------------------

In disciplined convex programming, a scalar expression is classified by
its *curvature*. There are four categories of curvature: *constant*,
*affine*, *convex*, and *concave*. For a function
:math:`f:\mathbf{R}^n\rightarrow\mathbf{R}` defined on all
:math:`\mathbf{R}^n`, the categories have the following meanings:

.. math::

	\begin{array}{lll}
		\text{constant} & f(\alpha x + (1-\alpha)y) = f(x)                             & \forall x,y\in\mathbf{R}^n,~\alpha\in\mathbf{R} \\
		\text{affine}   & f(\alpha x + (1-\alpha)y) = \alpha f(x) + (1-\alpha) f(y)    & \forall x,y\in\mathbf{R}^n,~\alpha\in\mathbf{R} \\
		\text{convex}   & f(\alpha x + (1-\alpha)y) \leq \alpha f(x) + (1-\alpha) f(y) & \forall x,y\in\mathbf{R}^n,~\alpha\in[0,1] \\
		\text{concave}  & f(\alpha x + (1-\alpha)y) \geq \alpha f(x) + (1-\alpha) f(y) & \forall x,y\in\mathbf{R}^n,~\alpha\in[0,1]
	\end{array}		

Of course, there is significant overlap in these categories. For
example, constant expressions are also affine, and (real) affine
expressions are both convex and concave.

Convex and concave expressions are real by definition. Complex constant
and affine expressions can be constructed, but their usage is more
limited; for example, they cannot appear as the left- or right-hand side
of an inequality constraint.

Top-level rules
---------------

CVX supports three different types of disciplined convex programs:

-  A *minimization problem*, consisting of a convex objective function
   and zero or more constraints.
-  A *maximization problem*, consisting of a concave objective function
   and zero or more constraints.
-  A *feasibility problem*, consisting of one or more constraints and no
   objective.

Constraints
-----------

Three types of constraints may be specified in disciplined convex
programs:

-  An *equality constraint*, constructed using ``==``, where both sides
   are affine.
-  A *less-than inequality constraint*, using ``<=``, where the left
   side is convex and the right side is concave.
-  A *greater-than inequality constraint*, using ``>=``, where the left
   side is concave and the right side is convex.

*Non*-equality constraints, constructed using ``~=``, are never allowed.
(Such constraints are not convex.)

One or both sides of an equality constraint may be complex; inequality
constraints, on the other hand, must be real. A complex equality
constraint is equivalent to two real equality constraints, one for the
real part and one for the imaginary part. An equality constraint with a
real side and a complex side has the effect of constraining the
imaginary part of the complex side to be zero.

As discussed in :ref:`sets`, CVX enforces set membership
constraints (*e.g.*, :math:`x\in S`) using equality constraints. The
rule that both sides of an equality constraint must be affine applies to
set membership constraints as well. In fact, the returned value of set
atoms like ``semidefinite()`` and ``lorentz()`` is affine, so it is
sufficient to simply verify the remaining portion of the set membership
constraint. For composite values like ``{ x, y }``, each element must be
affine.

.. _strict:

Strict inequalities
~~~~~~~~~~~~~~~~~~~

As mentioned in :ref:`constraints`, strict inequalities ``<``, ``>`` are interpreted 
in an identical fashion to nonstrict inequalities ``>=``, ``<=``. It is important to 
note that CVX cannot guarantee that an inequality will be strictly satisfied
at the solution it computes. This is not simply a choice we have made in CVX; it is
a natural consequence of both the underlying mathematics and the
design of convex optimization solvers.
For that reason, we *strongly* discourage the use of strict inequalities in CVX, 
and a future version may remove them altogether.

When a strict inequality is essential to your model, you may need to take additional
steps to ensure compliance. In some cases, this can be accomplished through 
*normalization*. For instance, consider a set of homogeneous equations and inequalities:

.. math::

	A x = 0, \quad C x \preceq 0, \quad x \succ 0
	
Except for the strict inequality, :math:`x=0` would be an acceptable solution; indeed
the need to avoid the origin is the very reason for the strict inequality. However, note
that if a given :math:`x` satisfies these constraints, then so does 
:math:`\alpha x` for all :math:`\alpha>0`. By eliminating this degree of freedom with
normalization, we can eliminate the strict inequality; for instance:

.. math::

	A x = 0, \quad C x \preceq 0, \quad x \succ 0, \quad \mathbf{1}^T x = 1
	
If normalization is not a valid approach for your model, you may simply need to convert
the strict inequality into a non-strict one by adding a small offset; *e.g.*, convert
``x > 0`` to, say, ``x >= 1e-4``. Note that the bound needs to be large enough so
that the underlying solver considers it numerically significant.

Finally, note that for some functions like ``log(x)`` and ``inv_pos(x)``, which have domains
defined by strict inequalities, the domain restriction is handled *by the function itself*.
You do not need to add an explicit constraint ``x > 0`` to your model to guarantee
that the solution is positive.

.. _expressions:

Expression rules
----------------

So far, the rules as stated are not particularly restrictive, in that
all convex programs (disciplined or otherwise) typically adhere to them.
What distinguishes disciplined convex programming from more general
convex programming are the rules governing the construction of the
expressions used in objective functions and constraints.

Disciplined convex programming determines the curvature of scalar
expressions by recursively applying the following rules. While this list
may seem long, it is for the most part an enumeration of basic rules of
convex analysis for combining convex, concave, and affine forms: sums,
multiplication by scalars, and so forth.

-  A valid constant expression is

   -  any well-formed Matlab expression that evaluates to a finite
      value.

-  A valid affine expression is

   -  a valid constant expression;
   -  a declared variable;
   -  a valid call to a function in the atom library with an affine
      result;
   -  the sum or difference of affine expressions;
   -  the product of an affine expression and a constant.

-  A valid convex expression is

   -  a valid constant or affine expression;
   -  a valid call to a function in the atom library with a convex
      result;
   -  an affine scalar raised to a constant power :math:`p\geq 1`,
      :math:`p\neq3,5,7,9,...`;
   -  a convex scalar quadratic form---see
      :ref:`quadforms`;
   -  the sum of two or more convex expressions;
   -  the difference between a convex expression and a concave
      expression;
   -  the product of a convex expression and a nonnegative constant;
   -  the product of a concave expression and a nonpositive constant;
   -  the negation of a concave expression.

-  A valid concave expression is

   -  a valid constant or affine expression;
   -  a valid call to a function in the atom library with a concave
      result;
   -  a concave scalar raised to a power :math:`p\in(0,1)`;
   -  a concave scalar quadratic form---see
      :ref:`quadforms`;
   -  the sum of two or more concave expressions;
   -  the difference between a concave expression and a convex
      expression;
   -  the product of a concave expression and a nonnegative constant;
   -  the product of a convex expression and a nonpositive constant;
   -  the negation of a convex expression.

If an expression cannot be categorized by this ruleset, it is rejected
by CVX. For matrix and array expressions, these rules are applied on
an elementwise basis. We note that the set of rules listed above is
redundant; there are much smaller, equivalent sets of rules.

Of particular note is that these expression rules generally forbid
*products* between nonconstant expressions, with the exception of scalar
quadratic forms. For
example, the expression ``x*sqrt(x)`` happens to be a convex function of
``x``, but its convexity cannot be verified using the CVX ruleset,
and so is rejected. (It can be expressed as ``pow_p(x,3/2)``, however.) 
We call this the *no-product rule*, and paying close attention to it will 
go a long way to insuring that the
expressions you construct are valid.

Functions
---------

In CVX, functions are categorized in two attributes: *curvature*
(*constant*, *affine*, *convex*, or *concave*) and *monotonicity*
(*nondecreasing*, *nonincreasing*, or *nonmonotonic*). Curvature
determines the conditions under which they can appear in expressions
according to the expression rules given above. Monotonicity determines
how they can be used in function compositions, as we shall see in the
next section.

For functions with only one argument, the categorization is
straightforward. Some examples are given in the table below.

.. tabularcolumns:: CCCC

===============   ====================   ===========  ===============
 Function          Meaning                Curvature    Monotonicity   
===============   ====================   ===========  ===============
 ``sum( x )``      :math:`\sum_i x_i`     affine       nondecreasing  
 ``abs( x )``      :math:`|x|`            convex       nonmonotonic   
 ``log( x )``      :math:`\log x`         concave      nondecreasing   
 ``sqrt( x )``     :math:`\sqrt x`        concave      nondecreasing  
===============   ====================   ===========  ===============

Following standard practice in convex analysis, convex functions are
interpreted as :math:`+\infty` when the argument is outside the domain
of the function, and concave functions are interpreted as
:math:`-\infty` when the argument is outside its domain. In other words,
convex and concave functions in CVX are interpreted as their
*extended-valued extensions*.

This has the effect of automatically constraining the argument of a
function to be in the function's domain. For example, if we form
``sqrt(x+1)`` in a CVX specification, where ``x`` is a variable,
then ``x`` will automatically be constrained to be larger than or equal
to :math:`-1`. There is no need to add a separate constraint, ``x>=-1``,
to enforce this.

Monotonicity of a function is determined in the extended sense, *i.e.*,
*including the values of the argument outside its domain*. For example,
``sqrt(x)`` is determined to be nondecreasing since its value is
constant (:math:`-\infty`) for negative values of its argument; then
jumps *up* to :math:`0` for argument zero, and increases for positive
values of its argument.

CVX does *not* consider a function to be convex or concave if it is
so only over a portion of its domain, even if the argument is
constrained to lie in one of these portions. As an example, consider the
function :math:`1/x`. This function is convex for :math:`x>0`, and
concave for :math:`x<0`. But you can never write ``1/x`` in CVX
(unless ``x`` is constant), even if you have imposed a constraint such
as ``x>=1``, which restricts ``x`` to lie in the convex portion of
function :math:`1/x`. You can use the CVX function ``inv_pos(x)``,
defined as :math:`1/x` for :math:`x>0` and :math:`\infty` otherwise, for
the convex portion of :math:`1/x`; CVX recognizes this function as
convex and nonincreasing. In CVX, you can express the concave
portion of :math:`1/x`, where :math:`x` is negative, using
``-inv_pos(-x)``, which will be correctly recognized as concave and
nonincreasing.

For functions with multiple arguments, curvature is always considered
*jointly*, but monotonicity can be considered on an
*argument-by-argument* basis. For example, the function
``quad_over_lin(x,y)``

.. math:: 

	f_{\text{quad\_over\_lin}}(x,y) = \begin{cases} |x|^2/y & y > 0 \\ +\infty & y\leq 0  \end{cases}

is jointly convex in both :math:`x` and :math:`y`, but it is monotonic
(nonincreasing) only in :math:`y`.

Some functions are convex, concave, or affine only for a *subset* of its
arguments. For example, the function ``norm(x,p)`` where ``p \geq 1`` is
convex only in its first argument. Whenever this function is used in a
CVX specification, then, the remaining arguments must be constant,
or CVX will issue an error message. Such arguments correspond to a
function's parameters in mathematical terminology; *e.g.*,

.. math:: 

	f_p(x):\mathbf{R}^n\rightarrow\mathbf{R}, \quad f_p(x) \triangleq \|x\|_p

So it seems fitting that we should refer to such arguments as
*parameters* in this context as well. Henceforth, whenever we speak of a
CVX function as being convex, concave, or affine, we will assume
that its parameters are known and have been given appropriate, constant
values.

Compositions
------------

A basic rule of convex analysis is that convexity is closed under
composition with an affine mapping. This is part of the DCP ruleset as
well:

-  A convex, concave, or affine function may accept an affine expression
   (of compatible size) as an argument. The result is convex, concave,
   or affine, respectively.

For example, consider the function ``square(x)``, which is provided in
the CVX atom library. This function squares its argument; *i.e.*, it
computes ``x.*x``. (For array arguments, it squares each element
independently.) It is in the CVX atom library, and known to be
convex, provided its argument is real. So if ``x`` is a real variable of
dimension :math:`n`, ``a`` is a constant :math:`n`-vector, and ``b`` is
a constant, the expression

::

    square( a' * x + b )

is accepted by CVX, which knows that it is convex.

The affine composition rule above is a special case of a more
sophisticated composition rule, which we describe now. We consider a
function, of known curvature and monotonicity, that accepts multiple
arguments. For *convex* functions, the rules are:

-  If the function is nondecreasing in an argument, that argument must
   be convex.
-  If the function is nonincreasing in an argument, that argument must
   be concave.
-  If the function is neither nondecreasing or nonincreasing in an
   argument, that argument must be affine.

If each argument of the function satisfies these rules, then the
expression is accepted by CVX, and is classified as convex. Recall
that a constant or affine expression is both convex and concave, so any
argument can be affine, including as a special case, constant.

The corresponding rules for a concave function are as follows:

-  If the function is nondecreasing in an argument, that argument must
   be concave.
-  If the function is nonincreasing in an argument, that argument must
   be convex.
-  If the function is neither nondecreasing or nonincreasing in an
   argument, that argument must be affine.

In this case, the expression is accepted by CVX, and classified as
concave.

For more background on these composition rules, see `Convex
Optimization <http://www.stanford.edu/~boyd/cvxbook>`_, Section 3.2.4.
In fact, with the exception of scalar quadratic expressions, the entire
DCP ruleset can be thought of as special cases of these six rules.

Let us examine some examples. The maximum function is convex and
nondecreasing in every argument, so it can accept any convex expressions
as arguments. For example, if ``x`` is a vector variable, then

::

    max( abs( x ) )

obeys the first of the six composition rules and is therefore accepted
by CVX, and classified as convex. As another example, consider the
sum function, which is both convex and concave (since it is affine), and
nondecreasing in each argument. Therefore the expressions

::

    sum( square( x ) )
    sum( sqrt( x ) )

are recognized as valid in CVX, and classified as convex and
concave, respectively. The first one follows from the first rule for
convex functions; and the second one follows from the first rule for
concave functions.

Most people who know basic convex analysis like to think of these
examples in terms of the more specific rules: a maximum of convex
functions is convex, and a sum of convex (concave) functions is convex
(concave). But these rules are just special cases of the general
composition rules above. Some other well known basic rules that follow
from the general composition rules are:

-  a nonnegative multiple of a convex (concave) function is convex
   (concave);
-  a nonpositive multiple of a convex (concave) function is concave
   (convex).

Now we consider a more complex example in depth. Suppose ``x`` is a
vector variable, and ``A``, ``b``, and ``f`` are constants with
appropriate dimensions. CVX recognizes the expression

::

    sqrt(f'*x) + min(4,1.3-norm(A*x-b))

as concave. Consider the term ``sqrt(f'*x)``. CVX recognizes that
``sqrt`` is concave and ``f'*x`` is affine, so it concludes that
``sqrt(f'*x)`` is concave. Now consider the second term
``min(4,1.3-norm(A*x-b))``. CVX recognizes that ``min`` is concave
and nondecreasing, so it can accept concave arguments. CVX
recognizes that ``1.3-norm(A*x-b)`` is concave, since it is the
difference of a constant and a convex function. So CVX concludes
that the second term is also concave. The whole expression is then
recognized as concave, since it is the sum of two concave functions.

The composition rules are sufficient but not necessary for the
classification to be correct, so some expressions which are in fact
convex or concave will fail to satisfy them, and so will be rejected by
CVX. For example, if ``x`` is a vector variable, the expression

::

        sqrt( sum( square( x ) ) )

is rejected by CVX, because there is no rule governing the
composition of a concave nondecreasing function with a convex function.
Of course, the workaround is simple in this case: use ``norm( x )``
instead, since ``norm`` is in the atom library and known by CVX to
be convex.

Monotonicity in nonlinear compositions
--------------------------------------

Monotonicity is a critical aspect of the rules for nonlinear
compositions. This has some consequences that are not so obvious, as we
shall demonstrate here by example. Consider the expression

::

    square( square( x ) + 1 )

where ``x`` is a scalar variable. This expression is in fact convex,
since :math:`(x^2+1)^2 = x^4+2x^2+1` is convex. But CVX will reject
the expression, because the outer ``square`` cannot accept a convex
argument. Indeed, the square of a convex function is not, in general,
convex: for example, :math:`(x^2-1)^2 = x^4-2x^2+1` is not convex.

There are several ways to modify the expression above to comply with the
ruleset. One way is to write it as ``x^4 + 2*x^2 + 1``, which CVX
recognizes as convex, since CVX allows positive even integer powers
using the ``^`` operator. (Note that the same technique, applied to the
function :math:`(x^2-1)^2`, will fail, since its second term is
concave.)

Another approach is to use the alternate outer function ``square_pos``,
included in the CVX library, which represents the function
:math:`(x_+)^2`, where :math:`x_+ = \max\{0,x\}`. Obviously, ``square``
and ``square_pos`` coincide when their arguments are nonnegative. But
``square_pos`` is nondecreasing, so it can accept a convex argument.
Thus, the expression

::

    square_pos( square( x ) + 1 )

is mathematically equivalent to the rejected version above (since the
argument to the outer function is always positive), but it satisfies the
DCP ruleset and is therefore accepted by CVX.

This is the reason several functions in the CVX atom library come in
two forms: the "natural" form, and one that is modified in such a way
that it is monotonic, and can therefore be used in compositions. Other
such "monotonic extensions" include ``sum_square_pos`` and
``quad_pos_over_lin``. If you are implementing a new function yourself,
you might wish to consider if a monotonic extension of that function
would also be useful.

.. _quadforms:

Scalar quadratic forms
----------------------

In its pure form, the DCP ruleset forbids even the use of simple
quadratic expressions such as ``x * x`` (assuming ``x`` is a scalar
variable). For practical reasons, we have chosen to make an exception to
the ruleset to allow for the recognition of certain specific quadratic
forms that map directly to certain convex quadratic functions (or their
concave negatives) in the CVX atom library:

=====================   =============================
``x .* x``              ``square( x )`` (real ``x``)
``conj( x ) .* x``      ``square_abs( x )``                
``y' * y``              ``sum_square_abs( y )``            
``(A*x-b)'*Q*(Ax-b)``   ``quad_form( A*x - b, Q )`` 
=====================   =============================

CVX detects the quadratic expressions such as those on the left
above, and determines whether or not they are convex or concave; and if
so, translates them to an equivalent function call, such as those on the
right above.

CVX examines each *single* product of affine expressions, and each
*single* squaring of an affine expression, checking for convexity; it
will not check, for example, sums of products of affine expressions. For
example, given scalar variables ``x`` and ``y``, the expression

::

    x ^ 2 + 2 * x * y + y ^2

will cause an error in CVX, because the second of the three terms
``2 * x * y``, is neither convex nor concave. But the equivalent
expressions

::

    ( x + y ) ^ 2
    ( x + y ) * ( x + y )

will be accepted. 

CVX actually completes the square when it comes
across a scalar quadratic form, so the form need not be symmetric. For
example, if ``z`` is a vector variable, ``a``, ``b`` are constants, and
``Q`` is positive definite, then

::

    ( z + a )' * Q * ( z + b )

will be recognized as convex. Once a quadratic form has been verified by
CVX, it can be freely used in any way that a normal convex or
concave expression can be, as described in :ref:`expressions`.

Quadratic forms should actually be used *less frequently* in disciplined
convex programming than in a more traditional mathematical programming
framework, where a quadratic form is often a smooth substitute for a
nonsmooth form that one truly wishes to use. In CVX, such
substitutions are rarely necessary, because of its support for nonsmooth
functions. For example, the constraint

::

    sum( ( A * x - b ) .^ 2 ) <= 1

is equivalently represented using the Euclidean norm:

::

    norm( A * x - b ) <= 1

With modern solvers, the second form is more naturally represented using
a second-order cone constraint---so the second form may actually be more
efficient. In fact, in our experience, the non-squared form will often
be handled more accurately. So we strongly encourage you to re-evaluate
the use of quadratic forms in your models, in light of the new
capabilities afforded by disciplined convex programming.

