.. _advanced:

===============
Advanced topics
===============

.. note::

	In this section we describe a number of the more advanced capabilities
	of CVX. We recommend that you *skip* this section at first, until
	you are comfortable with the basic capabilities described above.

.. _quad-forms:

Eliminating quadratic forms
---------------------------

One particular reformulation that we *strongly* encourage is to eliminate quadratic
forms---that is, functions like ``sum_square``, ``sum(square(.))`` or ``quad_form``---whenever
it is possible to construct equivalent models using ``norm`` instead.
Our experience tells us that quadratic forms often pose a numerical challenge for
the underlying solvers that CVX uses.

We acknowledge that this advice goes against conventional wisdom: quadratic forms 
are the prototypical smooth convex function, while norms are nonsmooth and therefore 
unwieldy. But with the *conic* solvers that CVX uses, this wisdom is *exactly backwards*. 
It is the *norm* that is best suited for conic formulation and solution. Quadratic forms
are handled by *converting* them to a conic form---using norms, in fact! This conversion 
process poses some interesting scaling challenges. It is better if the modeler can eliminate
the need to perform this conversion.

For a simple example of such a change, consider the objective

::

  minimize( sum_square( A * x - b ) )
  
In exact arithmetic, this is precisely equivalent to
	
::

  minimize( square_pos( norm( A * x - b ) ) )
  
But equivalence is also preserved if we eliminate the square altogether:  

::

  minimize( norm( A * x - b ) )

The optimal value of ``x`` is identical in all three cases, but this last version is
likely to produce more accurate results. Of course, if you *need* the value of the
squared norm, you can always recover it by squaring the norm after the fact.

Conversions using ``quad_form`` can sometimes be a bit more difficult. For instance, consider

::

	quad_form( A * x - b, Q ) <= 1
	
where ``Q`` is a positive definite matrix. The equivalent ``norm`` version is

::

	norm( Qsqrt * ( A * x - b ) ) <= 1

where ``Qsqrt`` is an appropriate matrix square root of ``Q``. One option is to compute
the symmetric square root ``Qsqrt = sqrtm(Q)``, but this computation destroys sparsity.
If ``Q`` is sparse, it is likely worth the effort to compute a sparse Cholesky-based
square root:

::

			[ Qsqrt, p, S  ] = chol( Q );
			Qsqrt = Qsqrt * S;
			
Sometimes an effective reformulation requires a practical understanding of what it
means for problems to be equivalent. For instance, suppose we wanted to add an
:math:`\ell_1` regularization term to the objective above, weighted by some fixed, 
positive ``lambda``:

::

  	minimize( sum_square( A * x - b ) + lambda * norm( x, 1 ) )
  	
In this case, we typically do not care about the *specific* values of ``lambda``; rather
we are varying it over a range to study the tradeoff between the residual of ``A*x-b``
and the 1-norm of ``x``. The same tradeoff can be studied by examining this modified model:

::

  	minimize( norm( A * x - b ) + lambda2 * norm( x, 1 ) )

This is not precisely the same model; setting ``lambda`` and ``lambda2`` to the same value
will not yield identical values of ``x``. But both models *do* trace the same tradeoff
curve---only the second form is likely to produce more accurate results.

.. _indexed-dual:

Indexed dual variables
----------------------

In some models, the *number* of constraints depends on the model
parameters---not just their sizes. It is straightforward to build such
models in CVX using, say, a Matlab ``for`` loop. In order to assign
each of these constraints a separate dual variable, we must find a way
to adjust the number of dual variables as well. For this reason, CVX
supports *indexed dual variables*. In reality, they are simply standard
Matlab cell arrays whose entries are CVX dual variable objects.

Let us illustrate by example how to declare and use indexed dual
variables. Consider the following semidefinite program from the
`SeDuMi <http://sedumi.ie.lehigh.edu>`_ examples:

.. math::

	\begin{array}{ll} 
		\text{minimize} & \sum_{i=1}^n (n-i) X_{ii} \\ 
		\text{subject to} & \sum_{i=1}^n X_{i,i+k} = b_k, ~ k = 1,2,\dots,n \\ 
		& X \succeq 0 
	\end{array}

This problem minimizes a weighted sum of the main diagonal of a positive
semidefinite matrix, while holding the sums along each diagonal
constant. The parameters of the problem are the elements of the vector
:math:`b\in\mathbf{R}^n`, and the optimization variable is a symmetric
matrix :math:`X\in\mathbf{R}^{n\times n}`. The CVX version of this
model is

::

    cvx_begin
        variable X( n, n ) symmetric
        minimize( ( n - 1 : -1 : 0 ) * diag( X ) );
        for k = 0 : n-1,
            sum( diag( X, k ) ) == b( k+1 );
        end
        X == semidefinite(n);
    cvx_end

If we wish to obtain dual information for the :math:`n` simple equality
constraints, we need a way to assign each constraint in the ``for`` loop
a separate dual variable. This is accomplished as follows:

::

    cvx_begin
        variable X( n, n ) symmetric
        dual variables y{n}
        minimize( ( n - 1 : -1 : 0 ) * diag( X ) );
        for k = 0 : n-1,
            sum( diag( X, k ) ) == b( k+1 ) : y{k+1};
        end
        X == semidefinite(n);
    cvx_end

The statement ``dual variables y{n}`` allocates a cell array of
:math:`n` dual variables, and stores the result in the Matlab variable
``Z``. The equality constraint in the ``for`` loop has been augmented
with a reference to ``y{k+1}``, so that each constraint is assigned a
separate dual variable. When the ``cvx_end`` command is issued, CVX
will compute the optimal values of these dual variables, and deposit
them into an :math:`n`-element cell array ``y``.

This example admittedly is a bit simplistic. With a bit of careful
arrangement, it is possible to rewrite this model so that the :math:`n`
equality constraints can be combined into a single vector constraint,
which in turn would require only a single vector dual variable. [3]_
For a more complex example that is not amenable to such a
simplification, see the file

::

    examples/cvxbook/Ch07_statistical_estim/cheb.m

in the CVX distribution. In that problem, each constraint in the
``for`` loop is a linear matrix inequality, not a scalar linear
equation; so the indexed dual variables are symmetric matrices, not
scalars.

.. _successive:

The successive approximation method
-----------------------------------

.. note::

  If you were referred to this web page by CVX's warning message: welcome! 
  Please read this section carefully to fully understand why using 
  functions like ``log``, ``exp``, etc. within CVX models requires special care.

Prior to version 1.2, the functions ``exp``, ``log``, ``log_det``,
and other functions from the exponential family could not be used within
CVX. Unfortunately, CVX utilizes symmetric primal/dual solvers that
simply cannot support those functions natively [4]_, and a variety of factors
prevent us from incorporating other types of solvers into CVX.

Nevertheless, support for these functions was requested quite frequently.
For this reason, we constructed a *successive approximation* heuristic that
allows the symmetric primal/dual solvers to support the exponential
family of functions. A precise description of the approach is beyond the
scope of this text, but roughly speaking, the method proceeds as follows:

1. Choose an initial approximation centerpoint :math:`x_c=0`.
2. Construct a polynomial approximation for each log/exp/entropy term 
   which is accurate in the neighborhood of :math:`x_c`.
3. Solve this approximate model to obtain its optimal point :math:`\bar{x}`. 
4. If :math:`\bar{x}` satisfies the optimality conditions for
   the *orignal* model to sufficient precision, exit.
5. Otherwise, shift :math:`x_c` towards :math:`\bar{x}`, and repeat steps 2-5.

Again, this is a highly simplified description of the
approach; for instance, we actually employ both the primal and dual
solutions to guide our judgements for shifting :math:`x_c` and
terminating.

This approach has proven surprisingly effective for many problems. 
*However, as with many heuristic approaches, it 
is not perfect.* It will sometimes fail to converge even for problems known to have solutions. 
Even when it does converge, it is several times slower than the standard solver,
due to its iterative approach. Therefore, it is best to use it sparingly and carefully.
Here are some specific usage tips:

- First, confirm that the log/exp/entropy terms are truly necessary for your model. In 
  many cases, an exactly equivalent model can be constructed without them, and that should
  always be preferred. For instance, the constraint
  
  ::
  
  	  sum_log(x) >= 10
  	  
  can be expressed in terms of the `geo_mean` function as
  
  ::
       
      geo_mean(x) >= log(10)^(1/length(x))
  	  
  Many determinant maximization problems are commonly written using `log_det`, but in 
  fact that is often unnecessary. For instance, consider the objective

  ::
  
      minimize( log_det(X) )
      
  CVX actually converts this internally to this:
  
  ::
  
      minimize( n*log(det_rootn(X)) )
      
  So what you can do instead is simply remove the logarithm, and solve this instead:
  
  ::
  
      minimize( det_rootn(X) )

  The value of ``log_det(X)`` can simply be computed after the model is completed.
  Unfortunately, this only 
  works if ``log_det`` is the only term in the objective; so, for instance, this
  function cannot, unfortunately, be converted in a similar fashion:
  
  ::
  
        minimize( log_det(X) + trace(C*X) )
   
- Second, try different solvers. For instance, SeDuMi and MOSEK
  tend to be more effective with the successive approximation method
  than SDPT3. So if the default solver choice fails to give a 
  solution to your model, try switching to one of these solvers.
  
- Third, try smaller instances of your problem. If they succeed where
  the larger instance fails, then at least you can confirm if the model
  is behaving as you hope before considering alternative options like
  a different solver.  
  
The bottom line, unfortunately, is that we cannot guarantee that 
the successive approximation approach will successfully handle your
specific models. If you encounter problems, you are invited to submit
a bug report, but we will not be able to promise a fix.

Suppressing the warning
~~~~~~~~~~~~~~~~~~~~~~~

Because of all of these caveats, we believe that it is necessary to
issue a warning when it is used so that users understand its
experimental nature. This warning appears the first time you 
attempt to specify a model in CVX that uses an function that
requires the successive approximation method. In fact, that warning
may very well have brought you to this section of the manual.

If you wish to suppress this warning in the future, simply issue
the command

::

    cvx_expert true

before you construct your model. If you wish to suppress this
message for all future sessions of MATLAB, follow this command
with the ``cvx_save_prefs`` command.

.. _powerfunc:

Power functions and p-norms
---------------------------

In order to implement the convex or concave branches of the power
function :math:`x^p` and :math:`p`-norms :math:`\|x\|_p` for general
values of :math:`p`, CVX uses an enhanced version of the SDP/SOCP
conversion method described by [AG00]_.
This approach is exact---as long as the exponent :math:`p` is rational.
To determine integral values :math:`p_n,p_d` such that
:math:`p_n/p_d=p`, CVX uses Matlab's ``rat`` function with its
default tolerance of :math:`10^{-6}`. There is currently no way to
change this tolerance. See the
`MATLAB documentation <http://www.mathworks.com/help/techdoc/ref/rat.html>`_ 
for the ``rat`` function for more details.

The complexity of the resulting model depends roughly on the size of the
values :math:`p_n` and :math:`p_d`. Let us introduce a more precise
measure of this complexity. For :math:`p=2`, a constraint
:math:`x^p\leq y` can be represented with exactly one :math:`2\times 2`
LMI:

.. math:: 

	x^2 \leq y \quad\Longleftrightarrow\quad \begin{bmatrix} y & x \\ x & 1 \end{bmatrix} \succeq 0.
	
For other values of :math:`p=p_n/p_d`, CVX generates a number of
:math:`2\times 2` LMIs that depends on both :math:`p_n` and :math:`p_d`;
we denote this number by :math:`k(p_n,p_d)`. (In some cases additional
linear constraints are also generated, but we ignore them for this
analysis.) For instance, for :math:`p=3/1`, we have

.. math::

   x^3\leq y, x\geq 0 \quad\Longleftrightarrow\quad \exists z ~ 
       \begin{bmatrix} z & x \\ x & 1 \end{bmatrix} \succeq 0. ~
       \begin{bmatrix} y & z \\ z & x \end{bmatrix} \succeq 0.

So :math:`k(3,1)=2`. An empirical study has shown that for
:math:`p=p_n/p_d>1`, we have

.. math:: 

	k(p_n,p_d)\leq\log_2 p_n+\alpha(p_n)

where the :math:`\alpha(p_n)` term grows very slowly compared to the
:math:`\log_2` term. Indeed, for :math:`p_n\leq 4096`, we have verified
that :math:`\alpha(p_n)` is usually 1 or 2, but occasionally 0 or 3.
Similar results are obtained for :math:`0 < p < 1` and :math:`p < 0`.

The cost of this SDP representation is relatively small for nearly all
useful values of :math:`p`. Nevertheless, CVX issues a warning
whenever :math:`k(p_n,p_d)>10` to insure that the user is not surprised
by any unexpected slowdown. In the event that this threshold does not
suit you, you may change it using the command
:samp:`cvx_power_warning({thresh})`, where :samp:`{thresh}` is the desired
cutoff value. Setting the threshold to ``Inf`` disables it completely.
As with the command ``cvx_precision``, you can place a call to
``cvx_power_warning`` within a model to change the threshold for a
single model; or outside of a model to make a global change. The command
always returns the *previous* value of the threshold, so you can save it
and restore it upon completion of your model, if you wish. You can query
the current value by calling ``cvx_power_warning`` with no arguments.

.. _overdetermined:

Overdetermined problems
-----------------------

The status message ``Overdetermined`` commonly occurs when structure
in a variable or set is not properly recognized. For example, consider
the problem of finding the smallest diagonal addition to a matrix
:math:`W\in\mathbf{R}^{n\times n}` to make it positive semidefinite:

.. math::

   \begin{array}{ll}
       \text{minimize}   & \operatorname*{\textrm{Tr}}(D) \\
       \text{subject to} & W + D \succeq 0 \\
                         & D ~ \text{diagonal}
   \end{array}

In CVX, this problem might be expressed as follows:

::

    n = size(W,1);
    cvx_begin
        variable D(n,n) diagonal;
        minimize( trace( D ) );
        subject to
            W + D == semidefinite(n);
    cvx_end

If we apply this specification to the matrix ``W=randn(5,5)``, a warning
is issued,

::

    Warning: Overdetermined equality constraints;
        problem is likely infeasible.

and the variable ``cvx_status`` is set to ``Overdetermined``.

What has happened here is that the unnamed variable returned by
statement ``semidefinite(n)`` is *symmetric*, but :math:`W` is fixed and
*unsymmetric*. Thus the problem, as stated, is infeasible. But there are
also :math:`n^2` equality constraints here, and only :math:`n+n*(n+1)/2`
unique degrees of freedom---thus the problem is overdetermined. We can
correct this problem by replacing the equality constraint with

::

            sym( W ) + D == semidefinite(n);

``sym`` is a function we have provided that extracts the symmetric part
of its argument; that is, ``sym(W)`` equals ``0.5 * ( W + W' )``.
	
.. _newfunc:

Adding new functions to the atom library
-----------------------------------------

CVX allows new convex and concave functions to be defined and added
to the atom library, in two ways, described in this section. The first
method is simple, and can (and should) be used by many users of CVX,
since it requires only a knowledge of the basic DCP ruleset. The second
method is very powerful, but a bit complicated, and should be considered
an advanced technique, to be attempted only by those who are truly
comfortable with convex analysis, disciplined convex programming, and
CVX in its current state.

Please let us know if you have implemented a convex or concave
function that you think would be useful to other users; we will be happy
to incorporate it in a future release.

New functions via the DCP ruleset
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The simplest way to construct a new function that works within CVX
is to construct it using expressions that fully conform to the DCP
ruleset. Consider, for instance, the deadzone function

.. math:: 

	f(x) = \max \{ |x|-1, 0 \} = \begin{cases} 0 & |x| \leq 1\\ x-1 & x > 1 \end{cases}

To implement this function in CVX, simply create a file
``deadzone.m`` containing

::

    function y = deadzone( x )
    y = max( abs( x ) - 1, 0 )

This function works just as you expect it would outside of
CVX --- that is, when its argument is numerical. But thanks to Matlab's
operator overloading capability, it will also work within CVX if
called with an affine argument. CVX will properly conclude that the
function is convex, because all of the operations carried out conform to
the rules of DCP: ``abs`` is recognized as a convex function; we can
subtract a constant from it, and we can take the maximum of the result
and ``0``, which yields a convex function. So we are free to use
``deadzone`` anywhere in a CVX specification that we might use
``abs``, for example, because CVX knows that it is a convex
function.

Let us emphasize that when defining a function this way, the expressions
you use *must* conform to the DCP ruleset, just as they would if they
had been inserted directly into a CVX model. For example, if we
replace ``max`` with ``min`` above; *e.g.*,

::

    function y = deadzone_bad( x )
    y = min( abs( x ) - 1, 0 )

then the modified function fails to satisfy the DCP ruleset. The function
will work *outside* of a CVX specification, happily computing the
value :math:`\min \{|x|-1,0\}` for a *numerical* argument :math:`x`. But
inside a CVX specification, invoked with a nonconstant argument, it
will generate an error.

.. _newfunc-psp:

New functions via partially specified problems
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A more advanced method for defining new functions in CVX relies on
the following basic result of convex analysis. Suppose that
:math:`S\subset\mathbf{R}^n\times\mathbf{R}^m` is a convex set and
:math:`g:(\mathbf{R}^n\times\mathbf{R}^m)\rightarrow(\mathbf{R}\cup+\infty)`
is a convex function. Then

.. math:: 

	f:\mathbf{R}^n\rightarrow(\mathbf{R}\cup+\infty), \quad f(x) \triangleq \inf\left\{\,g(x,y)\,~|~\,\exists y,~(x,y)\in S \,\right\}

is also a convex function. (This rule is sometimes called the *partial
minimization rule*.) We can think of the convex function :math:`f` as
the optimal value of a family of convex optimization problems, indexed
or parametrized by :math:`x`,

.. math::

   \begin{array}{ll}
       \mbox{minimize} & g(x,y) \\
       \mbox{subject to} & (x,y) \in S
   \end{array}

with optimization variable :math:`y`.

One special case should be very familiar: if :math:`m=1` and
:math:`g(x,y)\triangleq y`, then

.. math:: 

	f(x) \triangleq \inf\left\{\,y\,~|~\,\exists y,~(x,y)\in S\,\right\}

gives the classic *epigraph* representation of :math:`f`:

.. math:: 

	\operatorname{\textbf{epi}}f = S+ \left( \{ 0 \} \times \mathbf{R}_+ \right),

where :math:`0 \in \mathbf{R}^n`.

In CVX you can define a convex function in this very manner, that
is, as the optimal value of a parameterized family of disciplined convex
programs. We call the underlying convex program in such cases an
*incomplete specification*---so named because the parameters (that is,
the function inputs) are unknown when the specification is constructed.
The concept of incomplete specifications can at first seem a bit
complicated, but it is very powerful mechanism that allows CVX to
support a wide variety of functions.

Let us look at an example to see how this works. Consider the
unit-halfwidth Huber penalty function :math:`h(x)`:

.. math:: 

	h:\mathbf{R}\rightarrow\mathbf{R}, \quad h(x) \triangleq \begin{cases} x^2 & |x| \leq 1 \\ 2|x|-1 & |x| \geq 1 \end{cases}.

We can express the Huber function in terms of the following family of
convex QPs, parameterized by :math:`x`:

.. math::

   \begin{array}{ll}
       \text{minimize}   & 2 v + w^2 \\
       \text{subject to} & | x | \leq v + w \\
                         & w \leq 1, ~ v \geq 0
   \end{array}

with scalar variables :math:`v` and :math:`w`. The optimal value of this
simple QP is equal to the Huber penalty function of :math:`x`. We note
that the objective and constraint functions in this QP are (jointly)
convex in :math:`v`, :math:`w` *and* :math:`x`.

We can implement the Huber penalty function in CVX as follows:

::

    function cvx_optval = huber( x )
    cvx_begin
        variables w v;
        minimize( w^2 + 2 * v );
        subject to
            abs( x ) <= w + v;
            w <= 1; v >= 0;
    cvx_end

If ``huber`` is called with a numeric value of ``x``, then upon reaching
the ``cvx_end`` statement, CVX will find a complete specification,
and solve the problem to compute the result. CVX places the optimal
objective function value into the variable ``cvx_optval``, and function
returns that value as its output. Of course, it's very inefficient to
compute the Huber function of a numeric value :math:`x` by solving a QP.
But it does give the correct value (up to the core solver accuracy).

What is most important, however, is that if ``huber`` is used within a
CVX specification, with an affine CVX expression for its
argument, then CVX will do the right thing. In particular, CVX
will recognize the Huber function, called with affine argument, as a
valid convex expression. In this case, the function ``huber`` will
contain a special Matlab object that represents the function call in
constraints and objectives. Thus the function ``huber`` can be used
anywhere a traditional convex function can be used, in constraints or
objective functions, in accordance with the DCP ruleset.

There is a corresponding development for concave functions as well.
Given a convex set :math:`S` as above, and a concave function
:math:`g:(\mathbf{R}^n\times\mathbf{R}^m)\rightarrow(\mathbf{R}\cup-\infty)`,
the function

.. math:: 

	f:\mathbf{R}\rightarrow(\mathbf{R}\cup-\infty), \quad f(x) \triangleq \sup\left\{\,g(x,y)\,~|~\,\exists y,~(x,y)\in S \,\right\}

is concave. If :math:`g(x,y)\triangleq y`, then

.. math:: 

	f(x) \triangleq \sup\left\{\,y\,~|~\,\exists y,~(x,y)\in S\,\right\}

gives the *hypograph* representation of :math:`f`:

.. math:: 

	\operatorname{\textbf{hypo}}f = S - \mathbf{R}_+^n.

In CVX, a concave incomplete specification is simply one that uses a
``maximize`` objective instead of a ``minimize`` objective; and if
properly constructed, it can be used anywhere a traditional concave
function can be used within a CVX specification.

For an example of a concave incomplete specification, consider the
function

.. math:: 

	f:\mathbf{R}^{n\times n}\rightarrow\mathbf{R}, \quad f(X) = \lambda_{\min}(X+X^T)

Its hypograph can be represented using a single linear matrix
inequality:

.. math:: 

	\operatorname{\textbf{hypo}}f = \left\{\, (X,t) \,~|~\, f(X) \geq t \,\right\} = \left\{\, (X,t) \,~|~\, X + X^T - t I \succeq 0 \,\right\}

So we can implement this function in CVX as follows:

::

    function cvx_optval = lambda_min_symm( X )
    n = size( X, 1 );
    cvx_begin
        variable y;
        maximize( y );
        subject to
            X + X' - y * eye( n ) == semidefinite( n );
    cvx_end

If a numeric value of ``X`` is supplied, this function will return
``min(eig(X+X'))`` (to within numerical tolerances). However, this
function can also be used in CVX constraints and objectives, just
like any other concave function in the atom library.

There are two practical issues that arise when defining functions using
incomplete specifications, both of which we will illustrate using our
``huber`` example above. First of all, as written the function works
only with scalar values. To apply it (elementwise) to a vector requires
that we iterate through the elements in a ``for`` loop---a *very*
inefficient enterprise, particularly in CVX. A far better approach
is to extend the ``huber`` function to handle vector inputs. This is, in
fact, rather simple to do: we simply create a *multiobjective* version
of the problem:

::

    function cvx_optval = huber( x )
    sx = size( x );
    cvx_begin
        variables w( sx ) v( sx );
        minimize( w .^ 2 + 2 * v );
        subject to
            abs( x ) <= w + v;
            w <= 1; v >= 0;
    cvx_end

This version of ``huber`` will in effect create ``sx`` "instances" of
the problem in parallel; and when used in a CVX specification, will
be handled correctly.

The second issue is that if the input to ``huber`` is numeric, then
direct computation is a far more efficient way to compute the result
than solving a QP. (What is more, the multiobjective version cannot be
used with numeric inputs.) One solution is to place both versions in one
file, with an appropriate test to select the proper version to use:

::

    function cvx_optval = huber( x )
    if isnumeric( x ),
        xa   = abs( x );
        flag = xa < 1;
        cvx_optval = flag .* xa.^2 + (~flag) * (2*xa-1);
    else,
        sx = size( x );
        cvx_begin
            variables w( sx ) v( sx );
            minimize( w .^ 2 + 2 * v );
            subject to
                abs( x ) <= w + v;
                w <= 1; v >= 0;
        cvx_end
    end

Alternatively, you can create two separate versions of the function, one
for numeric input and one for CVX expressions, and place the CVX
version in a subdirectory called ``@cvx``. (Do not include this
directory in your Matlab ``path``; only include its parent.) Matlab will
automatically call the version in the ``@cvx`` directory when one of the
arguments is a CVX variable. This is the approach taken for the
version of ``huber`` found in the CVX atom library.

One good way to learn more about using incomplete specifications is to
examine some of the examples already in the CVX atom library. Good
choices include ``huber``, ``inv_pos``, ``lambda_min``, ``lambda_max``,
``matrix_frac``, ``quad_over_lin``, ``sum_largest``, and others. Some
are a bit difficult to read because of diagnostic or error-checking
code, but these are relatively simple.

.. [3]
   Indeed, a future version of CVX will support the use of the
   Matlab function ``spdiags``, which will reduce the entire for loop to
   the single constraint ``spdiags(X,0:n-1)==b``.
   
.. [4]
   Technically there are a couple of exceptions here. First of all, 
   SDPT3 does, in fact, support the existence of logarithms and ``log_det``
   terms in the objective function. However, it doesn't support such terms
   within constraints. Unfortunately, because CVX does not differentiate
   between objective terms and constraint terms internally, it is not able
   to utilize this capability of SDPT3. Secondly, this section was written
   before the inclusion of MOSEK support in CVX, and CVX does indeed provide
   support for smooth nonlinearities in its solver. But this capability
   is not easy to use in MATLAB.
   
   