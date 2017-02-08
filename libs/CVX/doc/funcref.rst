.. _funcref:

===============
Reference guide
===============

In this section we describe each operator, function, set, and command that you are 
likely to encounter in CVX. In some cases, limitations of the underlying solver
place certain restrictions or caveats on their use:

-  Functions marked with a dagger (†) are not supported natively by the
   solvers that CVX uses. They are handled using a successive
   approximation method which makes multiple calls to the underlying
   solver, achieving the same final precision. If you use one of these
   functions, you will be warned that successive approximation will be
   used. This technique is discussed further in
   :ref:`successive`. As this section discusses, this is an experimental
   approach that works well in many cases, but cannot be guaranteed.

-  Functions involving powers (*e.g.*, ``x^p``) and :math:`p`-norms
   (*e.g.*, ``norm(x,p)``) are marked with a double dagger (‡). CVX
   represents these functions exactly when :math:`p` is a rational
   number. For irrational values of ``p``, a nearby rational is selected
   instead. See :ref:`powerfunc` for details on
   how both cases are handled.

Arithmetic operators
--------------------

Matlab's standard arithmetic operations for addition ``+``, subtraction ``-``, 
multiplication, ``*`` ``.*``, division ``/`` ``./`` ``\`` ``.\``, and 
exponentiation ``^`` ``.^`` have been overloaded to work in
CVX whenever appropriate---that is, whenever their use is consistent
with both standard mathematical and Matlab conventions *and* the DCP
ruleset. For example:

-  Two CVX expressions can be added together if they are of the same
   dimension (or one is scalar) and have the same curvature (*i.e.*,
   both are convex, concave, or affine).

-  A CVX expression can be multiplied or divided by a scalar
   constant. If the constant is positive, the curvature is preserved; if
   it is negative, curvature is reversed.

-  An affine column vector CVX expression can be multiplied by a
   constant matrix of appropriate dimensions; or it can be left-divided
   by a non-singular constant matrix of appropriate dimension.

Numerous other combinations are possible, of course. The use of the exponentiation 
operators ``^`` ``.^`` are somewhat limited;
see the definitions of ``power`` in :ref:`nonlinear` below.

Matlab's basic matrix manipulation and arithmetic operations have been
extended to work with CVX expressions as well, including:

-  Concatenation: ``[ A, B ; C, D ]``
-  Indexing: ``x(n+1:end)``, ``X([3,4],:)``, *etc.*
-  Indexed assignment, including deletion: ``y(2:4) = 1``,
   ``Z(1:4,:) = []``, *etc.*
-  Transpose and conjugate transpose: ``Z.'``, ``y'``

.. _builtin:

Built-in functions
-------------------

Linear
~~~~~~

A number of Matlab's basic linear and bilinear functions either work automatically
with ``cvx`` expressions or have been extended to do so, including:
``conj``, ``conv``, ``cumsum``, ``diag``, ``dot``,
``find``, ``fliplr``, ``flipud``, ``flipdim``,
``horzcat``, ``hankel``, ``ipermute``, ``kron``, ``mean``,
``permute``, ``repmat``, ``reshape``, ``rot90``, 
``sparse``, ``sum``, ``trace``, ``tril``, ``triu``,   
``toeplitz``, ``vertcat``.

Most should behave identically with CVX expressions as they do with
numeric expressions. Those that perform some sort of summation, such as
``cumsum``, ``sum``, or multiplication, such as ``conv``, ``dot`` or
``kron``, can only be used in accordance with the disciplined convex
programming rules. For example, ``kron(X,Y)`` is valid only if either
``X`` or ``Y`` is constant; and ``trace(Z)`` is valid only if the
elements along the diagonal have the same curvature.

.. _nonlinear:

Nonlinear
~~~~~~~~~

``abs``
    absolute value for real and complex arrays. Convex.

† ``exp``
    exponential. Convex and nondecreasing.

† ``log``
    logarithm. Concave and nondecreasing.

``max``
    maximum. Convex and nondecreasing.

``min``
    minimum. Concave and nondecreasing.

``norm``
    norms for real and complex vectors and matrices. Convex. Thus
    function follows the Matlab conventions closely. Thus the
    one-argument version ``norm(x)`` computes the 2-norm for vectors,
    and the 2-norm (maximum singular value) for matrices. The
    two-argument version ``norm(x,p)`` is supported as follows:

    -  ‡ For vectors, all values :math:`p\geq 1` are accepted.
    -  For matrices, ``p`` must be ``1``, ``2``, ``Inf``, or ``'Fro'``.

``polyval``
    polynomial evaluation. ``polyval(p,x)``, where ``p`` is a vector of
    length ``n``, computes

    ::

            p(1) * x.^(n-1) + p(2) * x.^(n-2) + ... + p(n-1) * x + p(n)

    This function can be used in CVX in two ways:

    -  If ``p`` is a variable and ``x`` is a constant, then
       ``polyval(x,p)`` computes a linear combination of the elements of
       ``p``. The combination must satisfy the DCP rules for addition
       and scaling.
    -  If ``p`` is a constant and ``x`` is a variable, then
       ``polyval(x,p)`` constructs a polynomial function of the variable
       ``x``. The polynomial must be affine, convex, or concave, and
       ``x`` must be real and affine.
       
‡ ``power(x,p)``
    ``x^p`` and ``x.^p``, where ``x`` is a real variable and and ``p``
    is a real constant. For ``x^p``, both ``x`` and ``p`` must be
    scalars. Only those values of ``p`` which can reasonably and
    unambiguously interpreted as convex or concave are accepted:

    -  :math:`p=0`. Constant. ``x.^p`` is treated as identically 1.
    -  :math:`0 < p < 1`. Concave. The argument :math:`x` must be
       concave (or affine), and is implicitly constrained to be
       nonnegative.
    -  :math:`p = 1`. Affine. ``x.^p`` is simply ``x``.
    -  :math:`p \in \{2,4,6,8,...\}`. Convex. Argument :math:`x` must be
       affine.
    -  :math:`p > 1`, :math:`p\not\in\{2,3,4,5,...\}`. Convex. Argument
       :math:`x` must be affine, and is implicitly constrained to be
       nonnegative.

    Negative and odd integral values of :math:`p` are not permitted, but
    see the functions ``pow_p``, ``pow_pos``, and ``pow_abs`` in the
    next section for useful alternatives.

† ``power(p,x)``
    ``p.^x`` and ``p^x``, where ``p`` is a real constant and ``x`` is a
    real variable. For ``p^x``, both ``p`` and ``x`` must be scalars.
    Valid values of ``p`` include:

    -  :math:`p \in \{0,1\}`. Constant.
    -  :math:`0 < p < 1`. Convex and nonincreasing; ``x`` must be
       concave.
    -  :math:`p > 1`. Convex and nondecreasing; ``x`` must be convex.

    Negative values of ``p`` are not permitted.

``std``
    Standard deviation. Convex.

``sqrt``
    Square root. Implicitly constrains its argument to be nonnegative.
    Concave and nondecreasing.

``var``
    Variance. Convex.

.. _newfuncs:

New functions
--------------

Even though these functions were developed specifically for CVX,
they work outside of a CVX specification as well, when supplied with
numeric arguments.

``avg_abs_dev``
    The average absolute deviation about the mean :math:`\mu(x)` of :math:`x`. Convex.
	.. math::
	
		f_{\text{aad}}(x) = \frac{1}{n} \sum_{i=1}^n |x_i-\mu(x)| = \frac{1}{n} \sum_{i=1}^n \left| x_i - {\textstyle\frac{1}{n}\sum_i x_i}\right| = \frac{1}{n}\left\| (I-\tfrac{1}{n}\textbf{1}\textbf{1}^T)x \right\|_1.
		
``avg_abs_dev_med``
    The average absolute deviation about the median :math:`\mathop{\textrm{m}}(x)` of :math:`x`. Convex.
	.. math::
	
		f_{\text{aadm}}(x) = \frac{1}{n} \sum_{i=1}^n |x_i-\mathop{\textrm{m}}(x)| = \inf_y \frac{1}{n} \sum_{i=1}^n |x_i-y|
		
``berhu(x,M)``
    The reversed Huber function (hence, Berhu), defined as
	.. math:: 

		f_{\text{berhu}}(x,M) \triangleq \begin{cases} |x| & |x| \leq M \\ (|x|^2+M^2)/2M & |x| \geq M \end{cases}

    Convex. If :math:`M` is omitted, :math:`M=1` is assumed; but if supplied, it must be a positive constant.
    Also callable with three arguments as ``berhu(x,M,t)``, which computes ``t+t*berhu(x/t,M)``, 
    useful for concomitant scale estimation (see [Owen06]_).

``det_inv``
    determinant of inverse of a symmetric (or Hermitian) positive
    definite matrix, :math:`\det X^{-1}`, which is the same as the
    product of the inverses of the eigenvalues. When used inside a
    CVX specification, ``det_inv`` constrains the matrix to be
    symmetric (if real) or Hermitian (if complex) and positive
    semidefinite. When used with numerical arguments, ``det_inv``
    returns ``+Inf`` if these constraints are not met. Convex.

``det_rootn``
    :math:`n`-th root of the determinant of a semidefinite matrix,
    :math:`(\det X)^{1/n}`. When used inside a CVX specification,
    ``det_rootn`` constrains the matrix to be symmetric (if real) or
    Hermitian (if complex) and positive semidefinite. When used with
    numerical arguments, ``det_rootn`` returns ``-Inf`` if these
    constraints are not met. Concave.

``det_root2n``
    the :math:`2n`-th root of the determinant of a semidefinite matrix;
    *i.e.*, ``det_root2n(X)=sqrt(det_rootn(X))``. Concave. Maintained
    solely for back-compatibility purposes.

† ``entr``
    the elementwise entropy function: ``entr(x)=-x.*log(x)``. Concave.
    Returns ``-Inf`` when called with a constant argument that has a
    negative entry.

``geo_mean``
    the geometric mean of a vector,
    :math:`\left( \prod_{k=1}^n x_k \right)^{1/n}`. When used inside a
    CVX specification, ``geo_mean`` constrains the elements of the
    vector to be nonnegative. When used with numerical arguments,
    ``geo_mean`` returns ``-Inf`` if any element is negative. Concave
    and increasing.

``huber(x,M)``
    The Huber function, defined as
	.. math:: 

		f_{\text{huber}}(x,M) \triangleq \begin{cases} |x|^2 & |x| \leq M \\ 2M|x|-M^2 & |x| \geq M \end{cases}

    Convex. If $x$ is a vector or array, the function is applied on an elementwise basis. If $M$ is omitted, then $M=1$ is assumed; but if it supplied, it must be a positive constant. Also callable as ``huber(x,M,t)``, which computes ``t+t*huber(x/t,M)``, useful for concomitant scale estimation (see [Owen06]_).

``huber_circ(x,M)``
    The circularly symmetric Huber function, defined as
	.. math:: 

		f_{\text{huber\_circ}}(x,M) \triangleq \begin{cases} \|x\|_2^2 & \|x\|_2 \leq M \\ 2M\|x\|_2-M^2 & \|x\|_2 \geq M \end{cases}

    Convex. Same (and implemented) as ``huber_pos(norm(x),M)``.

``huber_pos(x,M)``
    The same as the Huber function for nonnegative ``x``; zero for
    negative ``x``. Convex and nondecreasing.

``inv_pos``
    The inverse of the positive portion, :math:`1/\max\{x,0\}`. Inside
    CVX specification, imposes constraint that its argument is
    positive. Outside CVX specification, returns :math:`+\infty` if
    :math:`x\leq 0`. Convex and decreasing.

† ``kl_div``
    Kullback-Leibler distance:
    
    .. math::
    
    	f_{\text{kl}}(x,y) \triangleq \begin{cases} x\log(x/y)-x+y & x,y>0 \\ 0 & x=y=0 \\ +\infty & \text{otherwise} \end{cases}
    	
    Convex. Outside CVX specification, returns :math:`+\infty` if arguments aren't in the
    domain.

``lambda_max``
    maximum eigenvalue of a real symmetric or complex Hermitian matrix.
    Inside CVX, imposes constraint that its argument is symmetric
    (if real) or Hermitian (if complex). Convex.

``lambda_min``
    minimum eigenvalue of a real symmetric or complex Hermitian matrix.
    Inside CVX, imposes constraint that its argument is symmetric
    (if real) or Hermitian (if complex). Concave.

``lambda_sum_largest(X,k)``
    sum of the largest :math:`k` values of a real symmetric or complex
    Hermitian matrix. Inside CVX, imposes constraint that its
    argument is symmetric (if real) or Hermitian (if complex). Convex.

``lambda_sum_smallest(X,k)``
    sum of the smallest :math:`k` values of a real symmetric or complex
    Hermitian matrix. Inside CVX, imposes constraint that its
    argument is symmetric (if real) or Hermitian (if complex). Concave.

``log_det``
    log of determinant of a positive definite matrix,
    :math:`\log \det(X)`. When used inside a CVX specification,
    ``log_det`` constrains its argument to be symmetric (if real) or
    Hermitian (if complex) and positive definite. With numerical
    argument, ``log_det`` returns ``-Inf`` if these constraints are not
    met. Concave.

‡ ``log_normcdf(x)``
    logarithm of cumulative distribution function of standard normal
    random variable. Concave and increasing. The current implementation
    is a fairly crude SDP-representable approximation, with modest
    accuracy over the interval :math:`[-4,4]`; we intend to replace it
    with a much better approximation at some point.
    
† ``log_prod(x)``
	:math:`\log\prod_i x_i` if when :math:`x` is positive; :math:`-\infty` otherwise. 
	Concave and nonincreasing. Equivalent to ``sum_log(x)``.

† ``log_sum_exp(x)``
    the logarithm of the sum of the elementwise exponentials of ``x``.
    Convex and nondecreasing.

``logsumexp_sdp``
    a polynomial approximation to the log-sum-exp function with global
    absolute accuracy. This can be used to estimate the log-sum-exp
    function without using the successive approximation method.

``matrix_frac(x,Y)``
    matrix fractional function, :math:`x^TY^{-1}x`. In CVX, imposes constraint 
    that :math:`Y` is symmetric (or Hermitian) and positive definite; outside CVX, 
    returns :math:`+\infty` unless :math:`Y=Y^T\succ 0`. Convex.

``norm_largest(x,k)``
    For real and complex vectors, returns the sum of the largest ``k``
    *magnitudes* in the vector ``x``. Convex.

``norm_nuc(X)``
    The sum of the singular values of a real or complex matrix ``X``.
    (This is the dual of the usual spectral matrix norm, *i.e.*, the
    largest singular value.) Convex.

‡ ``norms(x,p,dim)``, ``norms_largest(x,k,dim)``
    Computes *vector* norms along a specified dimension of a matrix or
    N-d array. Useful for sum-of-norms and max-of-norms problems.
    Convex.

``poly_env(p,x)``
    Computes the value of the *convex or concave envelope* of the
    polynomial described by ``p`` (in the ``polyval`` sense). ``p`` must
    be a real constant vector whose length ``n`` is 0, 1, 2, 3, or some
    other *odd* length; and ``x`` must be real and affine. The sign of
    the first nonzero element of ``p`` determines whether a convex
    (positive) or concave (negative) envelope is constructed. For
    example, consider the function :math:`p(x)\triangleq (x^2-1)^2=x^4-2x^2+1`, 
    depicted along with its convex envelope in the figure below.

    The two coincide when :math:`|x|\geq 1`, but deviate when
    :math:`|x|<1`. Attempting to call ``polyval([1,0,2,0,1],x)`` in a
    CVX model would yield an error, but a call to ``poly_env([1,0,2,0,1],x)`` 
    yields a valid representation of the envelope. For convex or concave 
    polynomials, this function produces the same result as ``polyval``.
    
    .. figure:: envelope.pdf

       The polynomial function :math:`p(x)=x^4-2x^2+1` and its convex envelope.
       
``pos(x)``
    :math:`\max\{x,0\}`, for real :math:`x`. Convex and increasing.

‡ ``pow_abs(x,p)``
    :math:`|x|^p` for :math:`x\in\mathbf{R}` or :math:`x\in\mathbf{C}`
    and :math:`p\geq 1`. Convex.

‡ ``pow_pos(x,p)``
    :math:`\max\{x,0\}^p` for :math:`x\in\mathbf{R}` and
    :math:`p\geq 1`. Convex and nondecreasing.

‡ ``pow_p(x,p)``
	for :math:`x\in\mathbf{R}` and real constant :math:`p`, computes nonnegative convex
	and concave branches of the power function:

	.. math::
		\begin{array}{ccl}
			p\leq 0 & f_p(x) \triangleq \begin{cases} x^p & x > 0 \\ +\infty & x \leq 0 \end{cases} & \text{convex, nonincreasing} \\
			0 < p \leq 1 & f_p(x) \triangleq \begin{cases} x^p & x \geq 0 \\ -\infty & x < 0 \end{cases} & \text{concave, nondecreasing} \\
			p \geq 1 & f_p(x) \triangleq \begin{cases} x^p & x \geq 0 \\ +\infty & x < 0 \end{cases} & \text{convex, nonmonotonic}
		\end{array}

``prod_inv(x)``
	:math:`\prod_i x_i^{-1}` when :math:`x` is positive; :math:`+\infty` otherwise. Convex
	and nonincreasing.

``quad_form(x,P)``
    :math:`x^TPx` for real :math:`x` and symmetric :math:`P`, and
    :math:`x^HPx` for complex :math:`x` and Hermitian :math:`P`. Convex
    in :math:`x` for :math:`P` constant and positive semidefinite;
    concave in :math:`x` for :math:`P` constant and negative
    semidefinite.
    
.. note::
	Quadratic functions such as ``quad_form``, ``sum_square`` can often be replaced
	by the ``norm`` function without sacrificing equivalence. For numerical reasons,
	this alternate formulation is *preferred*. Please see :ref:`quad-forms` for
	more information.

``quad_over_lin(x,y)``
    :math:`x^Tx/y` for :math:`x \in \mathbf{R}^n`, :math:`y >0`; for
    :math:`x \in \mathbf{C}^n`, :math:`y>0`, :math:`x^*x/y`. In CVX
    specification, adds constraint that :math:`y>0`. Outside CVX
    specification, returns :math:`+\infty` if :math:`y\leq 0`. Convex,
    and decreasing in :math:`y`.

``quad_pos_over_lin(x,y)``
    ``sum_square_pos( x )/y`` for :math:`x\in\mathbf{R}^n`, :math:`y>0`.
    Convex, increasing in :math:`x`, and decreasing in :math:`y`.

† ``rel_entr(x)``
    Scalar relative entropy; ``rel_entr(x,y)=x.*log(x/y)``. Convex.

``sigma_max``
    maximum singular value of real or complex matrix. Same as ``norm``.
    Convex.

``square``
    :math:`x^2` for :math:`x \in \mathbf{R}`. Convex.

``square_abs``
    :math:`|x|^2` for :math:`x\in\mathbf{R}` or :math:`x\in\mathbf{C}`.

``square_pos``
    :math:`\max\{x,0\}^2` for :math:`x\in\mathbf{R}`. Convex and
    increasing.

``sum_largest(x,k)``
    sum of the largest :math:`k` values, for real vector :math:`x`. Convex and nondecreasing.

† ``sum_log(x)``
    :math:`\sum_i\log(x_i)` when :math:`x` is positive; :math:`-\infty` otherwise.
    Concave and nondecreasing.
    
``sum_smallest(x,k)``
    sum of the smallest :math:`k` values, *i.e.*, equivalent to ``-sum_largest(-x,k)``. Concave and nondecreasing.

``sum_square``
    Equivalent to ``sum(square(x))``, but more efficient. Convex. Works only for real values.

``sum_square_abs``
    Equivalent to ``sum(square_abs(x))``, but more efficient. Convex.

``sum_square_pos``
    Equivalent to ``sum(square_pos(x))``, but more efficient. Works only for real values. 
    Convex and increasing.

``trace_inv(X)``
    trace of the inverse of an SPD matrix ``X``, which is the same as
    the sum of the inverses of the eigenvalues. Convex. Outside of
    CVX, returns ``+Inf`` if argument is not positive definite.

``trace_sqrtm(X)``
    trace of the matrix square root of a positive semidefinite matrix
    ``X``. which is the same as the sum of the squareroots of the
    eigenvalues. Concave. Outside of CVX, returns ``+Inf`` if
    argument is not positive semidefinite.
    
.. _sets-ref:    

Sets
----

CVX currently supports the following sets; in each case, ``n`` is a
positive integer constant.

``nonnegative(n)``
	.. math:: 
	
		R^n_+ \triangleq \left\{\,x\in\mathbf{R}^n\,~|~\,x_i\geq 0,~i=1,2,\dots,n\,\right\}

``simplex(n)``
    .. math:: 
    
    	R^n_{1+} \triangleq \left\{\,x\in\mathbf{R}^n\,~|~\,x_i\geq 0,~i=1,2,\dots,n,~\textstyle\sum_ix_i=1\,\right\}

``lorentz(n)``
    .. math:: 
    
    	\mathbf{Q}^n \triangleq \left\{\,(x,y)\in\mathbf{R}^n\times\mathbf{R}\,~|~\,\|x\|_2\leq y\,\right\}

``rotated_lorentz(n)``
    .. math:: 
    
    	\mathbf{Q}^n_r \triangleq \left\{\,(x,y,z)\in\mathbf{R}^n\times\mathbf{R}\times\mathbf{R}\,~|~\,\|x\|_2\leq \sqrt{yz},~y,z\geq 0\,\right\}
    	
``complex_lorentz(n)``
    .. math:: 
    
    	\mathbf{Q}^n_c \triangleq \left\{\,(x,y)\in\mathbf{C}^n\times\mathbf{R}\,~|~\,\|x\|_2\leq y\,\right\}

``rotated_complex_lorentz(n)``
    .. math:: 
    
    	\mathbf{Q}^n_{rc} \triangleq \left\{\,(x,y,z)\in\mathbf{C}^n\times\mathbf{R}\times\mathbf{R}\,~|~\,\|x\|_2\leq \sqrt{yz},~y,z\geq 0\,\right\}

``semidefinite(n)``
    .. math:: 
    
    	\mathbf{S}^n_+ \triangleq \left\{\,X\in\mathbf{R}^{n\times n}\,~|~\,X=X^T,~X\succeq 0\,\right\}

``hermitian_semidefinite(n)``
    .. math:: 
    
    	\mathbf{H}^n_+ \triangleq \left\{\,Z\in\mathbf{C}^{n\times n}\,~|~\,Z=Z^H,~Z\succeq 0\,\right\}

``nonneg_poly_coeffs(n)``
    The cone of all coefficients of nonnegative polynomials of degree :math:`n`; :math:`n` must be even: 
    
    .. math:: 
    
    	\mathbf{P}_{+,n} \triangleq \left\{\,p\in\mathbf{R}^n[n+1]\,~|~\,\sum_{i=0}^n p_{i+1} x^{n-i} \geq 0 ~ \forall x\in\mathbf{R}\,\right\}

``convex_poly_coeffs(n)``
    The cone of all coefficients of convex polynomials of degree :math:`n`; :math:`n` must be even:
    
    .. math:: 
    
    	\mathbf{P}_{+,n} \triangleq \left\{\,p\in\mathbf{R}^n[n+1]\,~|~\,\sum_{i=0}^{n-2} (n-i)(n-i-1) p_{i+1} x^{n-i-2} \geq 0 ~ \forall x\in\mathbf{R}\,\right\}

``exp_cone``
    .. math:: 
    
    	\mathbf{E} \triangleq \text{cl}\left\{\,(x,y,z)\in\mathbf{R}\times\mathbf{R}\times\mathbf{R}\,~|~\,y>0,~ye^{x/y}\leq z\,\right\}

``geo_mean_cone(n)``
    .. math:: 
    
    	\mathbf{G}_n \triangleq \text{cl}\left\{\,(x,y)\in\mathbf{R}^n\times\mathbf{R}^n\times\mathbf{R}^n\,~|~\,x\geq 0,~(\prod_{i=1}^n x_i)^{1/n} \geq y\,\right\}
    
Commands
---------

``cvx_begin``
	Begins a new CVX model. If a model is already in progress, it will issue a warning
	and clear it. See :ref:`begin-end` for a full description, including the modifying
	keywords that control solver output, SDP mode, GDP mode, etc.

``cvx_clear``
	Clears any model being constructed. Useful when an error has been made and it is
	necessary to start from the beginning. Whereas ``cvx_begin`` issues a warning if
	called with a model in progress, ``cvx_clear`` is silent.
	
``cvx_end``
	Signals the end of a CVX model. In typical use, this instructs CVX to begin the solution process.
	See :ref:`begin-end`.
	
``cvx_expert``
	Controls the issuance of warnings when models requiring the use of successive
	approximation are employed; see :ref:`successive` more details.

``cvx_power_warning``
	Controls if and when CVX issues warnings during the construction of models involving
	rational power functions (i.e., ``x^p``, where ``x`` is a variable and ``p`` is a constant);
	see :ref:`powerfunc`.

``cvx_precision``
	Controls solver precision; see :ref:`solver-precision`.

``cvx_quiet``
	Enables or disables screen output during the solution process; see :ref:`solver-output`.
	Also see :ref:`begin-end` for the newer, preferred syntax ``cvx_begin quiet``.

``cvx_save_prefs``
	Saves the current states for ``cvx_expert``, ``cvx_power_warning``, ``cvx_precision``, 
	and ``cvx_solver`` to disk, so that their values are retained when quitting and
	re-starting MATLAB. The file is saved in MATLAB's preference directory, which can
	be located by typing the ``prefdir`` command.

``cvx_setup``
	The setup script used to install and configure CVX; see :ref:`install`.

``cvx_solver``
	Selects the solver to be employed when solving CVX models; see :ref:`solver-selection`.
	
``cvx_solver_settings``
	Allows the user to deliver advanced, solver-specific settings to the solver that CVX
	does not otherwise support; see :ref:`solver-settings`.

``cvx_version``
	Prints information about the current versions of CVX, Matlab, and the operating system.
	When submitting bug reports, please include the output of this command.

``cvx_where``
	Returns the directory where CVX is installed.

``dual variable``, ``dual variables``
	Creates one or more dual variables to be connected to constraints in the current model;
	see :ref:`dual-variables`.

``expression``, ``expressions``
	Creates one or more expression holders; see :ref:`assignment`.

``maximise``, ``maximize``
	Specifies a maximization objective; see :ref:`objectives`.

``minimise``, ``minimize``
	Specifies a minimization objective; see :ref:`objectives`.

``variable``, ``variables``
	Creates one or more variables for use in the current CVX model; see :ref:`variables`.
