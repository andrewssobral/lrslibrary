.. _sdp-mode:

=============================
Semidefinite programming mode
=============================

Those who are familiar with *semidefinite programming* (SDP) know that
the constraints that utilize the set ``semidefinite(n)`` in the discussion
on :ref:`sets` above are, in practice, typically expressed using
*linear matrix inequality* (LMI) notation. For example, given
:math:`X=X^T\in\mathbf{R}^{n \times n}`, the constraint
:math:`X\succeq 0` denotes that :math:`X\in\mathbf{S}^n_+`; that is,
that :math:`X` is positive semidefinite.

CVX provides a special *SDP mode* that allows this LMI notation
to be employed inside CVX models using Matlab's standard inequality
operators ``>=``, ``<=``. In order to use it, one simply
begins a model with the statement ``cvx_begin sdp`` or ``cvx_begin SDP``
instead of simply ``cvx_begin``.

When SDP mode is engaged, CVX
interprets certain inequality constraints in a different manner. To be
specific:

-  Equality constraints are interpreted the same (*i.e.*, elementwise).
-  Inequality constraints involving vectors and scalars are interpreted
   the same; *i.e.*, elementwise.
-  Inequality constraints involving non-square matrices are
   *disallowed*; attempting to use them causes an error. If you wish to
   do true elementwise comparison of matrices ``X`` and ``Y``, use a
   vectorization operation ``X(:) <= Y(:)`` or ``vec( X ) <= vec( Y )``.
   (``vec`` is a function provided by CVX that is equivalent to the
   colon operation.)
-  Inequality constraints involving real, square matrices are
   interpreted as follows:

   ==========  =======   ============================
   ``X >= Y``  becomes   ``X - Y == semidefinite(n)``
   ``X <= Y``  becomes   ``Y - X == semidefinite(n)``
   ==========  =======   ============================

   If either side is complex, then the inequalities are interpreted as follows:

   +--------------+----------+------------------------------------------+
   | ``X >= Y``   | becomes  | ``X - Y == hermitian_semidefinite(n)``   |
   +--------------+----------+------------------------------------------+
   | ``X <= Y``   | becomes  | ``Y - X == hermitian_semidefinite(n)``   |
   +--------------+----------+------------------------------------------+

-  There is one additional restriction: both ``X`` and ``Y`` must be the
   same size, or one must be the scalar zero. For example, if ``X`` and
   ``Y`` are matrices of size ``n``,

   +----------------------+------+----------------------+-------------+
   | ``X >= 1``           | or   | ``1 >= Y``           | *illegal*   |
   +----------------------+------+----------------------+-------------+
   | ``X >= ones(n,n)``   | or   | ``ones(n,n) >= Y``   | *legal*     |
   +----------------------+------+----------------------+-------------+
   | ``X >= 0``           | or   | ``0 >= Y``           | *legal*     |
   +----------------------+------+----------------------+-------------+

   In effect, CVX enforces a stricter interpretation of the
   inequality operators for LMI constraints.
   
-  Note that LMI constraints enforce symmetry (real or Hermitian, as
   appropriate) on their inputs. Unlike
   `SDPSOL <http://www.stanford.edu/~boyd/old_software/SDPSOL.html>`_,
   CVX does not extract the symmetric part for you: you must take
   care to insure symmetry yourself. Since CVX supports the
   declaration of symmetric matrices, this is reasonably
   straightforward. If CVX cannot determine that an LMI is symmetric
   to within a reasonable numeric tolerance, a warning will be issued.
   We have provided a function ``sym(X)`` that extracts the symmetric
   part of a square matrix; that is, ``sym(X) = 0.5*(X+X')``.
   
-  A dual variable, if supplied, will be applied to the converted
   equality constraint. It will be given a positive semidefinite value
   if an optimal point is found.

So, for example, the CVX model found in the file
:file:`examples/closest_toeplitz_sdp.m`,

::

    cvx_begin
        variable Z(n,n) hermitian toeplitz
        dual variable Q
        minimize( norm( Z - P, 'fro' ) )
        Z == hermitian_semidefinite( n ) : Q;
    cvx_end

can also be written as follows:

::

    cvx_begin sdp
        variable Z(n,n) hermitian toeplitz
        dual variable Q
        minimize( norm( Z - P, 'fro' ) )
        Z >= 0 : Q;
    cvx_end

Many other examples in the CVX example library utilize semidefinite
constraints; and all of them use SDP mode. To find them, simply search
for the text ``cvx_begin sdp`` in the ``examples/`` subdirectory tree
using your favorite file search tool. One of these examples is
reproduced in :ref:`indexed-dual`.

Since semidefinite programming is popular, some may wonder why SDP mode
is not the default behavior. The reason for this is that we place a
strong emphasis on maintaining consistency between Matlab's native
behavior and that of CVX. Using the ``>=``, ``<=``, ``>``,
``<`` operators to create LMIs represents a deviation from that ideal.
For example, the expression ``Z >= 0`` in the example above constrains
the variable ``Z`` to be positive semidefinite. But after the model has
been solved and ``Z`` has been replaced with a numeric value, the
expression ``Z >= 0`` will test for the *elementwise* nonnegativity of
``Z``. To verify that the numeric value of ``Z`` is, in fact, positive
semidefinite, you must perform a test like ``min(eig(Z)) >= 0``.
