.. _support:

=======
Support
=======

The user base for CVX has grown to such an extent that full email-based
support is no longer feasible for our free version of CVX. Therefore, we 
have created several avenues for obtaining support.

For help on how to *use* CVX, this users' guide is your first line of support.
Please make sure that you have attempted to find an answer to your question
here before you pursue other avenues. We have placed this document
`online <http://cvxr.com/cvx/doc>`_ and made it searchable in order to 
help you find the answers to the questions you may have.

With a package like CVX that encapsulates such mathematical complexity, it can sometimes
be unclear if a problem with a model is due to model formulation or a bug in CVX. See
:ref:`What is a bug? <whatbug>` below to help you discern the difference,
and to determine the most appropriate channel for support.

The CVX Forum
-------------

If your answers cannot be found here, consider posting your question
to the `CVX Forum <http://ask.cvxr.com>`_. This is a community forum
that allows our users to submit questions and to answer other people's questions.
The forum uses the open-source `Askbot <http://www.askbot.com>`_ system, and its format
should be familiar to anyone who participates in `OR-Exchange <http://www.or-exchange.com>`_,
`Stack Overflow <http://stackoverflow.com>`_, or any one of the `Stack Exchange <http://stackexchange.com>`_
family of sites. 

We highly encourage our expert users who enjoy helping others to participate in
this forum. We hope that it will not only serve as a resource for diagnosing problems
and issues, but a clearinghouse for advanced usage tips and tricks.

Bug reports
-----------

If you believe you have found a *bug* in CVX or in one of the underlying solvers, 
then we encourage you to submit a bug report---either by email to
to ``cvx@cvxr.com`` or through our
`web-based support portal <http://support.cvxr.com>`_. Please include the following in your
bug report so we can fully reproduce the problem:

1. the CVX model and supporting data that caused the error. 
2. a copy of any error messages that it produced
3. the CVX version number and build number
4. the version number of Matlab that you are running
5. the name and version of the operating system you are using

The easiest way to supply items 3-5 is to type ``cvx_version`` at the command
prompt and copy its output into your email message.

Please note that we do not own all of Matlab's toolboxes. We cannot debug a model that
employs functions from a toolbox we do not own.

.. _whatbug:

What *is* a bug?
-----------------

Certain issues are unambiguously bugs, and you should feel free to report them 
immediately. In particular, CVX often attempts to catch unexpected errors in key
places---including ``cvx_setup``, ``cvx_end``, etc. It will report those errors and
instruct you to report them to us. If your model produces a MATLAB error that CVX 
did not itself generate, and you cannot readily tie it to a syntax error in your 
model, please report that as well.

That said, because disciplined convex programming is a new concept for many, we often 
receive support requests for problems with their models that are not, in fact, bugs.
A particularly common class of support requests are reports of models being rejected
due to ``Disciplined`` ``convex`` ``programming`` ``error`` messages. For instance,
the following code

::

	variable x(10)
	norm(x) == 1
	
will produce this error in CVX:

::

	Error using cvxprob/newcnstr (line 181)
	Disciplined convex programming error:
	  Invalid constraint: {convex} == {real constant}
   	
Disciplined convex programming errors indicate that the model fails
to adhere to the rules in the :ref:`DCP ruleset <dcp>`. In nearly all cases,
this is *not* a bug, and should not be reported as such. Rather,
the underlying issue falls into one of two categories:

1. The model is *not convex* (mixed-integer or otherwise). A model with the
   nonlinear equation above would fall squarely in this category. CVX simply
   cannot solve such problems. In some cases,
   it is possible to transform a problem that is non-convex into one that is convex (for
   example, :ref:`geometric programs <gp-mode>`). This has not been the case for any
   problem submitted by our users---so far.
   
2. The model *is* convex, but it is still written in a manner that violates the rules.
   For instance, given the same vector ``x`` above, the constraint

   ::

      sqrt( sum( square( x ) ) ) <= 1
	 	
   is convex, but it violates the ruleset---so it is rejected. However, the 
   mathematically equivalent form
   
   ::

      norm( x ) <= 1   
		
   is acceptable.
   If your error is of this type, you will need to find a way to express your 
   problem in a DCP-compliant manner. We have attempted to supply all of the commonly used 
   functions that the underlying solvers can support; so if you cannot easily rewrite 
   your problem using the functions supplied, it may not be possible. If you think this
   is a possibility, you may wish to see if the wizards on the
   `CVX Forum <http://ask.cvxr.com>`_ have suggestions for you.
   
In rare cases, users have discovered that certain models were rejected with
``Disciplined`` ``convex`` ``programming`` ``error`` messages
even though they satisfied the :ref:`DCP ruleset <dcp>`.
We have not received a bug report of this type in quite a long time, however, so we
suspect our users have helped us iron out these issues.

Handling numerical issues
--------------------------

No developer likes to tell their customers that their software may not work for them.
Alas, we have no choice. The fact is that we cannot guarantee that CVX will be able to
solve your problem, *even if* it is formulated properly, even if it avoids the use of
integer or binary variables, even if it is of reasonable size, even if it avoids the use
of our experimental :ref:`exponential cone support <successive>`.

We blame the solvers---but we must come to their defense, too. 
Even the best and most mature solvers will struggle with a particular problem that
seems straightforward. Another solver may have no difficulty with that one, but fail 
to find an accurate solution on another. While sometimes these challenges are due
to bugs in the solver's implementation, quite often it is simply due to limits imposed
by the nature of finite numerical precision computation. And different problems push
those limits to different degrees. So the fact is that no solver is perfect, but 
no solver *can* be.

When we consider :ref:`mixed-integer problems <what-is-midcp>`, the situation is even worse.
Solvers must perform what is effectively an exhaustive search among
the integer variables to determine the correct solution. Yes, there are some intelligent
and innovative ways to speed up that search, and the performance of mixed-integer
solvers has improved dramatically over the years. But there will always be models for
which the exhaustive search will simply take too long.

None of this is much comfort if it is *your* model the solver is struggling with. Here
are some practical tips if you encounter this problem:

*Try a different solver.* 
  Use the ``cvx_solver`` command for this. If you are using 
  Gurobi or MOSEK, don't hesitate to try one of the free solvers if they are compatible
  with your problem.
  
*Reduce the precision.* 
  Consider inserting ``cvx_precision medium`` or even 
  ``cvx_precision low`` into your problem to see if that allows the solver to exit
  successfully. Of course, if it does succeed, make sure to check the results to see
  if they are acceptable for your application. If they are not, consider some of the
  other advice here to see if the solvability of your model may be improved.
  
*Remove constraints.* 
  If you think that one or more of the constraints might
  not be active at the solution, try removing them. If the solver terminates, you can
  confirm that your guess was correct by examining the solution to the modified problem.
  
*Add constraints*. 
  Consider adding simple bounds to the constraints to reduce the size of
  the feasible set. This will sometimes improve the numerical conditioning of the problem.
  Make them as tight as you can without impinging on the optimal set. If this modified
  problem is successfully solved, check the solution to see if any of the added bounds
  are active. If they are, relax them, and try again.
  
*Watch for scaling issues.* 
  Scaling issues are the most vexing problems for numerical
  solvers to deal with. Solvers will often re-scale the problem to reduce the dynamic
  range of the numerical coefficients, but doing so sometimes leads to undesirable
  effects on the solution. It is better to avoid scaling issues during the modeling
  process. Don't mix values of wildly different magnitudes, such as ``1e-3``
  and ``1e20``. Even better, try to avoid any numerical values (both in fixed parameters
  and likely values of the variables) that exceed ``1e8`` in absolute value.
  
*Try equivalent reformulations.* 
  It is quite likely that your model can be expressed
  in a variety of different ways. Certainly, you should begin with the most obvious
  and natural formulation; but if you encounter numerical issues, a reformulation may
  often solve them. One reformulation we highly recommend is to eliminate quadratic
  forms; see :ref:`this section <quad-forms>` for more details.
  
*Reach out to the CVX Forum.* 
  Share your struggles with the 
  `larger CVX community <http://ask.cvxr.com>`_! Perhaps
  they will have concrete suggestions for improving the solvability of your model.

CVX Professional support
-------------------------

Paid CVX Professional users will receive support through a trouble-ticket support system,
so that they can have confidence that their issues are being addressed promptly. This
infrastructure is still under development; we will update this section and the Web site
with more information once it has been completed.

