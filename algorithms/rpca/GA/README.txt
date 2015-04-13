CODE FOR GRASSMANN AVERAGES
===========================

author:        Søren Hauberg <sohau@dtu.dk>
version:       0.3
release date:  June 17, 2014
license:       GPLv3+

INSTALLATION
============

The built-in Matlab routines for computing medians and trimmed means are
dreadfully slow as they rely on sorting. We therefore distribute faster
implementations of these routines. We distribute binaries for 64 bit Linux
and Mac. If you run a different platform, you need to compile these routines
by typing

  compile

from the Matlab prompt. If this fails you will only be able to run the
'grassmann_average' routine. Please contact me if this happens, then I
will do my best to help you out.

If you are compile the code on other platforms than those distributed,
please get in touch with me if you are interested in making your
compiled binaries available to others.

It is possible to compile these functions to run in parallel using OpenMP. This provides a substantial speed-up on multi-core machines.

EXAMPLE
=======

As a simple example, we will first generate samples from a random Gaussian
distribution, and then estimate its leading component.

  D = 2; % we consider a two-dimensional problem
  N = 100; % we will generate 100 observations

  %% Generate a random Covariance matrix
  tmp = randn(D);
  Sigma = tmp.' * tmp;

  %% Sample from the corresponding Gaussian
  X = mvnrnd(zeros(D, 1), Sigma, N);

  %% Estimate the leading component
  comp = grassmann_average(X, 1); % the second input is the number of component to estimate

  %% Plot the results
  plot(X(:, 1), X(:, 2), 'ko', 'markerfacecolor', [255,153,51]./255);
  axis equal
  hold on
  plot(3*[-comp(1), comp(1)], 3*[-comp(2), comp(2)], 'k', 'linewidth', 2)
  hold off
  axis off

The call to 'grassmann_average' can be replaced with a corresponding call to 'grassmann_median' or 'trimmed_grassmann_average'.

BUGS
====

I am currently unaware of any bugs in the code, but if you find any please
get in touch, such that I can update the code.

HACKING GUIDE
=============

The main algorithm to understand is the 'Grassmann Average' which is implemented in 'grassmann_average.m'. This is the best starting point in order to understand the code. The two other algorithms ('grassmann_median' and 'trimmed_grassmann_average') or just small modifications of the core algorithm.

CHANGELOG
=========

Version 0.1 introduced a bug during the release refactoring which prevented multiple components from being estimated correctly.

Version 0.3 fixes a bug in 'trimmed_mean' which halves the trimming parameter.

ACKNOWLEDGEMENTS
================

The Gram Schmidt implementation comes from NNLS-SDP:
https://github.com/MPF-Optimization-Laboratory

Søren Hauberg is funded in parts by both the Villum Foundation and the Danish Council for Independent Research (Natural Sciences).

REFERENCES
==========

If you find the code useful or interesting, please cite the following paper:

  Grassmann Averages for Scalable Robust PCA
  Hauberg, S., Feragen, A. and Black, M.J.
  In Proceedings IEEE Conf. on Computer Vision and Pattern Recognition (CVPR), IEEE, Piscataway, NJ, USA. June 2014.

