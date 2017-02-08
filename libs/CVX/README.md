## CVX: A system for disciplined convex programming

### Introduction

CVX is a Matlab package for convex optimizaton.
To learn more about what CVX is and how to use it, please visit our
[web site](http://cvxr.com/cvx), read the
[users' guide](http://cvxr.com/cvx/doc), and browse the
[example library](http://cvxr.com/cvx/examples). 

The best way to obtain CVX is to visit the
[download page](http://cvxr.com/cvx/download), which provides
pre-built archives containing standard and professional versions of CVX
tailored for specific operating systems.
This repository provides an alternate means of obtaining CVX 
for those who prefer to obtain their software via clone, fork, or
subrepo, or who simply like to browse source code.

### About this repository

This is a filtered mirror of the main branch of our internal 
development repository, with all administrative and non-redistributable
files removed. The differences between the files found here and
the packages offered on our [download](http://cvxr.com/cvx/download)
page are as follows:

* The functionality supporting the use of commercial solvers
  [Gurobi](http://gurobi.com) and [MOSEK](http://mosek.com) is not
  present. This functionality is exclusive to the non-redistributable
  Professional version, available for download on the web site.
* The [documentation](http://cvxr.com/cvx/doc) is not compiled. 
  The soruce code is provided in the `doc/` subdirectory, and
  requires the [Sphinx](http://sphinx-doc.org) Python documentation
  generator and a LaTeX system such as 
  [TeXLive](http://tug.org/texlive/) to generate it.
* The solvers [SDPT3](https://github.com/sqlp/sdpt3/) and 
  [SeDuMi](https://github.com/sqlp/sedumi/) are provided as 
  *submodules*. In other words, they are not actually *included*
  in the repository itself; instead, *links* to their separate
  GitHub repositories are included. Make sure you use the
  `--recursive` flag when cloning this repository to download
  the solvers along with CVX.

Needless to say, working from the raw source is not a straighforward
process. We know this first-hand! This is for hardcore GitHub users.

### Support

We intend to keep this repository complete and up to date, but
its use is completely unsupported. That said, if you are having
issues related specifically to the repository itself, please feel
free to [submit an report](https://github.com/cvxr/CVX/issues) to the
repository's issue tracker. Please do *not* submit other types
of issues (CVX bug reports, usage questions, etc.) to this tracker;
they will likely be ignored.

We cannot provide direct email support for CVX without
a paid contract. However, we have created and assembled a variety of
avenues for obtaining help with CVX in particular or optimization in general.
Please see the [Support section](http://cvxr.com/cvx/doc/support.html)
of the documentation for more details.

### License

Most of the files in this repository are governed by the terms of our 
[GPLv3](http://www.gnu.org/licenses/gpl-3.0.html)-based 
[CVX Standard License](http://cvxr.com/cvx/doc/license.html). Please
see the files 
[LICENSE.txt](https://github.com/cvxr/CVX/blob/master/README.txt) and 
[GPL.txt](https://github.com/cvxr/CVX/blob/master/GPL.txt), 
included with the distribution, for more details about the license.

The contents of the example library in the `examples/` subdirectory
are *public domain*. We do ask that if you use any of this content in
your own work, that you acknowledge the source and any specific authors
cited therein.

Thank you for your interest in CVX!    
Michael Grant and Stephen Boyd    
[CVX Research, Inc.](http://cvxr.com)    
(c) 2014. All rights reserved.
