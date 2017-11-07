YALL1: Your ALgorithms for L1

Version 1.4

COPYRIGHT (c) 2010 Yin Zhang, Junfeng Yang, and Wotao Yin.

YALL1 is distributed under the terms of the GNU General Public License 3.0.

http://www.gnu.org/copyleft/gpl.html

Permission to use, copy, modify, and distribute this software for
any purpose without fee is hereby granted, provided that this entire
notice is included in all copies of any software which is or includes
a copy or modification of this software and in all copies of the
supporting documentation for such software.
This software is being provided "as is", without any express or
implied warranty.  In particular, the authors do not make any
representation or warranty of any kind concerning the merchantability
of this software or its fitness for any particular purpose.

--------------------------------------------------------------------------------

First run "Run_me_1st" from this directory. Then, you can try
any of the demo files, but in order to run the 3rd-party codes
SPGl1 and l1_ls (provided you already installed them), you will
need to edit relevant demo files.

In order to run demo_hard.m, you need to download from the
YALL1 site the 2 data files: hard150.mat and hard8nz.mat, and
put them in the same folder or the search path.

In order to use the discrete Walsh-Hadamard transform, the script
Run_Me_1st.m will try to "mex" the file fastWHtrans.cpp in the
Utilities directory, which will require a relevant compiler
installed on your system.
