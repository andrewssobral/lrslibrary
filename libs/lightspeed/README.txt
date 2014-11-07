This is a library of efficient and useful matlab functions, with an
emphasis on statistics.
See Contents.m for a synopsis.

You can place the lightspeed directory anywhere.
To make sure lightspeed is always in your path, create a startup.m
file in your matlab directory, if you don't already have one, and add
a line like this:
  addpath(genpath('c:\matlab\lightspeed'))
Replace 'c:\matlab\lightspeed' with the location of the lightspeed directory.

There are some Matlab Extension (MEX) files that need to be compiled.
This can be done in matlab via:
  cd c:\matlab\lightspeed
  install_lightspeed

I recommend using Microsoft Visual C++ as the mex compiler, though this is not required.  You can set the mex compiler by typing 'mex -setup' in matlab.

To use Microsoft Visual C++ 2010 with Matlab 7.10 (R2010a), you will need to download this patch:
http://www.mathworks.com/support/solutions/en/data/1-D5W493/?solution=1-D5W493

To use Microsoft Visual C++ with Matlab 7.0 (R14), you will need to download 
R14 service pack 2, as described here:
http://www.mathworks.com/support/solutions/en/data/1-UMEKK/?solution=1-UMEKK

To compile mex files on a Snow Leopard upgrade: Go to mexopts.sh in your $HOME/.matlab/ directory, and change the line SDKROOT='/Developer/SDKs/MacOSX10.5.sdk' to SDKROOT='/Developer/SDKs/MacOSX10.6.sdk'. That line is not updated during updating mac OSX, so you need to do manually. This file also exists in the standard matlab bin, so if you run mex with the -v option it will tell you which mexopts.sh file it's looking at.  You may also need to change some lines in install_lightspeed.m.  By default, install_lightspeed.m is set up for 64-bit MacOSX 10.6 with gcc-4.0.  If you are using some other version of MacOSX or some other compiler, then you need to edit a few lines (see the comments in that file).

You can find timing tests in the tests/ subdirectory.  
The test_lightspeed.m script will run all tests, and is a good way to check 
that lightspeed installed properly.


Tom Minka
