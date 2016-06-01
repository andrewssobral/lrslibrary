/**************************************************************************************
  GCO_MATLAB - a Matlab wrapper for Olga Veksler's C++ graph-cut optimization library

  GCO_MATLAB Author(s):
     Andrew Delong <firstname.lastname@gmail.com>
     Anton Osokin <lastname.firstname@gmail.com>

     We're especially grateful to Lena Gorelick for helpful suggestions and for 
     tracking down so many bugs!

  GCoptimization Author(s):
     Olga Veksler <firstname@csd.uwo.ca>
  
  Description:
     This download provides a Matlab wrapper for the latest version of 'GCoptimization',
     Olga Veksler's multi-label optimization library written in C++. 
     
  Note: A wrapper for an earlier version GCoptimization was authored by Shai Bagon and
        is available at http://www.wisdom.weizmann.ac.il/~bagon/matlab.html
        
  Revision History:
     Oct 14, 2014; - Added GCO_ENERGYTERM and GCO_ENERGYTERMTYPE macros to make 
                     float/double terms easier. Please note that float/double energy terms
                     can cause expansion/swap steps to report a very small increase in energy 
                     due to accumulated rounding error inside the maxflow library.
                     If Inf or NaN values appear in the energy terms, behaviour is undefined.
     May 18, 2014; - Support Matlab R2014a by removing use of mxCreateReference
     Jan 15, 2014; - Compiles with gcc 4.6+ even without -fpermissive
     Apr 12, 2011; - Fixed bug when sparse data costs had a dense bucket (thanks Joseph Tighe!)
     Nov 25, 2010; - Detect MACI64 correctly (thanks Francis Lauzon)
     Aug 31, 2010; - Compiles with gcc 4.4.1 on opensuse 11.2 (thanks Wei Liu)
     Aug  7, 2010; - Fixed bug when data costs are computed in a callback (thanks Evan Herbst!)
                   - Fixed bug where setAllNeighbours didn't apply neighbourhood (Evan Herbst again)
     Jul 22, 2010; - Compiles with gcc 4.4.1 (thanks Julius Ziegler for the patch!)
     Jul  8, 2010; - Fixed crash in greedy code path when all labels get added (thanks Yangyan Li!)
     Apr 25, 2010; - Faster code path for sparse data costs; fixed related bug in higher-order labels
     Apr 21, 2010; - Added basic "verbose" mode (print cycle, energy, timings etc)
                   - Expansion cycles now focus on labels for which the energy decreased (faster)
     Apr 19, 2010; - Added sparse datacost support
                   - Allow GCO_SetLabelOrder to specify exact label order
     Apr 13, 2010; - Potts model is now the default if SetNeighbors is called without SetSmoothCost
                   - Fixed bug in higher-order label costs
                   - Expansion is now interruptable from MATLAB; temporary memory is freed
                   - Added GCO_ListHandles and allow GCO_Delete to accept multiple handles
                   - Better error message if a bad handle is passed to GCO_*
     Nov 17, 2009; - Fixed integer overflow in label-cost construction
                   - Fixed bug where greedy algorithm would sometimes skip a label
     Nov  6, 2009; - Fixed bug in re-setting label costs after Expansion
                   - Fixed bug in GCO_LoadLib
     Oct 27, 2009; - Removed support for arbitrary smoothcost (too slow, hard to maintain)
                   - Added support for re-setting data, smooth, and label costs after Expansion
                   - Added support for label subset costs
                   - Changed build process to directly use MEX command
                   - Added integer overflow checks into GCoptimization
                   - GCoptimization now uses maxflow-3.0 library
     Sep 12, 2009; - Added support for arbitrary smoothcost from Matlab via a function_handle
                   - Added GCO_UnitTest
                   - Build script now handles spaces in paths properly
     Aug 23, 2009; - First version for internal testing
     
***************************************************************************************/

0. System Requirements

- Matlab 7.4.0 (R2007a) or above for 32-bit.
  Matlab 7.6.0 (R2008) or above for 64-bit.

- Mex must be pre-configured to use a C++ compiler. (Run "mex -setup" if you have 
  not already.) The C++ code requires at least Visual C++ 2005 (VC8).
    

----------------------------------------------------------------------------------------
1. Installation

- The package should contain the following files:

     GCO_MATLAB files:
          gco\matlab\GCO_*.m            ; the Matlab commands that you can run
          gco\matlab\gco_matlab.cpp     ; the library that links Matlab to GCoptimization
          
     GCoptimization files:
          gco\*.{h,cpp}                 ; the GCoptimization C++ library

- Start Matlab, and make gco\matlab your working directory or add it to your path.

- To test your installation of GCO_MATLAB, run the GCO_UnitTest command.
  You should hopefully see output like below.
     >> GCO_UnitTest
     BuildLib PASSED
     LoadLib PASSED
     Create/Delete PASSED
     ...
     >>

- 

----------------------------------------------------------------------------------------
2. Getting Started -- A basic example, and important usage notes

Once GCO_UnitTest passes, you should be able run the example sequence of commands below.

     >> h = GCO_Create(4,3);             % Create new object with NumSites=4, NumLabels=3
     >> GCO_SetDataCost(h,[
        0 9 2 0;                         % Sites 1,4 prefer  label 1
        3 0 3 3;                         % Site  2   prefers label 2 (strongly)
        5 9 0 5;                         % Site  3   prefers label 3
        ]);
     >> GCO_SetSmoothCost(h,[
        0 1 2;      % 
        1 0 1;      % Linear (Total Variation) pairwise cost
        2 1 0;      %
        ]);
     >> GCO_SetNeighbors(h,[
        0 1 0 0;     % Sites 1 and 2 connected with weight 1
        0 0 1 0;     % Sites 2 and 3 connected with weight 1
        0 0 0 2;     % Sites 3 and 4 connected with weight 2
        0 0 0 0;
        ]);
     >> GCO_Expansion(h);                % Compute optimal labeling via alpha-expansion 
     >> GCO_GetLabeling(h)
     ans =
           1                             % Optimal labeling is (1,2,1,1)
           2
           1
           1
     >> [E D S] = GCO_ComputeEnergy(h)   % Energy = Data Energy + Smooth Energy
     E =
         4
     D = 
         2
     S = 
         2
     >> GCO_Delete(h);                   % Delete the GCoptimization object when finished


*** Before using the MATLAB wrapper, please note the following: ***
 
 - Sites and labels are identified with 1-based indices (i.e.  1..N  and *not*  0..N-1)

 - By default, all numeric costs should be int32, not single or double!
   To use single/double energy terms with the library, please type "help GCO_BuildLib"
   at the MATLAB command prompt.
   You'll receive a conversion warning if you pass in a large matrix of the wrong type.
   The only function that accepts double is GCO_SetNeighbors, because it needs a sparse matrix
   and sparse matrices only support double in MATLAB.
   ** The weights themselves must still be integer valued!! ** (1.0, 17.0, 42.0 etc)

----------------------------------------------------------------------------------------
3. GCO_MATLAB functions

Run 'help' in MATLAB to see the documentation of each function, e.g.
   >> help GCO_SetSmoothCost

Most of the GCO_MATLAB functions are one-to-one with the C++ methods in the 
GCoptimization library.
For more detailed documentation, please refer to the C++ library itself.
Relevant files are:
    GCO_README.TXT
    example.cpp
    GCoptimization.h
