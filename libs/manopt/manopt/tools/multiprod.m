function c = multiprod(a, b, idA, idB)
%MULTIPROD  Multiplying 1-D or 2-D subarrays contained in two N-D arrays.
%   C = MULTIPROD(A,B) is equivalent  to C = MULTIPROD(A,B,[1 2],[1 2])
%   C = MULTIPROD(A,B,[D1 D2]) is eq. to C = MULTIPROD(A,B,[D1 D2],[D1 D2])
%   C = MULTIPROD(A,B,D1) is equival. to C = MULTIPROD(A,B,D1,D1)
%
%   MULTIPROD performs multiple matrix products, with array expansion (AX)
%   enabled. Its first two arguments A and B are "block arrays" of any
%   size, containing one or more 1-D or 2-D subarrays, called "blocks" (*).
%   For instance, a 5×6×3 array may be viewed as an array containing five
%   6×3 blocks. In this case, its size is denoted by 5×(6×3). The 1 or 2
%   adjacent dimensions along which the blocks are contained are called the
%   "internal dimensions" (IDs) of the array (°).
%
%   1) 2-D by 2-D BLOCK(S) (*)
%         C = MULTIPROD(A, B, [DA1 DA2], [DB1 DB2]) contains the products
%         of the P×Q matrices in A by the R×S matrices in B. [DA1 DA2] are
%         the IDs of A; [DB1 DB2] are the IDs of B.
%
%   2) 2-D by 1-D BLOCK(S) (*)
%         C = MULTIPROD(A, B, [DA1 DA2], DB1) contains the products of the
%         P×Q matrices in A by the R-element vectors in B. The latter are
%         considered to be R×1 matrices. [DA1 DA2] are the IDs of A; DB1 is
%         the ID of B.
%
%   3) 1-D by 2-D BLOCK(S) (*)
%         C = MULTIPROD(A, B, DA1, [DB1 DB2]) contains the products of the 
%         Q-element vectors in A by the R×S matrices in B. The vectors in A
%         are considered to be 1×Q matrices. DA1 is the ID of A; [DB1 DB2]
%         are the IDs of B.
%
%   4) 1-D BY 1-D BLOCK(S) (*)
%      (a) If either SIZE(A, DA1) == 1 or SIZE(B, DB1) == 1, or both,
%             C = MULTIPROD(A, B, DA1, DB1) returns products of scalars by 
%             vectors, or vectors by scalars or scalars by scalars.
%      (b) If SIZE(A, DA1) == SIZE(B, DB1), 
%             C = MULTIPROD(A, B, [0 DA1], [DB1 0]) or 
%             C = MULTIPROD(A, B, DA1, DB1) virtually turns the vectors
%             contained in A and B into 1×P and P×1 matrices, respectively,
%             then returns their products, similar to scalar products.
%             Namely, C = DOT2(A, B, DA1, DB1) is equivalent to 
%             C = MULTIPROD(CONJ(A), B, [0 DA1], [DB1 0]).
%      (c) Without limitations on the length of the vectors in A and B,
%             C = MULTIPROD(A, B, [DA1 0], [0 DB1]) turns the vectors
%             contained in A and B into P×1 and 1×Q matrices, respectively,
%             then returns their products, similar to outer products.
%             Namely, C = OUTER(A, B, DA1, DB1) is equivalent to
%             C = MULTIPROD(CONJ(A), B, [DA1 0], [0 DB1]).
%
%   Common constraints for all syntaxes:
%      The external dimensions of A and B must either be identical or 
%      compatible with AX rules. The internal dimensions of each block
%      array must be adjacent (DA2 == DA1 + 1 and DB2 == DB1 + 1 are
%      required). DA1 and DB1 are allowed to be larger than NDIMS(A) and
%      NDIMS(B). In syntaxes 1, 2, and 3, Q == R is required, unless the
%      blocks in A or B are scalars. 
%
%   Array expansion (AX):
%      AX is a powerful generalization to N-D of the concept of scalar
%      expansion. Indeed, A and B may be scalars, vectors, matrices or
%      multi-dimensional arrays. Scalar expansion is the virtual
%      replication or annihilation of a scalar which allows you to combine
%      it, element by element, with an array X of any size (e.g. X+10,
%      X*10, or []-10). Similarly, in MULTIPROD, the purpose of AX is to
%      automatically match the size of the external dimensions (EDs) of A
%      and B, so that block-by-block products can be performed. ED matching
%      is achieved by means of a dimension shift followed by a singleton
%      expansion:
%      1) Dimension shift (see SHIFTDIM).
%            Whenever DA1 ~= DB1, a shift is applied to impose DA1 == DB1.
%            If DA1 > DB1, B is shifted to the right by DA1 - DB1 steps.
%            If DB1 > DA1, A is shifted to the right by DB1 - DA1 steps.
%      2) Singleton expansion (SX).
%            Whenever an ED of either A or B is singleton and the
%            corresponding ED of the other array is not, the mismatch is
%            fixed by virtually replicating the array (or diminishing it to
%            length 0) along that dimension.
% 
%   MULTIPROD is a generalization for N-D arrays of the matrix
%   multiplication function MTIMES, with AX enabled. Vector inner, outer,
%   and cross products generalized for N-D arrays and with AX enabled are
%   performed by DOT2, OUTER, and CROSS2 (MATLAB Central, file #8782).
%   Elementwise multiplications (see TIMES) and other elementwise binary
%   operations with AX enabled are performed by BAXFUN (MATLAB Central,
%   file #23084). Together, these functions make up the “ARRAYLAB toolbox”.
%
%   Input and output format:
%      The size of the EDs of C is determined by AX. Block size is
%      determined as follows, for each of the above-listed syntaxes:
%      1) C contains P×S matrices along IDs MAX([DA1 DA2], [DB1 DB2]).
%      2) Array     Block size     ID(s)
%         ----------------------------------------------------
%         A         P×Q  (2-D)     [DA1 DA2]
%         B         R    (1-D)     DB1
%         C (a)     P    (1-D)     MAX(DA1, DB1)
%         C (b)     P×Q  (2-D)     MAX([DA1 DA2], [DB1 DB1+1])
%         ----------------------------------------------------
%         (a) The 1-D blocks in B are not scalars (R > 1).
%         (b) The 1-D blocks in B are scalars (R = 1).
%      3) Array     Block size     ID(s)
%         ----------------------------------------------------
%         A           Q  (1-D)     DA1
%         B         R×S  (2-D)     [DB1 DB2]
%         C (a)       S  (1-D)     MAX(DA1, DB1)
%         C (b)     R×S  (2-D)     MAX([DA1 DA1+1], [DB1 DB2])
%         ----------------------------------------------------
%         (a) The 1-D blocks in A are not scalars (Q > 1).
%         (b) The 1-D blocks in A are scalars (Q = 1).
%      4)     Array     Block size         ID(s)
%         --------------------------------------------------------------
%         (a) A         P        (1-D)     DA1
%             B         Q        (1-D)     DB1
%             C         MAX(P,Q) (1-D)     MAX(DA1, DB1)
%         --------------------------------------------------------------
%         (b) A         P        (1-D)     DA1
%             B         P        (1-D)     DB1
%             C         1        (1-D)     MAX(DA1, DB1)
%         --------------------------------------------------------------
%         (c) A         P        (1-D)     DA1
%             B         Q        (1-D)     DB1
%             C         P×Q      (2-D)     MAX([DA1 DA1+1], [DB1 DB1+1])
%         --------------------------------------------------------------
%
%   Terminological notes:
%   (*) 1-D and 2-D blocks are generically referred to as "vectors" and 
%       "matrices", respectively. However, both may be also called
%       “scalars” if they have a single element. Moreover, matrices with a
%       single row or column (e.g. 1×3 or 3×1) may be also called “row
%       vectors” or “column vectors”.
%   (°) Not to be confused with the "inner dimensions" of the two matrices
%       involved in a product X * Y, defined as the 2nd dimension of X and
%       the 1st of Y (DA2 and DB1 in syntaxes 1, 2, 3).
%
%   Examples:
%    1) If  A is .................... a 5×(6×3)×2 array,
%       and B is .................... a 5×(3×4)×2 array,
%       C = MULTIPROD(A, B, [2 3]) is a 5×(6×4)×2 array.
%
%       A single matrix A pre-multiplies each matrix in B
%       If  A is ........................... a (1×3)    single matrix,
%       and B is ........................... a 10×(3×4) 3-D array,
%       C = MULTIPROD(A, B, [1 2], [3 4]) is a 10×(1×4) 3-D array.
%
%       Each matrix in A pre-multiplies each matrix in B (all possible
%       combinations)
%       If  A is .................... a (6×3)×5   array,
%       and B is .................... a (3×4)×1×2 array,
%       C = MULTIPROD(A, B, [1 2]) is a (6×4)×5×2 array.
%
%   2a) If  A is ........................... a 5×(6×3)×2 4-D array,
%       and B is ........................... a 5×(3)×2   3-D array,
%       C = MULTIPROD(A, B, [2 3], [2]) is   a 5×(6)×2   3-D array.
%
%   2b) If  A is ........................... a 5×(6×3)×2 4-D array,
%       and B is ........................... a 5×(1)×2   3-D array,
%       C = MULTIPROD(A, B, [2 3], [2]) is   a 5×(6×3)×2 4-D array.
%
%   4a) If both A and B are .................. 5×(6)×2   3-D arrays,
%       C = MULTIPROD(A, B, 2) is .......... a 5×(1)×2   3-D array, while
%   4b) C = MULTIPROD(A, B, [2 0], [0 2]) is a 5×(6×6)×2 4-D array
%
%   See also DOT2, OUTER, CROSS2, BAXFUN, MULTITRANSP, MULTITRACE, MULTISCALE.

% $ Version: 2.1 $
% CODE      by:            Paolo de Leva
%                          (Univ. of Rome, Foro Italico, IT)    2009 Jan 24
%           optimized by:  Paolo de Leva
%                          Jinhui Bai (Georgetown Univ., D.C.)  2009 Jan 24
% COMMENTS  by:            Paolo de Leva                        2009 Feb 24
% OUTPUT    tested by:     Paolo de Leva                        2009 Feb 24
% -------------------------------------------------------------------------

error( nargchk(2, 4, nargin) ); % Allow 2 to 4 input arguments
switch nargin % Setting IDA and/or IDB
    case 2, idA = [1 2]; idB = [1 2];
    case 3, idB = idA;
end

% ESC 1 - Special simple case (both A and B are 2D), solved using C = A * B

     if ndims(a)==2 && ndims(b)==2 && ...
         isequal(idA,[1 2]) && isequal(idB,[1 2])
         c = a * b; return
     end

% MAIN 0 - Checking and evaluating array size, block size, and IDs

     sizeA0 = size(a);
     sizeB0 = size(b);
     [sizeA, sizeB, shiftC, delC, sizeisnew, idA, idB, ...
     squashOK, sxtimesOK, timesOK, mtimesOK, sumOK] = ...
                                           sizeval(idA,idB, sizeA0,sizeB0);

% MAIN 1 - Applying dimension shift (first step of AX) and 
%          turning both A and B into arrays of either 1-D or 2-D blocks

     if sizeisnew(1), a = reshape(a, sizeA); end    
     if sizeisnew(2), b = reshape(b, sizeB); end

% MAIN 2 - Performing products with or without SX (second step of AX)

     if squashOK % SQUASH + MTIMES (fastest engine)
         c = squash2D_mtimes(a,b, idA,idB, sizeA,sizeB, squashOK); 
     elseif timesOK % TIMES (preferred w.r. to SX + TIMES)
         if sumOK, c = sum(a .* b, sumOK);
         else      c =     a .* b; end
     elseif sxtimesOK % SX + TIMES
         if sumOK, c = sum(bsxfun(@times, a, b), sumOK);
         else      c =     bsxfun(@times, a, b); end
     elseif mtimesOK % MTIMES (rarely used)
         c = a * b;
     end

% MAIN 3 - Reshaping C (by inserting or removing singleton dimensions)

     [sizeC sizeCisnew] = adjustsize(size(c), shiftC, false, delC, false);
     if sizeCisnew, c = reshape(c, sizeC); end


function c = squash2D_mtimes(a, b, idA, idB, sizeA, sizeB, squashOK)
% SQUASH2D_MTIMES  Multiproduct with single-block expansion (SBX).
%    Actually, no expansion is performed. The multi-block array is
%    rearranged from N-D to 2-D, then MTIMES is applied, and eventually the
%    result is rearranged back to N-D. No additional memory is required.
%    One and only one of the two arrays must be single-block, and its IDs
%    must be [1 2] (MAIN 1 removes leading singletons). Both arrays
%    must contain 2-D blocks (MAIN 1 expands 1-D blocks to 2-D).

    if squashOK == 1 % A is multi-block, B is single-block (squashing A)

        % STEP 1 - Moving IDA(2) to last dimension
        nd = length(sizeA);
        d2 = idA(2);    
        order = [1:(d2-1) (d2+1):nd d2]; % Partial shifting
        a = permute(a, order); % ...×Q

        % STEP 2 - Squashing A from N-D to 2-D  
        q = sizeB(1);
        s = sizeB(2);
        lengthorder = length(order);
        collapsedsize = sizeA(order(1:lengthorder-1)); 
        n = prod(collapsedsize);
        a = reshape(a, [n, q]); % N×Q    
        fullsize = [collapsedsize s]; % Size to reshape C back to N-D

    else % B is multi-block, A is single-block (squashing B)

        % STEP 1 - Moving IDB(1) to first dimension
        nd = length(sizeB);
        d1 = idB(1);    
        order = [d1 1:(d1-1) (d1+1):nd]; % Partial shifting
        b = permute(b, order); % Q×...

        % STEP 2 - Squashing B from N-D to 2-D  
        p = sizeA(1);
        q = sizeA(2);
        lengthorder = length(order);
        collapsedsize = sizeB(order(2:lengthorder)); 
        n = prod(collapsedsize);
        b = reshape(b, [q, n]); % Q×N
        fullsize = [p collapsedsize]; % Size to reshape C back to N-D

    end

    % FINAL STEPS - Multiplication, reshape to N-D, inverse permutation
    invorder(order) = 1 : lengthorder;
    c = permute (reshape(a*b, fullsize), invorder);


function [sizeA, sizeB, shiftC, delC, sizeisnew, idA, idB, ...
          squashOK, sxtimesOK, timesOK, mtimesOK, sumOK] = ...
                                          sizeval(idA0,idB0, sizeA0,sizeB0)
%SIZEVAL   Evaluation of array size, block size, and IDs
%    Possible values for IDA and IDB:
%        [DA1 DA2], [DB1 DB2]
%        [DA1 DA2], [DB1]
%        [DA1],     [DB1 DB2]
%        [DA1],     [DB1]
%        [DA1 0],   [0 DB1]
%        [0 DA1],   [DB1 0]
%
%    sizeA/B     Equal to sizeA0/B0 if RESHAPE is not needed in MAIN 1
%    shiftC, delC    Variables controlling MAIN 3.
%    sizeisnew   1x2 logical array; activates reshaping of A and B.
%    idA/B       May change only if squashOK ~= 0
%    squashOK    If only A or B is a multi-block array (M-B) and the other
%                is single-block (1-B), it will be rearranged from N-D to
%                2-D. If both A and B are 1-B or M-B arrays, squashOK = 0.
%                If only A (or B) is a M-B array, squashOK = 1 (or 2).
%    sxtimesOK, timesOK, mtimesOK    Flags controlling MAIN 2 (TRUE/FALSE).
%    sumOK       Dimension along which SUM is performed. If SUM is not
%                needed, sumOK = 0.

% Initializing output arguments

    idA = idA0;
    idB = idB0;
     squashOK = 0;
    sxtimesOK = false;
      timesOK = false;
     mtimesOK = false;
        sumOK = 0;
    shiftC = 0;
    delC = 0;

% Checking for gross input errors

    NidA = numel(idA);
    NidB = numel(idB);
    idA1 = idA(1);
    idB1 = idB(1);
    if  NidA>2 || NidB>2 || NidA==0 || NidB==0 || ...
           ~isreal(idA1) ||    ~isreal(idB1)   || ...
        ~isnumeric(idA1) || ~isnumeric(idB1)   || ...
                 0>idA1  ||          0>idB1    || ... % negative 
         idA1~=fix(idA1) ||  idB1~=fix(idB1)   || ... % non-integer
         ~isfinite(idA1) ||  ~isfinite(idB1) % Inf or NaN               
        error('MULTIPROD:InvalidDimensionArgument', ...
        ['Internal-dimension arguments (e.g., [IDA1 IDA2]) must\n', ...
         'contain only one or two non-negative finite integers']);
    end

% Checking Syntaxes containing zeros (4b/c)

    declared_outer = false;
    idA2 = idA(NidA); % It may be IDA1 = IDA2 (1-D block)
    idB2 = idB(NidB);

    if any(idA==0) || any(idB==0)
        
        % "Inner products": C = MULTIPROD(A, B, [0 DA1], [DB1 0])
        if idA1==0 && idA2>0 && idB1>0 && idB2==0
            idA1 = idA2;
            idB2 = idB1;
        % "Outer products": C = MULTIPROD(A, B, [DA1 0], [0 DB1]) 
        elseif idA1>0 && idA2==0 && idB1==0 && idB2>0
            declared_outer = true;
            idA2 = idA1;
            idB1 = idB2;
        else
            error('MULTIPROD:InvalidDimensionArgument', ...
            ['Misused zeros in the internal-dimension arguments\n', ...
            '(see help heads 4b and 4c)']);
        end
        NidA = 1; 
        NidB = 1;
        idA = idA1;
        idB = idB1;

    elseif (NidA==2 && idA2~=idA1+1) || ...  % Non-adjacent IDs
           (NidB==2 && idB2~=idB1+1)
        error('MULTIPROD:InvalidDimensionArgument', ...
        ['If an array contains 2-D blocks, its two internal dimensions', ... 
        'must be adjacent (e.g. IDA2 == IDA1+1)']);
    end

% ESC - Case for which no reshaping is needed (both A and B are scalars)

    scalarA = isequal(sizeA0, [1 1]);
    scalarB = isequal(sizeB0, [1 1]);
    if scalarA && scalarB
        sizeA = sizeA0;
        sizeB = sizeB0;
        sizeisnew = [false false];
        timesOK = true; return
    end

% Computing and checking adjusted sizes
% The lengths of ADJSIZEA and ADJSIZEB must be >= IDA(END) and IDB(END)

    NsA = idA2 - length(sizeA0); % Number of added trailing singletons
    NsB = idB2 - length(sizeB0);
    adjsizeA = [sizeA0 ones(1,NsA)];
    adjsizeB = [sizeB0 ones(1,NsB)];
    extsizeA = adjsizeA([1:idA1-1, idA2+1:end]); % Size of EDs
    extsizeB = adjsizeB([1:idB1-1, idB2+1:end]);
    p = adjsizeA(idA1);
    q = adjsizeA(idA2);
    r = adjsizeB(idB1);
    s = adjsizeB(idB2);    
    scalarsinA = (p==1 && q==1);
    scalarsinB = (r==1 && s==1);
    singleA = all(extsizeA==1);
    singleB = all(extsizeB==1);
    if q~=r && ~scalarsinA && ~scalarsinB && ~declared_outer
       error('MULTIPROD:InnerDimensionsMismatch', ...
             'Inner matrix dimensions must agree.');
    end

% STEP 1/3 - DIMENSION SHIFTING (FIRST STEP OF AX)
%   Pipeline 1 (using TIMES) never needs left, and may need right shifting.
%   Pipeline 2 (using MTIMES) may need left shifting of A and right of B.

    shiftA = 0;
    shiftB = 0;
    diffBA = idB1 - idA1;    
    if scalarA % Do nothing
    elseif singleA && ~scalarsinB, shiftA = -idA1 + 1; %  Left shifting A
    elseif idB1 > idA1,            shiftA = diffBA;    % Right shifting A        
    end    
    if scalarB % Do nothing
    elseif singleB && ~scalarsinA, shiftB = -idB1 + 1; %  Left shifting B
    elseif idA1 > idB1,            shiftB = -diffBA;   % Right shifting B
    end

% STEP 2/3 - SELECTION OF PROPER ENGINE AND BLOCK SIZE ADJUSTMENTS

    addA  = 0; addB  = 0;
    delA  = 0; delB  = 0;
    swapA = 0; swapB = 0;
    idC1 = max(idA1, idB1);
    idC2 = idC1 + 1;
    checktimes = false;

    if (singleA||singleB) &&~scalarsinA &&~scalarsinB % Engine using MTIMES

        if singleA && singleB 
            mtimesOK = true;
            shiftC=idC1-1; % Right shifting C
            idC1=1; idC2=2;
        elseif singleA
            squashOK = 2;
            idB = [idB1, idB1+1] + shiftB;
        else % singleB
            squashOK = 1;
            idA = [idA1, idA1+1] + shiftA;
        end

        if NidA==2 && NidB==2 % 1) 2-D BLOCKS BY 2-D BLOCKS
            % OK 
        elseif NidA==2        % 2) 2-D BLOCKS BY 1-D BLOCKS
            addB=idB1+1; delC=idC2;
        elseif NidB==2        % 3) 1-D BLOCKS BY 2-D BLOCKS
            addA=idA1; delC=idC1;
        else                  % 4) 1-D BLOCKS BY 1-D BLOCKS
            if declared_outer
                addA=idA1+1; addB=idB1;
            else
                addA=idA1; addB=idB1+1; delC=idC2;
            end
        end    

    else % Engine using TIMES (also used if SCALARA || SCALARB)
        
        sxtimesOK = true;

        if NidA==2 && NidB==2 % 1) 2-D BLOCKS BY 2-D BLOCKS

            if scalarA || scalarB
                timesOK=true;                
            elseif scalarsinA && scalarsinB % scal-by-scal
                checktimes=true;
            elseif scalarsinA || scalarsinB || ... % scal-by-mat
                (q==1 && r==1)  % vec-by-vec ("outer")
            elseif p==1 && s==1 % vec-by-vec ("inner")
                swapA=idA1; sumOK=idC1; checktimes=true;
            elseif s==1 % mat-by-vec
                swapB=idB1; sumOK=idC2;
            elseif p==1 % vec-by-mat
                swapA=idA1; sumOK=idC1;
            else % mat-by-mat
                addA=idA2+1; addB=idB1; sumOK=idC2; delC=idC2;
            end

        elseif NidA==2 % 2) 2-D BLOCKS BY 1-D BLOCKS

            if scalarA || scalarB
                timesOK=true;                
            elseif scalarsinA && scalarsinB % scal-by-scal
                addB=idB1; checktimes=true;
            elseif scalarsinA % scal-by-vec
                delA=idA1;
            elseif scalarsinB % mat-by-scal
                addB=idB1;
            elseif p==1 % vec-by-vec ("inner")
                delA=idA1; sumOK=idC1; checktimes=true;
            else % mat-by-vec
                addB=idB1; sumOK=idC2; delC=idC2;
            end

        elseif NidB==2 % 3) 1-D BLOCKS BY 2-D BLOCKS

            if scalarA || scalarB
                timesOK=true;                
            elseif scalarsinA && scalarsinB % scal-by-scal
                addA=idA1+1; checktimes=true;
            elseif scalarsinB % vec-by-scal
                delB=idB2;
            elseif scalarsinA % scal-by-mat
                addA=idA1+1;
            elseif s==1 % vec-by-vec ("inner")
                delB=idB2; sumOK=idC1; checktimes=true;
            else % vec-by-mat
                addA=idA1+1; sumOK=idC1; delC=idC1;
            end

        else % 4) 1-D BLOCKS BY 1-D BLOCKS

            if scalarA || scalarB
                timesOK=true;                
            elseif declared_outer % vec-by-vec ("outer")
                addA=idA1+1; addB=idB1;
            elseif scalarsinA && scalarsinB % scal-by-scal
                checktimes=true;
            elseif scalarsinA || scalarsinB % vec-by-scal
            else % vec-by-vec
                sumOK=idC1; checktimes=true;
            end
        end
    end

% STEP 3/3 - Adjusting the size of A and B. The size of C is adjusted
%            later, because it is not known yet.

    [sizeA, sizeisnew(1)] = adjustsize(sizeA0, shiftA, addA, delA, swapA);
    [sizeB, sizeisnew(2)] = adjustsize(sizeB0, shiftB, addB, delB, swapB);

    if checktimes % Faster than calling BBXFUN
        diff = length(sizeB) - length(sizeA);
        if isequal([sizeA ones(1,diff)], [sizeB ones(1,-diff)])
            timesOK = true;
        end
    end


function [sizeA, sizeisnew] = adjustsize(sizeA0, shiftA, addA, delA, swapA)
% ADJUSTSIZE  Adjusting size of a block array.

    % Dimension shifting (by adding or deleting trailing singleton dim.)
    if     shiftA>0, [sizeA,newA1] = addsing(sizeA0, 1, shiftA);
    elseif shiftA<0, [sizeA,newA1] = delsing(sizeA0, 1,-shiftA); 
    else   sizeA = sizeA0;  newA1  = false;
    end
    % Modifying block size (by adding, deleting, or moving singleton dim.)
    if      addA, [sizeA,newA2] = addsing(sizeA, addA+shiftA, 1); % 1D-->2D 
    elseif  delA, [sizeA,newA2] = delsing(sizeA, delA+shiftA, 1); % 2D-->1D
    elseif swapA, [sizeA,newA2] = swapdim(sizeA,swapA+shiftA); % ID Swapping
    else                 newA2  = false;
    end
    sizeisnew = newA1 || newA2;


function [newsize, flag] = addsing(size0, dim, ns)
%ADDSING   Adding NS singleton dimensions to the size of an array.
%   Warning: NS is assumed to be a positive integer.
%   Example: If the size of A is ..... SIZE0 = [5 9 3]
%            NEWSIZE = ADDSING(SIZE0, 3, 2) is [5 9 1 1 3]

    if dim > length(size0)
        newsize = size0;
        flag = false;
    else 
        newsize = [size0(1:dim-1), ones(1,ns), size0(dim:end)];
        flag = true;
    end


function [newsize, flag] = delsing(size0, dim, ns)
%DELSING   Removing NS singleton dimensions from the size of an array.
%   Warning: Trailing singletons are not removed
%   Example: If the size of A is SIZE0 = [1 1 1 5 9 3]
%            NEWSIZE = DELSING(SIZE, 1, 3) is  [5 9 3]

    if dim > length(size0)-ns % Trailing singletons are not removed
        newsize = size0;
        flag = false;
    else % Trailing singl. added, so NEWSIZE is guaranteed to be 2D or more
        newsize = size0([1:dim-1, dim+ns:end, dim]);
        flag = true;
    end


function [newsize, flag] = swapdim(size0, dim)
%SWAPDIM   Swapping two adjacent dimensions of an array (DIM and DIM+1).
%   Used only when both A and B are multi-block arrays with 2-D blocks.
%   Example: If the size of A is .......... 5×(6×3)
%            NEWSIZE = SWAPIDS(SIZE0, 2) is 5×(3×6)

    newsize = [size0 1]; % Guarantees that dimension DIM+1 exists.
    newsize = newsize([1:dim-1, dim+1, dim, dim+2:end]);
    flag = true;
