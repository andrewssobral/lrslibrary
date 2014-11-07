% Tensor Toolbox (Sandia National Labs)
% Version 2.5 01-FEB-2012
%
% Tensor Toolbox for dense, sparse, and decomposed n-way arrays.
% 
%   cp_als         - Compute a CP decomposition of any type of tensor.
%   cp_apr         - Compute nonnegative CP with alternating Poisson regression.
%   cp_nmu         - Compute nonnegative CP with multiplicative updates.
%   cp_opt         - Fits a CP model to a tensor via optimization.
%   cp_wopt        - Fits a weighted CP model to a tensor via optimization.
%   create_guess   - Creates initial guess for CP or Tucker fitting.
%   create_problem - Create test problems for tensor factorizations.
%   export_data    - Export tensor-related data to a file.
%   import_data    - Import tensor-related data to a file.
%   khatrirao      - Khatri-Rao product of matrices.
%   parafac_als    - Deprecated. Use CP_ALS instead.
%   sptendiag      - Creates a sparse tensor with v on the diagonal.
%   sptenrand      - Sparse uniformly distributed random tensor.
%   sshopm         - Shifted power method for finding a real eigenpair of a real tensor.
%   sshopmc        - Shifted power method for real/complex eigenpair of a real tensor.
%   tendiag        - Creates a tensor with v on the diagonal.
%   teneye         - Create identity tensor of specified size.
%   tenones        - Ones tensor.
%   tenrand        - Uniformly distributed pseudo-random tensor.
%   tenzeros       - Create zeros tensor.
%   tucker_als     - Higher-order orthogonal iteration.
%
% @TENSOR
%   and         - Logical AND (&) for tensors.
%   collapse    - Collapse tensor along specified dimensions.
%   contract    - Contract tensor along two dimensions (array trace).
%   ctranspose  - is not defined for tensors.
%   disp        - Command window display of a tensor.
%   display     - Command window display of a tensor.
%   double      - Convert tensor to double array.
%   end         - Last index of indexing expression for tensor.
%   eq          - Equal (==) for tensors.
%   find        - Find subscripts of nonzero elements in a tensor.
%   full        - Convert to a (dense) tensor.
%   ge          - Greater than or equal (>=) for tensors.
%   gt          - Greater than (>) for tensors.
%   innerprod   - Efficient inner product with a tensor.
%   isequal     - for tensors.
%   issymmetric - Verify that a tensor X is symmetric in specified modes.
%   ldivide     - Left array divide for tensor.
%   le          - Less than or equal (<=) for tensor.
%   lt          - Less than (<) for tensor.
%   minus       - Binary subtraction (-) for tensors.
%   mldivide    - Slash left division for tensors.
%   mrdivide    - Slash right division for tensors.
%   mtimes      - tensor-scalar multiplication.
%   mttkrp      - Matricized tensor times Khatri-Rao product for tensor.
%   ndims       - Return the number of dimensions of a tensor.
%   ne          - Not equal (~=) for tensors.
%   nnz         - Number of nonzeros for tensors. 
%   norm        - Frobenius norm of a tensor.
%   not         - Logical NOT (~) for tensors.
%   nvecs       - Compute the leading mode-n vectors for a tensor.
%   or          - Logical OR (|) for tensors.
%   permute     - Permute tensor dimensions.
%   plus        - Binary addition (+) for tensors. 
%   power       - Elementwise power (.^) operator for a tensor.
%   rdivide     - Right array divide for tensors.
%   reshape     - Change tensor size.
%   scale       - Scale along specified dimensions of tensor.
%   size        - Tensor dimensions.
%   squeeze     - Remove singleton dimensions from a tensor.
%   subsasgn    - Subscripted assignment for a tensor.
%   subsref     - Subscripted reference for tensors.
%   symmetrize  - Symmetrize a tensor X in specified modes.
%   tenfun      - Apply a function to each element in a tensor.
%   tensor      - Create tensor.
%   times       - Array multiplication for tensors.
%   transpose   - is not defined on tensors.
%   ttm         - Tensor times matrix.
%   ttsv        - Tensor times same vector in multiple modes.
%   ttt         - Tensor mulitplication (tensor times tensor).
%   ttv         - Tensor times vector.
%   uminus      - Unary minus (-) for tensors.
%   uplus       - Unary plus (+) for tensors.
%   xor         - Logical EXCLUSIVE OR for tensors.
%
% @SPTENSOR
%   and        - Logical AND (&) for sptensors.
%   collapse   - Collapse sparse tensor along specified dimensions.
%   contract   - Contract sparse tensor along two dimensions (array trace).
%   ctranspose - is not defined for sparse tensors.
%   disp       - Command window display of a sparse tensor.
%   display    - Command window display of a sparse tensor.
%   divide     - Divide an SPTENSOR by a nonnegative KTENSOR.
%   double     - Converts a sparse tensor to a dense multidimensional array.
%   elemfun    - Manipulate the nonzero elements of a sparse tensor.
%   end        - Last index of indexing expression for sparse tensor.
%   eq         - Equal (==) for sptensors.
%   find       - Find subscripts of nonzero elements in a sparse tensor.
%   full       - Convert a sparse tensor to a (dense) tensor.
%   ge         - Greater than or equal for sptensors.
%   gt         - Greater than for sptensors.
%   innerprod  - Efficient inner product with a sparse tensor.
%   isequal    - for sptensors.
%   ldivide    - Array right division for sparse tensors.
%   le         - Less than or equal for sptensors.
%   lt         - Less than for sptensors.
%   minus      - Binary subtraction for sparse tensors. 
%   mldivide   - Slash left division for sparse tensors.
%   mrdivide   - Slash right division for sparse tensors.
%   mtimes     - sptensor-scalar multiplication.
%   mttkrp     - Matricized tensor times Khatri-Rao product for sparse tensor.
%   ndims      - Number of dimensions of a sparse tensor.
%   ne         - Not equal (~=) for sptensors.
%   nnz        - Number of nonzeros in sparse tensor.
%   norm       - Frobenius norm of a sparse tensor.
%   not        - Logical NOT (~) for sptensors.
%   nvecs      - Compute the leading mode-n vectors for a sparse tensor.
%   ones       - Replace nonzero elements of sparse tensor with ones.
%   or         - Logical OR (|) for sptensors.
%   permute    - Rearrange the dimensions of a sparse tensor.
%   plus       - Binary addition for sparse tensors. 
%   rdivide    - Array right division for sparse tensors.
%   reshape    - Reshape sparse tensor.
%   scale      - Scale along specified dimensions for sparse tensors.
%   size       - Sparse tensor dimensions.
%   spmatrix   - Converts a two-way sparse tensor to sparse matrix.
%   sptensor   - Create a sparse tensor.
%   squeeze    - Remove singleton dimensions from a sparse tensor.
%   subsasgn   - Subscripted assignment for sparse tensor.
%   subsref    - Subscripted reference for a sparse tensor.
%   times      - Array multiplication for sparse tensors.
%   transpose  - is not defined on sparse tensors.
%   ttm        - Sparse tensor times matrix.
%   ttt        - Sparse tensor times sparse tensor.
%   ttv        - Sparse tensor times vector.
%   uminus     - Unary minus (-) for sptensor.
%   uplus      - Unary plus (+) for sptensor.
%   xor        - Logical XOR for sptensors.
%
% @TTENSOR
%   disp      - Command window display of a ttensor.
%   display   - Command window display of a ttensor.
%   double    - Convert ttensor to double array.
%   end       - Last index of indexing expression for ttensor.
%   full      - Convert a ttensor to a (dense) tensor.
%   innerprod - Efficient inner product with a ttensor.
%   isequal   - True if each component of two ttensor's is numerically equal.
%   mtimes    - Implement scalar multiplication for a ttensor.
%   mttkrp    - Matricized tensor times Khatri-Rao product for ttensor.
%   ndims     - Return the number of dimensions for a ttensor.
%   norm      - Norm of a ttensor.
%   nvecs     - Compute the leading mode-n vectors for a ttensor.
%   permute   - Permute dimensions for a ttensor.
%   size      - Size of a ttensor.
%   subsasgn  - Subscripted reference for a ttensor.
%   subsref   - Subscripted reference for a ttensor.
%   ttensor   - Tensor stored as a Tucker operator (decomposed).
%   ttm       - Tensor times matrix for ttensor.
%   ttv       - Tensor times vector for ttensor.
%   uminus    - Unary minus for ttensor.
%   uplus     - Unary plus for ttensor.
%
% @KTENSOR
%   arrange      - Arranges the rank-1 components of a ktensor.
%   datadisp     - Special display of a ktensor.
%   disp         - Command window display for a ktensor.
%   display      - Command window display for a ktensor.
%   double       - Convert a ktensor to a double array.
%   end          - Last index of indexing expression for ktensor.
%   extract      - Creates a new ktensor with only the specified components.
%   fixsigns     - Fix sign ambiguity of a ktensor.
%   full         - Convert a ktensor to a (dense) tensor.
%   innerprod    - Efficient inner product with a ktensor.
%   isequal      - True if each component of two ktensor's is numerically equal.
%   ktensor      - Tensor stored as a Kruskal operator (decomposed).
%   minus        - Binary subtraction for ktensor.  
%   mtimes       - Implement A*B (scalar multiply) for ktensor.
%   mttkrp       - Matricized tensor times Khatri-Rao product for ktensor.
%   ncomponents  - Number of components for a ktensor.
%   ndims        - Number of dimensions for a ktensor.
%   norm         - Frobenius norm of a ktensor.
%   normalize    - Normalizes the columns of the factor matrices.
%   nvecs        - Compute the leading mode-n vectors for a ktensor.
%   permute      - Permute dimensions of a ktensor.
%   plus         - Binary addition for ktensor.
%   redistribute - Distribute lambda values to a specified mode. 
%   score        - Checks if two ktensors match except for permutation.
%   size         - Size of ktensor.
%   subsasgn     - Subscripted assignement for ktensor.
%   subsref      - Subscripted reference for a ktensor.
%   times        - Element-wise multiplication for ktensor.
%   tocell       - Convert X to a cell array.
%   ttm          - Tensor times matrix for ktensor.
%   ttv          - Tensor times vector for ktensor.
%   uminus       - Unary minus for ktensor. 
%   uplus        - Unary plus for a ktensor. 
%
% @TENMAT
%   ctranspose - Complex conjugate transpose for tenmat.
%   disp       - Command window display of a matricized tensor (tenmat).
%   display    - Command window display of a tenmat.
%   double     - Convert tenmat to double array.
%   end        - Last index of indexing expression for tenmat.
%   minus      - Binary subtraction (-) for tenmat.
%   mtimes     - Multiplies two tenmat objects.
%   norm       - Frobenius norm of a tenmat.
%   plus       - Binary addition (+) for tenmat. 
%   size       - Size of tenmat.
%   subsasgn   - Subscripted assignment for tenmat.  
%   subsref    - Subscripted reference for tenmat.
%   tenmat     - Create a matricized tensor.
%   tsize      - Tensor size of tenmat.
%   uminus     - Unary minus (-) for tenmat.
%   uplus      - Unary plus (+) for tenmat.
%
% @SPTENMAT
%   aatx     - Implicitly compute A * A' * x for sptenmat.
%   disp     - Command window display of a sptenmat.
%   display  - Command window display of a sptenmat.
%   double   - Convert a sptenmat to a sparse matrix.
%   end      - Last index of indexing expression for sptenmat.
%   full     - Convert a sptenmat to a (dense) tenmat.
%   nnz      - Return number of nonzeros in a sptenmat.
%   norm     - Frobenius norm of a sptenmat.
%   size     - Return size of sptenmat.
%   sptenmat - Matricized sparse tensor stored as a sparse 2D array.
%   subsasgn - Subscripted assignment for sptenmat.  
%   subsref  - Subscripted reference for a sptenmat.
%   tsize    - Tensor size of sptenmat.
%   uminus   - Unary minus (-) for sptenmat.
%   uplus    - Unary plus (+) for sptenmat.
%
% MET
%   ttm_me           - Memory-efficient sptensor times matrix.
%   ttm_me_mem       - Estimates intermediate memory comsumption for ttm_me.
%   ttm_me_partition - Finds best order for ttm_me.
%   tucker_me        - Memory-efficient Tucker higher-order orthogonal iteration.
%   tucker_me_test   - Very simple tests of tucker_me.
