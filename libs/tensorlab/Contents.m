% TENSORLAB
% Version 2.01, 2014-02-05
%
% BLOCK TERM DECOMPOSITION
% Algorithms
%    btd_minf     - BTD by unconstrained nonlinear optimization.
%    btd_nls      - BTD by nonlinear least squares.
% Initialization
%    btd_rnd      - Pseudorandom initialization for BTD.
% Utilities
%    btdgen       - Generate full tensor given a BTD.
%    btdres       - Residual of a BTD.
%
% CANONICAL POLYADIC DECOMPOSITION
% Algorithms
%    cpd          - Canonical polyadic decomposition.
%    cpd_als      - CPD by alternating least squares.
%    cpd_minf     - CPD by unconstrained nonlinear optimization.
%    cpd_nls      - CPD by nonlinear least squares.
%    cpd3_sd      - CPD by simultaneous diagonalization.
%    cpd3_sgsd    - CPD by simultaneous generalized Schur decomposition.
% Initialization
%    cpd_gevd     - CPD by a generalized eigenvalue decomposition.
%    cpd_rnd      - Pseudorandom initialization for CPD.
% Line and plane search
%    cpd_aels     - CPD approximate enhanced line search.
%    cpd_els      - CPD exact line search.
%    cpd_eps      - CPD exact plane search.
%    cpd_lsb      - CPD line search by Bro.
% Utilities
%    cpderr       - Errors between factor matrices in a CPD.
%    cpdgen       - Generate full tensor given a polyadic decomposition.
%    cpdres       - Residual of a polyadic decomposition.
%    rankest      - Estimate rank.
%
% COMPLEX OPTIMIZATION
% Nonlinear least squares
%    nls_gncgs    - Nonlinear least squares by Gauss-Newton with CG-Steihaug.
%    nls_gndl     - Nonlinear least squares by Gauss-Newton with dogleg trust region.
%    nls_lm       - Nonlinear least squares by Levenberg-Marquardt.
%    nlsb_gndl    - Bound-constrained NLS by projected Gauss-Newton dogleg TR.
% Unconstrained nonlinear optimization
%    minf_lbfgs   - Minimize a function by L-BFGS with line search.
%    minf_lbfgsdl - Minimize a function by L-BFGS with dogleg trust region.
%    minf_ncg     - Minimize a function by nonlinear conjugate gradient.
% Utilities
%    deriv        - Approximate gradient and Jacobian.
%    ls_mt        - Strong Wolfe line search by More-Thuente.
%    mpcg         - Modified preconditioned conjugate gradients method.
%
% LOW MULTILINEAR RANK APPROXIMATION
% Algorithms
%    lmlra        - Low multilinear rank approximation.
%    lmlra_hooi   - LMLRA by higher-order orthogonal iteration.
%    lmlra_minf   - LMLRA by unconstrained nonlinear optimization.
%    lmlra_nls    - LMLRA by nonlinear least squares.
%    lmlra3_dgn   - LMLRA by a differential-geometric Newton method.
%    lmlra3_rtr   - LMLRA by a Riemannian trust region method.
%    mlsvd        - (Truncated) multilinear singular value decomposition.
% Initialization
%    lmlra_aca    - LMLRA by adaptive cross-approximation.
%    lmlra_rnd    - Pseudorandom initialization for LMLRA.
% Utilities
%    lmlraerr     - Errors between factor matrices in a LMLRA.
%    lmlragen     - Generate full tensor given a core tensor and factor matrices.
%    lmlrares     - Residual of a LMLRA.
%    mlrank       - Multilinear rank.
%    mlrankest    - Estimate multilinear rank.
%
% STRUCTURED DATA FUSION
% Algorithms
%   sdf_minf      - Structured data fusion by unconstrained nonlinear optimization.
%   sdf_nls       - Structured data fusion by nonlinear least squares.
% Structure
%   struct_abs        - Absolute value.
%   struct_band       - Band matrix.
%   struct_cell2mat   - Convert the contents of a cell array into a matrix.
%   struct_conj       - Complex conjugate.
%   struct_ctranspose - Complex conjugate transpose.
%   struct_diag       - Diagonal matrix.
%   struct_gram       - Gramian matrix.
%   struct_hankel     - Hankel matrix.
%   struct_inv        - Matrix inverse.
%   struct_invsqrtm   - Matrix inverse square root.
%   struct_invtransp  - Matrix inverse transpose.
%   struct_LL1        - Structure of third factor matrix in a rank-(Lr,Lr,1) BTD.
%   struct_log        - Natural logarithm.
%   struct_matvec     - Matrix-vector and matrix-matrix product.
%   struct_nonneg     - Nonnegative array.
%   struct_normalize  - Normalize columns to unit norm.
%   struct_orth       - Rectangular matrix with orthonormal columns.
%   struct_plus       - Plus.
%   struct_poly       - Matrix with columns as polynomials.
%   struct_power      - Array power.
%   struct_rational   - Matrix with columns as rational functions.
%   struct_rbf        - Matrix with columns as sums of Gaussian RBF kernels.
%   struct_sigmoid    - Constrain array elements to an interval.
%   struct_sqrt       - Square root.
%   struct_sum        - Sum of elements.
%   struct_times      - Times.
%   struct_toeplitz   - Toeplitz matrix.
%   struct_transpose  - Transpose.
%   struct_tridiag    - Tridiagonal matrix.
%   struct_tril       - Lower triangular matrix.
%   struct_triu       - Upper triangular matrix.
%   struct_vander     - Vandermonde matrix.
%
% UTILITIES
% Clustering
%    gap          - Optimal clustering based on the gap statistic.
%    kmeans       - Cluster multivariate data using the k-means++ algorithm.
% Polynomials
%    polymin      - Minimize a polynomial.
%    polymin2     - Minimize bivariate and real polyanalytic polynomials.
%    polyval2     - Evaluate bivariate and univariate polyanalytic polynomials.
%    polysol2     - Solve a system of two bivariate polynomials.
%    ratmin       - Minimize a rational function.
%    ratmin2      - Minimize bivariate and real polyanalytic rational functions.
% Statistics
%    cum4         - Fourth-order cumulant tensor.
%    scov         - Shifted covariance matrices.
% Tensors
%    dotk         - Dot product in K-fold precision.
%    fmt          - Format data set.
%    frob         - Frobenius norm.
%    ful          - Convert formatted data set to an array.
%    kr           - Khatri-Rao product.
%    kron         - Kronecker product.
%    mat2tens     - Tensorize a matrix.
%    mtkronprod   - Compute a matricized tensor Kronecker product.
%    mtkrprod     - Compute a matricized tensor Khatri-Rao product.
%    noisy        - Generate a noisy version of a given array.
%    sumk         - Summation in K-fold precision.
%    tens2mat     - Matricize a tensor.
%    tens2vec     - Vectorize a tensor.
%    tmprod       - Mode-n tensor-matrix product.
%    vec2tens     - Tensorize a vector.
% Visualization
%    slice3       - Visualize a third-order tensor with slices.
%    spy3         - Visualize a third-order tensor's sparsity pattern.
%    surf3        - Visualize a third-order tensor with surfaces.
%    voxel3       - Visualize a third-order tensor with voxels.
