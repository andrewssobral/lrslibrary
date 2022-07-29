[![View LRSLibrary on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/48404-lrslibrary)

Last Page Update: **29/07/2022**, Previous Page Update: **07/03/2020**

Latest Library Version: **1.0.11** (see Release Notes for more info)

LRSLibrary
----------
*Low-Rank and Sparse* tools for Background Modeling and Subtraction in Videos.

The *LRSLibrary* provides a collection of **low-rank and sparse decomposition** algorithms in MATLAB. The library was designed for moving object detection in videos, but it can be also used for other computer vision and machine learning problems (for more information, please see [here](https://en.wikipedia.org/wiki/Robust_principal_component_analysis#Applications) and [here](https://en.wikipedia.org/wiki/Low-rank_approximation#Applications)). Currently the LRSLibrary offers more than **100** algorithms based on *matrix* and *tensor* methods. The LRSLibrary was tested successfully in several MATLAB versions (e.g. R2014, R2015, R2016, R2017, on both x86 and x64 versions). It requires minimum **R2014b**.

<p align="center"><img src="https://raw.githubusercontent.com/andrewssobral/lrslibrary/master/figs/lrs_results2.png" /></p>
<p align="center"><img src="https://raw.githubusercontent.com/andrewssobral/lrslibrary/master/figs/lrs-opt.gif" /></p>

See also:
```
Presentation about Matrix and Tensor Tools for Computer Vision
http://www.slideshare.net/andrewssobral/matrix-and-tensor-tools-for-computer-vision

MTT: Matlab Tensor Tools for Computer Vision
https://github.com/andrewssobral/mtt

IMTSL: Incremental and Multi-feature Tensor Subspace Learning
https://github.com/andrewssobral/imtsl
```

Citation
---------
If you use this library for your publications, please cite it as:
```
@incollection{lrslibrary2015,
author    = {Sobral, Andrews and Bouwmans, Thierry and Zahzah, El-hadi},
title     = {LRSLibrary: Low-Rank and Sparse tools for Background Modeling and Subtraction in Videos},
booktitle = {Robust Low-Rank and Sparse Matrix Decomposition: Applications in Image and Video Processing},
publisher = {CRC Press, Taylor and Francis Group.}
year      = {2015}
}
```
Additional reference:
```
@article{bouwmans2015,
author    = {Bouwmans, Thierry and Sobral, Andrews and Javed, Sajid and Jung, Soon Ki and Zahzah, El-hadi},
title     = {Decomposition into Low-rank plus Additive Matrices for Background/Foreground Separation: {A} Review for a Comparative Evaluation with a Large-Scale Dataset},
journal   = {CoRR},
volume    = {abs/1511.01245}
year      = {2015},
url       = {http://arxiv.org/abs/1511.01245}
}
```

## Stargazers over time

[![Stargazers over time](https://starchart.cc/andrewssobral/lrslibrary.svg)](https://starchart.cc/andrewssobral/lrslibrary)

Install
---
Just do the following steps:

* First, clone the repository: 

```$ git clone --recursive https://github.com/andrewssobral/lrslibrary.git```

* Then, open your MATLAB and run the following setup script: 

```>> lrs_setup```

That's all!

GUI
---
The *LRSLibrary* provides an easy-to-use graphical user interface (GUI) for background modeling and subtraction in videos. First, run the setup script **lrs_setup** (or **run('C:/lrslibrary/lrs_setup')**), then run **lrs_gui**, and enjoy it!

<p align="center">(Click in the image to see the video)</p>
<p align="center">
<a href="https://www.youtube.com/watch?v=zziJ7-WnvV8" target="_blank">
<img src="https://raw.githubusercontent.com/andrewssobral/lrslibrary/master/figs/lrslibrary_gui2.png" width="500" border="0" />
</a>
</p>

Each algorithm is classified by its cpu time consumption with the following icons:
<p align="center"><img src="https://raw.githubusercontent.com/andrewssobral/lrslibrary/master/figs/time_legend.png" width="300" /></p>

The algorithms were grouped in eight categories: **RPCA** for Robust PCA, **ST** for Subspace Tracking, **MC** for Matrix Completion, **TTD** for Three-Term Decomposition, **LRR** for Low-Rank Representation, **NMF** for Non-negative Matrix Factorization, **NTF** for Non-negative Tensor Factorization, or **TD** for standard Tensor Decomposition.

List of the algorithms available in LRSLibrary
----------------------------------------------
* RPCA: Robust PCA (44)
* * RPCA: Robust Principal Component Analysis [(De la Torre and Black, 2001)](http://users.salleurl.edu/~ftorre/papers/rpca/rpca.pdf) [website](http://users.salleurl.edu/~ftorre/papers/rpca2.html)

* * PCP: Principal Component Pursuit [(Candes et al. 2009)](http://arxiv.org/abs/0912.3599)

* * FPCP: Fast PCP [(Rodriguez and Wohlberg, 2013)](http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=6738015)

* * R2PCP: Riemannian Robust Principal Component Pursuit [(Hintermüller and Wu, 2014)](http://link.springer.com/article/10.1007/s10851-014-0527-y)

* * AS-RPCA: Active Subspace: Towards Scalable Low-Rank Learning [(Liu and Yan, 2012)](http://dl.acm.org/citation.cfm?id=2421487)

* * ALM: Augmented Lagrange Multiplier [(Tang and Nehorai 2011)](http://dx.doi.org/10.1109/CISS.2011.5766144)

* * EALM: Exact ALM [(Lin et al. 2009)](http://arxiv.org/abs/1009.5055) [website](http://perception.csl.illinois.edu/matrix-rank/sample_code.html)

* * IALM: Inexact ALM [(Lin et al. 2009)](http://arxiv.org/abs/1009.5055)  [website](http://perception.csl.illinois.edu/matrix-rank/sample_code.html)

* * IALM_LMSVDS: IALM with LMSVDS [(Liu et al. 2012)](http://epubs.siam.org/doi/abs/10.1137/120871328)

* * IALM_BLWS: IALM with BLWS [(Lin and Wei, 2010)](http://arxiv.org/abs/1012.0365)

* * APG_PARTIAL: Partial Accelerated Proximal Gradient [(Lin et al. 2009)](http://arxiv.org/abs/1009.5055)  [website](http://perception.csl.illinois.edu/matrix-rank/sample_code.html)

* * APG: Accelerated Proximal Gradient [(Lin et al. 2009)](http://arxiv.org/abs/1009.5055)  [website](http://perception.csl.illinois.edu/matrix-rank/sample_code.html)

* * DUAL: Dual RPCA [(Lin et al. 2009)](http://arxiv.org/abs/1009.5055)  [website](http://perception.csl.illinois.edu/matrix-rank/sample_code.html)

* * SVT: Singular Value Thresholding [(Cai et al. 2008)](http://arxiv.org/abs/0810.3286)  [website](http://perception.csl.illinois.edu/matrix-rank/sample_code.html)

* * ADM: Alternating Direction Method [(Yuan and Yang, 2009)](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.400.8797) 

* * LSADM: LSADM [(Goldfarb et al. 2010)](http://arxiv.org/abs/0912.4571)

* * L1F: L1 Filtering [(Liu et al. 2011)](http://arxiv.org/abs/1108.5359)

* * DECOLOR: Contiguous Outliers in the Low-Rank Representation [(Zhou et al. 2011)](http://arxiv.org/abs/1109.0882) [website1](https://sites.google.com/site/eeyangc/software/decolor) [website2](https://fling.seas.upenn.edu/~xiaowz/dynamic/wordpress/?p=144)

* * RegL1-ALM: Low-Rank Matrix Approximation under Robust L1-Norm [(Zheng et al. 2012)](https://sites.google.com/site/yinqiangzheng/home/zheng_CVPR12_robust%20L1-norm%20low-rank%20matrix%20factorization.pdf) [website](https://sites.google.com/site/yinqiangzheng/)

* * GA: Grassmann Average [(Hauberg et al. 2014)](http://files.is.tue.mpg.de/black/papers/RGA2014.pdf) [website](http://ps.is.tuebingen.mpg.de/project/Robust_PCA)

* * GM: Grassmann Median [(Hauberg et al. 2014)](http://files.is.tue.mpg.de/black/papers/RGA2014.pdf) [website](http://ps.is.tuebingen.mpg.de/project/Robust_PCA)

* * TGA: Trimmed Grassmann Average [(Hauberg et al. 2014)](http://files.is.tue.mpg.de/black/papers/RGA2014.pdf) [website](http://ps.is.tuebingen.mpg.de/project/Robust_PCA)

* * STOC-RPCA: Online Robust PCA via Stochastic Optimization [(Feng et al. 2013)](http://guppy.mpe.nus.edu.sg/~mpexuh/papers/Stochastic_online_pca.pdf) [website](https://sites.google.com/site/jshfeng/)

* * MoG-RPCA: Mixture of Gaussians RPCA [(Zhao et al. 2014)](http://jmlr.org/proceedings/papers/v32/zhao14.pdf) [website](http://www.cs.cmu.edu/~deyum/index.htm)

* * noncvxRPCA: Robust PCA via Nonconvex Rank Approximation [(Kang et al. 2015)](http://arxiv.org/abs/1511.05261)

* * NSA1: Non-Smooth Augmented Lagrangian v1 [(Aybat et al. 2011)](http://arxiv.org/abs/1105.2126)

* * NSA2: Non-Smooth Augmented Lagrangian v2 [(Aybat et al. 2011)](http://arxiv.org/abs/1105.2126)

* * PSPG: Partially Smooth Proximal Gradient [(Aybat et al. 2012)](http://arxiv.org/abs/1309.6976)

* * flip-SPCP-sum-SPG: Flip-Flop version of Stable PCP-sum solved by Spectral Projected Gradient [(Aravkin et al. 2014)](https://github.com/stephenbeckr/fastRPCA)

* * flip-SPCP-max-QN: Flip-Flop version of Stable PCP-max solved by Quasi-Newton [(Aravkin et al. 2014)](https://github.com/stephenbeckr/fastRPCA)

* * Lag-SPCP-SPG: Lagrangian SPCP solved by Spectral Projected Gradient [(Aravkin et al. 2014)](https://github.com/stephenbeckr/fastRPCA)

* * Lag-SPCP-QN: Lagrangian SPCP solved by Quasi-Newton [(Aravkin et al. 2014)](https://github.com/stephenbeckr/fastRPCA)

* * FW-T: SPCP solved by Frank-Wolfe method [(Mu et al. 2014)](http://arxiv.org/abs/1403.7588) [website](https://sites.google.com/site/mucun1988/publi)

* * BRPCA-MD: Bayesian Robust PCA with Markov Dependency [(Ding et al. 2011)](http://people.ee.duke.edu/~lcarin/LRS_09.pdf) [website](http://people.ee.duke.edu/~lcarin/BCS.html)

* * BRPCA-MD-NSS: BRPCA-MD with Non-Stationary Noise [(Ding et al. 2011)](http://people.ee.duke.edu/~lcarin/LRS_09.pdf) [website](http://people.ee.duke.edu/~lcarin/BCS.html)

* * VBRPCA: Variational Bayesian RPCA [(Babacan et al. 2011)](http://arxiv.org/abs/1102.5288)

* * PRMF: Probabilistic Robust Matrix Factorization [(Wang et al. 2012)](http://winsty.net/papers/prmf.pdf) [website](http://winsty.net/prmf.html)

* * OPRMF: Online PRMF [(Wang et al. 2012)](http://winsty.net/papers/prmf.pdf) [website](http://winsty.net/prmf.html)

* * MBRMF: Markov BRMF [(Wang and Yeung, 2013)](http://winsty.net/papers/brmf.pdf) [website](http://winsty.net/brmf.html)

* * TFOCS-EC: TFOCS with equality constraints [(Becker et al. 2011)](https://github.com/cvxr/TFOCS/raw/gh-pages/TFOCS.pdf) [website](http://cvxr.com/tfocs/demos/rpca/)

* * TFOCS-IC: TFOCS with inequality constraints [(Becker et al. 2011)](https://github.com/cvxr/TFOCS/raw/gh-pages/TFOCS.pdf) [website](http://cvxr.com/tfocs/demos/rpca/)

* * GoDec: Go Decomposition [(Zhou and Tao, 2011)](http://www.icml-2011.org/papers/41_icmlpaper.pdf) [website](https://sites.google.com/site/godecomposition/home)

* * SSGoDec: Semi-Soft GoDec [(Zhou and Tao, 2011)](http://www.icml-2011.org/papers/41_icmlpaper.pdf) [website](https://sites.google.com/site/godecomposition/home)

* * GreGoDec: Greedy Semi-Soft GoDec Algotithm [(Zhou and Tao, 2013)](http://jmlr.org/proceedings/papers/v31/zhou13b.pdf) [website](https://sites.google.com/site/godecomposition/home)

* ST: Subspace Tracking (3)
* * GRASTA: Grassmannian Robust Adaptive Subspace Tracking Algorithm [(He et al. 2012)](http://www.citeulike.org/user/lambertch/article/12543964) [website](https://sites.google.com/site/hejunzz/grasta)

* * GOSUS: Grassmannian Online Subspace Updates with Structured-sparsity [(Xu et al. 2013)](http://pages.cs.wisc.edu/~jiaxu/projects/gosus/gosus-iccv2013.pdf) [website](http://pages.cs.wisc.edu/~jiaxu/projects/gosus/)

* * pROST: Robust PCA and subspace tracking from incomplete observations using L0-surrogates [(Hage and Kleinsteuber, 2013)](http://arxiv.org/abs/1210.0805) [website](http://www.gol.ei.tum.de/index.php?id=37&L=1)

* * ReProCS: Provable Dynamic Robust PCA or Robust Subspace Tracking [(Narayanamurthy and Vaswani, 2017a)](https://arxiv.org/abs/1705.08948)

* * MEDRoP: Memory Efficient Dynamic Robust PCA [(Narayanamurthy and Vaswani, 2017b)](https://arxiv.org/abs/1712.06061)

* MC: Matrix Completion (15)
* * PG-RMC: Nearly Optimal Robust matrix Completion [(Cherapanamjeri et al. 2016)](https://arxiv.org/abs/1606.07315)

* * FPC: Fixed point and Bregman iterative methods for matrix rank minimization [(Ma et al. 2008)](http://arxiv.org/pdf/0905.1643.pdf) [website](http://www1.se.cuhk.edu.hk/~sqma/FPCA.html)

* * GROUSE: Grassmannian Rank-One Update Subspace Estimation [(Balzano et al. 2010)](http://arxiv.org/pdf/1006.4046.pdf) [website](http://sunbeam.ece.wisc.edu/grouse/)

* * IALM-MC: Inexact ALM for Matrix Completion [(Lin et al. 2009)](https://arxiv.org/abs/1009.5055) [website](http://perception.csl.illinois.edu/matrix-rank/sample_code.html)

* * LMaFit: Low-Rank Matrix Fitting [(Wen et al. 2012)](http://link.springer.com/article/10.1007%2Fs12532-012-0044-1) [website](http://lmafit.blogs.rice.edu/)

* * LRGeomCG: Low-rank matrix completion by Riemannian optimization [(Bart Vandereycken, 2013)](http://web.math.princeton.edu/~bartv/papers/84576.pdf) [website1](http://web.math.princeton.edu/~bartv/matrix_completion.html) [website2](http://www.manopt.org/reference/examples/low_rank_matrix_completion.html)

* * MC_logdet: Top-N Recommender System via Matrix Completion [(Kang et al. 2016)](https://arxiv.org/abs/1601.04800)

* * MC-NMF: Nonnegative Matrix Completion [(Xu et al. 2011)](https://arxiv.org/abs/1103.1168)

* * OP-RPCA: Robust PCA via Outlier Pursuit [(Xu et al. 2012)](http://guppy.mpe.nus.edu.sg/~mpexuh/papers/OutlierPursuit-TIT.pdf) [website](http://guppy.mpe.nus.edu.sg/~mpexuh/publication.html)

* * OptSpace: Matrix Completion from Noisy Entries  [(Keshavan et al. 2009)](http://arxiv.org/pdf/0906.2027v1.pdf) [website](http://web.engr.illinois.edu/~swoh/software/optspace/code.html)

* * OR1MP: Orthogonal rank-one matrix pursuit for low rank matrix completion [(Wang et al. 2015)](https://arxiv.org/abs/1404.1377)

* * RPCA-GD: Robust PCA via Gradient Descent [(Yi et al. 2016)](https://arxiv.org/abs/1605.07784)

* * ScGrassMC: Scaled Gradients on Grassmann Manifolds for Matrix Completion [(Ngo and Saad, 2012)](https://papers.nips.cc/paper/4713-scaled-gradients-on-grassmann-manifolds-for-matrix-completion.pdf)

* * SVP: Guaranteed Rank Minimization via Singular Value Projection [(Meka et al. 2009)](https://arxiv.org/abs/0909.5457)

* * SVT: A singular value thresholding algorithm for matrix completion [(Cai et al. 2008)](http://arxiv.org/pdf/0810.3286.pdf) [website](http://svt.stanford.edu/)

* LRR: Low Rank Recovery (6)
* * EALM: Exact ALM [(Lin et al. 2009)](http://arxiv.org/abs/1009.5055)

* * IALM: Inexact ALM [(Lin et al. 2009)](http://arxiv.org/abs/1009.5055)

* * ADM: Alternating Direction Method [(Lin et al. 2011)](http://arxiv.org/abs/1109.0367)

* * LADMAP: Linearized ADM with Adaptive Penalty [(Lin et al. 2011)](http://arxiv.org/abs/1109.0367)

* * FastLADMAP: Fast LADMAP [(Lin et al. 2011)](http://arxiv.org/abs/1109.0367)

* * ROSL: Robust Orthonormal Subspace Learning [(Shu et al. 2014)](https://dl.dropboxusercontent.com/u/10893363/Homepage/CVPR2014_ROSL.pdf) [website](https://sites.google.com/site/xianbiaoshu/)

* TTD: Three-Term Decomposition (4)
* * 3WD: 3-Way-Decomposition [(Oreifej et al. 2012)](http://www.cs.ucf.edu/~oreifej/papers/3-Way-Decomposition.pdf) [website](http://vision.eecs.ucf.edu/projects/Turbulence/)

* * MAMR: Motion-Assisted Matrix Restoration [(Ye et al. 2015)](http://projects.medialab-tju.org/bf_separation/download/2015_TCSVT.pdf) [website](http://projects.medialab-tju.org/bf_separation/)

* * RMAMR: Robust Motion-Assisted Matrix Restoration [(Ye et al. 2015)](http://projects.medialab-tju.org/bf_separation/download/2015_TCSVT.pdf) [website](http://projects.medialab-tju.org/bf_separation/)

* * ADMM: Alternating Direction Method of Multipliers [(Parikh and Boyd, 2014)](http://projects.medialab-tju.org/bf_separation/download/2015_TCSVT.pdf) [website1](http://stanford.edu/~boyd/admm.html) [website2](http://web.stanford.edu/~boyd/papers/prox_algs.html)

* NMF: Non-Negative Matrix Factorization (14)
* * NMF-MU: NMF solved by Multiplicative Updates

* * NMF-PG: NMF solved by Projected Gradient

* * NMF-ALS: NMF solved by Alternating Least Squares

* * NMF-ALS-OBS: NMF solved by Alternating Least Squares with Optimal Brain Surgeon

* * PNMF: Probabilistic Non-negative Matrix Factorization

* * ManhNMF: Manhattan NMF [(Guan et al. 2013)](http://arxiv.org/abs/1207.3438) [website](https://sites.google.com/site/nmfsolvers/)

* * NeNMF: NMF via Nesterovs Optimal Gradient Method [(Guan et al. 2012)](http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=6166359) [website](https://sites.google.com/site/nmfsolvers/)

* * LNMF: Spatially Localized NMF [(Li et al. 2001)](http://dx.doi.org/10.1109/CVPR.2001.990477)

* * ENMF: Exact NMF [(Gillis and Glineur, 2012)](http://arxiv.org/abs/1009.0880) [website](https://sites.google.com/site/nicolasgillis/code)

* * nmfLS2: Non-negative Matrix Factorization with sparse matrix [(Ji and Eisenstein, 2013)](http://www.cc.gatech.edu/~jeisenst/papers/ji-emnlp-2013.pdf) [website](https://github.com/jiyfeng/tfkld)

* * Semi-NMF: Semi Non-negative Matrix Factorization

* * Deep-Semi-NMF: Deep Semi Non-negative Matrix Factorization [(Trigeorgis et al. 2014)](http://trigeorgis.com/uploads/downloads/file/1/cameraready.pdf) [website](http://trigeorgis.com/papers/deepseminmfmodel-2014)

* * iNMF: Incremental Subspace Learning via NMF [(Bucak and Gunsel, 2009)](http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=4298684) [website](http://www.cse.msu.edu/~bucakser/inmf_genel.html)

* * DRMF: Direct Robust Matrix Factorization [(Xiong et al. 2011)](http://ieeexplore.ieee.org/xpl/freeabs_all.jsp?arnumber=6137289) [website](http://www.autonlab.org/autonweb/20605.html)

* NTF: Non-Negative Tensor Factorization (6)
* * betaNTF: Simple beta-NTF implementation [(Antoine Liutkus, 2012)](http://www.mathworks.com/matlabcentral/fileexchange/38109-nonnegative-matrix-and-tensor-factorization--nmf--ntf--with-any-beta-divergence)

* * bcuNTD: Non-negative Tucker Decomposition by block-coordinate update (Xu and Yin, 2012) [website](http://www.math.ucla.edu/~wotaoyin/papers/bcu/ntd/index.html)

* * bcuNCP: Non-negative CP Decomposition by block-coordinate update (Xu and Yin, 2012) [website](http://www.math.ucla.edu/~wotaoyin/papers/bcu/ncp/index.html)

* * NTD-MU: Non-negative Tucker Decomposition solved by Multiplicative Updates [(Zhou et al. 2012)](http://dx.doi.org/10.1109/TSP.2012.2190410)

* * NTD-APG: Non-negative Tucker Decomposition solved by Accelerated Proximal Gradient [(Zhou et al. 2012)](http://dx.doi.org/10.1109/TSP.2012.2190410)

* * NTD-HALS: Non-negative Tucker Decomposition solved by Hierarchical ALS [(Zhou et al. 2012)](http://dx.doi.org/10.1109/TSP.2012.2190410)

* TD: Tensor Decomposition (12)
* * HoSVD: Higher-order Singular Value Decomposition (Tucker Decomposition)

* * HoRPCA-IALM: HoRPCA solved by IALM [(Goldfarb and Qin, 2013)](http://arxiv.org/abs/1311.6182) [website](https://sites.google.com/site/tonyqin/research)

* * HoRPCA-S: HoRPCA with Singleton model solved by ADAL [(Goldfarb and Qin, 2013)](http://arxiv.org/abs/1311.6182) [website](https://sites.google.com/site/tonyqin/research)

* * HoRPCA-S-NCX: HoRPCA with Singleton model solved by ADAL (non-convex) [(Goldfarb and Qin, 2013)](http://arxiv.org/abs/1311.6182) [website](https://sites.google.com/site/tonyqin/research)

* * Tucker-ADAL: Tucker Decomposition solved by ADAL [(Goldfarb and Qin, 2013)](http://arxiv.org/abs/1311.6182) [website](https://sites.google.com/site/tonyqin/research)

* * Tucker-ALS: Tucker Decomposition solved by ALS

* * CP-ALS: PARAFAC/CP decomposition solved by ALS

* * CP-APR: PARAFAC/CP decomposition solved by Alternating Poisson Regression [(Chi et al. 2011)](http://arxiv.org/abs/1112.2414)

* * CP2: PARAFAC2 decomposition solved by ALS [(Bro et al. 1999)](http://www.mathworks.com/matlabcentral/fileexchange/1089-parafac2)

* * RSTD: Rank Sparsity Tensor Decomposition [(Yin Li, 2010)](www.pami.sjtu.edu.cn/demo/RSTD.pdf) [website](http://yinli.cvpr.net/)

* * t-SVD: Tensor SVD in Fourrier Domain [(Zhang et al. 2013)](http://arxiv.org/abs/1307.0805)

* * OSTD: Online Stochastic Tensor Decomposition [(Sobral et al. 2015)](http://www.cv-foundation.org/openaccess/content_iccv_2015_workshops/w24/papers/Sobral_Online_Stochastic_Tensor_ICCV_2015_paper.pdf) [website](https://github.com/andrewssobral/ostd)

* **Some remarks**:
* * The FW-T algorithm of Mu et al. (2014) works only with [CVX library](http://cvxr.com/cvx/). Download and install it in: *lrslibrary/libs/cvx/*.

* * The DECOLOR algorithm of Zhou et al. (2011) don't works in MATLAB R2014a(x64), but works successfully in MATLAB R2013b(x64) and both R2014a(x86) and R2013b(x86).

Usage example
----------------------------
For complete details and examples, please see the **demo.m** file.
```Matlab

%% First run the setup script
lrs_setup; % or run('C:/lrslibrary/lrs_setup')

%% Load configuration
lrs_load_conf;

%% Load video
input_avi = fullfile(lrs_conf.lrs_dir,'dataset','demo.avi');
output_avi = fullfile(lrs_conf.lrs_dir,'output','output.avi');

%% Processing videos
%
% Robust PCA
process_video('RPCA', 'FPCP', input_avi, output_avi);
% Subspace Tracking
process_video('ST', 'GRASTA', input_avi, output_avi);
% Matrix Completion
process_video('MC', 'GROUSE', input_avi, output_avi);
% Low Rank Recovery
process_video('LRR', 'FastLADMAP', input_avi, output_avi);
% Three-Term Decomposition
process_video('TTD', '3WD', input_avi, output_avi);
% Non-Negative Matrix Factorization
process_video('NMF', 'ManhNMF', input_avi, output_avi);
% Non-Negative Tensor Factorization
process_video('NTF', 'bcuNCP', input_avi, output_avi);
% Tensor Decomposition
process_video('TD', 'Tucker-ALS', input_avi, output_avi);

%% Processing matrices and tensors
%
load('dataset/trafficdb/traffic_patches.mat');
V = im2double(imgdb{100});
show_3dvideo(V);

%% Matrix-based algorithms
%
[M,m,n,p] = convert_video3d_to_2d(V);
show_2dvideo(M,m,n);

% Robust PCA
out = run_algorithm('RPCA', 'FPCP', M);
% Subspace Tracking
out = run_algorithm('ST', 'GRASTA', M);
% Matrix Completion
obs = 0.5; % Percentual of observed entries [0...1]
[params.Idx, params.Omega] = subsampling(M, obs);
out = run_algorithm('MC', 'GROUSE', M, params);
% Low Rank Recovery
out = run_algorithm('LRR', 'FastLADMAP', M);
% Three-Term Decomposition
out = run_algorithm('TTD', '3WD', M);
% Non-Negative Matrix Factorization
out = run_algorithm('NMF', 'ManhNMF', M);

% Show results
show_results(M,out.L,out.S,out.O,p,m,n);

%% Tensor-based algorithms
T = tensor(V);

% Non-Negative Tensor Factorization
out = run_algorithm('NTF', 'bcuNCP', T);
% Tensor Decomposition
out = run_algorithm('TD', 'Tucker-ALS', T);

% Show results
show_3dtensors(T,out.L,out.S,out.O);
```
<p align="center"><img src="https://raw.githubusercontent.com/andrewssobral/lrslibrary/master/figs/lrs_results.png" width="650" /></p>

CPU time consumption
--------------------
The figure below shows the average CPU time consumption and the speed classification of each algorithm to decompose a *2304x51* matrix or *48x48x51* tensor data. Both matrix and tensor data were built from *dataset/demo.avi* file. The experiments were performed in a Intel Core i7-3740QM CPU 2.70GHz with 16Gb of RAM running MATLAB R2013b and Windows 7 Professional SP1 64 bits.

***A complete review evaluating the algorithms in many specific criterias will be published in a paper journal soon***

<p align="center"><img src="https://raw.githubusercontent.com/andrewssobral/lrslibrary/master/figs/algorithms_by_speed.png" /></p>

About LRSLibrary
----------------
The *LRSLibrary* has been developed by [Andrews Sobral](https://sites.google.com/site/andrewssobral) thanks especially to [Thierry Bouwmans](https://sites.google.com/site/thierrybouwmans) for his continued support and for collaborating on this important initiative. I'm grateful to all authors who have contributed in some way to the success of the LRSLibrary.

How to contribute with LRSLibrary project
-----------------------------------------
Everyone is invited to cooperate with the LRSLibrary project by sending to us any implementation of low-rank and sparse decomposition algorithms.

Option 1: email it to me (andrewssobral **at** gmail **dot** com).

Option 2: fork the library on GitHub, push your changes, then send me a pull request.

For the Option 2 you can follow the instructions below:

1) Create the algorithm folder in **lrslibrary/algorithms/XYZ/** where "**XYZ**" is one of the following methods: **RPCA** (for Robust PCA), **ST** (for Subspace Tracking), **MC** (for Matrix Completion), **LRR** (for Low Rank Recovery), **TTD** (for Three-Term Decomposition), **NMF** (for Non-Negative Matrix Factorization), **NTF** (for Non-Negative Tensor Factorization), and **TD** (for Tensor Decomposition).

2) Create a **run_alg.m** script file to execute your algorithm, all **run_alg.m** files are called automatically by the lrslibrary and each of them receives a full rank matrix/tensor named "**M**" or "**T**" for a matrix-based or a tensor-based algorithm, respectively. For the matrix **M**, each column is a vectorized version of a gray-scale frame. The number of columns in **M** is the number of frames given a video, and the number of rows represents the number of pixels in the frame. For the tensor **T** each frontal slice represents a frame given a video.

2.1) For the matrix/tensor completion case you need to define a random sub-sampling matrix/tensor named "**Idx**" that indicates which elements from the input matrix/tensor (**M** or **T**) are observed (setting **1**'s and **0**'s for observed and unobserved elements, respectively).

3) The output of your **run_alg.m** script file are a matrix/tensor named "**L**" that represents the low-rank component, and a matrix/tensor named "**S**" that represents the sparse components.

4) The last step is to add your algorithm details in the **lrs_algorithms.dat** file (open it as text file):
https://github.com/andrewssobral/lrslibrary/blob/master/lrs_algorithms.dat
For example:
**XYZ | ABB | NAME_REFERENCE | SPEED_CLASS**
where:
* **XYZ** is one of the following options: **RPCA**, **ST**, **MC**, **LRR**, **TTD**, **NMF**, **NTF**, or **TD**.
* **ABB** abbreviation name for your method, i.e: **APG** for Accelerated Proximal Gradient.
* **NAME_REFERENCE** the full name of your method and reference, i.e: **Accelerated Proximal Gradient (Lin et al. 2009)**.
* **SPEED_CLASS** represents the speed classification of your algorithm for the demo case. For this, you need to compare the CPU time consumption of your algorithm from all others following this section:
https://github.com/andrewssobral/lrslibrary#cpu-time-consumption

i.e.:
**RPCA | FPCP | Fast PCP (Rodriguez and Wohlberg, 2013) | 1**

That's all ;-)

License
-------
The LRSLibrary is free and open source for academic/research purposes (non-commercial)¹.

***¹ Some algorithms of the LRSLibrary are free for commercial purposes and others not. First you need to contact the authors of your desired algorithm and check with them the appropriate license.***

Problems or Questions
---------------------
If you have any problems or questions, please contact the author: Andrews Sobral (andrewssobral **at** gmail **dot** com)

Release Notes:
--------------
* Version 1.0.11: Fixed IALM with LMSVDS (Liu et al. 2012).

* Version 1.0.10: Added ReProCS (Narayanamurthy and Vaswani, 2017a) and MEDRoP (Narayanamurthy and Vaswani, 2017b). Thanks to Praneeth Narayanamurthy.

* Version 1.0.9: Fixed some issues for matrix completion methods.

* Version 1.0.8: Added PG-RMC (Cherapanamjeri et al. 2016) and fixed some small issues.

* Version 1.0.7: Code refactoring: *process\_matrix()*, *process\_tensor()*, *run\_algorithm\_###()* were excluded. A standard interface called *run\_algorithm* was created. For each algorithm, there is a *run\_alg.m* script for execution. Added 10 new algorithms: OSTD (Sobral et al. 2015), Nonconvex RPCA (Kang et al. 2015), MC_logdet (Kang et al. 2016), RPCA-GD (Yi et al. 2016), LMaFit (Wen et al. 2012), MC-NMF (Xu et al. 2011), ScGrassMC (Ngo and Saad, 2012), SVP (Meka et al. 2009), OR1MP (Wang et al. 2015), IALM-MC (Lin et al. 2009). OP-RPCA was moved from RPCA to MC category.

* Version 1.0.6: Added three new algorithms: STOC-RPCA: Online Robust PCA via Stochastic Optimization of Feng et al. (2013), MoG-RPCA: Mixture of Gaussians RPCA of Zhao et al. (2014), and OP-RPCA: Robust PCA via Outlier Pursuit of Xu et al. (2012).

* Version 1.0.5: Added three new method categories ST (Subspace Tracking), MC (Matrix Completion), and  TTD (Three-Term Decomposition). Added fifteen new algorithms: 3WD - 3-Way-Decomposition of Oreifej et al. (2012), MAMR and Robust MAMR - Motion-Assisted Matrix Restoration of Ye et al. (2015), ADMM - Alternating Direction Method of Multipliers of Parikh and Boyd (2014), GreGoDec - Greedy Semi-Soft GoDec Algotithm of Zhou and Tao (2013), GRASTA (Grassmannian Robust Adaptive Subspace Tracking Algorithm) of He et al. (2012), GOSUS (Grassmannian Rank-One Update Subspace Estimation) of Balzano et al. (2010), OptSpace - A Matrix Completion Algorithm of Keshavan et al. (2009), FPC - Fixed point and Bregman iterative methods for matrix rank minimization of Ma et al. (2008), SVT - A singular value thresholding algorithm for matrix completion of Cai et al. (2008), LRGeomCG - Low-rank matrix completion by Riemannian optimization of Bart Vandereycken (2013), RPCA - Robust Principal Component Analysis of De la Torre and Black (2001), GA - Grassmann Average, GM - Grassmann Median, and TGA - Trimmed Grassmann Average of Hauberg et al. (2014). Also fixed some errors.

* Version 1.0.4: Added a setup script, and configuration file. Also fixed some errors.

* Version 1.0.3: Added three new algorithms: FW-T (Mu et al. 2014), iNMF (Bucak and Gunsel, 2009) and DRMF (Xiong et al. 2011).

* Version 1.0.2: Added four new algorithms: GOSUS (Xu et al. 2013), pROST (Hage and Kleinsteuber, 2013), RegL1-ALM (Zheng et al. 2012) and ROSL (Shu et al. 2014).

* Version 1.0.1: Added RPCA-SPCP algorithms of Aravkin et al. (2014), thanks to Professor [Stephen Becker](http://amath.colorado.edu/faculty/becker/).

* Version 1.0.0: First version with 64 algorithms.
