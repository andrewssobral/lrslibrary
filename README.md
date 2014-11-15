Last update: **16/11/2014**

Library Version: **1.0.1**

LRSLibrary
----------
*Low-Rank and Sparse* tools for Background Modeling and Subtraction in Videos.

The *LRSLibrary* provides a collection of **low-rank and sparse decomposition** algorithms in MATLAB. The library was designed for motion segmentation in videos, but it can be also used or adapted for other computer vision problems, please see this [page](http://perception.csl.illinois.edu/matrix-rank/applications.html). Currently the LRSLibrary contains a total of **68** *matrix-based* and *tensor-based* algorithms. The LRSLibrary was tested successfully in MATLAB R2013b both x86 and x64 versions.

<p align="center"><img src="https://sites.google.com/site/andrewssobral/lrs_results2.png" /></p>

```
See also:

Presentation about Matrix and Tensor Tools for Computer Vision 
http://www.slideshare.net/andrewssobral/matrix-and-tensor-tools-for-computer-vision

MTT: Matlab Tensor Tools for Computer Vision
https://github.com/andrewssobral/mtt

IMTSL: Incremental and Multi-feature Tensor Subspace Learning
https://github.com/andrewssobral/imtsl
```

Citation
---------
A paper about the LRSLibrary will be published in a journal, but if you use this code for your publications, please cite it currently as:
```
@inproceedings{asobral2014,
    author       = "Sobral, A. and Baker, C. G. and Bouwmans, T. and Zahzah, E.",
    title        = "Incremental and Multi-feature Tensor Subspace Learning applied for Background Modeling and Subtraction",
    booktitle    = "International Conference on Image Analysis and Recognition (ICIAR'14)",
    year         = "2014",
    month        = "October",
    publisher    = "Lecture Notes in Computer Science (Springer LNCS)",
    url          = "https://github.com/andrewssobral/imtsl"
}
@article{bouwmans2014,
  author         = "Thierry Bouwmans and El Hadi Zahzah",
  title          = "Robust \{PCA\} via Principal Component Pursuit: A review for a comparative evaluation in video surveillance",
  journal        = "Computer Vision and Image Understanding (CVIU)",
  volume         = "122",
  pages          = "22--34",
  year           = "2014"
}
```

GUI
---
The *LRSLibrary* provides an easy-to-use graphical user interface (GUI) for background modeling and subtraction in videos. Just execute **run.m** and enjoy it ;)

<p align="center"><img src="https://sites.google.com/site/andrewssobral/lrslibrary_gui2.png" /></p>

Each algorithm is classified by its cpu time consumption with the following icons:
<p align="center"><img src="https://sites.google.com/site/andrewssobral/time_legend.png?width=300" /></p>


List of the algorithms available in LRSLibrary
----------------------------------------------
* RPCA: Robust PCA (34)
* * PCP: Principal Component Pursuit [(Candes et al. 2009)](http://arxiv.org/abs/0912.3599) 

* * FPCP: Fast PCP [(Rodriguez and Wohlberg 2013)](http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=6738015) 

* * AS-RPCA: Active Subspace: Towards Scalable Low-Rank Learning [(Liu and Yan, 2012)](http://dl.acm.org/citation.cfm?id=2421487)

* * R2PCP: Riemannian Robust Principal Component Pursuit [(Hintermüller and Wu, 2014)](http://link.springer.com/article/10.1007/s10851-014-0527-y) 

* * ALM: Augmented Lagrange Multiplier [(Tang and Nehorai 2011)](http://dx.doi.org/10.1109/CISS.2011.5766144) 

* * EALM: Exact ALM [(Lin et al. 2009)](http://arxiv.org/abs/1009.5055) [website](http://perception.csl.illinois.edu/matrix-rank/sample_code.html)

* * IALM: Inexact ALM [(Lin et al. 2009)](http://arxiv.org/abs/1009.5055)  [website](http://perception.csl.illinois.edu/matrix-rank/sample_code.html)

* * IALM_LMSVDS: IALM with LMSVDS [(Liu et al. 2012)](http://epubs.siam.org/doi/abs/10.1137/120871328) 

* * IALM_BLWS: IALM with BLWS [(Lin and Wei 2010)](http://arxiv.org/abs/1012.0365) 

* * APG_PARTIAL: Partial Accelerated Proximal Gradient [(Lin et al. 2009)](http://arxiv.org/abs/1009.5055)  [website](http://perception.csl.illinois.edu/matrix-rank/sample_code.html)

* * APG: Accelerated Proximal Gradient [(Lin et al. 2009)](http://arxiv.org/abs/1009.5055)  [website](http://perception.csl.illinois.edu/matrix-rank/sample_code.html)

* * DUAL: Dual RPCA [(Lin et al. 2009)](http://arxiv.org/abs/1009.5055)  [website](http://perception.csl.illinois.edu/matrix-rank/sample_code.html)

* * SVT: Singular Value Thresholding [(Cai et al. 2008)](http://arxiv.org/abs/0810.3286)  [website](http://perception.csl.illinois.edu/matrix-rank/sample_code.html)

* * ADM: Alternating Direction Method [(Yuan and Yang 2009)](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.400.8797) 

* * LSADM: LSADM [(Goldfarb et al. 2010)](http://arxiv.org/abs/0912.4571) 

* * L1F: L1 Filtering [(Liu et al. 2011)](http://arxiv.org/abs/1108.5359) 

* * DECOLOR: Contiguous Outliers in the Low-Rank Representation [(Zhou et al. 2011)](http://arxiv.org/abs/1109.0882) 

* * NSA1: Non-Smooth Augmented Lagrangian v1 [(Aybat et al. 2011)](http://arxiv.org/abs/1105.2126) 

* * NSA2: Non-Smooth Augmented Lagrangian v2 [(Aybat et al. 2011)](http://arxiv.org/abs/1105.2126) 

* * PSPG: Partially Smooth Proximal Gradient [(Aybat et al. 2012)](http://arxiv.org/abs/1309.6976) 

* * flip-SPCP-sum-SPG: Flip-Flop version of Stable PCP-sum solved by Spectral Projected Gradient [(Aravkin et al. 2014)](https://github.com/stephenbeckr/fastRPCA)

* * flip-SPCP-max-QN: Flip-Flop version of Stable PCP-max solved by Quasi-Newton [(Aravkin et al. 2014)](https://github.com/stephenbeckr/fastRPCA)

* * Lag-SPCP-SPG: Lagrangian SPCP solved by Spectral Projected Gradient [(Aravkin et al. 2014)](https://github.com/stephenbeckr/fastRPCA)

* * Lag-SPCP-QN: Lagrangian SPCP solved by Quasi-Newton [(Aravkin et al. 2014)](https://github.com/stephenbeckr/fastRPCA)

* * BRPCA-MD: Bayesian Robust PCA with Markov Dependency [(Ding et al. 2011)](http://people.ee.duke.edu/~lcarin/LRS_09.pdf) [website](http://people.ee.duke.edu/~lcarin/BCS.html)

* * BRPCA-MD-NSS: BRPCA-MD with Non-Stationary Noise [(Ding et al. 2011)](http://people.ee.duke.edu/~lcarin/LRS_09.pdf) [website](http://people.ee.duke.edu/~lcarin/BCS.html) 

* * VBRPCA: Variational Bayesian RPCA [(Babacan et al. 2011)](http://arxiv.org/abs/1102.5288) 

* * PRMF: Probabilistic Robust Matrix Factorization [(Wang et al. 2012)](http://winsty.net/papers/prmf.pdf) [website](http://winsty.net/prmf.html)

* * OPRMF: Online PRMF [(Wang et al. 2012)](http://winsty.net/papers/prmf.pdf) [website](http://winsty.net/prmf.html)

* * MBRMF: Markov BRMF [(Wang and Yeung 2013)](http://winsty.net/papers/brmf.pdf) [website](http://winsty.net/brmf.html) 

* * TFOCS-EC: TFOCS with equality constraints [(Becker et al. 2011)](https://github.com/cvxr/TFOCS/raw/gh-pages/TFOCS.pdf) [website](http://cvxr.com/tfocs/demos/rpca/)

* * TFOCS-IC: TFOCS with inequality constraints [(Becker et al. 2011)](https://github.com/cvxr/TFOCS/raw/gh-pages/TFOCS.pdf) [website](http://cvxr.com/tfocs/demos/rpca/)

* * GoDec: Go Decomposition [(Zhou and Tao 2011)](http://www.icml-2011.org/papers/41_icmlpaper.pdf) [website](https://sites.google.com/site/godecomposition/home)

* * SSGoDec: Semi-Soft GoDec [(Zhou and Tao 2011)](http://www.icml-2011.org/papers/41_icmlpaper.pdf) [website](https://sites.google.com/site/godecomposition/home)

* LRR: Low Rank Recovery (5)
* * EALM: Exact ALM [(Lin et al. 2009)](http://arxiv.org/abs/1009.5055) 

* * IALM: Inexact ALM [(Lin et al. 2009)](http://arxiv.org/abs/1009.5055) 

* * ADM: Alternating Direction Method [(Lin et al. 2011)](http://arxiv.org/abs/1109.0367)

* * LADMAP: Linearized ADM with Adaptive Penalty [(Lin et al. 2011)](http://arxiv.org/abs/1109.0367)

* * FastLADMAP: Fast LADMAP [(Lin et al. 2011)](http://arxiv.org/abs/1109.0367) 

* NMF: Non-Negative Matrix Factorization (12)
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

* NTF: Non-Negative Tensor Factorization (6)
* * betaNTF: Simple beta-NTF implementation [(Antoine Liutkus, 2012)](http://www.mathworks.com/matlabcentral/fileexchange/38109-nonnegative-matrix-and-tensor-factorization--nmf--ntf--with-any-beta-divergence) 

* * bcuNTD: Non-negative Tucker Decomposition by block-coordinate update (Xu and Yin, 2012) [website](http://www.math.ucla.edu/~wotaoyin/papers/bcu/ntd/index.html)

* * bcuNCP: Non-negative CP Decomposition by block-coordinate update (Xu and Yin, 2012) [website](http://www.math.ucla.edu/~wotaoyin/papers/bcu/ncp/index.html)

* * NTD-MU: Non-negative Tucker Decomposition solved by Multiplicative Updates [(Zhou et al. 2012)](http://dx.doi.org/10.1109/TSP.2012.2190410) 

* * NTD-APG: Non-negative Tucker Decomposition solved by Accelerated Proximal Gradient [(Zhou et al. 2012)](http://dx.doi.org/10.1109/TSP.2012.2190410) 

* * NTD-HALS: Non-negative Tucker Decomposition solved by Hierarchical ALS [(Zhou et al. 2012)](http://dx.doi.org/10.1109/TSP.2012.2190410)  

* TD: Tensor Decomposition (11)
* * HoSVD: Higher-order Singular Value Decomposition (Tucker Decomposition) 

* * HoRPCA-IALM: HoRPCA solved by IALM [(Goldfarb and Qin, 2013)](http://arxiv.org/abs/1311.6182) [website](https://sites.google.com/site/tonyqin/research)

* * HoRPCA-S: HoRPCA with Singleton model solved by ADAL [(Goldfarb and Qin, 2013)](http://arxiv.org/abs/1311.6182) [website](https://sites.google.com/site/tonyqin/research)

* * HoRPCA-S-NCX: HoRPCA with Singleton model solved by ADAL (non-convex) [(Goldfarb and Qin, 2013)](http://arxiv.org/abs/1311.6182) [website](https://sites.google.com/site/tonyqin/research)

* * Tucker-ADAL: Tucker Decomposition solved by ADAL [(Goldfarb and Qin, 2013)](http://arxiv.org/abs/1311.6182) [website](https://sites.google.com/site/tonyqin/research)

* * Tucker-ALS: Tucker Decomposition solved by ALS 

* * CP-ALS: PARAFAC/CP decomposition solved by ALS 

* * CP-APR: PARAFAC/CP decomposition solved by Alternating Poisson Regression [(Chi et al. 2011)](http://arxiv.org/abs/1112.2414) 

* * CP2: PARAFAC2 decomposition solved by ALS [(Bro et al. 1999)](http://www.mathworks.com/matlabcentral/fileexchange/1089-parafac2) 

* * RSTD: Rank Sparsity Tensor Decomposition [(Yin Li 2010)](www.pami.sjtu.edu.cn/demo/RSTD.pdf) [website](http://yinli.cvpr.net/) 

* * t-SVD: Tensor SVD in Fourrier Domain [(Zhang et al. 2013)](http://arxiv.org/abs/1307.0805) 

* **Some remarks**:
* * The ADM algorithm of Yuan and Yang (2009) works only on win32 platforms (mexsvd.mexw32).

* * The DECOLOR algorithm of Zhou et al. (2011) don't works in MATLAB R2014a(x64), but works successfully in MATLAB R2013b(x64) and both R2014a(x86) and R2013b(x86).

Usage example
----------------------------
For complete details and examples, please see the **demo.m** file.
```Matlab

%% Processing videos
%
% Robust PCA
process_video('RPCA', 'FPCP', 'dataset/demo.avi', 'output/demo_FPCP.avi');
% Low Rank Recovery
process_video('LRR', 'FastLADMAP', 'dataset/demo.avi', 'output/demo_LRR-FastLADMAP.avi');
% Non-Negative Matrix Factorization
process_video('NMF', 'ManhNMF', 'dataset/demo.avi', 'output/demo_ManhNMF.avi');
% Non-Negative Tensor Factorization
process_video('NTF', 'bcuNCP', 'dataset/demo.avi', 'output/demo_bcuNCP.avi');
% Tensor Decomposition
process_video('TD', 'Tucker-ALS', 'dataset/demo.avi', 'output/demo_Tucker-ALS.avi');

%% Processing matrices and tensors
%
load('dataset/trafficdb/traffic_patches.mat');
V = im2double(imgdb{100});
show_3dvideo(V);

%% Matrix-based algorithms
%
[M, m, n, p] = convert_video3d_to_2d(V);
show_2dvideo(M,m,n);
% Robust PCA
results = process_matrix('RPCA', 'FPCP', M, []);
show_results(M,results.L,results.S,results.O,p,m,n);
% Low Rank Recovery
results = process_matrix('LRR', 'FastLADMAP', M, []);
show_results(M,results.L,results.S,results.O,p,m,n);
% Non-Negative Matrix Factorization
results = process_matrix('NMF', 'ManhNMF', M, []);
show_results(M,results.L,results.S,results.O,p,m,n);

%% Tensor-based algorithms
%
add_tensor_libs;
T = tensor(V);
% Non-Negative Tensor Factorization
results = process_tensor('NTF', 'bcuNCP', T);
show_3dtensors(T,results.L,results.S,results.O);
% Tensor Decomposition
results = process_tensor('TD', 'Tucker-ALS', T);
show_3dtensors(T,results.L,results.S,results.O);
rem_tensor_libs;
```
<p align="center"><img src="https://sites.google.com/site/andrewssobral/lrs_results.png?width=650" /></p>

CPU time consumption
--------------------
The figure below shows the average CPU time consumption and the speed classification of each algorithm to decompose a *2304x51* matrix or *48x48x51* tensor data. Both matrix and tensor data were built from *dataset/demo.avi* file. The experiments were performed in a Intel Core i7-3740QM CPU 2.70GHz with 16Gb of RAM running MATLAB R2013b and Windows 7 Professional SP1 64 bits.
<p align="center"><img src="https://sites.google.com/site/andrewssobral/algorithms_by_speed.png" /></p>

About LRSLibrary
----------------
The *LRSLibrary* has been developed by [Andrews Sobral](https://sites.google.com/site/andrewssobral) thanks especially to [Thierry Bouwmans](https://sites.google.com/site/thierrybouwmans) for his continued support and for collaborating on this important initiative. I'm grateful to all authors who have contributed in some way to the success of the LRSLibrary.

How to contribute with LRSLibrary project
-----------------------------------------
Everyone is invited to cooperate with the LRSLibrary project by sending to us any implementation of low-rank and sparse decomposition algorithms.

License
-------
The source code is available only for academic/research purposes (non-commercial).

Problems or Questions
---------------------
If you have any problems or questions, please contact the author: Andrews Sobral (andrewssobral@gmail.com)

Release Notes:
--------------
* Version 1.0.1: Added RPCA-SPCP algorithms from Aravkin et al. (2014), thanks to Professor [Stephen Becker](http://amath.colorado.edu/faculty/becker/)

* Version 1.0.0: First version with 64 algorithms.
