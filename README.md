Last update: **08/11/2014**

Library Version: **1.0.0**

LRSLibrary
----------
*Low-Rank and Sparse* tools for Background Modeling and Subtraction in Videos.

The *LRSLibrary* provides a collection of **low-rank and sparse decomposition** algorithms in MATLAB for video motion segmentation. Currently the LRSLibrary contains a total of **64** *matrix-based* and *tensor-based* algorithms. The LRSLibrary was tested successfully in MATLAB R2013b both x86 and x64 versions.

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
* RPCA: Robust PCA (30)
* * PCP: Principal Component Pursuit (Candes et al. 2009) 

* * FPCP: Fast PCP (Rodriguez and Wohlberg 2013) 

* * AS-RPCA: Active Subspace:  Towards Scalable Low-Rank Learning (Liu and Yan, 2012)

* * R2PCP: Riemannian Robust Principal Component Pursuit (Hintermüller and Wu, 2014) 

* * ALM: Augmented Lagrange Multiplier (Tang and Nehorai 2011) 

* * EALM: Exact ALM (Lin et al. 2009) 

* * IALM: Inexact ALM (Lin et al. 2009) 

* * IALM_LMSVDS: IALM with LMSVDS (Liu et al. 2012) 

* * IALM_BLWS: IALM with BLWS (Lin and Wei 2010) 

* * APG_PARTIAL: Partial Accelerated Proximal Gradient (Lin et al. 2009) 

* * APG: Accelerated Proximal Gradient (Lin et al. 2009) 

* * DUAL: Dual RPCA (Lin et al. 2009) 

* * SVT: Singular Value Thresholding (Cai et al. 2008) 

* * ADM: Alternating Direction Method (Yuan and Yang 2009) 

* * LSADM: LSADM (Goldfarb et al. 2010) 

* * L1F: L1 Filtering (Liu et al. 2011) 

* * DECOLOR: Contiguous Outliers in the Low-Rank Representation (Zhou et al. 2011) 

* * NSA1: Non-Smooth Augmented Lagrangian v1 (Aybat et al. 2011) 

* * NSA2: Non-Smooth Augmented Lagrangian v2 (Aybat et al. 2011) 

* * PSPG: Partially Smooth Proximal Gradient (Aybat et al. 2012) 

* * BRPCA-MD: Bayesian Robust PCA with Markov Dependency (Ding et al. 2011) 

* * BRPCA-MD-NSS: BRPCA-MD with Non-Stationary Noise (Ding et al. 2011) 

* * VBRPCA: Variational Bayesian RPCA (Babacan et al. 2011) 

* * PRMF: Probabilistic Robust Matrix Factorization (Wang et al. 2012) 

* * OPRMF: Online PRMF (Wang et al. 2012) 

* * MBRMF: Markov BRMF (Wang and Yeung 2013) 

* * TFOCS-EC: TFOCS with equality constraints (Becker et al. 2011) 

* * TFOCS-IC: TFOCS with inequality constraints (Becker et al. 2011) 

* * GoDec: Go Decomposition (Zhou and Tao 2011) 

* * SSGoDec: Semi-Soft GoDec (Zhou and Tao 2011) 

* LRR: Low Rank Recovery (5)
* * EALM: Exact ALM (Lin et al. 2009) 

* * IALM: Inexact ALM (Lin et al. 2009) 

* * ADM: Alternating Direction Method (Lin et al. 2011) 

* * LADMAP: Linearized ADM with Adaptive Penalty (Lin et al. 2011) 

* * FastLADMAP: Fast LADMAP (Lin et al. 2011) 

* NMF: Non-Negative Matrix Factorization (12)
* * NMF-MU: NMF solved by Multiplicative Updates 

* * NMF-PG: NMF solved by Projected Gradient 

* * NMF-ALS: NMF solved by Alternating Least Squares 

* * NMF-ALS-OBS: NMF solved by Alternating Least Squares with Optimal Brain Surgeon 

* * PNMF: Probabilistic Non-negative Matrix Factorization 

* * ManhNMF: Manhattan NMF (Guan et al. 2013) 

* * NeNMF: NMF via Nesterovs Optimal Gradient Method (Guan et al. 2012) 

* * LNMF: Spatially Localized NMF (Li et al. 2001) 

* * ENMF: Exact NMF (Gillis and Glineur, 2012) 

* * nmfLS2: Non-negative Matrix Factorization with sparse matrix (Ji and Eisenstein, 2013) 

* * Semi-NMF: Semi Non-negative Matrix Factorization (Trigeorgis et al. 2014) 

* * Deep-Semi-NMF: Deep Semi Non-negative Matrix Factorization (Trigeorgis et al. 2014) 

* NTF: Non-Negative Tensor Factorization (6)
* * betaNTF: Simple beta-NTF implementation (Antoine Liutkus, 2012) 

* * bcuNTD: Non-negative Tucker Decomposition by block-coordinate update (Xu and Yin, 2012) 

* * bcuNCP: Non-negative CP Decomposition by block-coordinate update (Xu and Yin, 2012) 

* * NTD-MU: Non-negative Tucker Decomposition solved by Multiplicative Updates (Zhou et al. 2012) 

* * NTD-APG: Non-negative Tucker Decomposition solved by Accelerated Proximal Gradient (Zhou et al. 2012) 

* * NTD-HALS: Non-negative Tucker Decomposition solved by Hierarchical ALS (Zhou et al. 2012) 

* TD: Tensor Decomposition (11)
* * HoSVD: Higher-order Singular Value Decomposition (Tucker Decomposition) 

* * HoRPCA-IALM: HoRPCA solved by IALM (Goldfarb and Qin, 2013) 

* * HoRPCA-S: HoRPCA with Singleton model solved by ADAL (Goldfarb and Qin, 2013) 

* * HoRPCA-S-NCX: HoRPCA with Singleton model solved by ADAL (non-convex) (Goldfarb and Qin, 2013) 

* * Tucker-ADAL: Tucker Decomposition solved by ADAL (Goldfarb and Qin, 2013) 

* * Tucker-ALS: Tucker Decomposition solved by ALS 

* * CP-ALS: PARAFAC/CP decomposition solved by ALS 

* * CP-APR: PARAFAC/CP decomposition solved by Alternating Poisson Regression (Chi et al. 2011) 

* * CP2: PARAFAC2 decomposition solved by ALS (Bro et al. 1999) 

* * RSTD: Rank Sparsity Tensor Decomposition (Yin Li 2010) 

* * t-SVD: Tensor SVD in Fourrier Domain (Zhang et al. 2013) 

* **Some remarks**:
* * The ADM algorithm of Yuan and Yang (2009) works only on win32 platforms (mexsvd.mexw32).

* * The DECOLOR algorithm of Zhou et al. (2011) don't works in MATLAB R2014a(x64), but works successfully in MATLAB R2013b(x64) and both R2014a(x86) and R2013b(x86).

Usage example
----------------------------
For complete details and examples, please see the **demo.m** file.
```Matlab
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
```

CPU time consumption
--------------------
The figure below shows the average CPU time consumption and the speed classification of each algorithm to decompose a *2304x51* matrix or *48x48x51* tensor data. Both matrix and tensor data were built from *dataset/demo.avi* file. The experiments were performed in a Intel Core i7-3740QM CPU 2.70GHz with 16Gb of RAM running MATLAB R2013b and Windows 7 Professional SP1 64 bits.
<p align="center"><img src="https://sites.google.com/site/andrewssobral/algorithms_by_speed.png" /></p>

About LRSLibrary
----------------
The *LRSLibrary* has been developed by [Andrews Sobral](https://sites.google.com/site/andrewssobral) thanks especially to [Thierry Bouwmans](https://sites.google.com/site/thierrybouwmans). We are grateful to all authors who have contributed in some way to the success of the LRSLibrary.

How to contribute with LRSLibrary project
-----------------------------------------
Everyone is invited to cooperate with the LRSLibrary project by sending any implementation of low-rank and sparse decomposition algorithms.

License
-------
The source code is available only for academic/research purposes (non-commercial).

Problems or Questions
---------------------
If you have any problems or questions, please contact the author: Andrews Sobral (andrewssobral@gmail.com)

Release Notes:
--------------
* Version 1.0.0: First version with 64 algorithms.
