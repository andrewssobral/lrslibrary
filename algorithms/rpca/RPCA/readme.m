 This Software computes the Robust Principal Component Analysis or Robust Singular Value Decomposition
 explained in the references:

 REFERENCES:
 De la Torre, F. and Black, M. J., Robust principal component analysis for computer vision. Int. Conf. on Computer Vision 2001, Vancouver. 
 De la Torre, F. and Black, M. J., A Framework for Robust Subspace Learning. In preparation.
 
This is EXPERIMENTAL software, so with high probability is possible to have bugs. If you find any bug, please
send an e-mail to  ftorre@salleURL.edu.
Use at your own risk. No warranty is implied by this distribution. 

There are three Matlab(it has been tested in version 5.3) programs:

rpca: An example program of how to run RPCA
rob_pca: Computes Robust PCA
weighted_pca: Computes Weighted PCA

data_2_days is a file containing the data gathered during 2 days with a static camera.
Beware! you need at least a 800M RAM memory to run the algorithm in Matlab with this data!
It takes several hours to converge!

If you do not have this memory space, just compute the algorithm with half of the data
(e.g. add this line after loading the data in rpca.m Data=Data(:,1:end/2);)

Link to download original source code + data: https://gofile.io/?c=QIsyrz
