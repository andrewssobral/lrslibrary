% Example of using ExactNMF on the slack matrix of the regular 9-gon with
% 1000 attemps

load SlackMatrices;  
[V,W]=ExactNMF(S9gon, 7, 1000);