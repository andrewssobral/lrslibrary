function [r]=RMSE(err)
% r = RMSE(err)
% author: Liang Xiong (lxiong@cs.cmu.edu)

r = sqrt(mean(err(:).^2));
