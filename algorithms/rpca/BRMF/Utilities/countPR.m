function [p r] = countPR(truth, est)
    correct = sum(sum(double(truth) + double(est) == 2));
    p = correct / sum(sum(truth));
    r = correct / sum(sum(est));
end