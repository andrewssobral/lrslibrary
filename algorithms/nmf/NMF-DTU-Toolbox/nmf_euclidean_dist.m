function err = nmf_eucl_dist(X,Y)

err = sum(sum((X-Y).^2));