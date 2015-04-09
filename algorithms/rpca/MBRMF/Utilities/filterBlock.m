function ret = filterBlock(in)
    in = in + 1e-8;
%     idx = [2 4 5 6 8];
%     in = in(idx, :);
    mid = in(5, :);
%     ret = sum(abs(bsxfun(@minus, log(mid), log(in)))) / 8;
    temp = abs(bsxfun(@minus, log(mid), log(in)));
    ret = sum((temp)) / 8;
end