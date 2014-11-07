function D = EuDist2(fea_a,fea_b,bSqrt)
% Euclidean Distance matrix
%   D = EuDist(fea_a,fea_b)
%   fea_a:    nSample_a * nFeature
%   fea_b:    nSample_b * nFeature
%   D:      nSample_a * nSample_a
%       or  nSample_a * nSample_b


if ~exist('bSqrt','var')
    bSqrt = 1;
end


if (~exist('fea_b','var')) | isempty(fea_b)
    [nSmp, nFea] = size(fea_a);

    aa = sum(fea_a.*fea_a,2);
    ab = fea_a*fea_a';
    
    aa = full(aa);
    ab = full(ab);

    if bSqrt
        D = sqrt(repmat(aa, 1, nSmp) + repmat(aa', nSmp, 1) - 2*ab);
        D = real(D);
    else
        D = repmat(aa, 1, nSmp) + repmat(aa', nSmp, 1) - 2*ab;
    end
    
    D = max(D,D');
    D = D - diag(diag(D));
    D = abs(D);
else
    [nSmp_a, nFea] = size(fea_a);
    [nSmp_b, nFea] = size(fea_b);
    
    aa = sum(fea_a.*fea_a,2);
    bb = sum(fea_b.*fea_b,2);
    ab = fea_a*fea_b';

    aa = full(aa);
    bb = full(bb);
    ab = full(ab);

    if bSqrt
        D = sqrt(repmat(aa, 1, nSmp_b) + repmat(bb', nSmp_a, 1) - 2*ab);
        D = real(D);
    else
        D = repmat(aa, 1, nSmp_b) + repmat(bb', nSmp_a, 1) - 2*ab;
    end
    
    D = abs(D);
end

