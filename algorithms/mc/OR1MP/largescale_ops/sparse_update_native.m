function S = sparse_update_native ( S, v )
idx = find(S~=0);

for ii = 1: length(idx)
    S(idx(ii)) = v(ii);
end

end