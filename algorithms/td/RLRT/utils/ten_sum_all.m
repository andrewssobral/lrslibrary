function result = ten_sum_all( Xs )
% sums up a cell array of tensors

result = tenzeros( size(Xs{1}) );
for i = 1:length(Xs)
    result = result + Xs{i};
end

end