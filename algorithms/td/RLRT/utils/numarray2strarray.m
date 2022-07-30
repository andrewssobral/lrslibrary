function strcell = numarray2strarray( vec )

N = length(vec);
strcell = cell( 1, N );
for i = 1:N
    strcell{i} = num2str( vec(i) );
end

end