function bytesize(in, fid)
% BYTESIZE writes the memory usage of the provide variable to the given file
% identifier. Output is written to screen if fid is 1, empty or not provided.

if nargin == 1 || isempty(fid)
    fid = 1;
end

s = whos('in');
fprintf(fid,[Bytes2str(s.bytes) '\n']);
end

function str = Bytes2str(NumBytes)
% BYTES2STR Private function to take integer bytes and convert it to
% scale-appropriate size.

scale = floor(log(NumBytes)/log(1024));
switch scale
    case 0
        str = [sprintf('%.0f',NumBytes) ' b'];
    case 1
        str = [sprintf('%.2f',NumBytes/(1024)) ' kb'];
    case 2
        str = [sprintf('%.2f',NumBytes/(1024^2)) ' Mb'];
    case 3
        str = [sprintf('%.2f',NumBytes/(1024^3)) ' Gb'];
    case 4
        str = [sprintf('%.2f',NumBytes/(1024^4)) ' Tb'];
    case -inf
        % Size occasionally returned as zero (eg some Java objects).
        str = 'Not Available';
    otherwise
       str = 'Over a petabyte!!!';
end
end