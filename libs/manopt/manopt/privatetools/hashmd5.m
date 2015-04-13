function h = hashmd5(inp)
% Computes the MD5 hash of input data.
%
% function h = hashmd5(inp)
% 
% Returns a string containing the MD5 hash of the input variable. The input
% variable may be of any class that can be typecast to uint8 format, which
% is fairly non-restrictive.

% This file is part of Manopt: www.manopt.org.
% This code is a stripped version of more general hashing code by
% Michael Kleder, Nov 2005.
% Change log: 
%   Aug. 8, 2013 (NB) : Made x a static (persistent) variable, in the hope
%                       it will speed it up. Furthermore, the function is
%                       now Octave compatible.

    is_octave = exist('OCTAVE_VERSION', 'builtin');
        
    persistent x;
    if isempty(x) && ~is_octave
        x = java.security.MessageDigest.getInstance('MD5');
    end

    inp=inp(:);
    % Convert strings and logicals into uint8 format
    if ischar(inp) || islogical(inp)
        inp=uint8(inp);
    else % Convert everything else into uint8 format without loss of data
        inp=typecast(inp,'uint8');
    end
    
    % Create hash
    if ~is_octave
        x.update(inp);
        h = typecast(x.digest, 'uint8');
        h = dec2hex(h)';
        % Remote possibility: all hash bytes < 128, so pad:
        if(size(h,1))==1
            h = [repmat('0',[1 size(h,2)]);h];
        end
        h = lower(h(:)');
    else
        h = md5sum(char(inp'), true);
    end
	
end
