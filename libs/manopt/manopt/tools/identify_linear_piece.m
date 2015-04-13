function [range, poly] = identify_linear_piece(x, y, window_length)
% Identify a segment of the curve (x, y) that appears to be linear.
%
% function [range poly] = identify_linear_piece(x, y, window_length)
%
% This function attempts to identify a contiguous segment of the curve
% defined by the vectors x and y that appears to be linear. A line is fit
% through the data over all windows of length window_length and the best
% fit is retained. The output specifies the range of indices such that
% x(range) is the portion over which (x, y) is the most linear and the
% output poly specifies a first order polynomial that best fits (x, y) over
% that range, following the usual matlab convention for polynomials
% (highest degree coefficients first).
%
% See also: checkdiff checkgradient checkhessian

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, July 8, 2013.
% Contributors: 
% Change log: 

    residues = zeros(length(x)-window_length, 1);
    polys = zeros(2, length(residues));
    for i = 1 : length(residues)
        range = i:i+window_length;
        [poly, meta] = polyfit(x(range), y(range), 1);
        residues(i) = meta.normr;
        polys(:, i) = poly';
    end
    [unused, best] = min(residues); %#ok<ASGLU>
    range = best:best+window_length;
    poly = polys(:, best)';

end
