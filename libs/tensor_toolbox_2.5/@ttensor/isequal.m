function [tf, tf_core, tf_U] = isequal(A,B)
%ISEQUAL True if each component of two ttensor's is numerically equal.

tf = false;
tf_core = false;
tf_U = false;

if ~isa(B,'ttensor')
    return;
end

if ndims(A) ~= ndims(B)
    return;
end

tf_core = isequal(A.core, B.core);
tf_U = cellfun(@isequal, A.u, B.u);
tf = tf_core & all(tf_U);

