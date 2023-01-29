function [is_matched] = row_match(A,x)

%row_match: Function which finds rows in a matrix that match a given
%vector.
%
% © 2005-2013 Michael Cardiff, All Rights Reserved.
%
%[is_matched] = row_match(A,x)
%where:
%   -is_matched is a vector where is_matched(i) = 1 if and only if row i of
%   A is exactly the same as the vector x.
%   -A is the input matrix being checked
%   -x is the vector being matched.

[m,n] = size(A);
[o,p] = size(x);

if o > p
    x = x';
    [o,p] = deal(p,o);
end

if o ~= 1
    error('x must be a vector');
end

if n ~= p
    error('Number of columns in A must be equal to the length of x');
end

xmat = repmat(x,m,1);

is_matched = prod(1.*(xmat == A),2);
