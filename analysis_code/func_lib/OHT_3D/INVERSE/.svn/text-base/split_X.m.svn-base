function [X_split] = split_X(X,s,col_split,cutoff)

%split_X: Simple function to split an assignment X-matrix (n x p) given a
%threshhold value within that zone. Useful for interactive zonation of
%geostatistical solutions. Can be applied multiple times in order to refine
%solutions.
%   X_split = split_X(X,s,col_split,cutoff)
%   where:
%       X is the current assignment matrix (n x p)
%       s is the current best parameter estimate (n x 1)
%       col_split is the zone (column) of X which should be thresholded
%       (scalar)
%       cutoff is the threshold value (scalar)

p = size(X,2);

min_incol = min(s(find(X(:,col_split)==1)));
max_incol = max(s(find(X(:,col_split)==1)));

if cutoff < min_incol || cutoff > max_incol
    X_split = X;
else
    X_split = zeros(size(X,1),p+1);
    split_found = 0;
    for i=1:p
        if i ~= col_split
            X_split(:,i+split_found) = X(:,i);
        else
            X_split(:,i:i+1) = [X(:,i).*(s<cutoff) X(:,i).*(s>=cutoff)];
            split_found = 1;
        end
    end
end