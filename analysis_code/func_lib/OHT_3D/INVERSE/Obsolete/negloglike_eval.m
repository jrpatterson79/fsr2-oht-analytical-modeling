function [NLAP] = NLAP_eval(y,X,s,beta,Q,R,h_function)

%negloglike_eval:Function which calculates the negative log-likelihood
%using a general model (h_function), and given either a Q matrix or a
%function that calculates a Q matrix-vector product
%
% [negloglike] = negloglike_eval(y,X,s,beta,Q,R,h_function)
% where:
%   OUTPUTS:
%       -negloglike is the neg

%TODO - Make arguments to minres optional
%TODO - Add options for methods for performing the inversions step of Q

resid = y - h_function(s);

NLAP_resid_part = 0.5.*resid'*(R\resid);
xi = s-X*beta;

if isa(Q,'function_handle')
    NLAP_geostat_part = 0.5.*xi'*minres(Q,xi,1e-4,10000);
else
    NLAP_geostat_part = 0.5.*xi'*(Q\xi);
end
    
NLAP = NLAP_resid_part + NLAP_geostat_part;
