function [s_hat] = lin_tikh_inv(y,C_dd_inv_pf,L,H)

%lin_tikh_inv.m: Function used to compute the linearized solution to a
%roughened tikhonov inverse problem
%   [s_hat] = lin_invert_solve(y,C_dd_inv_pf,L,H)
%   where:
%       -y is the data. If solving a nonlinear case, y should be passed as
%       y-h(s_tilde) + H*s_tilde
%       -C_dd_inv_pf is the square-root of the inverse of the data
%       covariance matrix (i.e., 1/ the standard deviations)
%       -L is the Tikhonov roughening matrix
%       -H is the forward model (or sensitivity matrix, in the nonlinear
%       case)
%

num_reqout = 1;

m = size(y,1);
n = size(H,2);
r = size(L,1);


b = zeros((m+r),1);

A_fun = @(s,transparg) tikh_mult_func(s,transparg,H,C_dd_inv_pf,L);

s_hat = lsqr(A_fun,b);

end