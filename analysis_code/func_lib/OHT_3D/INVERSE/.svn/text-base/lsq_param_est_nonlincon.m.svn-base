function [s_best, s_cov] = lsq_param_est_nonlincon(d,h_func,R,s_init,varargin)

%lsq_param_est_nonlin: Function which performs least squares parameter estimation
%for over-determined inverse problems. Can estimate best values of the
%parameters either with or without prior distribution information.
%[s_best, s_cov] =
%lsq_param_est_nonlin(d,h_func,R,s_init,{s_step},{s_prior},{s_covar})
%where:
%   -d is the data (n x 1)
%   -h_func is the forward mapping from parameters to data (R^m --> R^n),
%   supplied as an anonymous function with a single (m x 1) vector input
%   and a single (n x 1) vector output
%   -R is the data error covariance matrix (n x n)
%   -s_init is an initial guess at the parameter values (m x 1)
%   -{s_step} (optional) is a stepsize to be used for the parameters during
%   computation of the approximate Jacobian with finite difference methods.
%   s_step(i) is the step size that will be used when evaluating the ith
%   column of the Jacobian. Default is 1e-3 for all variables.
%   -{s_prior} (optional) is a prior estimate of parameter values (m x 1)
%   -{s_covar} (optional) is a prior covariance matrix (m x m) of the
%   parameters. If s_prior is supplied but s_covar is empty or not
%   supplied, then s_covar is assumed to be eye(m,m)
%   -{A},{b},{Aeq},{beq},{lb},{ub},{nonlcon} are constraints supplied on
%   the parameter values. See fmincon for more info on method for supplying
%   these. NOTE: Even if constraints are used, h_func should still return
%   values if parameter values are supplied that do not meet these
%   constraints.

num_params = size(s_init,1);

s_step = 1e-3*ones(num_params,1);
s_prior = [];
s_covar = [];
A = []; b = []; Aeq = []; beq = []; lb = []; ub = []; nonlcon = [];

nreqin = 4;
if nargin > nreqin
    [s_step, s_prior, s_covar, A, b, Aeq, beq, lb, ub, nonlcon] = ...
        process_extra_args(varargin,s_step, s_prior, s_covar, A, b, ...
        Aeq, beq, lb, ub, nonlcon);
end

R_inv = inv(R);
R_npf = chol(R_inv);

if ~isempty(s_prior)
    Q_inv = inv(s_covar);
    neglogAP_func = @(pvec) .5*sum(((d-h_func(pvec))'*R_npf).^2) ...
        + .5*(pvec-s_prior)'*Q_inv*(pvec-s_prior);
else
    neglogAP_func = @(pvec) .5*sum(((d-h_func(pvec))'*R_npf).^2);
end

paramopt_options = optimset('Display','iter','FunValCheck','on');
s_best = fmincon(neglogAP_func,s_init,A,b,Aeq,beq,lb,ub,nonlcon,paramopt_options);

H_hat = local_sens_fd(h_func,s_best,s_step);

if ~isempty(s_prior)
    s_cov = inv(H_hat'*R_inv*H_hat + Q_inv);
else
    s_cov = inv(H_hat'*R_inv*H_hat);
end
