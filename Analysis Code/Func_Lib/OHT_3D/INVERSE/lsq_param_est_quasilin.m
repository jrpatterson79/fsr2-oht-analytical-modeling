function [s_best, s_cov] = lsq_param_est_quasilin(d,h_func,R,s_init,varargin)

%lsq_param_est_quasilin: Function which performs least squares parameter estimation
%for over-determined inverse problems. Can estimate best values of the
%parameters either with or without prior distribution information.
%[s_best, s_cov] =
%lsq_param_est_quasilin(d,h_func,R,s_init,[s_step],[s_prior],[s_covar])
%where:
%   -d is the data (n x 1)
%   -h_func is the forward mapping from parameters to data (R^m --> R^n),
%   supplied as an anonymous function with a single (m x 1) vector input
%   and a single (n x 1) vector output
%   -R is the data error covariance matrix (n x n)
%   -s_init is an initial guess at the parameter values (m x 1)
%   -[s_step] (optional) is a stepsize to be used for the parameters during
%   computation of the approximate Jacobian with finite difference methods.
%   s_step(i) is the step size that will be used when evaluating the ith
%   column of the Jacobian.
%   -[s_prior] (optional) is a prior estimate of parameter values (m x 1)
%   -[s_covar] (optional) is a prior covariance matrix (m x m) of the
%   parameters. If s_prior is supplied but s_covar is empty or not
%   supplied, then s_covar is assumed to be eye(m,m)

m = numel(s_init);
s_step = .01.*s_init;
s_prior = [];
s_covar = eye(m);
max_linsearch = 40;
tol_linsearch = 1e-3;

%Assign additional arguments, if supplied by user
[s_step,s_prior,s_covar] = process_extra_args(varargin,s_step,s_prior,s_covar);

R_inv = inv(R);
R_npf = chol(R_inv);

s_hat = s_init;
% h_of_s_hat = h_func(s_hat);
if isempty(s_prior)
    obj_func = @(s) .5*sum(((d - h_func(s))'*R_npf).^2);
else
    Q_inv = inv(s_covar);
    obj_func = @(s) .5*sum(((d - h_func(s))'*R_npf).^2) + .5*(s - s_prior)'*Q_inv*(s - s_prior);
end
obj_hat = obj_func(s_hat);

rel_obj_thresh = .001;
rel_s_thresh = .001;

rel_obj_change = 1;
rel_s_change = 1;

while (rel_obj_change > rel_obj_thresh) || (rel_s_change > rel_s_thresh)
    s_tilde = s_hat;
    obj_tilde = obj_hat;
    
    H_tilde = local_sens_fd(h_func,s_tilde,s_step);
    dp = d - h_func(s_tilde) + H_tilde*s_tilde;
    if isempty(s_prior)
        s_hat = (H_tilde'*R_inv*H_tilde)\(H_tilde'*R_inv*dp);
    else
        s_hat = s_prior + (H_tilde'*R_inv*H_tilde + Q_inv)\(H_tilde'*R_inv*(dp - H_tilde*s_prior));
    end
    h_of_s_hat = h_func(s_hat);
    obj_hat = obj_func(s_hat);
    
    clear linsearch_func
    linsearch_func = @(x) obj_func(s_tilde + (s_hat - s_tilde).*x);
    if obj_hat < obj_tilde
        disp('Performing linesearch, starting at s_hat');
        lin_options = optimset('MaxIter', max_linsearch,'Display','iter','TolFun',tol_linsearch*obj_hat);
        step = fminsearch(linsearch_func,1,lin_options)
    else
        disp('Performing linesearch, starting at s_tilde');
        lin_options = optimset('MaxIter', max_linsearch,'Display','iter','TolFun',tol_linsearch*obj_tilde);
        step = fminsearch(linsearch_func,0,lin_options)
    end
    s_hat = s_tilde + step.*(s_hat - s_tilde);
    obj_hat = obj_func(s_hat);

    rel_obj_change = (obj_tilde - obj_hat)./obj_tilde;
    rel_s_change = max(abs((s_tilde - s_hat)./s_tilde));
    
    obj_best = min(obj_hat, obj_tilde);
    disp(['Best objective = ', num2str(obj_best)]);
    
end

if obj_tilde < obj_hat
    s_best = s_tilde;
    H_best = H_tilde;
else
    s_best = s_hat;
    H_best = local_sens_fd(h_func,s_best,s_step);
end

if isempty(s_prior)
    s_cov = inv(H_best'*R_inv*H_best);
else
    %For gaussian distributions, variance is equal to inverse Hessian of L
    s_cov = inv(H_best'*R_inv*H_best + Q_inv);
end
