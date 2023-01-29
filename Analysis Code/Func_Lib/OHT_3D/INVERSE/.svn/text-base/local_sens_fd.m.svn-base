function H_lin = local_sens_fd(h_func,params,p_step)

%local_sens_fd: Simple function for calculating the sensitivity matrix via
%finite differences using any forward model.
%
%[H_lin] = local_sens_fd(h_func,params,p_step)
%
%where:
%   -h_func is an anonymous function that takes (mx1) parameters as input
%   and outputs (nx1) observations 
%   -params is a (mx1) vector of parameter values at which the sensitivity
%   is being evaluated
%   -p_step is the step size used when calculating the finite difference
%   approximation. it may either be a scalar or an (mx1) vector of values
%   to be used.

m = numel(params);
d_init = h_func(params);
n = size(d_init,1);
p = numel(p_step);

H_lin = zeros(n,m);
parfor i = 1:m
    tic;
    if p ~= 1
        this_p_step = p_step(i);
    else
        this_p_step = p_step;
    end
    try_params = params;
    try_params(i) = try_params(i) + this_p_step;
    d_new = h_func(try_params);
    H_lin(:,i) = (d_new - d_init)./this_p_step;
    %Uncomment to show derivative calculation progress
    disp(['Unknown # ', num2str(i), ' of ', num2str(m)]);
    toc
end