function [H_partial] = local_sens_partial(h_func, s_curr, steps, filter)

%local_sens_partial: Function for calculating the sensitivity matrix for a
%specified set of parameters (out of a long vector of parameters)
%
%[H_partial] = local_sens_partial(h_func,s_curr,steps,filter)
%where:
%   -H_partial is the partially-developed sensitivity matrix. Columns of
%   H_partial will be 0 if the sensitivity is not calculated with respect
%   to these parameters
%   -h_func is the forward model function which takes a (m x 1) vector of
%   parameters and produces an (n x 1) vector of measurements
%   -s_curr is the local estimate of the parameters around which the
%   sensitivity matrix should be evaluated, an (m x 1) vector
%   -steps is the vector (or scalar) of step sizes that should be used when
%   calculating the sensitivities using the finite difference method. It
%   may either be an (m x 1) vector of values, or a scalar.
%   -filter is a boolean (m x 1) vector that tells whether sensitivities
%   should be calculated for a specific parameter. If filter(i) ~= 0, the
%   sensitivity with respect to parameter i is calculated, otherwise it is
%   not.

h_curr = h_func(s_curr);

n = size(h_curr,1);
m = size(s_curr,1);
p = numel(steps);
H_partial = zeros(n,m);
H_partial_unsorted = zeros(n,m);

eval_list = find(filter ~= 0);
num_evaled = numel(eval_list);

parfor i = 1:1:num_evaled
    if p ~= 1
        this_step = steps(i);
    else
        this_step = steps;
    end
    param_ed = eval_list(i);
    s_new = s_curr;
    s_new(param_ed) = s_new(param_ed) + this_step;
    h_new = h_func(s_new);
    H_partial_unsorted(:,i) = (h_new - h_curr)./this_step;
    disp(['Unknown #', num2str(i), ' of ', num2str(num_evaled)]);
end
    
for i = 1:1:num_evaled
    H_partial(:,eval_list(i)) = H_partial_unsorted(:,i);
end

    