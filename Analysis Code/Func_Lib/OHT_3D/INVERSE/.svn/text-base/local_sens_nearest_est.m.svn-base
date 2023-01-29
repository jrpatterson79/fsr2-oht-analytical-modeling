function [H] = local_sens_nearest_est(h_func,s_curr,steps,filter,s_locs)

%local_sens_nearest_est: Function which estimates a full sensitivity matrix
%by evaluating the sensitivity for a number of parameters and then
%interpolating the values for the rest of the field using nearest neighbor
%interpolation. 
%
%[H] =
%local_sens_nearest_est(h_func,s_curr,steps,filter,s_locs)
%
%where:
%   -H is the output approximate sensitivity matrix (n x m)
%   -h_func is the forward model function which takes an (m x 1) vector of
%   parameters and returns an (n x 1) vector of measurements
%   -s_curr is the current value of the parameters (m x 1 vector) at which
%   the sensitivity is being evaluated
%   -steps is the step size for finite difference estimation of the
%   sensitivity. may either be an (m x 1) vector or a scalar.
%   -filter is an (m x 1) vector indicating which parameters should have
%   actual sensitivities computed. If a non-zero value is found in filter,
%   sensitivity for that parameter is computed via finite differences,
%   otherwise it is interpolated.
%   -s_locs is the spatial location of each parameter value (m x d), for a
%   d-dimensional problem

params_evaled = (filter ~= 0);
params_evaled_ind = find(filter ~= 0);
params_notevaled = (filter == 0);
params_notevaled_ind = find(filter == 0);
m = size(s_curr,1);

%First, fill in the values at the locations where we're actually evaluating
[H_partial] = local_sens_partial(h_func,s_curr,steps,filter);
n = size(H_partial,1);
H = H_partial;

%Then, assign the rest based on nearest neighbors (use of dsearchn)
known_locs = s_locs(params_evaled_ind,:);
unknown_locs = s_locs(params_notevaled_ind,:);

nearest_array = dsearchn(known_locs,unknown_locs);
for i = 1:1:n
    H(i,params_notevaled_ind) = H(i,params_evaled_ind(nearest_array));
end
