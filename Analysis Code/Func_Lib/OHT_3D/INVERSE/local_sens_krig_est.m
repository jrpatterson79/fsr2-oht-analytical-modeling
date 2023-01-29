function [H] = local_sens_krig_est(h_func,s_curr,steps,filter,s_locs,GCF_param_func, GCF_params, varargin)

%local_sens_krig_est: Function which estimates a full sensitivity matrix by
%evaluating the sensitivity for a number of parameters and then kriging the
%values for the rest of the field. It is assumed that the sensitivity field
%can be adequately modeled as a random field with constant (unknown) mean.
%
%[H] =
%local_sens_krig_est(h_func,s_curr,steps,filter,s_locs,GCF_param_func,GCF_p
%arams, [opt_GCF_params])
%
%where:
%   -H is the output sensitivity matrix estimate (n x m)
%   -h_func is the forward model function which takes an (m x 1) vector of
%   parameters and returns an (n x 1) vector of measurements
%   -s_curr is the current value of the parameters (m x 1 vector) at which
%   the sensitivity is being evaluated
%   -steps is the step size for finite difference estimation of the
%   sensitivity. may either be an (m x 1) vector or a scalar.
%   -filter is an (m x 1) vector indicating which parameters should have
%   actual sensitivities computed. If a non-zero value is found in filter,
%   sensitivity for that parameter is computed via finite differences,
%   otherwise it is kriged.
%   -s_locs is the spatial location of each parameter value
%   -GCF_param_func is the generalized covariance function that is being
%   used for kriging. Must accept two inputs 1) a distance matrix of any
%   size and 2) a set of GCF structural parameter values.
%   -GCF_params is a set of structural parameters to be used in the
%   GCF_param_func. Can either be supplied as a single (1 x vp) vector of
%   structural parameters or a (n x vp) matrix, where a different
%   structural parameter set is used for each dependent observation (row in
%   the H matrix)
%   -[opt_GCF_params] optional. If set to 1, the values in GCF_params are
%   optimized before the kriging is carried out. 0 by default. The values
%   given in GCF_params are used as initial guesses.

opt_GCF_params = 0;
nreqin = 7;
if nargin > nreqin
    opt_GCF_params = varargin{1};
end

params_evaled = find(filter ~= 0);
params_notevaled = find(filter == 0);
m = size(s_curr,1);

%First, fill in the values at the locations where we're actually evaluating
[H_partial] = local_sens_partial(h_func,s_curr,steps,filter);
n = size(H_partial,1);
H = H_partial;

%Then, krig the rest
known_locs = s_locs(params_evaled,:);
unknown_locs = s_locs(params_notevaled,:);

for i = 1:1:n
    known_vals = (H_partial(i,params_evaled))';
    if size(GCF_params,1) == 1
        obs_GCF_params = GCF_params;
    else
        obs_GCF_params = GCF_params(i,:);
    end
    if opt_GCF_params == true
        [GCF_best_params] = krig_struct_est(known_locs,known_vals,GCF_param_func,obs_GCF_params);
    else
        GCF_best_params = obs_GCF_params;
    end
    GCF_func = @(dist) GCF_param_func(dist,GCF_best_params);
    [kriged_vals] = simple_krig(unknown_locs,known_locs,known_vals,GCF_func);
    clear GCF_func
    H(i,params_notevaled) = kriged_vals;
end