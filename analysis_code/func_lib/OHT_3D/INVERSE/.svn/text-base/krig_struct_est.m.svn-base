function [struct_params] = krig_struct_est(data_loc,data_val,GCF_param_func,init_params,varargin)

%krig_struct_est: Function which estimates the structural parameters of the
%variogram using general nonlinear estimation routines (fmincon). For the
%linear variogram, exact solutions exist and this function may be
%inefficient. For other variograms, this function represents a numerical
%method for solving the structural parameter estimation problem. 
%
%[struct_params] =
%krig_struct_est(data_loc,data_val,GCF_param_func,init_params,[drift_func],
%[theta_mins],[theta_maxs],[optim_options])
%where:
%   -struct_params (vp x 1) are the optimized structural parameters of the variogram
%   -data_loc (n x dim) is the list of measurement locations
%   -data_val (n x 1) is the list of measured values at those locations
%   -GCF_param_func is a function which takes two arguments, a set of
%   distances (any sized matrix) and the variogram parameters (vp x 1), and
%   returns the generalized covariance matrix.
%   -init_params is an initial guess at the variogram parameters (vp x 1)
%   -drift_func (optional) is a function which takes the unknown locations
%   and calculates the drift function values for each point. If not
%   supplied, a constant mean is assumed
%   -theta_mins (optional) is a (vp x 1) set of minimum values for the
%   covariance parameters. If not supplied, a value close to 0 is assumed
%   -theta_maxs (optional) is a (vp x 1) set of maximum values for the
%   covariance parameters.
%   -optim_options (optional) is a option set to be used during
%   optimization. The optimization routine utilized is fmincon, so any
%   options supplied should apply to this routine.

drift_func = [];
theta_mins = []; 
theta_maxs = [];
optim_options = [];

[drift_func theta_mins theta_maxs optim_options] = process_extra_args(varargin,...
    drift_func,theta_mins,theta_maxs,optim_options);

num_data = size(data_loc,1);
vp = numel(init_params);
if ~isempty(drift_func)
    X = drift_func(data_loc);
else
    X = ones(num_data,1);
end
p = size(X,2);

%Generally, covaraince parameters must be positive
if isempty(theta_mins);
    theta_mins = 1e-6*ones(vp,1);
end
%Theta_maxs are not set, since they can generally have any value above 0

%Transformation matrix for creating generalized increments
T = null(X')';
z = T*data_val;

data_dist = dimdist(data_loc,[],1);

nll_func = @(params) 0.5*(sum(log(eig(T*GCF_param_func(data_dist,params)*T')))) + 0.5*z'*((T*GCF_param_func(data_dist,params)*T')\z);

if isempty(optim_options)
    optim_options = optimset('Display','iter');
end

struct_params = fmincon(nll_func,init_params,[],[],[],[],theta_mins,theta_maxs,[],optim_options);
