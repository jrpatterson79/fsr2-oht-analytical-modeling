function [theta_opt,varargout] = RML_data_weighting(d,forward_func,params_init,theta_init,d_sd_func,varargin)

%RML_data_weighting: Function which calculates optimal data weights using
%linearized restricted maximum likelihood (RML) theory. Useable for "bulk"
%parameter forward models, where there are generally fewer parameters than
%unknowns.
%
%[theta_opt,{R_mat}] =
%RML_data_weighting(d,forward_func,params_init,theta_init,d_sd_func,{A},{b}
%,{Aeq},{beq},{lb},{ub},{nonlcon},{param_step},{theta_lb})
%
%where:
%   -theta_opt is the optimal structural parameter values (data weights)
%   obtained after optimization.
%   -R_mat (optional output) is the best estimate of the data covariance
%   matrix generated using the best theta_opt.
%   -d is the data (n x 1) vector.
%   -forward_func is the forward model (map from Rm --> Rn).
%   -params_init is an initial guess at the model's parameter values (m x
%   1) vector.
%   -theta_init is an initial guess at the structural (data weighting)
%   parameter values (p x 1) vector.
%   -d_sd_func is a function (Rp --> Rn) which converts the structural
%   parameter values (thetas) into a vector of expected data standard
%   deviations. NOTE: it is assumed that data errors are uncorrelated and
%   thus the data error matrix can simply be described by its diagonal
%   elements. (the square of the data standard deviation vector)
%   -A,b,Aeq,beq,lb,ub,nonlcon are constraints that must be satisfied on
%   the parameter values, and are used during parameter optimization. See
%   fmincon for more info.
%   -param_step is a (m x 1) vector containing parameter step sizes to be
%   used when approximating the H sensitivity matrix. By default, steps are
%   1e-3.
%   -theta_lb is a (p x 1) vector containing lower bounds on the theta
%   values to be used. By default, these are all 1e-4.

theta_best = theta_init;
theta_curr = theta_init - 1; 
params_best = params_init;

num_reqin = 5; num_constrargin = 7;
num_unknowns = length(params_init);

param_step = 1e-3*ones(num_unknowns,1);
theta_lb = 1e-4*ones(size(theta_init));

A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
nonlcon = [];

if nargin > num_reqin
    A = varargin{1};
    b = varargin{2};
    Aeq = varargin{3};
    beq = varargin{4};
    lb = varargin{5};
    ub = varargin{6};
    nonlcon = varargin{7};
end

if nargin > (num_reqin+num_constrargin)
    if ~isempty(varargin{num_constrargin+1})
        param_step = varargin{num_constrargin+1};
    end
end
if nargin > (num_reqin+num_constrargin+1)
    if ~isempty(varargin{num_constrargin+2})
        theta_lb = varargin{num_constrargin+2};
    end
end

while any(theta_curr ~= theta_best)
    
    d_sd_vec = d_sd_func(theta_init);
    negloglike_func = @(pvec) .5*sum(((d-forward_func(pvec))./d_sd_vec).^2);

    display('Finding best parameters, for linearized theory');
    paramopt_options = optimset('Display','iter','FunValCheck','on');
    params_best = fmincon(negloglike_func,params_init,A,b,Aeq,beq,lb,ub,nonlcon,paramopt_options);

    %Estimate linearized uncertainty around best estimate
    H = local_sens_fd(forward_func,params_best,param_step);

    d_p = d - forward_func(params_best) + H*params_best;
    T = (null(H'))';
    z = T*d_p;

    R_func = @(theta) diag(d_sd_func(theta).^2);
    theta_NLL = @(theta) .5*sum(log(eig(T*R_func(theta)*T'))) + .5*z'/(T*R_func(theta)*T')*z;
    display('Optimizing theta values');
    thetaopt_options = optimset('Display','iter');
    theta_curr = theta_best;
    theta_best = fmincon(theta_NLL,theta_init,[],[],[],[],theta_lb,[],[],thetaopt_options);
    
end

theta_opt = theta_best;

if nargout > 1
    varargout{1} = R_func(theta_opt);
end