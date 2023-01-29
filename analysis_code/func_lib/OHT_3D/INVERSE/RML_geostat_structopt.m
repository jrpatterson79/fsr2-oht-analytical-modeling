function [theta_opt1, theta_opt2] = RML_geostat_structopt(y,H,X,R_func,Q_func,theta_init1,theta_init2,varargin)

theta_min1 = [];
theta_min2 = [];
theta_max1 = [];
theta_max2 = [];
optim_options = [];

%Sm_value is used as a default lower bound for the theta parameters
sm_value = 1e-6;
%Lg_value is the negative log-likelihood value returned if a value outside
%of the bounds is tried for the theta parameters.
lg_value = 1000;

[theta_min1, theta_min2, theta_max1, theta_max2, optim_options] = ...
    process_extra_args(varargin, theta_min1, theta_min2, theta_max1, ...
    theta_max2, optim_options);

if size(theta_init1,1) == 1
    theta_init1 = theta_init1';
end
if size(theta_init2,1) == 1
    theta_init2 = theta_init2';
end
theta_init = [theta_init1; theta_init2];
num_theta = numel(theta_init)

if size(theta_min1,1) == 1
    theta_min1 = theta_min1';
end
if size(theta_min2,1) == 1
    theta_min2 = theta_min2';
end

theta_mins = [theta_min1; theta_min2];

if size(theta_max1,1) == 1
    theta_max1 = theta_max1';
end
if size(theta_max2,1) == 1
    theta_max2 = theta_max2';
end

theta_maxs = [theta_max1; theta_max2];

if isempty(theta_mins)
    theta_mins = sm_value*ones(num_theta,1);
end

PHI = H*X;
T = null(PHI')';
num_Rth = numel(theta_init1);
num_Qth = numel(theta_init2);

PSI_func = @(theta) H*Q_func(theta((num_Rth+1):end))*H' + R_func(theta(1:num_Rth));

nll_func = @(theta) .5*sum(log(eig(T*PSI_func(theta)*T'))) ...
    + .5*y'*T'*((T*PSI_func(theta)*T')\(T*y));

constr_nll_func = @(theta) model_constrainer(nll_func,theta,lg_value,[],[],[],[],theta_mins,theta_maxs,[]);

if isempty(optim_options)
    optim_options = optimset('Display','iter');
end

% [theta_opt] = fminsearch(nll_func,theta_init,optim_options);
[theta_opt] = fminunc(nll_func,theta_init,optim_options);
% [theta_opt] = fmincon(constr_nll_func,theta_init,[],[],[],[],theta_mins,theta_maxs,[],optim_options);

theta_opt1 = theta_opt(1:num_Rth);
theta_opt2 = theta_opt((num_Rth+1):end);
