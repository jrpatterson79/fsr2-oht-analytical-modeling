function varargout = Lev_Marq(test_list, s, y, R_inv, lambda_init, delta, soln, varargin)
    
% Lev_Marq: This function conducts non-linear gradient inversion using the Levenberg-Marquardt algorithm. 
% The inversion uses the fully-confined or leaky-confined analytical models developed by Rasmussen et al. (2003) as the forward model.

% Inputs:
%   test_list (num_obs x 4) - Matrix where each row contains the necessary information of a single oscillatory flow test to generate phasor coefficients for a given set of                        
%                             parameters. The columns are [Period(s) omega(1/s) Q_max(m^3/s) Radial distance(m)]
%   s (num_param x 1) - column vector containing the initial parameter guesses. 
%       Confined - [ln(T); ln(S)]
%       Leaky - [ln(T); ln(S); ln(L)]
%   y (num_freq*2 x 1) - column vector containing phasor coefficients representing the calibration data. If more than one frequency is used during inversion, the vector order is
%   R_inv (num_data x num_data) - Data error covariance matrix
%   lambda_init (scalar) - Initial Levenberg-Marquardt stabilization parameter
%   delta (num_param x 1) - Column vector of small increment to perturb parameters to numerically approximate the Jacobian matrix
%   varargin:
%       1: max_iter (scalar) - Maximum gradient step iterations
%       2: s_close (scalar) - Parameter change convergence criteria
%       3: obj_close (scalar) - Objective function value change convergence criteria

%Outputs:
%   s_hat (num_param x 1) - column vector containing optimal parameter estimates
%   varargout
%       1: s_update(num_iter x num_param) - Matrix of parameter values at each gradient step.
%       2: out_flag: Inversion output flag. 1 if converged, 0 if maximum iterations exceeded before convergence

% Code developed by Jeremy Patterson
% Created Dec 2018; Updated May 2022

%Closure Criteria
max_linevals = 100;
if nargin == 7
    max_iter = 75;
    obj_close = 1e-6;
    s_close = 1e-6;
elseif nargin == 8
    max_iter = varargin{1};
    obj_close = 1e-6;
    s_close = 1e-6;
elseif nargin == 9
    max_iter = varargin{1};
    s_close = varargin{2};
    obj_close = 1e-6;
elseif nargin == 10
    max_iter = varargin{1};
    s_close = varargin{2};
    obj_close = varargin{3};
end

% Forward Model
h = @(s) RasSoln(test_list, s, soln);
% Objective Function
obj_func = @(s) (1/2) * (y - h(s))' * R_inv * (y - h(s));

% Initiate parameter estimation
s_curr = s;
lambda = lambda_init;

num_func_evals = 0;
iter = 0

num_data = numel(y);
num_param = numel(s_curr);

while iter < max_iter
    
    % Calculate the modeled data with current and new parameters
    y_curr = h(s_curr);    
    
    %Calculate the jacob of the forward model
    [J] = jacob(s_curr, delta, test_list, soln);
    
    % Calculate the Gauss-Newton or Levenberg-Marquardt step
    step = -(J' * J + lambda * eye(num_param,num_param))\(-J'*(y_curr - y));
    
    % Do the line search along the Gauss-newton step direction
    step_obj = @(alpha) obj_func(s_curr + (alpha .* step));
    options = optimset('Display', 'iter', 'MaxFunEvals', max_linevals);
    alpha_best = fminsearch(step_obj, 0.5, options);    
    
    % Calculate relative parameter change
    s_new = s_curr + (alpha_best .* step);
    s_change = max(abs((s_new - s_curr) ./ s_curr))
    
    % Calculate relative objective function change
    obj_curr = log10(obj_func(s_curr));
    obj_new = log10(obj_func(s_new));
    obj_change = abs((obj_curr - obj_new)./obj_curr)
    
    if (obj_change <= obj_close) && (s_change <= s_close)
        s_hat = s_curr;
        out_flag = 1;
        varargout = {s_hat, s_update, out_flag};
        return;
    else
               
        if obj_new < obj_curr
            s_update(iter+1,:) = s_curr;
            
            s_curr = s_new;
            obj_curr = obj_new;
            
            if lambda <= 1e-12
                lambda = 1e-12
            else
                lambda = lambda * 1e-1
            end
            iter = iter + 1
            
        else
            s_update(iter+1,:) = s_curr;
            
            if lambda >= 1e10
                lambda = 1e10
            else
                lambda = lambda * 1e1
            end
            iter = iter + 1
        end
    end
end
s_hat = s_curr;
out_flag = 0;
varargout = {s_hat, s_update, out_flag};
warning('Max iterations exceeded')