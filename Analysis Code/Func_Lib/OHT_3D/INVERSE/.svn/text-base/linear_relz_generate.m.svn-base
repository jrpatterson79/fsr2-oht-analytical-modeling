function [param_relz] = linear_relz_generate(mean,cov_mat,num_relz,varargin)

%linear_relz_generate: Function for generating conditional realizations of
%a parameter vector given its mean and covariance matrix (must be positive
%definite!). Realizations can be constrained to meet linear and nonlinear
%constraints as well.
%
%[param_relz] = linear_relz_generate(mean,cov_mat,num_relz,{A},{b},{Aeq}, 
%{beq},{lb},{ub},{nonlcon})
%
%where:
%   -param_relz is a matrix (m x num_relz) where each column is a parameter
%   realization.
%   -mean is the mean (m x 1) vector for the parameters
%   -cov_mat is the covariance matrix (m x m) of the parameters
%   -num_relz is the number of realizations that are desired (scalar)
%   -A,b,Aeq, beq,lb,ub,nonlcon (optional) are linear and non-linear
%   constraints that will be met by conditional realizations, if required.

num_params = size(mean,1);
num_reqin = 3;
constraints_supplied = false;

if nargin > num_reqin
    A = varargin{1};
    b = varargin{2};
    Aeq = varargin{3};
    beq = varargin{4};
    lb = varargin{5};
    ub = varargin{6};
    nonlcon = varargin{7};
    constraints_supplied = true;
end

Qpf = (chol(cov_mat))';

param_relz = zeros(num_params,num_relz);
if constraints_supplied == false
    meanmat = repmat(mean,1,num_relz);
    randmat = randn(num_params,num_relz);
    param_relz = meanmat + Qpf*randmat;
else
    i = 1;
    while i <= num_relz
        rand_vec = randn(num_params,1);
        trial_param_relz = mean + Qpf*rand_vec;
        func_calculator = @(params) params;
        %The below lines make sure that the currently generated realization
        %meets all constraints
        outside_value = zeros(num_params,1);
        test_vec = model_constrainer(func_calculator,trial_param_relz,outside_value,A,b,Aeq,beq,lb,ub,nonlcon);
        if ~all(test_vec == outside_value)
            param_relz(:,i) = trial_param_relz;
            i = i + 1;
        end
    end
end

    

