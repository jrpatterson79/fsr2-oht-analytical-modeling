% FSR^2 Oscillatory Flow Testing Analysis
% Step 6 - Single Frequency Parameter Estimation

% This code represents the final data processing step associated with field
% experiments and analysis described in:
% Patterson, Jeremy R. & Cardiff, Michael (2023). Do Simple Analytical
% Models Capture Complex Fractured Bedrock Hydraulics? Oscillatory Flow
% Tests Suggest Not. Groundwater. https://doi.org/

% This code loads the test list .mat file from the 2019 oscillatory flow
% testing. It then conducts gradient-based Levenberg-Marquardt inversion
% under uncertainty to determine the optimal flow parameters that matches
% the collected data. The fully confined analytical solution developed by
% Rasmussen et al. (2003) is the forward model used during inversion. 

% Code developed by Jeremy Patterson
% Created Dec 2019; Updated May 2022

%% Clean Environment
close all; clear; clc

%% Specify Directory
file_dir = '/.../.../'; % Specify directory location containing the test list .mat file specified on line 17 (created in step 5)
save_dir = '/.../.../'; % Specify directory to save parameter estimation results as .mat file for visualization
print_dir = '/.../.../'; % Specify preferred directory to print pressure signal figures for each test
addpath('/.../.../') % Add directory location of the folder Func_Lib which contains necessary function files to run the parameter estimation

%% Load TestList
load ([file_dir 'TestList_19.mat'])
num_test = numel(test_list(:,1));
soln = 'confined';

%% Single Frequency Parameter Estimation
% Gradient Inversion Setup
max_iter = 100;
s_close = 1e-4;   % Parameter change convergence criteria
obj_close = 1e-4; % Objective function value change convergence criteria
lambda = 1e1;     % Levenberg-Marquardt stabilization parameter

if strcmp(soln, 'confined') == 1
    s_init = [-12; -12]; % Initial parameter guess
    delta = [0.1; 0.1];  % Parameter step for sensitivity (Jacobian) calculation
elseif strcmp(soln, 'leaky') == 1
    error('Leaky solution requires a multi-frequency approach')
end

s_opt_19 = zeros(num_test,3); 
s_unc_19 = zeros(num_test,3);

for i = 1 : num_test
    R_inv = inv(reshape(cov_y(i,:), 2, 2)); % Data error covariance matrix
    
    % LM Inversion and Modeled Data
    [pest, ~, out_flag(i,:)] = Lev_Marq(test_list(i,1:4), s_init, y_obs(i,:)', R_inv, lambda, delta, soln);   
    if pest(1) < -15 || pest(2) > 0
        [pest, ~, out_flag(i,:)] = Lev_Marq(test_list(i,1:4), [-10; -12.5], y_obs(i,:)', R_inv, lambda, delta, soln);
        s_opt_19(i,:) = [pest(1) pest(2) pest(1)-pest(2)];
    else
        s_opt_19(i,:) = [pest(1) pest(2) pest(1)-pest(2)];            % Optimal parameters [T S D]
    end
    y_opt_19(i,:) = RasSoln(test_list(i,1:4), s_opt_19(i,1:2), soln); % Modeled phasors with optimal parameters
        
    % Linearized Parameter Uncertainty Analysis
    J = jacob(pest, delta, test_list(i,1:4), soln);      % Jacobian
    s_cov = inv(J' * R_inv * J);                         % Parameter error covariance
    param_unc = 1.96 * sqrt(diag(s_cov));                % 95% confidence interval (T & S)
    D_unc = 1.96 * sqrt([1 -1] * s_cov * [1; -1]);       % Diffusivity confidence interval
    s_unc_19(i,:) = [param_unc(1) param_unc(2) D_unc]; 
end

%% Save Results (Uncomment the following lines to save output as .mat file)
% save_file = '2019PESTResults.mat';
% save([save_dir save_file])
