% FSR^2 Oscillatory Flow Testing Analysis
% Step 6 - Multi-Frequency Parameter Estimation

% This code loads the test list .mat file from the 2019 oscillatory flow
% testing. It then conducts a multi-frequency gradient-based
% Levenberg-Marquardt inversion under uncertainty to determine the optimal
% flow parameters that matches the collected data. The user chooses between
% the fully confined and leaky confined analytical solutions developed by
% Rasmussen et al. (2003) as the forward model to be used during inversion.

% This code reproduces Figures 7 & 11 seen in:
% Patterson, Jeremy R. & Cardiff, Michael (2023). Do Simple Analytical
% Models Capture Complex Fractured Bedrock Hydraulics? Oscillatory Flow
% Tests Suggest Not. Groundwater. https://doi.org/

% Code developed by Jeremy Patterson
% Created Dec 2019; Updated May 2022
%% Clean Environment
close all; clear; clc

%% Specify Directory
file_dir = '/.../.../'; % Specify directory location containing the test list .mat file specified on line 17 (created in step 5)
addpath('/.../.../') % Add directory location of the folder Func_Lib which contains necessary function files to run the parameter estimation

%% Load Data
load ([file_dir 'TestList_19.mat'])
num_test = numel(test_list(:,1));

%% General Inversion Parameters
% Analytical Forward Models
soln = {'confined', 'leaky'};

% Gradient Inversion Setup
max_iter = 100;
s_close = 1e-4;   % Parameter change convergence criteria
obj_close = 1e-4; % Objective function value change convergence criteria

s_init = {[-12; -12], [-10; -12; -21]}; % Initial parameter guesses 
delta = {0.1*ones(2,1), 0.1*ones(3,1)}; % Parameter step for sensitivity (Jacobian) calculation 
lambda = 1e1;                           % Levenberg-Marquardt stabilization parameter

% Radius, Period, Flow rate list
rad_list = unique(test_list(:,4));
p_list = unique(test_list(:,1));
q_list = unique(test_list(:,3));

%% Confined Multi-Frequency Inversion - Heterogenous Flow Analysis
r_idx = find(test_list(:,4) == rad_list(end));
p_idx = find(test_list(:,1) == p_list(end));

% Multi-frequency inversion of all tests with r ~= 16 m
[y_mod_rad, rmse_rad, data_err_rad, ~, s_opt_rad, s_unc_rad] = MultiFreqPEST(test_list, y_obs, cov_y, s_init{1}, lambda, delta{1}, r_idx, soln{1});

% Multi-frequency inversion of all tests with P >= 300 s
[y_mod_per, rmse_per, data_err_per, ~, s_opt_per, s_unc_per] = MultiFreqPEST(test_list, y_obs, cov_y, s_init{1}, lambda, delta{1}, p_idx, soln{1});

% Figure 7
figure(7)
clf
subplot(1,2,1)
ax = gca;
plot([-0.1 0.1], [-0.1 0.1], 'k-', 'LineWidth', 2)
hold on
plot(y_obs(r_idx,1), y_mod_rad(1:2:end-1), 'ko', 'MarkerSize', 12,...
     'MarkerFaceColor', [0 0.4470 0.7410])
plot(y_obs(r_idx,2), y_mod_rad(2:2:end), 'ko', 'MarkerSize', 12,...
        'MarkerFaceColor', [0 0.4470 0.7410])
grid on
axis([-0.1 0.1 -0.1 0.1])
xlabel('\Phi_{obs}')
ylabel('\Phi_{mod}')
ax.FontSize = 34;
ax.XTick = [-0.1:5e-2:0.1];
ax.YTick = [-0.1:5e-2:0.1];
text(0.09, -0.0925, 'A', 'FontSize', 28)

subplot(1,2,2)
ax = gca;
plot([-0.25 0.25], [-0.25 0.25], 'k-', 'LineWidth', 2)
hold on
plot(y_obs(p_idx,1), y_mod_per(1:2:end-1), 'kd', 'MarkerSize', 12,...
     'MarkerFaceColor', [0.8500 0.3250 0.0980])
plot(y_obs(p_idx,2), y_mod_per(2:2:end), 'kd', 'MarkerSize', 12,...
        'MarkerFaceColor', [0.8500 0.3250 0.0980])
grid on
axis([-0.25 0.25 -0.25 0.25])
xlabel('\Phi_{obs}')
ylabel('\Phi_{mod}')
ax.FontSize = 34;
ax.XTick = [-0.25:0.1:0.25];
ax.YTick = [-0.25:0.1:0.25];
text(0.225, -0.23, 'B', 'FontSize', 28)
set(gcf, 'Position', [100 100 2025 2025/2.6667])

%% Multi-Frequency Inversion - Leaky Fracture Analysis 
lky_idx = find(test_list(:,1) >= 300);
[y_mod_lky, rmse_lky, data_err_lky, ~, s_opt_lky, s_unc_lky] = MultiFreqPEST(test_list, y_obs, cov_y, s_init{2}, lambda, delta{2}, lky_idx, soln{2});

% Figure 11
figure(11)
clf
ax = gca;
plot([-0.25 0.25], [-0.25 0.25], 'k-', 'LineWidth', 2)
hold on
plot(y_obs(lky_idx,1), y_mod_lky(1:2:end-1), 'kd', 'MarkerSize', 12,...
     'MarkerFaceColor', [0.8500 0.3250 0.0980])
plot(y_obs(lky_idx,2), y_mod_lky(2:2:end), 'kd', 'MarkerSize', 12,...
        'MarkerFaceColor', [0.8500 0.3250 0.0980])
grid on
axis([-0.25 0.25 -0.25 0.25])
xlabel('\Phi_{obs}')
ylabel('\Phi_{mod}')
ax.FontSize = 34;
ax.XTick = [-0.25:0.1:0.25];
ax.YTick = [-0.25:0.1:0.25];
set(gcf, 'Position', [100 100 975 975/1.3333])