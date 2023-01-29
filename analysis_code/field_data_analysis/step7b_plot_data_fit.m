% FSR^2 Oscillatory Flow Testing
% Step 7b - Plot modeled time series fit

% This code loads the single frequency parameter estimation results and
% reproduces Figure 6 as seen in:
% Patterson, Jeremy R. & Cardiff, Michael (2023). Do Simple Analytical
% Models Capture Complex Fractured Bedrock Hydraulics? Oscillatory Flow
% Tests Suggest Not. Groundwater. https://doi.org/

% Code developed by Jeremy Patterson
% Created Dec 2019; Updated May 2022

%% Clean Environment
close all; clear; clc

%% Specify Directory
file_dir = '/.../.../MAT_Files/'; % Add directory location to 2019PESTResults.mat file with parameter estimation results
addpath(genpath('/.../.../Func_Lib/')) % Add directory location to Func_Lib which contains all the needed function files to run this analysis

%% Load Data 
load ([file_dir '2019PESTResults.mat'])
srch = find(test_idx == 11 | test_idx == 15 | test_idx == 18 | test_idx == 19);
idx = srch(1:4:end);

input_dir = '/.../.../MAT_Files/Trimmed/'; % Provide the directory location to the trimmed mat files from the field data
filename = {'Test.11.mat',...
            'Test.15.mat',...
            'Test.18.mat',...
            'Test.19.mat'};

for j = 1 : numel(filename)
    file = [input_dir filename{j}];
    load(file)

    time{j} = test_time;
    obs_detrend{j} = detrend(obs_trim(:,1), 'linear');
    obs_mod{j} = y_opt_19(idx(j),1) .* cos(test_list(idx(j),2) .* test_time) +...
                -y_opt_19(idx(j),2) .* sin(test_list(idx(j),2) .* test_time);
end

% Figure 6 - Modeled time series fit with measured data
figure(6)
clf
for j = 1 : numel(filename)
    subplot(numel(filename),1,j)
    ax = gca;
    plot(time{j}, obs_detrend{j}, '.', 'Color', [0.7529 0 0], 'MarkerSize', 8)
    hold on
    plot(time{j}, obs_mod{j}, 'k-', 'LineWidth', 3)
    xlim([time{j}(1) time{j}(end)])
    ax.FontSize = 28;
    if j == 3
        ylabel('Head Change (m)')
    elseif j == 4
        xlabel('Time (s)')
    end
    title(['P = ' num2str(test_list(idx(j),1)) 's'], 'FontSize', 28)
end
set(gcf, 'Position', [100 100 2025 2025/2.25])
