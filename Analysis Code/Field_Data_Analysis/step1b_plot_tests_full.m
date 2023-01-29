% FSR^2 Oscillatory Flow Testing Analysis
% Step 1b - Open and plot each test file individually

% This code represents an intermediate data processing step associated with
% field experiments and analysis described in:
% Patterson, Jeremy R. & Cardiff, Michael (2023). Do Simple Analytical
% Models Capture Complex Fractured Bedrock Hydraulics? Oscillatory Flow
% Tests Suggest Not. Groundwater. https://doi.org/

% This file opens and plots the stimulation pressure signal for each group
% of tests as a quality control check and to determine trim points for
% follow on steady-periodic analysis.

% Code written by: Jeremy Patterson
% Created Aug 2019; Modified May 2022

%% Clean Environment
close all; clear; clc

%% Specify Directory
input_dir = '/.../.../';    % Specify the directory where raw test data is located
sys_dir = '/.../.../';      % If input directory has spaces in the directory path, specify system directory that Linux can read 

%% Load Test Database

% Create a vector of .mat files to plot each pressure indivdual pressure
% signal as an initial quality control check
dir = ['ls -1 ' sys_dir];
[~, tests_chr] = system(dir);

mat_files = textscan(tests_chr, '%s', 'Delimiter', '/r/n');
mat_list = mat_files{1};
num_matfiles = numel(mat_list);

for i = 1 : num_matfiles
    
    % Load individual test file
    input_file = [input_dir mat_list{i}];
    load (input_file)
    disp(mat_list{i})
    
    % Create time vector
    time_datenum = datenum(strcat(atm{1}, {' '}, atm{2}));    
    stim_datenum = datenum(strcat(stim{1}, {' '}, stim{2}));
    stim_time = (stim_datenum(1:end) - time_datenum(1)) * 86400;
    
    figure
    ax = gca;
    plot(stim_time, stim{3}, '-', 'LineWidth', 2)
    xlim([0 stim_time(end)])
    xlabel('Time (s)')
    ylabel('Pressure Change (kPa)')
    ax.FontSize = 28;
    set(gcf, 'Position', [100 100 1900 600])
    pause        
end
