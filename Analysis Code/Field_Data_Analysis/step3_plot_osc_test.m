% FSR^2 Oscillatory Flow Testing Analysis
% Step 3 - Plot Trimmed Oscillatory Flow Signals

% This code represents an intermediate data processing step associated with field
% experiments and analysis described in:
% Patterson, Jeremy R. & Cardiff, Michael (2023). Do Simple Analytical
% Models Capture Complex Fractured Bedrock Hydraulics? Oscillatory Flow
% Tests Suggest Not. Groundwater. https://doi.org/

% This script loads and plots the trimmed stimulation and pressure signals
% for each individual oscillatory flow test as a quality control check to
% ensure the steady-periodic state of the trimmed signal.

% Code developed by Jeremy Patterson
% Created Sep 2019; Modified May 2022

%% Clean Environment
close all; clear; clc

%% Specify Directory
database_dir = '/.../.../'; % Specify directory location for Excel testing database on line 20
input_dir = '/.../.../'; % Specify directory location containing trimmed .mat files
sys_dir = '/.../.../';    % If input directory has spaces in the directory path, specify system directory that Linux can read 
print_dir = '/.../.../'; % Specify preferred directory to print pressure signal figures for each test

%% Load Data
% Load test database as table
database_filename = '2019 Test Database.xlsx';
database = [database_dir database_filename];
TestDatabase = readtable(database);
num_tests = numel(TestDatabase.TestNumber);

for i = 1 : num_tests
    % Load trimmed .mat file
    input_file = [input_dir 'Test.' num2str(i) '.mat'];
    load (input_file)
    % Specify figure printing filename
    test_name = ['Test_' num2str(i)];
    print_file = [print_dir test_name];
    
    % Plot stimulation and observation pressure signals
    figure
    clf
    subplot(5,1,1)
    ax = gca;
    plot(test_time, stim_trim, 'r.', 'MarkerSize', 6)
    ylabel('\Delta h (m)')
    xlim([0 test_time(end)])
    ax.FontSize = 18;
    title('Stim Head ', 'FontSize', 16)

    subplot(5,1,2)
    ax = gca;
    plot(test_time, obs_trim(:,1), 'b.', 'MarkerSize', 6)
    ylabel('\Delta h (m)')
    xlim([0 test_time(end)])
    ax.FontSize = 18;
    title('Obs Head 1', 'FontSize', 16)
    
    subplot(5,1,3)
    ax = gca;
    plot(test_time, obs_trim(:,2), 'b.', 'MarkerSize', 6)
    ylabel('\Delta h (m)')
    xlim([0 test_time(end)])
    ax.FontSize = 18;
    title('Obs Head 2', 'FontSize', 16)
    
    subplot(5,1,4)
    ax = gca;
    plot(test_time, obs_trim(:,3), 'b.', 'MarkerSize', 6)
    ylabel('\Delta h (m)')
    xlim([0 test_time(end)])
    ax.FontSize = 18;
    title('Obs Head 3', 'FontSize', 16)
    
    subplot(5,1,5)
    ax = gca;
    plot(test_time, obs_trim(:,4), 'b.', 'MarkerSize', 6)
    xlabel('Time (s)')
    ylabel('\Delta h (m)')
    xlim([0 test_time(end)])
    ax.FontSize = 18;
    title('Obs Head 4', 'FontSize', 16)    
    set(gcf, 'Position', [0 0 1920 1080])
%     print(print_file, '-dpng')
%     close
end
