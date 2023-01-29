% FSR^2 Oscillatory Flow Testing Analysis Step 1 - Open each file
% containing a group of oscillatory flow tests, parse into individual
% oscillatory flow tests, and save individual test data into separate .mat
% files.

% This code represents the first data processing step associated with field
% experiments and analysis described in:
% Patterson, Jeremy R. & Cardiff, Michael (2023). Do Simple Analytical
% Models Capture Complex Fractured Bedrock Hydraulics? Oscillatory Flow
% Tests Suggest Not. Groundwater. https://doi.org/

% For this set of oscillatory flow experiments, multiple tests were
% conducted and saved as a group in a single output folder. This code opens
% the .txt file for each DAQ module and saves all of the output for a
% single group of tests as a .mat file for further processing. 

% Code written by: Jeremy Patterson 
% Created Aug 2019; Modified May 2022

%% Clean Environment
close all; clear; clc

%% Specify Directory
input_dir = '/.../.../';  % Specify the directory where raw test data is located
output_dir = '/.../.../'; % Specify the directory where parsed data should be saved .mat format
sys_dir = '/.../.../';    % If input directory has spaces in the directory path, specify system directory that Linux can read 

%% Import Text Files

% Folders with raw test data separated by field day
field_day = {'20.Aug';...
             '21.Aug';...
             '22.Aug';...
             '23.Aug';...
             '26.Aug';...
             '27.Aug';...
             '28.Aug'};

% Each signal is measured and recorded as a separate .txt file
filenames = {'14ZP1002-1.txt';...
             '14ZP1002-2.txt';...
             '14ZP1003-2.txt';
             '14ZP1004-2.txt';
             '14ZP1005-2.txt';
             '14ZP1006-2.txt'};

num_days = numel(field_day);
num_files = numel(filenames);

%% Parse and save data
for i = 1 : num_days
    % Create a vector of test group files for each test date
    dir = ['ls -1 ' sys_dir field_day{i} '/'];
    [~, folder_chr] = system(dir);
    
    folder_str = textscan(folder_chr, '%s', 'Delimiter', '/r/n');
    folderlist = folder_str{1};
    numf = numel(folderlist);
    
    % Open each test group folder individually
    for f = 1 : numf
        test_dir = [input_dir folderlist{f} '/'];
        
        % Open each pressure signal .txt file and parse into separate
        % tests to save as .mat files for further processing
        for j = 1 : num_files
            input_file = [input_dir field_day{i} '/' folderlist{f} '/' filenames{j}];
            FileIn = fopen(input_file);
            
            for k = 1 : 37
                metadata{j,k} = fgetl(FileIn);
            end
            
            FormatSpec = '%s %s %f';
            data = textscan(FileIn, FormatSpec);
            if j == 1
                atm = data;
            elseif j == 2
                stim = data;
            else
                obs{j-2} = data;
            end
        end
%         Uncomment the following lines to save individual test data as .mat file
%         output_file = [folderlist{f} '.Full.mat'];
%         save([output_dir output_file], 'atm', 'stim', 'obs')
        fclose(FileIn);
    end
end
