% FSR^2 Oscillatory Flow Testing Analysis
% Step 2 - Trim Oscillatory Flow Signals 

% This code represents the second data processing step associated with field
% experiments and analysis described in:
% Patterson, Jeremy R. & Cardiff, Michael (2023). Do Simple Analytical
% Models Capture Complex Fractured Bedrock Hydraulics? Oscillatory Flow
% Tests Suggest Not. Groundwater. https://doi.org/

% This file loads the full .mat file for each individual test. The data is
% resampled to a constant time step and then trimmed to the identified
% steady-periodic portion of the signal. The output is saved as a separate
% .mat file for follow-on analysis.

% Code created by Jeremy Patterson
% Created Aug 2019; Modified May 2022

%% Clean Environment
close all; clear; clc

%% Specify Directory
database_dir = '/.../.../'; % Specify directory location for Excel testing database on line 21
input_dir = '/.../.../'; % Specify directory location containing trimmed .mat files
output_dir = '/.../.../'; % Specify preferred directory to save .mat file with Fourier coefficients for each test

%% Load Test Database
database_filename = '2019 Test Database.xlsx';
database = [database_dir database_filename];
TestDatabase = readtable(database);

test_num = TestDatabase.TestNumber;
num_tests = numel(TestDatabase.TestNumber);

%% Resample / Trim Data
samp_freq = 125; %Sampling Frequency [Hz]
dt = 1 / samp_freq; % Uniform time step [s]
pad = dt*10; 

% Specify pressure signal trimming time stamps
time_start = TestDatabase.StartTime;
time_end = TestDatabase.EndTime;

for i = 1 : num_tests
    tic
    input_file = [input_dir TestDatabase.File{i}];
    if i == 1
        load (input_file)
        
        % Create time vector
        time_datenum = datenum(strcat(atm{1}, {' '}, atm{2}));
        time = (time_datenum(1:end) - time_datenum(1)) * 86400;
        % Resampled time vector
        time_resamp = [0 : dt : time(end)]'; 
        
        stim_datenum = datenum(strcat(stim{1}, {' '}, stim{2}));
        stim_time = (stim_datenum(1:end) - time_datenum(1)) * 86400;
        
        % Resample stimulation fluid and headspace pressure signals
        atm_resamp = interp1(time, atm{3}, time_resamp);
        stim_resamp = interp1(stim_time, stim{3}, time_resamp);
        
        % Resample observation pressure signals
        for k = 1 : 4
            obs_datenum = datenum(strcat(obs{k}{1}, {' '}, obs{k}{2}));
            obs_time{k} = (obs_datenum(1:end) - time_datenum(1)) * 86400;
            if obs_time{k}(1) > 0 
                obs_time{k}(1) = 0;
            end
            obs_resamp{k} = interp1(obs_time{k}, obs{k}{3}, time_resamp);
        end
        
        % Trimmed time vector
        trim_time = [time_resamp >= time_start(i)-pad &...
            time_resamp <= time_end(i)+pad];
        
        data_time = time_resamp(trim_time);
        test_time = data_time(1:end) - data_time(1);
        
        % Trimmed stimulation head signal - corrected for air pressure
        atm_trim = atm_resamp(trim_time);
        stim_trim = (stim_resamp(trim_time) - atm_trim) ./ 9.804139;
        
        % Trimmed observation head signal
        for j = 1 : 4
            obs_trim(:,j) = obs_resamp{j}(trim_time) ./ 9.804139;
        end
        
        % Save Trimmed Data - Uncomment the following lines to save the
        % output as .mat file
%         TestDatabase.Trim(i) = 'Y';
%         writetable(TestDatabase, database)
%         output_file = ['Test.' num2str(i) '.mat'];
%         save([output_dir output_file], 'test_time', 'stim_trim', 'obs_trim')
        clear atm_trim stim_trim obs_trim data_time test_time
    
    % Continue trimming pressure signals if the filename doesn't change
    elseif strcmp(TestDatabase.File(i), TestDatabase.File(i-1)) == 1
        
        % Trimmed time vector
        trim_time = [time_resamp >= time_start(i)-pad &...
            time_resamp <= time_end(i)+pad];
        
        data_time = time_resamp(trim_time);
        test_time = data_time(1:end) - data_time(1);
        
        % Trimmed stimulation head signal - corrected for air pressure 
        atm_trim = atm_resamp(trim_time);
        stim_trim = stim_resamp(trim_time);
        stim_trim = (stim_trim - atm_trim) ./ 9.804139;
        
        % Trimmed observation head signals
        for j = 1 : 4
            obs_trim(:,j) = obs_resamp{j}(trim_time) ./ 9.804139;
        end
        
        % Save Trimmed Data - Uncomment the following lines to save the
        % output as .mat file
%         TestDatabase.Trim(i) = 'Y';
%         writetable(TestDatabase, database)
%         output_file = ['Test.' num2str(i) '.mat'];
%         save([output_dir output_file], 'test_time', 'stim_trim', 'obs_trim')
        clear atm_trim stim_trim obs_trim data_time test_time
        
    % If filename changes, load the new file
    else
        load(input_file)
        % Create time vector
        time_datenum = datenum(strcat(atm{1}, {' '}, atm{2}));
        time = (time_datenum(1:end) - time_datenum(1)) * 86400;
        % Resampled time vector
        time_resamp = [0 : dt : time(end)]'; 
        
        stim_datenum = datenum(strcat(stim{1}, {' '}, stim{2}));
        stim_time = (stim_datenum(1:end) - time_datenum(1)) * 86400;
        
        % Resample stimulation fluid and air pressure signals
        atm_resamp = interp1(time, atm{3}, time_resamp);
        stim_resamp = interp1(stim_time, stim{3}, time_resamp);
        
        % Resample observation pressure signals
        for k = 1 : 4
            obs_datenum = datenum(strcat(obs{k}{1}, {' '}, obs{k}{2}));
            obs_time{k} = (obs_datenum(1:end) - time_datenum(1)) * 86400;
             if obs_time{k}(1) > 0 
                obs_time{k}(1) = 0;
            end
            obs_resamp{k} = interp1(obs_time{k}, obs{k}{3}, time_resamp);
        end
        
        % Trimmed time signal
        trim_time = [time_resamp >= time_start(i)-pad &...
            time_resamp <= time_end(i)+pad];
        
        data_time = time_resamp(trim_time);
        test_time = data_time(1:end) - data_time(1);
        
        % Trimmed stimulation head signal - corrected for air pressure
        atm_trim = atm_resamp(trim_time);
        stim_trim = stim_resamp(trim_time);
        stim_trim = (stim_trim - atm_trim) ./ 9.804139;
        
        % Trimmed observation head signals
        for j = 1 : 4
            obs_trim(:,j) = obs_resamp{j}(trim_time) ./ 9.804139;
        end
        
        % Save Trimmed Data - Uncomment the following lines to save output
        % as .mat file
%         TestDatabase.Trim(i) = 'Y';
%         writetable(TestDatabase, database)
%         output_file = ['Test.' num2str(i) '.mat'];
%         save([output_dir output_file], 'test_time', 'stim_trim', 'obs_trim')
        clear atm_trim stim_trim obs_trim data_time test_time
        
    end
    toc
end
