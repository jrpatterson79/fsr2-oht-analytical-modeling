% FSR^2 Oscillatory Flow Testing
% Step 5 - Generate Test List

% This code represents the an intermediate data processing step associated with field
% experiments and analysis described in:
% Patterson, Jeremy R. & Cardiff, Michael (2023). Do Simple Analytical
% Models Capture Complex Fractured Bedrock Hydraulics? Oscillatory Flow
% Tests Suggest Not. Groundwater. https://doi.org/

% This script opens the .mat files with extracted Fourier coefficents for all
% oscillatory flow tests conducted in the isolated fracture at 105' depth
% only. The loaded data is saved as a test list in matrix format and saved
% as a .mat file used for follow-on single and multi-frequency parameter estimation. 

% Code developed by Jeremy Patterson
% Created Dec 2019; Updated May 2022

%% Clean Environment
close all; clear; clc

%% Specify Directory
database_dir = '/.../.../'; % Specify directory location for Excel testing database on line 21
file_dir = '/.../.../'; % Specify directory location containing .mat files with extracted phasor coefficients 
save_dir = '/.../.../'; % Specify preferred directory to save test list as .mat file 
well_dir = '/.../.../'; % Specify directory with well survey file specified on line 36

%% Load Data
% Load testing database 
database_filename = '2019 Test Database.xlsx';
database = [database_dir database_filename];
TestDatabase = readtable(database);

InputAmp = TestDatabase.Amplitude_kPa_;
kept_test = [3:92]'; % Tests in fracture at 105' bls
num_tests = numel(kept_test);

%% Model Geometry

% Borehole Geometry
r_pipe = 0.0254;           % Riser pipe radius
pipe_area = pi * r_pipe^2; % Riser pipe cross-sectional area

% Surveyed Well Locations
well_locs = readtable([well_dir 'FSR2_WellSurvey.csv']); 

% Well locations in UTM coordinates with elevation to top of casing
utm_x = well_locs.utm_x;   
utm_y = well_locs.utm_y;   
elev = well_locs.Altitude; 

% Specify well A1 as origin
origin_idx = find(strcmp(well_locs.Well, 'A1') == 1);
origin = [utm_x(origin_idx) utm_y(origin_idx)];

% Translate UTM to site coordinates
local_x = utm_x - origin(1);
local_y = utm_y - origin(2);
depth = [elev(1) - TestDatabase.A1_Depth*0.3048...
         elev(2) - TestDatabase.B1_Depth*0.3048...
         elev(3) - TestDatabase.B2_Depth*0.3048...
         elev(4) - TestDatabase.B3_Depth*0.3048...
         elev(5) - TestDatabase.B4_Depth*0.3048];

%% Generate Testlist
test_idx = [];
test_list = [];
cov_y = [];
y_obs = [];

start = tic;
for i = 1 : num_tests
    % Load test .mat file
    disp(['Test Number ' num2str(TestDatabase.TestNumber(kept_test(i)))])
    input_file = ['Test_' num2str(kept_test(i)) '_Fourier.mat'];
    file = [file_dir input_file];
    load(file)
    
    P = TestDatabase.Period_s_(kept_test(i)); % Oscillation period (s)
    omega = (2 * pi) / P;                     % Angular frequency (1/s)
    Q_max = abs(Q_phasor);                    % Peak pumping rate (m^3/s)
    
    % Determine the stimulation well coordinates
    stim = find(strcmp(well_locs.Well, TestDatabase.StimulationWell(kept_test(i))) == 1);
    stim_loc = [utm_x(stim) utm_y(stim)];
    
    % Determine observation well coordinates and calculate radial distances
    if strcmp(TestDatabase.StimulationWell(kept_test(i)), 'A1') == 1
        obs = [find(strcmp(well_locs.Well, 'B1') == 1);...
            find(strcmp(well_locs.Well, 'B2') == 1);...
            find(strcmp(well_locs.Well, 'B3') == 1);...
            find(strcmp(well_locs.Well, 'B4') == 1)
            ];
        
        obs_loc = [utm_x(obs) utm_y(obs)];
        r = sqrt((stim_loc(1) - obs_loc(:,1)).^2 + (stim_loc(2) - obs_loc(:,2)).^2 + (depth(kept_test(i),stim) - depth(kept_test(i),obs)).^2');
        
    % Determine observation well coordinates and calculate radial distances
    elseif strcmp(TestDatabase.StimulationWell(kept_test(i)), 'B1') == 1
        obs = [find(strcmp(well_locs.Well, 'A1') == 1);...
            find(strcmp(well_locs.Well, 'B2') == 1);...
            find(strcmp(well_locs.Well, 'B3') == 1);...
            find(strcmp(well_locs.Well, 'B4') == 1)
            ];
        
        obs_loc = [utm_x(obs) utm_y(obs)];
        r = sqrt((stim_loc(1) - obs_loc(:,1)).^2 + (stim_loc(2) - obs_loc(:,2)).^2 + (depth(kept_test(i),stim) - depth(kept_test(i),obs)).^2');
        
    % Determine observation well coordinates and calculate radial distances
    elseif strcmp(TestDatabase.StimulationWell(kept_test(i)), 'B2') == 1
        obs = [find(strcmp(well_locs.Well, 'B1') == 1);...
            find(strcmp(well_locs.Well, 'A1') == 1);...
            find(strcmp(well_locs.Well, 'B3') == 1);...
            find(strcmp(well_locs.Well, 'B4') == 1)
            ];
        
        obs_loc = [utm_x(obs) utm_y(obs)];
        r = sqrt((stim_loc(1) - obs_loc(:,1)).^2 + (stim_loc(2) - obs_loc(:,2)).^2 + (depth(kept_test(i),stim) - depth(kept_test(i),obs)).^2');
        
    % Determine observation well coordinates and calculate radial distances
    elseif strcmp(TestDatabase.StimulationWell(kept_test(i)), 'B3') == 1
        obs = [find(strcmp(well_locs.Well, 'B1') == 1);...
            find(strcmp(well_locs.Well, 'B2') == 1);...
            find(strcmp(well_locs.Well, 'A1') == 1);...
            find(strcmp(well_locs.Well, 'B4') == 1)
            ];
        
        obs_loc = [utm_x(obs) utm_y(obs)];
        r = sqrt((stim_loc(1) - obs_loc(:,1)).^2 + (stim_loc(2) - obs_loc(:,2)).^2 + (depth(kept_test(i),stim) - depth(kept_test(i),obs)).^2');
        
    % Determine observation well coordinates and calculate radial distances
    elseif strcmp(TestDatabase.StimulationWell(kept_test(i)), 'B4') == 1
        obs = [find(strcmp(well_locs.Well, 'B1') == 1);...
            find(strcmp(well_locs.Well, 'B2') == 1);...
            find(strcmp(well_locs.Well, 'B3') == 1);...
            find(strcmp(well_locs.Well, 'A1') == 1)
            ];
        
        obs_loc = [utm_x(obs) utm_y(obs)];
        r = sqrt((stim_loc(1) - obs_loc(:,1)).^2 + (stim_loc(2) - obs_loc(:,2)).^2 + (depth(kept_test(i),stim) - depth(kept_test(i),obs)).^2');
        
    end
    % Stimulation signal amplitude [m]
    stim_amp = abs(stim_phasor);
    for j = 1 : numel(r)
        
        obs_amp = abs(obs_phasors(j,:)); % Observation signal amplitude [m]
        dh = stim_amp - obs_amp;         % Head amplitude difference between stimulation and observation signals [m]
        
    test_list = [test_list;
                 P omega Q_max r(j) data_err(j,:) dh];
    y_obs = [y_obs; real(obs_phasors(j,:)) imag(obs_phasors(j,:))]; % Observation signal phasor
    cov_y = [cov_y; reshape(data_cov{j}, 1, [])];                   % Data error covariance matrix
    test_idx = [test_idx; kept_test(i)];                            % Test number from test database
    end
end

% Uncomment the following lines to save the test list as .mat file
% save_file = 'TestList_19.mat';
% save([save_dir save_file])
