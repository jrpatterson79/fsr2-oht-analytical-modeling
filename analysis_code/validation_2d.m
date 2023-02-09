% FSR^2 Oscillatory Flow Testing Analysis
% 2D Numerical Modeling Validation

% This code validates the 2D numercial model used to explore the impacts of
% borehole storage on oscillatory flow signals. The numerical model is
% validated using the fully confined analytical model developed by
% Rasmussen et al. (2003).

% Code developed by Jeremy Patterson
% Created Dec 2018; Updated May 2022

%% Clean Environment
close all; clear; clc

%% Specify Directory
well_dir = '/.../.../'; % Specify directory with well survey file specified on line 36
addpath(genpath('/.../.../Func_Lib/')) % Specify the directory location of the folder Func_Lib, which contains the needed function files to execute this code

%% Forward Model Setup
% Specify domain geometry
xmin = -1500; xmax = -xmin; dx = 2;
ymin = xmin; ymax = xmax; dy = dx;

domain = struct('x',[],'y',[],'z',[]);
domain.x = [xmin : dx : xmax];
domain.y = [ymin : dy : ymax];
domain.z = [0 1];

numx = numel(domain.x) - 1;
numy = numel(domain.y) - 1;
num_cells = numx * numy;

% Grid Domain
[coords, cgrid] = plaid_cellcenter_coord(domain);

% Specify boundary types and boundary values (x / y constant head
% boundaries, no flux in z)
bdry_types = [1; 1; 1; 1; 0; 0];
bdry_vals = zeros(6,1);
bdry_L = zeros(6,1);
bdrys = struct('types',bdry_types,'vals',bdry_vals,'leaks',bdry_L);

%% Load Well Data
% Surveyed Well Locations
wells = readtable([well_dir 'FSR2_WellSurvey.csv']); 
utm_x = wells.utm_x; %UTM-X well locations
utm_y = wells.utm_y; %UTM-Y well locations
alt = wells.Altitude;

% Specify central well as origin in site coordinates
origin_idx = find(strcmp(wells.Well, 'A1') == 1);
origin = [utm_x(origin_idx) utm_y(origin_idx)];

% Translate well coordinates to site coordinates
local_x = utm_x - origin(1);
local_y = utm_y - origin(2);
well_locs = [local_x local_y];

% Calculate radial distance between well pairs
r = sqrt((well_locs(:,1) - well_locs(1,1)).^2 + (well_locs(:,2) - well_locs(1,2)).^2);
num_wells = numel(r);

%% Create Test List
P = logspace(1, 2.6990, 10); % Stimulation periods (s)
Q_max = 7e-5;                % Peak flow rate (m^3/s)

test_list = [];
for i = 1 : numel(P)
    for j = 2 : num_wells
    test_list = [test_list; ...
                 (2*pi)/P(i) 1 Q_max j; ...      
                 ];
    end
end

%% Calculate Fracture Flow Properties
% Water and aquifer properties
rho = 998.23;           % Fluid density (kg/m^3)
g = 9.81;               % Acceleration due to gravity (m/s^2)
mu = 1.0016e-3;         % Fluid dynamic viscosity (Pa s)
RockComp = 1.41e-8;     % Rock compressibility (1 / Pa)
WaterComp = 1 / 2.1e9;  % Fluid compressibility (1 / Pa)
eta_r = 0.15;           % Host rock porosity (-)

% Fracture Aperture Properties
ln_aper = log(3e-4);  % Fracture aperture (m)
lnT = log((rho*g)/(12*mu) * exp(ln_aper).^3);  % Fracture transmissivity (m^2/s)
Ss_rock = rho * g * (RockComp + (eta_r * WaterComp)); % Host rock specific storage  (1/s) 
Ss_frac = rho * g * WaterComp; % Fracture specific storage (1/s)

V_frac = Ss_frac * dx * dy * 0.4;
V_rock = Ss_rock * dx * dy * 0.6;
Ss_eff = (V_frac + V_rock) / (dx*dy);
lnS = log(Ss_eff * exp(ln_aper)); % Fracture storativity (-)

%% Create the "official" input files needed by the steady periodic model.
%This step creates the input files needed to run all of the steady periodic models.
[inputs] = OHT_create_inputs(well_locs,test_list,domain);
y_fxn = @(params) OHT_run_distribKSs(params, domain, bdrys, inputs, 1);

num_omegas = size(inputs,1);
num_obs = size(test_list,1);

y_oht = y_fxn([lnT*ones(num_cells,1); lnS*ones(num_cells,1)]);

A_syn = y_oht(1:num_obs);
B_syn = y_oht(num_obs+1:2*num_obs);

amp_oht = sqrt(A_syn.^2 + B_syn.^2);
phase_oht = atan2(-B_syn, A_syn);

%% Fully confined analytical solution (Rasmussen et al. (2003))
soln = 'confined';

syn_data(:,1) = (2*pi) ./ test_list(:,1); % Stimulation period (s)
syn_data(:,2) = test_list(:,1);           % Angular Frequency (rad/s)
syn_data(:,3) = test_list(:,3);           % Max Pumping Rate (m^3/s)
syn_data(:,4) = r(test_list(:,4));        % Inter-well spacing (m)

for k = 1 : numel(syn_data(:,1))
    y_ana(k,:) = RasSoln(syn_data(k,:), [lnT; lnS], soln);
end

amp_ana = sqrt(y_ana(:,1).^2 + y_ana(:,2).^2);
phase_ana = atan2(-y_ana(:,2), y_ana(:,1));

%% Calculate Numerical approximation error
amp_rmse = sqrt(mean((amp_oht - amp_ana).^2));
phase_rmse = sqrt(mean((phase_oht - phase_ana).^2));

%% Figures
figure
clf
subplot(1,2,1)
ax = gca;
hold on
for j = 2 : numel(r)
    idx = find(syn_data(:,4) == r(j));
    plot(P, amp_oht(idx), 'k-', 'LineWidth', 2)
    plot(P, amp_ana(idx), 'k^', 'MarkerSize', 8)
end
xlim([0 P(end)])
xlabel('Period (s)')
ylabel('Amplitude (m)')
ax.XScale = 'log';
ax.FontSize = 18;

subplot(1,2,2)
ax = gca;
hold on
for j = 2 : numel(r)
    idx = find(syn_data(:,4) == r(j));
    plot(P, phase_oht(idx), 'k-', 'LineWidth', 2)
    plot(P, phase_ana(idx), 'k^', 'MarkerSize', 8)
end
xlim([0 P(end)]);
xlabel('Period (s)')
ylabel('Phase (rad)')
ax.XScale = 'log';
ax.FontSize = 18;
set(gcf, 'Position', [0 0 2025 2025/2.75])