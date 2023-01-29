% FSR^2 Oscillatory Flow Testing Analysis

% Synthetic modeling experiment to test borehole storage as a potential mechanism
% contributing to period dependent fracture flow parameters under
% oscillatory flow conditions. The code uses OHT3D
% (https://github.com/wischydro-cardiff/oscillatory-tomography) to generate
% synthetic data.

% This code reproduces Figure 12 as seen in:
% Patterson, Jeremy R. & Cardiff, Michael (2023). Do Simple Analytical
% Models Capture Complex Fractured Bedrock Hydraulics? Oscillatory Flow
% Tests Suggest Not. Groundwater. https://doi.org/

% Code developed by Jeremy Patterson
% Created Dec 2020; Updated May 2022
%% Clean Environment
clear; close all; clc

%% Specify Directory
% Specify the directory location of the folder Func_Lib, which contains the needed function files to execute this code
addpath('/.../.../') 

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

% Specify well locations
well_locs = [0 0;...
            -1 -3];

% Calculate radial distances between wells
r = sqrt((well_locs(:,1) - well_locs(1,1)).^2 + (well_locs(:,2) - well_locs(1,2)).^2);
num_wells = numel(well_locs(:,1));

%% Create Test List
P = logspace(1, 2.6990, 10); % Stimulation periods [s]
Q_max = 7e-5;                % Peak flow rate [m^3/s]

test_list = [];
for i = 1:1:numel(P)
    for j = 2 : num_wells
        test_list = [test_list; ...
            (2*pi)/P(i) 1 Q_max j];
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

Ss_rock = rho * g * (RockComp + (eta_r * WaterComp)); 
Ss_frac = rho * g * WaterComp;

V_frac = Ss_frac * dx * dy * 0.4;
V_rock = Ss_rock * dx * dy * 0.6;
Ss_eff = (V_frac+V_rock) / (dx*dy);

% Fracture Aperture Properties
ln_aper = log(3e-4);                          % Fracture aperture [m]
lnT = log((rho*g)/(12*mu) * exp(ln_aper).^3); % Fracture transmissivity [m]
lnS = log(Ss_eff * exp(ln_aper));             % Fracture storativity [-]

lnT_true = lnT * ones(num_cells,1);
lnS_true = lnS * ones(num_cells,1);

V_well = pi * r_pipe^2 * 1;                  % Water volume change in well per unit head change
V_aq = exp(lnS) * ((dx*dy) - (pi*r_pipe^2)); % Water volume change in aquifer per unit head change
S_eff = log((V_well+V_aq) / (dx*dy));        % Effective storativity 

% Find observation well location to assign effective storativity value.
for b = 2 : numel(well_locs(:,1))
    idx(b-1) = find(well_locs(b,1) == coords(:,1) & well_locs(b,2) == coords(:,2));
end
lnS_true(idx) = S_eff;

%% Create the "official" input files needed by the steady periodic model.
%This step creates the input files needed to run all of the steady periodic models. 
[inputs] = OHT_create_inputs(well_locs,test_list,domain);
y_fxn = @(params) OHT_run_distribKSs(params, domain, bdrys, inputs, 1);

num_omegas = size(inputs,1);
numobs = size(test_list,1);

% Generate the true phasor data
y_oht = y_fxn([lnT_true; lnS_true]);
y_obs = [y_oht(1:numobs) y_oht(numobs+1:2*numobs)];

%% Homogeneous Confined Inversion
soln = 'confined';

syn_data(:,1) = (2*pi) ./ test_list(:,1); %Stimulation period [s] 
syn_data(:,2) = test_list(:,1);           %Angular Frequency [rad/s]
syn_data(:,3) = test_list(:,3);           %Max Pumping Rate [m^3/s]
syn_data(:,4) = r(test_list(:,4));        %Radial Distance [m]

s_init = [lnT; mean(lnS_true)]; % Initial parameter guesses
lambda = 1e1;                   % Initial LM stabilization parameter
delta = [0.1; 0.1];             % Small parameter increment to numerically approximate Jacobian
R_inv = 1e-30 * eye(2);         % Data error covariance matrix
soln = 'confined';              % Analytical forward model for inversion

for i = 1 : numobs
    % Optimal parameters [T S]
   [s_opt(i,:), ~] = Lev_Marq(syn_data(i,:), s_init, y_obs(i,:)', R_inv, lambda, delta, soln);
   % Modeled phasors using optimal parameters
    y_mod(i,:) = RasSoln(syn_data(i,:), s_opt(i,:)', soln);
end
s_opt(:,3) = s_opt(:,1) - s_opt(:,2); % Hydraulic diffusivity

%% Figure 12
figure(12)
clf
subplot(2,3,[1:2,4:5])
ax = gca;
plot(syn_data(:,1), exp(s_opt(:,3)), 'ko', 'MarkerSize', 12,...
    'MarkerFaceColor', [0 0.4470 0.7410])
hold off
grid on
axis([1e1 5e2 1e-1 1e1])
ax.XScale = 'log';
ax.YScale = 'log';
xlabel('Oscillation Period (s)')
ylabel('D_{est} (m^2/s)')
ax.FontSize = 30;
text(400, 8, 'A', 'FontSize', 40, 'FontWeight', 'bold')

subplot(2,3,3)
ax = gca;
hold on
plot(syn_data(:,1), exp(s_opt(:,1)), 'ko', 'MarkerSize', 12,...
    'MarkerFaceColor', [0 0.4470 0.7410])

axis([1e1 5e2 1e-6 1e-4])
hold off
grid on
ylabel('T_{est} (m^2/s)')
ax.XScale = 'log';
ax.YScale = 'log';
ax.YAxisLocation = 'right';
text(350, 7e-5, 'B', 'FontSize', 40, 'FontWeight', 'bold')
ax.FontSize = 30;

subplot(2,3,6)
ax = gca;
hold on
plot(syn_data(:,1), exp(s_opt(:,2)), 'ko', 'MarkerSize', 12,...
    'MarkerFaceColor', [0 0.4470 0.7410])
axis([1e1 5e2 9e-6 1e-4])
ax.XScale = 'log';
ax.YScale = 'log';
hold off
grid on
xlabel('Oscillation Period (s)')    
ylabel('S_{est} (-)')
ax.YAxisLocation = 'right';
ax.FontSize = 30;
text(350, 8e-5, 'C', 'FontSize', 40, 'FontWeight', 'bold')
set(gcf, 'Position', [0 0 2025 2025/1.8])
