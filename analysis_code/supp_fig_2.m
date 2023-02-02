%% Clean Environment
close all; clear; clc

%% Specify Directory
addpath(genpath('/.../.../Func_Lib/')) % Add directory location of the folder Func_Lib which contains necessary function files to run the parameter estimation

%% Specify Model Run Type
% Choose homogeneous ('homog') or heterogeneous ('heterog') model run
model = 'heterog';
% Specify geostatistical model ('exponential', 'gaussian', 'linear')
geostat_model = 'exponential';

% Set seed for random number generator
seed = 20;
randn('state', seed)

%% Load Well Data
% Well locations
well_locs = [0 0;...
    0 -4;...
    -6 0;...
    0 8;...
    10 0];
num_wells = numel(well_locs(:,1));
%% Forward Model Setup
% Model domain
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
%% Random field parameters

% Fracture geostatistical properties
ln_aper_mean = log(3e-4);          % Mean aperture (m)
if strcmp(model, 'heterog') == 1
    var_b = 0.3;          % Aperture variance
elseif strcmp(model, 'homog') == 1
    var_b = 0;            % Aperture variance
end
L_x = 20; L_y = L_x/4;  % Correlation lengths (m)

% Fluid and Rock Properties
rho = 998.23;
g = 9.81;
mu = 1.0016e-3;
RockComp = 1.41e-8;
WaterComp = 1 / 2.1e9;
eta_r = 0.15;

Ss_rock = rho * g * (RockComp + (eta_r * WaterComp));
Ss_frac = rho * g * WaterComp;

V_frac = Ss_frac * dx(1) * dy(1) * 0.4;
V_rock = Ss_rock * dx(1) * dy(1) * 0.6;
Ss_eff = (V_frac + V_rock) / (dx(1)*dy(1));

% Mean fracture flow parameters
lnT_mean = log((rho*g*exp(ln_aper_mean).^3) ./ (12*mu));
lnS_mean = log(Ss_eff .* exp(ln_aper_mean));
lnD_mean = lnT_mean - lnS_mean;

%% Create Test List
P = logspace(1, 3, 20);  % Stimulation periods [s]
Q_max = 7e-5;

%Test list for OHT input
test_list = [];
r = [];
for i = 1:numel(P)
    for j = 1 : num_wells
        for k = j+1:num_wells
            r  = [r;
                sqrt((well_locs(j,1)-well_locs(k,1)).^2 + (well_locs(j,2)-well_locs(k,2)).^2)];
            test_list = [test_list; ...
                (2*pi)/P(i) j Q_max k];
        end
    end
end

%% Generate Data / Conduct Inversion
% Synthetic Data for parameter estimation
syn_data(:,1) = (2*pi) ./ test_list(:,1);
syn_data(:,2) = test_list(:,1);
syn_data(:,3) = test_list(:,3);
syn_data(:,4) = r;

% Generate data covariance matrix
data_err = sqrt(1e-20);
R_inv = inv(data_err*eye(2));

% Initiate parameter estimation
soln = 'confined';
s_init = [lnT_mean; lnS_mean];
delta = [0.1; 0.1];
lambda = 1e1;
max_iter = 100;

% Generate parameter covariance row
distmat_row = dimdist(coords(1,:),coords);
max_dist = max(max((distmat_row(:,:,1).^2 + ...
    distmat_row(:,:,2).^2).^.5));

if strcmp(geostat_model, 'exponential') == 1
    corr_row = exp(-(...
        (distmat_row(:,:,1)./L_x).^2 + ...
        (distmat_row(:,:,2)./L_y).^2).^.5);
    [corr_relz] = toepmat_vector_math(corr_row,'r',[],2,[numy numx]);
    ln_aper_true = corr_relz(:,1).*var_b.^.5 + ln_aper_mean;

elseif strcmp(geostat_model, 'gaussian') == 1
    corr_row = exp(-(...
        (distmat_row(:,:,1)./L_x).^2 + ...
        (distmat_row(:,:,2)./L_y).^2));
    [corr_relz] = toepmat_vector_math(corr_row,'r',[],2,[numy numx]);
    ln_aper_true = corr_relz(:,1).*var_b.^.5 + ln_aper_mean;

elseif strcmp(geostat_model, 'linear') == 1
    corr_row = (1/max_dist) .*...
        (max_dist - ((distmat_row(:,:,1).^2 + ...
        distmat_row(:,:,2).^2).^.5));
    [corr_relz] = toepmat_vector_math(corr_row,'r',[],2,[numy numx]);
    ln_aper_true = corr_relz(:,1).*(1.8622*var_b).^.5 + ln_aper_mean;

else
    error('Pick a valid variogram model')
end
p = pcolor(cgrid{1}, cgrid{2}, reshape(ln_aper_true, numy, numx));
p.LineStyle = 'none';
hold on
plot(well_locs(:,1), well_locs(:,2), 'ko', 'MarkerFaceColor', 'r', 'MarkerSize', 12)
axis([-30 30 -30 30])
pause 

% Convert aperture to flow parameters
lnT_true = log((rho*g*exp(ln_aper_true).^3) ./ (12*mu));
lnS_true = log(Ss_eff .* exp(ln_aper_true));

% Generate OHT inputs
[inputs] = OHT_create_inputs(well_locs,test_list,domain);
num_omegas = size(inputs,1);
num_obs = size(test_list,1);

% Generate phasor data
y_fxn = @(params) OHT_run_distribKSs(params, domain, bdrys, inputs, 1);
y_obs = y_fxn([lnT_true; lnS_true]);

for i = 1 : num_obs
    [pest, ~, out_flag] = Lev_Marq(syn_data(i,:), s_init, [y_obs(i); y_obs(num_obs+i)], R_inv, lambda, delta, soln, max_iter);
    J = jacob(pest, delta, syn_data(i,:), soln);
    y_opt(i,:) = RasSoln(syn_data(i,:), pest, soln);
    s_hat(i,:) = pest;
end

%% 
% Figure S2 - Aperture realization and parameter estimates
figure
clf
subplot(2,2,1)
ax1 = gca;
p = pcolor(cgrid{1}, cgrid{2}, reshape(log10(exp(ln_aper_true)), numy, numx));
p.LineStyle = 'none';
hold on
plot(well_locs(:,1), well_locs(:,2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', [0.7529 0 0])
axis([-40 40 -40 40])
axis square
xlabel('X distance (m)')
ylabel('Y distance (m)')
ax1.FontSize = 30;
c = colorbar;
caxis([-4.5 -2.5])
c.Label.String = 'log_{10}(Aperture [m])';
c.FontSize = 24;
text(35, 37, 'A', 'FontSize', 30, 'FontWeight', 'bold')

subplot(2,2,2)
ax2 = gca;
scatter(syn_data(:,1), exp(s_hat(:,1)-s_hat(:,2)), 125, syn_data(:,4),...
        's', 'filled', 'MarkerEdgeColor', 'k')
grid on
axis([P(1) P(end) 1e1 1e4])
ax2.XScale = 'log';
ax2.YScale = 'log';
xlabel('Oscillation Period (s)')
ylabel('D_{est} (m^2/s)')
ax2.FontSize = 30;
c = colorbar;
c.Label.String = 'Inter-well Spacing (m)';
c.FontSize = 24;
caxis([4 16])
text(750, 6.5e3, 'B', 'FontSize', 30, 'FontWeight', 'bold')

subplot(2,2,3)
ax3 = gca;
scatter(syn_data(:,1), exp(s_hat(:,1)), 125, syn_data(:,4),...
        's', 'filled', 'MarkerEdgeColor', 'k')
hold on
plot([P(1) P(end)], [exp(lnT_mean) exp(lnT_mean)], 'k--', 'LineWidth', 3)
grid on
axis([P(1) P(end) 1e-6 1e-3])
ax3.XScale = 'log';
ax3.YScale = 'log';
xlabel('Oscillation Period (s)')
ylabel('T_{est} (m^2/s)')
ax3.FontSize = 30;
ax2.FontSize = 30;
c = colorbar;
c.Label.String = 'Inter-well Spacing (m)';
c.FontSize = 24;
caxis([4 16])
text(750, 6.5e-4, 'C', 'FontSize', 30,'FontWeight', 'bold')

subplot(2,2,4)
ax4 = gca;
scatter(syn_data(:,1), exp(s_hat(:,2)), 125, syn_data(:,4),...
        's', 'filled', 'MarkerEdgeColor', 'k')
hold on
plot([P(1) P(end)], [exp(lnS_mean) exp(lnS_mean)], 'k--', 'LineWidth', 3)
grid on
axis([1e1 1e3 1e-9 1e-6])
ax4.XScale = 'log';
ax4.YScale = 'log';
xlabel('Oscillation Period (s)')
ylabel('S_{est} (-)')
ax4.FontSize = 30;
ax2.FontSize = 30;
c = colorbar;
c.Label.String = 'Inter-well Spacing (m)';
c.FontSize = 24;
caxis([4 16])
text(800, 6.5e-7, 'D', 'FontSize', 30,'FontWeight', 'bold')
set(gcf, 'Position', [100 100 2025/1.1 2025])