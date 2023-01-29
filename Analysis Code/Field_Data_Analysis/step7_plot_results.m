% FSR^2 Oscillatory Flow Testing
% Step 7 - Plot characterization results

% This code loads the single frequency parameter estimation results and
% reproduces Figures 4, 5, 8, and 9 as seen in:
% Patterson, Jeremy R. & Cardiff, Michael (2023). Do Simple Analytical
% Models Capture Complex Fractured Bedrock Hydraulics? Oscillatory Flow
% Tests Suggest Not. Groundwater. https://doi.org/

% Code developed by Jeremy Patterson
% Created Dec 2019; Updated May 2022
%% Clean Environment
close all; clear; clc

%% Specify Directory
file_dir = '/.../.../'; % Add directory location to 2019PESTResults.mat file with parameter estimation results indicated on line 87
addpath(genpath('/./Func_Lib/')) % Add directory location to Func_Lib which contains all the needed function files to run this analysis

%% Figure 4 Analysis
% Fourire
%Number of components allowed in the represented signal
numpeaks_req = 2;
%Relative power of a peak (relative to largest peak) required
relpower_req = 0.001;
%Average power of a peak required relative to the power of it's immediate
%neighbors
neighborrelpower_req = 2;
%FFT will use 2^FFT_incpower multiplier when creating list of frequencies /
%periods
FFT_incpower = 7;

input_dir ='/.../.../'; % Specify directory location to filename specified on line 29
filename = 'Test.100.mat';
file = [input_dir filename];
load(file)

dt = test_time(2) - test_time(1);
obs_detrend = detrend(obs_trim(:,2), 'linear');

C = FFT_incpower; %Signal extension factor
Fs = 1./dt; % sampling frequency
L = length(obs_detrend);

NFFT = 2^(nextpow2(L)+C);
Y = fft(obs_detrend,NFFT)/L;

f = Fs/2*linspace(0,1,NFFT/2+1); % frequency list
P = 1./f; %Period list

coeffs = 2*Y(1:NFFT/2+1);
power = abs(coeffs);
max_power = max(power);

[locs_pks_good,pkspow_good] = findpeaks_nolobes(power(end:-1:2), ...
    P(end:-1:2),relpower_req,neighborrelpower_req);

%Filter down to only the number of peaks desired
num_locs = numel(locs_pks_good);
if num_locs < numpeaks_req
    numpeaks_req = num_locs;
end

locs_pks_good = locs_pks_good(1:numpeaks_req);
pkspow_good = pkspow_good(1:numpeaks_req);

obs_mod = zeros(numel(test_time),1);
for k = 1 : numpeaks_req
    ind = find(P == locs_pks_good(k));
    obs_mod = obs_mod + ...
        real(coeffs(ind)) .* cos(((2*pi) / P(ind)) .* test_time) +...
       -imag(coeffs(ind)) .* sin(((2*pi) / P(ind)) .* test_time);
end

% Figure 4
figure(4)
clf
ax = gca;
plot(test_time, obs_mod, 'r-', 'LineWidth', 2)
hold on
plot(test_time, obs_recon, 'k-', 'LineWidth', 3)
axis([test_time(1) test_time(end) -1e-2 1e-2])
xlabel('Time (s)')
ylabel('Head Change (m)')
ax.FontSize = 34;
set(gcf, 'Position', [100 100 2025 2025/3])

%% Load Field Data
load ([file_dir '2019PESTResults.mat']) % Provide the directory location to this file on line 15 
idx = find(s_unc_19(:,1) < 1 & test_list(:,1) > 4);

%% Figure 5 - Plot Parameter Estimates Full
fit_idx = [33; 49; 61; 65];

figure(5)
clf
subplot(2,3,[1:2,4:5])
ax = gca;
scatter(test_list(idx,1), exp(s_opt_19(idx,3)), 125, test_list(idx,4),...
        'filled', 'MarkerEdgeColor', 'k')
hold on
plot(test_list(fit_idx,1), exp(s_opt_19(fit_idx,3)), 's', 'Color', [0.7529 0 0], 'LineWidth', 3,...
     'MarkerSize', 18)
hold on 
axis([10 310 9e-2 3e1])
ax.XScale = 'log';
ax.YScale = 'log';
hold off
grid on
xlabel('Oscillation Period (s)')
ylabel('D_{est} (m^2/s)')
ax.FontSize = 30;
c = colorbar;
caxis([4 16])
c.Label.String = 'Radial Distance (m)';
c.FontSize = 24;
text(275,25, 'A', 'FontSize', 40, 'FontWeight', 'bold')

subplot(2,3,3)
ax = gca;
scatter(test_list(idx,1), exp(s_opt_19(idx,1)), 125, test_list(idx,4),...
        'filled', 'MarkerEdgeColor', 'k')
hold on
plot(test_list(fit_idx,1), exp(s_opt_19(fit_idx,1)), 's', 'Color', [0.7529 0 0], 'LineWidth', 3,...
     'MarkerSize', 18)
ax.YTick = [1e-6 1e-5 1e-4];
axis([10 310 7e-7 3e-4])
ax.XScale = 'log';
ax.YScale = 'log';
hold off
grid on
ylabel('T_{est} (m^2/s)')
ax.YAxisLocation = 'right';
ax.FontSize = 30;
text(225, 2e-4, 'B', 'FontSize', 40, 'FontWeight' ,'bold')

subplot(2,3,6)
ax = gca;
scatter(test_list(idx,1), exp(s_opt_19(idx,2)), 125, test_list(idx,4),...
        'filled', 'MarkerEdgeColor', 'k')
hold on
plot(test_list(fit_idx,1), exp(s_opt_19(fit_idx,2)), 's', 'Color', [0.7529 0 0], 'LineWidth', 3,...
     'MarkerSize', 18)
ax.YTick = [1e-6 1e-5 1e-4];
axis([10 310 7e-7 3e-4])
ax.XScale = 'log';
ax.YScale = 'log';
hold off
grid on
xlabel('Oscillation Period (s)')    
ylabel('S_{est} (-)')
ax.YAxisLocation = 'right';
ax.FontSize = 30;
text(225, 2e-4, 'C', 'FontSize', 40, 'FontWeight', 'bold')
set(gcf, 'Position', [100 100 2025 2025/1.8])

%% Darcy Flow Analysis
% Figure 8 - Reynolds number analysis
p_idx = find(test_list(:,1) > 4);
% Fluid properties
rho = 998.23;      % Fluid density (kg / m^3)
mu = 1.0016e-3;    % Fluid dynamic viscosity (Pa s)
g = 9.81;         
r_pipe =  0.0254;  % Acceleration du to gravity (m/s^2)

% Hydraulic aperture and Reynolds Number analysis
v = test_list(:,3) ./ (2 * pi * r_pipe); % Fluid velocity (m/s)
lnb = (s_opt_19(:,1) - log((rho * g)/(12 * mu))) ./ 3;
lnK = s_opt_19(:,1) - lnb;
Re = (rho .* v .* 2 .* exp(lnb)) ./ mu;

% Figure 8 - Reynolds Numbers
figure(8)
clf
subplot(1,2,1)
ax = gca;
plot(test_list(p_idx,3), exp(lnb(p_idx))*1e3, 'kd', 'MarkerSize', 12,...
     'MarkerFaceColor', [0.8500 0.3250 0.0980]) 
grid on
ylim([0 0.7])
xlabel('Q_{max} (m^3/s)')
ylabel('Aperture_{est} (mm)')
ax.FontSize = 30;
text(1.85e-4, 0.67, 'A','FontSize', 40, 'FontWeight', 'bold')

subplot(1,2,2)
ax = gca;
plot(test_list(p_idx,3), Re(p_idx), 'ks', 'MarkerSize', 12,...
     'MarkerFaceColor', [0.9290 0.6940 0.1250])
grid on
xlabel('Q_{max} (m^3/s)')
ylabel('Reynolds Number (-)')
text(1.85e-4, 0.67, 'B', 'FontSize', 40, 'FontWeight', 'bold')
ax.FontSize = 30;
set(gcf, 'Position', [100 100 2025 2025/2.75])


% Figure 9 - Flow Linearity Analysis
col = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250],...
       [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], [0.6350 0.0780 0.1840]};
mkr = {'v', 'd', '^', 's', '>', 'o'};

% Flow linearity test
q_idx = {find(test_idx >= 62 & test_idx <= 64),...
         find(test_idx >= 34 & test_idx <= 37),...
         find(test_idx >= 53 & test_idx <= 55)};
num_q = numel(q_idx);

figure(9)
clf
for j = 1 : num_q
    P_title = test_list(q_idx{j}(1),1);
    r_list = unique(test_list(q_idx{j},4));
    r_idx = []; q = []; delta_h = []; r_check = [];

    subplot(1, num_q, j)
    ax = gca;
    hold on
    for i = 1 : numel(r_list)
        r_idx = find(test_list(q_idx{j},4) == r_list(i));
        q(:,i) = test_list(q_idx{j}(r_idx),3);
        delta_h(:,i) = sqrt(y_obs(q_idx{j}(r_idx),1).^2 + y_obs(q_idx{j}(r_idx),2).^2);
        r_check(:,i) = test_list(q_idx{j}(r_idx),4);

        A = [0; q(:,i)];
        m = A \ [0; delta_h(:,i)];
        dh_mod = A * m;

        txt = [num2str(round(r_list(i))) ' m'];
        plot(q(:,i), delta_h(:,i), mkr{i}, 'MarkerSize', 12, 'MarkerFaceColor', col{i},...
            'MarkerEdgeColor', 'k', 'DisplayName', txt)
        plot(A, dh_mod, 'Color', col{i}, 'LineWidth', 2, 'HandleVisibility', 'off')
    end
    grid on
    xlabel('Q_{max} (m^3/s)')
    ylabel('|\Delta h| (m)')
    ax.FontSize = 30;
    l = legend;
    l.Location = 'NorthWest';
    l.FontSize = 24;
    title(['P = ' num2str(P_title) 's'], 'FontSize', 28)
    if j == 1
        axis([0 5e-5 0 0.03])
        text(3.85e-5, 0.029, 'A','FontSize', 40, 'FontWeight', 'bold')
    elseif j == 2
        axis([0 1.2e-4 0 1.4e-1])
        ax.XTick = [0:4e-5:1.2e-4];
        text(1e-4, 0.136, 'B','FontSize', 40, 'FontWeight', 'bold')
    elseif j == 3
        axis([0 9e-5 0 0.15])
        text(7.7e-5, 0.146, 'C','FontSize', 40, 'FontWeight', 'bold')
    end
end
set(gcf, 'Position', [100 100 2025 2025/3.1])