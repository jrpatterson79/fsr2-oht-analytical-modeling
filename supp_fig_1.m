% This code generates Figures seen in supplementary information for:
% Patterson, Jeremy R. & Cardiff, Michael (2023). Do Simple Analytical
% Models Capture Complex Fractured Bedrock Hydraulics? Oscillatory Flow
% Tests Suggest Not. Groundwater. https://doi.org/

% Code by Jeremy Patterson
% Created Dec. 2019; Modified May 2022

%% Clean Environment
close all; clear; clc

%% Specify Directory
pub_dir = '/Users/jpatt/My Drive/PhD/Drafts/Patterson_Cardiff_Groundwater_2022/Figures/';
supp_dir = '/Users/jpatt/My Drive/PhD/Drafts/Patterson_Cardiff_Groundwater_2023/Supplemental_Materials/';
file_dir = '/Users/jpatt/My Drive/PhD/Groundwater_FieldStudies/MAT_Files/';

%% Load Data
load ([file_dir '2019PESTResults.mat'])

%% Flow Anisotropy Analysis
a1_idx = [find(test_idx >= 7 & test_idx <= 19 & round(test_list(:,4)) == 10);...
          find(test_idx >= 59 & test_idx <= 74 & round(test_list(:,4)) == 10)];

b4_idx = [find(test_idx >= 20 & test_idx <= 42 & round(test_list(:,4)) == 10);...
          find(test_idx >= 75 & test_idx <= 92 & round(test_list(:,4)) == 10)];

% Figure S1      
figure
clf
subplot(2,3,[1:2,4:5])
ax = gca;
plot(test_list(a1_idx,1), exp(s_opt_19(a1_idx,3)), 'ko', 'MarkerSize', 12,...
     'MarkerFaceColor', [0 0.4470 0.7410])
hold on
plot(test_list(b4_idx,1), exp(s_opt_19(b4_idx,3)), 'kd', 'MarkerSize', 12,...
     'MarkerFaceColor', [0.9290 0.6940 0.1250])
ax.XScale = 'log';
ax.YScale = 'log';
grid on
axis([10 500 0.1 10])
xlabel('Oscillation Period (s)')
ylabel('D_{est} (m^2/s)')
ax.FontSize = 30;
text(400, 8, 'A', 'FontSize', 40, 'FontWeight', 'bold')
l = legend('A1-B2', 'B1-B4');
l.FontSize = 24;
l.Location = 'NorthWest';

subplot(2,3,3)
ax = gca;
plot(test_list(a1_idx,1), exp(s_opt_19(a1_idx,1)), 'ko', 'MarkerSize', 12,...
     'MarkerFaceColor', [0 0.4470 0.7410])
hold on
plot(test_list(b4_idx,1), exp(s_opt_19(b4_idx,1)), 'kd', 'MarkerSize', 12,...
     'MarkerFaceColor', [0.9290 0.6940 0.1250])
ax.XScale = 'log';
ax.YScale = 'log';
grid on
axis([10 500 1e-6 1e-4])
ylabel('T_{est} (m^2/s)')
ax.YAxisLocation = 'right';
ax.FontSize = 30;
text(325, 7e-5, 'B', 'FontSize', 40, 'FontWeight', 'bold')

subplot(2,3,6)
ax = gca;
plot(test_list(a1_idx,1), exp(s_opt_19(a1_idx,2)), 'ko', 'MarkerSize', 12,...
     'MarkerFaceColor', [0 0.4470 0.7410])
hold on
plot(test_list(b4_idx,1), exp(s_opt_19(b4_idx,2)), 'kd', 'MarkerSize', 12,...
     'MarkerFaceColor', [0.9290 0.6940 0.1250])
ax.XScale = 'log';
ax.YScale = 'log';
grid on
axis([10 500 1e-6 1e-4])
xlabel('Oscillation Period (s)')    
ylabel('S_{est} (-)')
ax.YAxisLocation = 'right';
ax.FontSize = 30;
text(325, 7e-5, 'C', 'FontSize', 40, 'FontWeight', 'bold')
set(gcf, 'Position', [100 100 2025 2025/1.8])
print([supp_dir 'Figure_S1'], '-dpng', '-r300')