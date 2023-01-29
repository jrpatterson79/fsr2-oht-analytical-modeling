% FSR^2 Oscillatory Flow Testing Analysis
% Step 4 - Extract Fourier Coefficients

% This code represents the third data processing step associated with field
% experiments and analysis described in:
% Patterson, Jeremy R. & Cardiff, Michael (2023). Do Simple Analytical
% Models Capture Complex Fractured Bedrock Hydraulics? Oscillatory Flow
% Tests Suggest Not. Groundwater. https://doi.org/

% This file loads individual trimmed .mat files and uses the linear
% regression approach from Bakhos et al. (2014) to extract the Fourier
% coefficients for the stimulation pressure signal, observation pumping
% signal, and the pumping signals. The output for each test is saved as a
% separate .mat file for further analysis / parameter estimation.

% Code created by Jeremy Patterson
% Created: Dec 2019; Modified: May 2022 

%% Clean Environment
close all; clear; clc

%% Specify Directory
database_dir = '/.../.../'; % Specify directory location for Excel file '2019 Test Database.xlsx'
file_dir = '/.../.../'; % Specify directory location containing trimmed .mat files
save_dir = '/.../.../'; % Specify preferred directory to save .mat file with extracted phasor coefficients for each test

%% Load Test Database
database_filename = '2019 Test Database.xlsx';
database = [database_dir database_filename];
TestDatabase = readtable(database);
num_tests = numel(TestDatabase.TestNumber);

%% Extract Fourier Coefficients

for i = 1 : num_tests % Set i = 11 to reproduce Figure 3
    disp(['Test #' num2str(i)])

    input_file = ['Test.' num2str(i) '.mat'];
    file = [file_dir input_file];
    load(file)
    
    P_in = TestDatabase.Period_s_(i);
    omega = (2*pi) / P_in;
    L = length(stim_trim);
    stim_detrend = []; obs_detrend = [];
    
    % Detrend stimulation signal and extract Fourier coefficients
    stim_detrend = detrend(stim_trim, 'linear');
    [~, stim_phasor, ~] = periodic_LS_fit(test_time, stim_detrend, P_in);
    stim_recon = real(stim_phasor) .* cos(omega .* test_time) +...
                -imag(stim_phasor) .* sin(omega .* test_time);
    
    % Riser pipe geometry - Used to calculate flow rate from pressure signals
    r_pipe =  0.0254;          % Riser pipe radius [m]
    pipe_area = pi * r_pipe^2; % Riser pipe cross-sectional area [m^2]
    % Generate pumping signal (Eqn 1 in manuscript) and extract Fourier coefficients
    Q = (-omega .*...
        (real(stim_phasor) .* sin(omega .* test_time) +...
         imag(stim_phasor) .* cos(omega .* test_time))) .* pipe_area;
    [~, Q_phasor, ~] = periodic_LS_fit(test_time, Q, P_in);
    
    % Detrend observation signals and extract Fourier coefficients
    num_obs = numel(obs_trim(1,:));
    obs_recon = zeros(numel(test_time),num_obs);
    for j = 1 : num_obs
        obs_detrend(:,j) = detrend(obs_trim(:,j), 'linear');
        [data_cov{j}, obs_phasors(j,:), data_err(j,:)] = periodic_LS_fit(test_time, obs_detrend(:,j), P_in);
        obs_recon(:,j) = real(obs_phasors(j,:)) .* cos(omega .* test_time) +...
                        -imag(obs_phasors(j,:)) .* sin(omega .* test_time);
    end
    
    if strcmp(TestDatabase.StimulationWell(i), 'A1') == 1
        lbl = {'A1', 'B1', 'B2', 'B3', 'B4'};
    elseif strcmp(TestDatabase.StimulationWell(i), 'B1') == 1
        lbl = {'B1', 'A1', 'B2', 'B3', 'B4'};
    elseif strcmp(TestDatabase.StimulationWell(i), 'B2') == 1
        lbl = {'B2', 'B1', 'A1', 'B3', 'B4'};
    elseif strcmp(TestDatabase.StimulationWell(i), 'B3') == 1
        lbl = {'B3', 'B1', 'B2', 'A1', 'B4'};
    elseif strcmp(TestDatabase.StimulationWell(i), 'B4') == 1
        lbl = {'B4', 'B1', 'B2', 'B3', 'A1'};
    end

    % Uncomment lines 81-115 to plot Figure 3 
    % Uncomment lines 121-177 to plot measured and reconstructed signals
    % for each test as a quality control check on extracted Fourier
    % coefficients
%     figure(3)
%     clf
%     subplot(2,1,1)
%     ax = gca;
%     yyaxis left
%     plot(test_time, stim_detrend, '-', 'LineWidth', 4)
%     grid on
%     ax.GridAlpha = 0.4;
%     ax.GridColor = 'k';
%     xlim([test_time(1) test_time(end)])
%     ylabel('Head change (m)')
% 
%     yyaxis right
%     plot(test_time, Q, '-', 'LineWidth', 4)
%     xlim([test_time(1) test_time(end)])
%     xlabel('Test time (s)')
%     ylabel('Q (m^3/s)')
%     ax.FontSize = 30;
%     text(220, 2e-4, 'A', 'FontSize', 40, 'FontWeight', 'bold')
% 
%     subplot(2,1,2)
%     lbls = {'B1', 'B2','B3', 'B4'};
%     ax = gca;
%     hold on
%     plot(test_time, obs_detrend, '-','LineWidth', 3)
%     grid on
%     xlim([test_time(1) test_time(end)])
%     xlabel('Test time (s)')
%     ylabel('Head change (m)')
%     ax.FontSize = 30;
%     l = legend(lbls);
%     l.Location = 'NorthWest';
%     l.FontSize = 26;
%     text(220, 0.1, 'B', 'FontSize', 40, 'FontWeight', 'bold')
%     set(gcf, 'Position', [100 100 2025 2025/2.25])

%     Uncomment for 5x1 subplot that includes stimulation and observation
%     signals for each test
%     figure
%     clf
%     for idx = 1 : 5
%         subplot(5,1,idx)
%         ax = gca;
%         if idx == 1
%             plot(test_time, stim_detrend, 'k-', 'LineWidth', 3)
%             hold on
%             plot(test_time, stim_recon, 'r-', 'LineWidth', 3)
%             xlim([0 test_time(end)])
%             xlabel('Time (s)')
%             ylabel('\Delta h (m)')
%             title(lbl{idx})
%             ax.FontSize = 20;
%             l = legend('Observed', 'Reconstructed');
%             l.FontSize = 18;
%         if idx == 2
%             subplot(5,1,idx)
%             ax = gca;
%             plot(test_time, obs_detrend(:,idx), 'r', 'Color', [0.6350 0.0780 0.1840], 'LineWidth',3)% 'MarkerSize', 10)
%             hold on
%             plot(test_time,  obs_recon(:,idx), 'k-', 'LineWidth', 3)
%             xlim([0 test_time(end)])
%             ylim([-1e-2 1e-2])
%             ylabel('\Delta h (m)')
%             ax.FontSize = 28;
%             title(lbl{idx+1}, 'FontSize', 24, 'Color', [0.6350 0.0780 0.1840])
%             
%         elseif idx == 4
%             subplot(5,1,idx)
%             ax = gca;
%             plot(test_time,  obs_recon(:,idx), 'k-', 'LineWidth', 3)
%             hold on
%             plot(test_time, obs_detrend(:,idx), '-', 'Color', [0.6350 0.0780 0.1840], 'LineWidth',3)% 'MarkerSize', 10)
%             xlim([0 test_time(end)])
%             ylim([-2e-1 2e-1])
%             xlabel('Time (s)', 'FontSize', 28)
%             ylabel('\Delta h (m)', 'FontSize', 28)
%             ax.FontSize = 28;
%             title(lbl{idx+1}, 'FontSize', 24, 'Color', [0.6350 0.0780 0.1840])
%             
%         else
%             subplot(5,1,idx)
%             ax = gca;
%             plot(test_time,  obs_recon(:,idx), 'k-', 'LineWidth', 3)
%             hold on
%             plot(test_time, obs_detrend(:,idx), '-', 'Color', [0.6350 0.0780 0.1840], 'LineWidth',3)% 'MarkerSize', 10)
%             xlim([0 test_time(end)])
%             ylim([-2e-1 2e-1])
%             ylabel('\Delta h (m)')
%             ax.FontSize = 28;
%             title(lbl{idx+1}, 'FontSize', 24, 'Color', [0.6350 0.0780 0.1840])
%             set(gcf, 'Position', [100 100 975 975/1.2667])          
%         end
%     end 
%     
    
%     save_filename = ['Test_' num2str(i) '_Fourier.mat'];
%     save_file = [save_dir save_filename];
%     save(save_file, 'stim_phasor', 'Q_phasor', 'obs_phasors', 'data_err', 'data_cov', 'test_time')
end

