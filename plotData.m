%% Authors: Sven Schellenberger, Kilin Shi
% Copyright (C) 2020  Sven Schellenberger, Kilin Shi

%%%%%% IMPORTANT %%%%%%
% Copy the data from figshare (link to be generated) to a subfolder "datasets" inside the folder of
% this script or change datadir to according dataset location

% Set the subject-ID(s) and scneraio(s) you want to view in the
% 'Options'-section and start the script
%%%%%%%%%%%%%%%%%%%%%%%

%% Plot data
% This script generates exemplary plots of selected subject-ID(s) and
% scenario(s)
% Data is loaded, processed and then plotted

%% Init
close all;
clc; 
clear ;

addpath(genpath('utils'))

%% Options
%%%% Choose Subject-ID(s) 
%IDrange = 16; % GDN00XX - simply choose numbers or ranges from 01 to 30

IDrange = 7;
%%%% Choose scnerio(s) 
% possible scenarios are {'Resting' 'Valsalva' 'Apnea' 'TiltUp' 'TiltDown'}

scenarios = {'Resting'};

%scenarios = {'Resting' 'Valsalva' 'Apnea' 'TiltUp' 'TiltDown'};

%%%% Set path to datasets 
% Datasets can be found on figshare
path = 'project_data'; % In this case "datasets" folder is in this scripts folder 

scrsz = get(groot,'ScreenSize'); % For plotting
%% Extract and plot the data

% Manually selected ECG channel for each subject (Position in ECG_CHANNEL vector corresponds to ID)
ECG_CHANNEL = [2 2 2 2 2 1 2 2 2 2 2 2 2 2 1 2 2 2 2 2 1 1 2 2 2 2 2 2 2 2];

for indx = 1:length(IDrange)
    %% Iterate all subject IDs
    ID = sprintf('GDN%04d',IDrange(indx));
    fprintf('----------- Loading %s ------------\n', ID);

    for sz = 1:length(scenarios)
        %% Iterate all existing subject scenarios
        scenario = scenarios{sz};
        fprintf('---- Scenario %s\n', scenario);
        
        % search file
        path_id = [path,'\',ID];
        files_synced_mat = dir([path_id,'\*.mat']);
        found = [];
        for j = 1:length(files_synced_mat)
            found = strfind(files_synced_mat(j).name,scenario);
            if ~isempty(found)
                % load file
                load([path_id,'\',files_synced_mat(j).name]);
                break;
            end
        end
        if isempty(found)
            fprintf('---- skipped\n');
            continue
        end
        
        %% Prepare data
        output.(ID).(scenario) = struct;

        % old start
        % [~,~,~,radar_dist] = elreko(radar_i,radar_q,measurement_info{1},1); % Ellipse fitting, compensation and distance reconstruction
        % old end 
        % Usage of elreko
        % [radar_i_compensated,radar_q_compensated,phase_compensated,radar_dist] = elreko(radar_i,radar_q,measurement_info{1}(timestamp of dataset),0(Plot flag -> 1: plots on, 0: plots off));
        
    
        [radar_i_compensated,radar_q_compensated,phase_compensated,radar_dist] = elreko(radar_i,radar_q,measurement_info{1},1);


        %%
        % Ellipse reconstruction 


        [radar_i_fitted,radar_q_fitted] = ellipseCreator(radar_i,radar_q,1); % 1 - plot the ellipses, 0 - no plotting

        
        %signalplotter(ThreeBandFilter(phase_compensated,fs_radar))
        %heartRateSignal, FsHeart, lowBandSignal, FsLow, highBandSignal, FsHigh =
        %continue
        %%  30 sec WINDOW 
        window_length_s = 30; % 30 seconds as per the paper 
        samples_per_window = round(window_length_s * fs_radar);
        overlap = 0;          % You can add overlap (e.g., 50%) for smoother results, but starting with 0 is simpler.

       
        
        
        

%%  distance calculating && plotting : 
        lambda = 0.0125 ; % 0.0125 in meters, 12.5 in mm 

        radar_distance_from_phase_compensated = phase_compensated .* lambda ./(4*pi) ;
        %phase reversed by distance
        phase_from_radar_dist =  radar_dist.*(4*pi)./lambda;

        %phase generated after compensation by us:
        phase_calc_after_fit = atan( radar_q_fitted./radar_i_fitted) ; 
        radar_distance_calculated_compensated = phase_calc_after_fit .* lambda ./(4*pi) ;


        %cutoff = abs(max(fft(phase_from_radar_dist)))/fs_radar;
        
        %phase_simple_calc_filtered = butterworth_low_pass_filter(phase_simple_calc,6,cutoff,fs_radar);
        %phase_compensated_filtered = butterworth_low_pass_filter(phase_compensated,6,cutoff,fs_radar); 

        %radar_distance_by_me_filtered = phase_simple_calc_filtered .* lambda ./(4*pi) ;
  
    
        %% Bandpass filter
        f_cutoff_low = 0.05;
        f_cutoff_high = 0.5;
        order=4 ;
        d = designfilt('bandpassiir', ...
    'FilterOrder', order, ...
    'HalfPowerFrequency1', f_cutoff_low, ...
    'HalfPowerFrequency2', f_cutoff_high, ...
    'SampleRate', fs_radar, ...
    'DesignMethod', 'butter');

    radar_dist_after_BPF = filtfilt(d, radar_distance_calculated_compensated);




        %% LPF 
        f_cut = 0.5;            
        Wn_low = f_cut / (fs_radar/2); 
        [b_low, a_low] = butter(order, Wn_low, 'low');     
        %%

        %radar_dist_after_BPF = filtfilt(b_band, a_band, radar_distance_calculated_compensated);
        radar_dist_after_LPF = filtfilt(b_low, a_low, radar_distance_calculated_compensated);
       % radar_dist_from_phase_compensated_filtered = ( filtfilt(b_resp, a_resp,radar_distance_from_phase_compensated ) ) .*1000 ;
       % radar_dist_calc_mm = radar_distance_calculated_compensated.*1000 ; %unfiltered 
       %radar_dist_resp_filtered_mm = radar_dist_resp_meters .* 1000;
       %num_zero_crossings = sum(abs(diff(sign(radar_dist_after_BPF))) > 0);
       % A full respiratory cycle (one breath) corresponds to 2 zero-crossings.
       %num_breaths = num_zero_crossings / 2.0;          




       [radar_respiration, radar_pulse, radar_heartsound, tfm_respiration] = getVitalSigns(radar_dist, fs_radar, tfm_z0, fs_z0);
        
        % TFM
        if ECG_CHANNEL(IDrange(indx)) == 1
            tfm_ecg = fillmissing(tfm_ecg1,'constant',0); % Sometimes ECG is NaN -> set all occurrences to 0
        else
            tfm_ecg = fillmissing(tfm_ecg2,'constant',0); % Sometimes ECG is NaN -> set all occurrences to 0
        end
        tfm_ecg = filtButter(tfm_ecg,fs_ecg,4,[1 20],'bandpass');
        
        % Resample radar respiration to match tfm respiration
        radar_respiration_re = resample(radar_respiration,fs_z0,fs_radar);
        radar_dist_re = resample(radar_dist,fs_z0,fs_radar);
        if length(radar_respiration_re) > length(tfm_respiration)
            radar_respiration_re = radar_respiration_re(1:length(tfm_respiration));
            radar_dist_re = radar_dist_re(1:length(tfm_respiration));
        elseif length(radar_respiration_re) < length(tfm_respiration)
            tfm_respiration = tfm_respiration(1:length(radar_respiration_re));
        end
        sc_sec = getInterventions(scenario,tfm_intervention,fs_intervention);
        
        %% Evaluate data
        
        % Heartbeat detection
        radar_hs_states = getHsmmStates(radar_heartsound, fs_radar);
        radar_hs_locsR = statesToLocs(radar_hs_states, 1);
        [~,tfm_ecg_locsR] = twaveend(tfm_ecg(1:length(radar_hs_states)), fs_ecg,32*(fs_ecg/200),'p');
        
        
        %% Time vectors
        time_radar = 1/fs_radar:1/fs_radar:length(radar_dist)/fs_radar;
        time_ecg = 1/fs_ecg:1/fs_ecg:length(tfm_ecg)/fs_ecg;
        time_icg = 1/fs_icg:1/fs_icg:length(tfm_icg)/fs_icg;
        time_bp = 1/fs_bp:1/fs_bp:length(tfm_bp)/fs_bp;
        time_z0 = 1/fs_z0:1/fs_z0:length(tfm_z0)/fs_z0;
        time_respiration = 1/fs_z0:1/fs_z0:length(radar_respiration_re)/fs_z0;

        %% resampling
        %BPF
        % Resample radar respiration to match tfm respiration
        respiration_test_BPF = resample(radar_dist_after_BPF,fs_z0,fs_radar);
        radar_dist_re_BPF = resample(radar_distance_calculated_compensated,fs_z0,fs_radar);
        if length(respiration_test_BPF) > length(tfm_respiration)
            respiration_test_BPF = respiration_test_BPF(1:length(tfm_respiration));
            radar_dist_re_BPF = radar_dist_re_BPF(1:length(tfm_respiration));
        elseif length(respiration_test_BPF) < length(tfm_respiration)
            tfm_respiration = tfm_respiration(1:length(respiration_test_BPF));
        end

        %{
        T_sec = time_respiration(end) - time_respiration(1);
        T_min = T_sec / 60.0;

        if T_min > 0
          respiration_rate_BPM = num_breaths / T_min;
       else
          respiration_rate_BPM = 0.0;
       end
        %}

        %%
%% --- NEW SEGMENTED RESPIRATION RATE CALCULATION (Zero Crossing) ---
% This implements the paper's method of splitting signals into 30s windows.
window_length_s = 30; 
% --- Radar Rate Calculation ---
samples_per_window_radar = round(window_length_s * fs_radar);
overlap_samples = 0; 

radar_rates_BPM = [];
radar_rate_time_points = [];
total_samples_radar = length(radar_dist_after_BPF);
start_sample = 1;

while start_sample + samples_per_window_radar - 1 <= total_samples_radar
    end_sample = start_sample + samples_per_window_radar - 1;
    current_window = radar_dist_after_BPF(start_sample:end_sample);
    num_zero_crossings = sum(abs(diff(sign(current_window))) > 0);
    rate_BPM = (num_zero_crossings / 2.0) / (window_length_s / 60.0);
    
    radar_rates_BPM = [radar_rates_BPM, rate_BPM];
    time_s = (start_sample + end_sample) / 2 / fs_radar;
    radar_rate_time_points = [radar_rate_time_points, time_s];
    
    start_sample = start_sample + samples_per_window_radar - overlap_samples;
end

% Calculate overall mean for summary (Table 1 comparison)
if ~isempty(radar_rates_BPM)
     respiration_rate_BPM = mean(radar_rates_BPM);
else
     respiration_rate_BPM = 0.0;
end

% --- TFM Rate Calculation (Applying the same logic) ---
samples_per_window_tfm = round(window_length_s * fs_z0);

tfm_rates_BPM = [];
tfm_rate_time_points = [];
total_samples_tfm = length(tfm_respiration); % Use the TFM signal
start_sample = 1;

while start_sample + samples_per_window_tfm - 1 <= total_samples_tfm
    end_sample = start_sample + samples_per_window_tfm - 1;
    current_window = tfm_respiration(start_sample:end_sample);
    
    % The TFM signal (Impedance) is filtered, so ZC is applicable here as well.
    num_zero_crossings = sum(abs(diff(sign(current_window))) > 0);
    rate_BPM = (num_zero_crossings / 2.0) / (window_length_s / 60.0);
    
    tfm_rates_BPM = [tfm_rates_BPM, rate_BPM];
    time_s = (start_sample + end_sample) / 2 / fs_z0;
    tfm_rate_time_points = [tfm_rate_time_points, time_s];
    
    start_sample = start_sample + samples_per_window_tfm - overlap_samples;
end

%% --- NEW PLOT: Radar vs. TFM Aggregated Respiration Rate ---

figure('Position',[1 1 scrsz(3) scrsz(4)-80]);
sgtitle(['Aggregated Respiration Rate Comparison (30s Windows) - Subject ' ID ' - ' scenario]);

ax_agg = subplot(1,1,1);
hold on;

% 1. Plot Radar Respiration Rate
plot(radar_rate_time_points, radar_rates_BPM, 'b-o', 'DisplayName', 'Radar Rate (ZC, 30s)');

% 2. Plot TFM Respiration Rate
plot(tfm_rate_time_points, tfm_rates_BPM, 'r-x', 'DisplayName', 'TFM Rate (ZC, 30s Reference)');

title('Radar vs. TFM Respiration Rate over Time (30s Segments)');
ylabel('Respiration Rate (Breaths Per Minute - BrPM)');
xlabel('Time (s)');
legend('show');
grid on;
if exist('time_radar','var')
    xlim([time_radar(1) time_radar(end)]);
end



        %%
        %%RMSE CALC
        % RMSE Calculation Between Two Vectors
        rmse = sqrt(mean((radar_rates_BPM - tfm_rates_BPM).^2));
        fprintf('RMSE = %.4f\n', rmse);

        






        %LPF
        % Resample radar respiration to match tfm respiration
        respiration_test_LPF = resample(radar_dist_after_LPF,fs_z0,fs_radar);
        radar_dist_re_LPF = resample(radar_distance_calculated_compensated,fs_z0,fs_radar);
        if length(respiration_test_LPF) > length(tfm_respiration)
            respiration_test_LPF = respiration_test_LPF(1:length(tfm_respiration));
            radar_dist_re_LPF = radar_dist_re_LPF(1:length(tfm_respiration));
        elseif length(respiration_test_LPF) < length(tfm_respiration)
            tfm_respiration = tfm_respiration(1:length(respiration_test_LPF));
        end


        %%
        %respiration plotting
  
        figure('Position',[1 1 scrsz(3) scrsz(4)-80]);
        ax2(1) = subplot(1,1,1);
        hold on;

        plot(time_respiration, radar_respiration_re.*1000, 'k-', 'DisplayName', 'Radar Respiration');
        plot(time_respiration, tfm_respiration, 'r-', 'DisplayName', 'TFM Respiration');
        plot(time_respiration, respiration_test_BPF.*1000, 'b-', 'DisplayName', 'test_BPF Respiration');
        plot(time_respiration, respiration_test_LPF.*1000, 'g-', 'DisplayName', 'test_LPF Respiration');
        
        hold off;
        title('Compare respiration');
        ylabel('Rel. Distance(mm)');
        xlabel('Time(s)');
        legend('show');
        grid on;

        %% Plot all data
        figure('Position',[1 1 scrsz(3) scrsz(4)-80]);
        ax(1) = subplot(4,1,1);
        plot(time_radar,radar_dist.*1000,'k-');
        title('Radar')
        ylabel('Rel. Distance(mm)');
        % xlabel('Time(s)')

        ax(2) = subplot(4,1,2);
        plot(time_z0,tfm_z0,'k-')
        title('Impedance')
        ylabel('Voltage (mV)');
        % xlabel('Time(s)')

        ax(3) = subplot(4,1,3);
        plot(time_ecg,tfm_ecg1+4,'k-');
        hold on;
        plot(time_ecg,tfm_ecg2+2,'b-');
        plot(time_icg,tfm_icg,'r-')
        legend('Lead 1','Lead 2', 'ICG')
        title('Electrocardiogram & Impedancecardiogram')
        ylabel('norm. Amplitude');
        % xlabel('Time(s)')

        ax(4) = subplot(4,1,4);
        plot(time_bp,tfm_bp,'k-')
        title('Blood pressure')
        ylabel('norm. Amplitude');
        xlabel('Time(s)')
 
        linkaxes(ax,'x');
        xlim([time_radar(1) time_radar(end)]);

        %% Compare radar vital signs to reference raw signals
        figure('Position',[1 1 scrsz(3) scrsz(4)-80]);

        ax2(1) = subplot(2,1,1);
        plot(time_respiration,radar_respiration_re.*1000,'k-');
        title('Compare respiration')
        ylabel('Rel. Distance(mm)');
        yyaxis right;
        plot(time_respiration,tfm_respiration);
        ylabel('Impedance')
        xlabel('Time(s)')

        ax2(2) = subplot(2,1,2);
        plot(time_radar,radar_heartsound.*10^6,'k-')
        hold on;
        set(vline(radar_hs_locsR(1)/fs_radar,'b-'),'HandleVisibility','on'); % Turn the legend on for vline
        vline(radar_hs_locsR./fs_radar,'b-');
        set(vline(tfm_ecg_locsR(1)/fs_ecg,'r-'),'HandleVisibility','on'); % Turn the legend on for vline
        vline(tfm_ecg_locsR./fs_radar,'r-');
        title('Compare heartbeat detection')
        ylabel('Rel. Distance(um)');
        xlabel('Time(s)')
        legend('Radar heart sound','Radar S1','ECG R-peak')
        legend('Radar heart sound');
        
        linkaxes(ax2,'x');
        xlim([time_radar(1) time_radar(end)]);
        
        figure('Position',[1 1 scrsz(3) scrsz(4)-80]);
        ax(1) = subplot(5,1,1);
        plot(time_radar,radar_i_compensated,'b');
        hold on;
        plot(time_radar,radar_i_fitted,'r');
        legend('compensated','regular');
        title('radar_i compensated and regular');
        ylabel('Amplitude');
        % xlabel('Time(s)')

        ax(2) = subplot(5,1,2);
        plot(time_radar,radar_q_compensated,'b');
        hold on;
        plot(time_radar,radar_q_fitted,'r');
        legend('compensated','regular');
        title('radar_q compensated and regular');
        ylabel('Amplitude');
        % xlabel('Time(s)')

        ax(3) = subplot(5,1,3);
        plot(time_ecg,tfm_ecg1+4,'k-');
        hold on;
        plot(time_ecg,tfm_ecg2+2,'b-');
        plot(time_icg,tfm_icg,'r-')
        legend('Lead 1','Lead 2', 'ICG')
        title('Electrocardiogram & Impedancecardiogram')
        ylabel('norm. Amplitude');
        % xlabel('Time(s)')

        ax(4) = subplot(5,1,4);
        plot(time_radar,phase_compensated,'k-')
        title('phase compensated')
        ylabel('radians');
        xlabel('Time(s)');
        linkaxes(ax,'x');
        xlim([time_radar(1) time_radar(end)]);





       %{     
        %% DISTANCES :
        
        figure('Position',[1 1 scrsz(3) scrsz(4)-80]);
        ax(1) = subplot(5,1,1);
        plot(time_radar,radar_dist.*1000,'k-');
        title('Radar distance from the paper');
        ylabel('Rel. Distance(mm)');
        % xlabel('Time(s)')

     %   ax(2) = subplot(5,1,2);
     %   plot(time_radar,radar_dist_resp_filtered_mm,'b');
     %   title('Radar distance compensation and filtered by me');
     %   ylabel('Rel. Distance(mm)');
        % xlabel('Time(s)')

        ax(3) = subplot(5,1,3);
        plot(time_radar,radar_dist_calc_mm,'r');
        title('distance compensation unfiltered by me');
        ylabel('Rel. Distance(mm)');
        xlabel('Time(s)')

        ax(4) = subplot(5,1,4);
        plot(time_radar,radar_dist_from_phase_compensated_filtered,'b');
        title('Radar distance filtered by me(from phase compensated');
        ylabel('Rel. Distance(mm)');
        xlabel('Time(s)')

        ax(5) = subplot(5,1,5);
        plot(time_radar,radar_distance_from_phase_compensated.*1000,'r');
        title('Radar distance unfiltered by me(from phase compensated');
        ylabel('Rel. Distance(mm)');



        %% PHASE: 


        figure('Position',[1 1 scrsz(3) scrsz(4)-80]);
        ax(1) = subplot(3,1,1);
        plot(time_radar,phase_from_radar_dist,'k-')
        title('phase from inverse calc');
        ylabel('radians');
        xlabel('Time(s)');

        ax(2) = subplot(3,1,2);
        plot(time_radar,phase_calc_after_fit,'b')
        title('phase compensated by calculation')
        ylabel('radians');
        xlabel('Time(s)')

        ax(3) = subplot(3,1,3);
        plot(time_radar,phase_compensated,'r')
        title('phase compensated')
        ylabel('radians');
        xlabel('Time(s)');


      %  ax(4) = subplot(5,1,4);
      %  plot(time_radar,phase_simple_calc_filtered,'r');
      %  title('phase regular after butterworth LPF');
      %  ylabel('radians');
      %  xlabel('Time(s)');

       % ax(5) = subplot(5,1,5);
      %  plot(time_radar,fft(phase_simple_calc),'r')
      %  title('FREQUENCY RESPONSE OF PHASE');
       % ylabel('Amplitude');
      %  xlabel('omega');
       %}


%%plotting radar_distance and heartrate       
        figure('Name', 'Correlation Visualizations');
        sgtitle(' radar distance heartrate');
        subplot(1, 1, 1);
        plot(time_radar,radar_dist_after_LPF.*1000,'r')
        hold on
        plot(time_radar,radar_heartsound.*10^6,'b')
        title(['radar_ distance and heartsound ']);
        xlabel('radar dist');
        ylabel('radar heartsound');
        grid on;
        legend('radar_dist_after_LPF' ,'radar_heartsound' );



        %{
        %%%%%%%%%%% correlation between phase and ECG leads/ICG %%%%%%%%%%%%

            we want to check:
              ECG VS  PHASE
              ECG1 VS PHASE
              ECG2 VS PHASE
              ICG  VS PHASE 
       
        disp([' fs_icg :', num2str(fs_icg)]);
        disp([' fs_ecg :', num2str(fs_ecg)]);
        disp([' fs_radar :', num2str(fs_radar)]);
        %%ecg vs phase
        
        pearson_correlation_corr_phase_ecg = corr(phase_from_radar_dist, tfm_ecg); 
        R_ecg = corrcoef(phase_from_radar_dist, tfm_ecg);
        pearson_correlation_R_ecg = R_ecg(1, 2);

        
        %%ecg1 vs phase
        pearson_correlation_corr_phase_ecg1 = corr(phase_from_radar_dist, tfm_ecg1); 
        R_ecg1 = corrcoef(phase_from_radar_dist, tfm_ecg1);
        pearson_correlation_R_ecg1 = R_ecg1(1, 2);
       

        %%ecg2 vs phase
        pearson_correlation_corr_phase_ecg2 = corr(phase_from_radar_dist, tfm_ecg2); 
        R_ecg2 = corrcoef(phase_from_radar_dist, tfm_ecg2);
        pearson_correlation_R_ecg2 = R_ecg2(1, 2);

        %%icg vs phase
        N_phase= length(phase_from_radar_dist);
        N__icg = length(tfm_icg);

        x_new= linspace(1,N_phase,N__icg);
        x_target = (1:N_phase)';
        tfm_icg_interpolated = interp1(x_new, tfm_icg, x_target, 'linear');
        R_interpolated = corr(phase_from_radar_dist, tfm_icg_interpolated);
        R_icg = corrcoef(phase_from_radar_dist, tfm_icg_interpolated);
        pearson_correlation_R_icg = R_icg(1, 2);

%displays:

        disp(['Size of phase: ', num2str(length(phase_from_radar_dist))]);
        disp(['Size of new_array_interpolated: ', num2str(length(tfm_icg_interpolated))]);
        disp([' ecg total vs phase :  Correlation using corr: ', num2str(pearson_correlation_corr_phase_ecg)]);
        disp([' ecg_lead1 vs phase :  Correlation using corr: ', num2str(pearson_correlation_corr_phase_ecg1)]);
        disp([' ecg_lead2 vs phase :  Correlation using corr: ', num2str(pearson_correlation_corr_phase_ecg2)]);
        disp(['Correlation of ICG vs phase: R = ', num2str(R_interpolated)]);


        disp([' ecg total vs phase :  Correlation using corr: ', num2str(pearson_correlation_R_ecg)]);
        disp([' ecg_lead1 vs phase :  Correlation using corr: ', num2str(pearson_correlation_R_ecg1)]);
        disp([' ecg_lead2 vs phase :  Correlation using corr: ', num2str(pearson_correlation_R_ecg2)]);
        disp(['Correlation of ICG vs phase: R = ', num2str(pearson_correlation_R_icg)]);


        %%plotting        
        figure('Name', 'Correlation Visualizations');
        sgtitle('Correlation of Phase with ECG Arrays');
        
        % --- Plot 1: Phase vs. ECG1 ---
        subplot(1, 4, 1);
        scatter(phase_from_radar_dist, tfm_ecg, 5, 'filled'); % 5 is the marker size
        title(['Phase vs. ECG1 (R = ', num2str(pearson_correlation_corr_phase_ecg, '%.3f'), ')']);
        xlabel('Phase');
        ylabel('ECG');
        grid on;
        
        % --- Plot 2: Phase vs. ECG2 ---
        subplot(1, 4, 2);
        scatter(phase_from_radar_dist, tfm_ecg1, 5, 'filled');
        title(['Phase vs. ECG2 (R = ', num2str(pearson_correlation_corr_phase_ecg1, '%.3f'), ')']);
        xlabel('Phase');
        ylabel('ECG lead 1');
        grid on;
        
        % --- Plot 3: Phase vs. ECG3 ---
        subplot(1, 4, 3);
        scatter(phase_from_radar_dist, tfm_ecg2, 5, 'filled');
        title(['Phase vs. ECG3 (R = ', num2str(pearson_correlation_corr_phase_ecg2, '%.3f'), ')']);
        xlabel('Phase');
        ylabel('ECG lead 2');
        grid on;

                % --- Plot 3: Phase vs. ICG ---
        subplot(1, 4, 4);
        scatter(phase_from_radar_dist, tfm_icg_interpolated, 5, 'filled');
        title(['Phase vs. ICG (R = ', num2str(R_interpolated, '%.3f'), ')']);
        xlabel('Phase');
        ylabel('ICG');
        grid on;

 %} 
        if IDrange(indx) > 4
        
            %% Compare Radar to TFM aggregated data
            figure('Position',[1 1 scrsz(3) scrsz(4)-80]);
            ax3(1) = subplot(2,1,1);
            plot(time_radar,radar_dist.*1000);
            title('Radar distance')
            ylabel('Rel. Distance(mm)');
            ax3(2) = subplot(2,1,2);
            plot(radar_hs_locsR(1:end-1)./fs_radar,60./(diff(radar_hs_locsR)./fs_radar));
            hold on;
            plot(tfm_param_time,tfm_param.HR);
            legend('Radar','TFM','Location','south east');
            title('Heart rate comparison');
            xlabel('Time(s)');
            ylabel('Heart rate(BPM)');
            linkaxes(ax3,'x')
            xlim([time_radar(1) time_radar(end)]);
    %         xlim([110 140]); %GDN0023 Apnea

%%
    %MSE_Heartrate = mse(60./(diff(radar_hs_locsR)./fs_radar),tfm_param.HR);
    %disp(['MSE HR IS:',num2str(MSE_Heartrate)]);

            %% Show some TFM aggregated data
            figure('Position',[1 1 scrsz(3) scrsz(4)-80]);
            ax4(1) = subplot(4,1,1);
            plot(tfm_param_time,tfm_param.HR);
            ylabel('HR(BPM)')
            yyaxis right;
            plot(tfm_param_time,tfm_param.LVET);
            title('HR & LVET')
            ylabel('LVET(ms)')

            ax4(2) = subplot(4,1,2);
            plot(tfm_param_time,tfm_param.SV);
            ylabel('SV(ml)')
            yyaxis right;
            plot(tfm_param_time,tfm_param.HZV);
            title('SV & HZV')
            ylabel('HZV(l/min)')

            ax4(3) = subplot(4,1,3);
            plot(tfm_param_time,tfm_param.TPR);
            ylabel('TPR(dyne*s/cm^5)')
            yyaxis right;
            plot(tfm_param_time,tfm_param.TFC);
            title('TPR & TFC')
            ylabel('TFC(1/Ohm)')

            ax4(4) = subplot(4,1,4);
            plot(tfm_param_time,tfm_param.sBP);
            ylabel('sBP(mmHg)')
            yyaxis right;
            plot(tfm_param_time,tfm_param.dBP);
            title('sBP & dBP')
            ylabel('dBP(mmHg)')
            xlabel('Time(s)');

            linkaxes(ax4,'x');
            xlim([tfm_param_time(1) tfm_param_time(end)]);
            
        end
    end
end