% this is our main file for the project, from here 

close all;
clc; 
clear ;

%% options 
IDrange = 1;   % GDN00XX - simply choose numbers or ranges from 01 to 30
scenarios = {'Resting'}; %scenarios = {'Resting' 'Valsalva' 'Apnea' 'TiltUp' 'TiltDown'};
ECG_CHANNEL = [2 2 2 2 2 1 2 2 2 2 2 2 2 2 1 2 2 2 2 2 1 1 2 2 2 2 2 2 2 2];
path = 'project_data'; % In this case "datasets" folder is in this scripts folder 
b_USE_PAPER_DATA=false;




%% constants
lambda = 0.0125 ;


%% original code  
scrsz = get(groot,'ScreenSize');
addpath(genpath('utils'))



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
        [radar_i_compensated,radar_q_compensated,phase_compensated,radar_dist] = elreko(radar_i,radar_q,measurement_info{1},1);
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

        %% Time vectors
        time_radar = 1/fs_radar:1/fs_radar:length(radar_dist)/fs_radar;
        time_ecg = 1/fs_ecg:1/fs_ecg:length(tfm_ecg)/fs_ecg;
        time_icg = 1/fs_icg:1/fs_icg:length(tfm_icg)/fs_icg;
        time_bp = 1/fs_bp:1/fs_bp:length(tfm_bp)/fs_bp;
        time_z0 = 1/fs_z0:1/fs_z0:length(tfm_z0)/fs_z0;
        time_respiration = 1/fs_z0:1/fs_z0:length(radar_respiration_re)/fs_z0;





        %% Ellipse reconstruction 
        %% TODO FIX ellipse reconstruction
       
        [radar_i_fitted,radar_q_fitted] = ellipseCreator(radar_i,radar_q,1); % 1 - plot the ellipses, 0 - no plotting
        phase_calc_after_fit = atan( radar_q_fitted./radar_i_fitted) ; 
        radar_distance_calculated_compensated = phase_calc_after_fit .* lambda ./(4*pi) ;
        %TODO : wrap in one function




%% filters  
%%TODO finalize signal processing route, segmentaion to functions 
%TODO add route for Heart_rate calculation

%% Bandpass filter
        f_cutoff_low = 0.05;
        f_cutoff_high = 0.5;
        order=4 ;
        BPF_005_05 = designfilt('bandpassiir', ...
    'FilterOrder', order, ...
    'HalfPowerFrequency1', f_cutoff_low, ...
    'HalfPowerFrequency2', f_cutoff_high, ...
    'SampleRate', fs_radar, ...
    'DesignMethod', 'butter');

     if  b_USE_PAPER_DATA
             radar_dist_after_BPF = filtfilt(BPF_005_05,radar_dist);
     else 
            radar_dist_after_BPF = filtfilt(BPF_005_05,radar_distance_calculated_compensated);
     end

  %% LPF 
        f_cut = 0.5;            
        Wn_low = f_cut / (fs_radar/2); 
        [b_low, a_low] = butter(order, Wn_low, 'low'); 
       if  b_USE_PAPER_DATA
             radar_dist_after_LPF = filtfilt(b_low, a_low, radar_dist);
       else 
            radar_dist_after_LPF = filtfilt(b_low, a_low,radar_distance_calculated_compensated);
       end



       %% RESAMPLING TODO WRITE AS FUNCTION 
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



        %% Respiration Rate calculation
        %% TODO go over understand the process 

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





        %% TODO as function validation 
           rmse = sqrt(mean((radar_rates_BPM - tfm_rates_BPM).^2));
           fprintf('RMSE = %.4f\n', rmse);





        %% Plots
        %% TODO segment to functions


            %respiration plotting

            figure('Position',[1 1 scrsz(3) scrsz(4)-80]);
            ax2(1) = subplot(1,1,1);
            hold on;

            plot(time_respiration, radar_respiration_re.*1000, 'k-', 'DisplayName', 'Radar Respiration');
            plot(time_respiration, tfm_respiration, 'r-', 'DisplayName', 'TFM Respiration');
            plot(time_respiration, respiration_test_BPF.*1000, 'b-', 'DisplayName', 'test BPF Respiration');
            plot(time_respiration, respiration_test_LPF.*1000, 'g-', 'DisplayName', 'test LPF Respiration');

            hold off;
            title('Compare respiration');
            ylabel('Rel. Distance(mm)');
            xlabel('Time(s)');
            legend('show');
            grid on;
    end
end    
>>>>>>> e86559d9519dfda361a3eb8397ac0f963fae3af5
