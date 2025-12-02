%% this is a main sequence distributed to functions
%%Steps:
% 1. parameters and options
% 2. initialization
% 3. original code for reference
% 4. initial Radar processing- I-Q into a distance vector
% 5. frequency domain processing
% 6. time analysis to gather data
% 7. print our results
% 8. save results
% 9. Clean up (close figures)

% --- STEP 1: Global Initialization (Run only once) ---
b_CLEAN_START = true; 

if b_CLEAN_START
    clc; 
    clear all; 
    close all; 
end

IDrange = 1; %11:12;   
scenarios = {'Resting'}; 
ECG_CHANNEL = [2 2 2 2 2 1 2 2 2 2 2 2 2 2 1 2 2 2 2 2 1 1 2 2 2 2 2 2 2 2];
path = 'project_data'; 
b_USE_PAPER_DATA=1;
resampleFS=100; 
b_plot_ALL=1; 
scrsz = get(groot,'ScreenSize');
addpath(genpath('utils'))
windowSeconds=15; 
windowStep=1; 
saveBaseDir = 'SavedAnalysisFigures'; 
lambda = 0.0125 ;

%% 2. initialization - Loop Start
for indx = 1:length(IDrange)
    ID = sprintf('GDN%04d',IDrange(indx));
    fprintf('----------- Loading %s ------------\n', ID);
    
    for sz = 1:length(scenarios)
        scenario = scenarios{sz};
        fprintf('---- Scenario %s\n', scenario);
        
        % --- FIX 1: Initialize as an Object Array ---
        % 'gobjects(0)' creates an empty array specifically for Graphics Objects.
        % This prevents it from accidentally becoming a numeric array (doubles).
        current_figures = gobjects(0); 
        
        path_id = [path,'\',ID];
        files_synced_mat = dir([path_id,'\*.mat']);
        found = [];
        for j = 1:length(files_synced_mat)
            found = strfind(files_synced_mat(j).name,scenario);
            if ~isempty(found)
                load([path_id,'\',files_synced_mat(j).name]);
                break;
            end
        end
        if isempty(found)
            fprintf('---- skipped\n');
            continue
        end
        
%% 3. original code for reference
        output.(ID).(scenario) = struct;
        [radar_i_compensated,radar_q_compensated,phase_compensated,radar_dist] = elreko(radar_i,radar_q,measurement_info{1},1);
        [radar_respiration, radar_pulse, radar_heartsound, tfm_respiration] = getVitalSigns(radar_dist, fs_radar, tfm_z0, fs_z0);
        
        if ECG_CHANNEL(IDrange(indx)) == 1
            tfm_ecg = fillmissing(tfm_ecg1,'constant',0); 
        else
            tfm_ecg = fillmissing(tfm_ecg2,'constant',0); 
        end
        tfm_ecg = filtButter(tfm_ecg,fs_ecg,4,[1 20],'bandpass');
        
        time_respiration = 1/fs_z0:1/fs_z0:length(tfm_respiration)/fs_z0;
        time_ecg = 1/fs_ecg:1/fs_ecg:length(tfm_ecg)/fs_ecg;
        
%% 4. initial Radar processing
%% 5. frequency domain processing
        vTimeFs= 1/fs_radar:1/fs_radar:length(radar_dist)/fs_radar;
        radar_dist_RsFs=decimate(radar_dist,fs_radar/resampleFS);
        vTimeResample=decimate(vTimeFs,fs_radar/resampleFS);
        
        tic;
        vHeartSignal=HPF_05(radar_dist_RsFs,resampleFS);
        HPFtoc=toc;
        tic;
        vHeartSignalBand=HRfir(radar_dist_RsFs,resampleFS);
        
        % Capture HRfir figure if it exists
        HRfir_h = gcf;
        % FIX 2: Use isgraphics instead of isvalid
        if isgraphics(HRfir_h) && strcmp(get(HRfir_h, 'Name'), 'HRfir_Internal_Plot')
            current_figures(end+1) = HRfir_h;
        end

        BPFtoc=toc;
        tic;
        vRrSignal=BPF_005_05(radar_dist_RsFs,resampleFS);
        RRtoc=toc;

   
            %old HR filtered_Signals
           %{
            h = figure(); 
            set(h,'Position',[1 1 scrsz(3) scrsz(4)-80], 'Name', 'HeartRate_Filtered_Signals');
            current_figures(end+1) = h; 
            
            ax2(1) = subplot(1,1,1);
            hold on;
            plot(vTimeResample, vHeartSignal.*1000, 'g-', 'DisplayName', 'Radar after HPF for HR');
            plot(vTimeResample, vHeartSignalBand.*1000, 'b-', 'DisplayName', 'Radar after BPF for HR');
            plot(time_ecg, tfm_ecg, 'r-', 'DisplayName', 'TFM ecg');
            hold off;
            title(sprintf('Heart Rate Filtered Signals - ID: %s, Scenario: %s', ID, scenario));
            ylabel('Rel. Distance(mm)');
            xlabel('Time(s)');
            legend('show');
            grid on;
           %}
        
%% 6. time analysis
        tic;
        [vHrPeaks,vTimeHrPeaks] = movingWindowHR(vHeartSignalBand,resampleFS,windowSeconds,windowStep,"peaks");
        toc;
        tic;
        [vHrCorr,vTimeHrCorr] = movingWindowHR(vHeartSignalBand,resampleFS,windowSeconds,windowStep,"corr");
        toc;
        tic;
        [vHrGT,vTimeHrGT] = movingWindowHR(tfm_ecg,fs_ecg,windowSeconds,windowStep,"qrs");
        toc;
        tic;
        [vHrCorrMax,vTimeHrCorrMax] = movingWindowHR(vHeartSignalBand,resampleFS,windowSeconds,windowStep,"corrmax");
        toc;
        tic;
        [vHrFft,vTimeHrFft] = movingWindowHR(vHeartSignalBand,resampleFS,windowSeconds,windowStep,"fft");
        toc;
        tic;
        [qrs_amp_raw_ref,qrs_i_raw_ref,delay_ref] = pan_tompkin(tfm_ecg,fs_ecg,0); 
        HR_pan_tompkin_reference = (fs_ecg./(diff(qrs_i_raw_ref))) * 60;
        toc;
        tic;
        thresholdHRbnad= mean(abs((vHeartSignalBand)))*0.05; 
        [pksR,locsR,widthsR,promsR] = findpeaks(vHeartSignalBand, "MinPeakHeight",thresholdHRbnad,'MinPeakDistance',0.33*resampleFS);
        toc;
        vHrFromPeaks = 60*resampleFS./diff(locsR);
        vHrGtPeaks = 60*fs_ecg./diff(qrs_i_raw_ref);

%% 7. print our results
        if(b_plot_ALL)   
        h = figure(); 
        set(h,'Position',[1 1 scrsz(3) scrsz(4)-80], 'Name', 'Respiration_Comparison');
        current_figures(end+1) = h; 
        ax2(1) = subplot(1,1,1);
        hold on;
        plot(vTimeResample, vRrSignal.*1000, 'g-', 'DisplayName', 'Radar Respiration');
        plot(time_respiration, tfm_respiration, 'r-', 'DisplayName', 'TFM Respiration');         
        hold off;
        title(sprintf('Respiration Signals Comparison - ID: %s, Scenario: %s', ID, scenario));
        ylabel('Rel. Distance(mm)');
        xlabel('Time(s)');
        legend('show');
        grid on;
        h = figure(); 
        set(h, 'Name', 'Radar_ECG_Peak_Finders');
        current_figures(end+1) = h; 
        subplot(2,1,1)
        title(sprintf('Radar HR Filters & Peak Finder - ID: %s, Scenario: %s (Radar Signal)', ID, scenario));
        hold on;
        plot(vTimeResample,vHeartSignalBand*1e4,'DisplayName','filtered signal radar');
        plot(locsR/resampleFS,vHeartSignalBand(locsR)*1e4,'*', 'DisplayName','peaks');
        yline(thresholdHRbnad*1e4,'DisplayName','detector threshold');
        hold off;
        subplot(2,1,2)
        title(sprintf('ECG Signal & Peak Finder - ID: %s, Scenario: %s (ECG Signal)', ID, scenario))
        hold on;
        plot(time_ecg,tfm_ecg,'DisplayName','ECG signal');
        plot(qrs_i_raw_ref/fs_ecg,tfm_ecg(qrs_i_raw_ref),'*','DisplayName','ECG peaks');
        hold off;
        linkaxes(get(gcf,'Children'),'x');
        ylabel('Amp');
        xlabel('Time(s)');
        legend('show');
        grid on;

        h = figure(); 
        set(h, 'Name', 'Estimated_HR_IBI');
        current_figures(end+1) = h; 
        hold on;
        title(sprintf('Estimated HR using IBI (Instantaneous) - ID: %s, Scenario: %s', ID, scenario));
        plot(vHrFromPeaks,'DisplayName','RADAR HEART RATE');
        plot(vHrGtPeaks,'DisplayName','ECG HEART RATE');
        ylabel('Heart Rate (BPM)');
        xlabel('Sample Index of IBI');
        legend('show');
        grid on;
        hold off;

        vRMseCorrVsGt=sqrt(mean((vHrCorr- vHrGT).^2));
        vRMseCorrMaxVsGt=sqrt(mean((vHrCorrMax- vHrGT).^2));
        vRMsePeaksVsGt=sqrt(mean((vHrPeaks- vHrGT).^2));
        vMaeCorrVsGt=abs(mean((vHrCorr- vHrGT).^1));
        vMaeCorrMaxVsGt=abs(mean((vHrCorrMax- vHrGT).^1));
        vMaePeaksVsGt=abs(mean((vHrPeaks- vHrGT).^1));

        h = figure(); 
        set(h, 'Name', 'HR_Estimation_Windowed_Comparison');
        current_figures(end+1) = h; 
        hold on;
        plot(vTimeHrPeaks,vHrPeaks, 'g-', 'DisplayName', sprintf('HR-Peaks (RMSE: %.3f)',vRMsePeaksVsGt));
        plot(vTimeHrCorr,vHrCorr, 'b-', 'DisplayName', sprintf('HR-xCorr (RMSE: %.3f)',vRMseCorrVsGt));
        plot(vTimeHrCorrMax,vHrCorrMax, 'y-', 'DisplayName', sprintf('HR-xCorrMax (RMSE: %.3f)',vRMseCorrMaxVsGt));
        plot(vTimeHrGT, vHrGT, 'r-', 'DisplayName', 'TFM ecg HR Ground Truth');
        hold off;
        title(sprintf('Heart Rate Estimation Comparison (Windowed) - ID: %s, Scenario: %s', ID, scenario));
        ylabel('Heart Rate (BPM)');
        xlabel('Time(s)');
        legend('show');
        grid on;

        h = figure(); 
        set(h, 'Name', 'Radar_BPF_Signal_Segment_Peaks');
        current_figures(end+1) = h; 
        hold on;
        plot(vTimeResample, vHeartSignalBand*1e4);
        plot(locsR/resampleFS,vHeartSignalBand(locsR)*1e4,'*'); 
        ylabel('Rel. Distance(10^{-4} mm)');
        xlabel('Time(s)'); 
        title(sprintf('Radar BPF Signal Segment with Peaks - ID: %s, Scenario: %s', ID, scenario));
        hold off;

        [BAmeans,BAdiffs,BAmeanDiff,BACR,BAlinFit]=BlandAltman(vHrGT,vHrPeaks,2,0);

        BA_Hndl = gcf; 
        set(BA_Hndl, 'Name', 'Bland_Altman_Analysis'); 
        title(sprintf('Bland-Altman Analysis (HR Peaks vs GT) - ID: %s, Scenario: %s', ID, scenario));
        current_figures(end+1) = BA_Hndl;
       

%%% HELLO ADD HERE!
        h = figure(); 
        set(h, 'Name', 'Summary_Dashboard_Vertical_Linked', 'Position', [100 50 1000 1200]); 
        current_figures(end+1) = h; 
        
        % Initialize array to store axes handles for linking
        ax_link = [];

        % 1. Radar Heart signal WITH PEAKS
        ax_link(1) = subplot(4,1,1);
        hold on;
        plot(vTimeResample, vHeartSignalBand, 'b', 'DisplayName', 'Radar Signal');
        plot(locsR/resampleFS, vHeartSignalBand(locsR), 'r*', 'MarkerSize', 8, 'DisplayName', 'Radar Peaks');
        hold off;
        title(sprintf('Radar Heart Signal (BPF) - ID: %s, Scenario: %s', ID, scenario));
        xlabel('Time (s)'); ylabel('Amp');
        legend('show'); grid on; axis tight;
        
        % 2. ECG reference WITH PEAKS
        ax_link(2) = subplot(4,1,2);
        hold on;
        plot(time_ecg, tfm_ecg, 'r', 'DisplayName', 'ECG Signal');
        plot(qrs_i_raw_ref/fs_ecg, tfm_ecg(qrs_i_raw_ref), 'k*', 'MarkerSize', 8, 'DisplayName', 'QRS Peaks');
        hold off;
        title(sprintf('ECG Reference Signal - ID: %s, Scenario: %s', ID, scenario));
        xlabel('Time (s)'); ylabel('Amp');
        legend('show'); grid on; axis tight;
        
        % 3. HeartRate (Reference vs Calculated)
        ax_link(3) = subplot(4,1,3);
        hold on;
        plot(vTimeHrGT, vHrGT, 'r', 'LineWidth', 1.5, 'DisplayName', 'ECG Ground Truth');
        plot(vTimeHrPeaks, vHrPeaks, 'b--', 'LineWidth', 1.2, 'DisplayName', 'Radar (Peaks)');
        hold off;
        title(sprintf('Heart Rate Comparison - ID: %s, Scenario: %s', ID, scenario));
        xlabel('Time (s)'); ylabel('BPM');
        legend('Location', 'best');
        grid on; axis tight;
        
        % 4. Bland-Altman (NOT LINKED - Different X-axis domain)
        subplot(4,1,4);
        BlandAltman(vHrGT, vHrPeaks, 2, 0);
        title(sprintf('Bland-Altman Analysis - ID: %s, Scenario: %s', ID, scenario));

        % Link the time-domain axes (1, 2, and 3)
        linkaxes(ax_link, 'x');
        
        end

%% 8. save results
        saveFigures(current_figures, ID, scenario, saveBaseDir);

% --- STEP 9: Clean up ---
        for h = current_figures
            % FIX 3: Use isgraphics in the cleanup loop too
            if isgraphics(h)
                close(h);
            end
        end
    end
end