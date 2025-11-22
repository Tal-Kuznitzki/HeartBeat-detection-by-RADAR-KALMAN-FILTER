%% this is a main sequence destributed to functions
%%Steps:
% parameters and options
% initialization
% original code for reference
%then, our signal processing
% initial Radar processing- I-Q into a distance vector
% frequency domain processing
% time analysis to gather data
% print our results
% save results
%% 1. parameters options and constants:
IDrange = 1;   % GDN00XX - simply choose numbers or ranges from 01 to 30
scenarios = {'Resting'}; %scenarios = {'Resting' 'Valsalva' 'Apnea' 'TiltUp' 'TiltDown'};
ECG_CHANNEL = [2 2 2 2 2 1 2 2 2 2 2 2 2 2 1 2 2 2 2 2 1 1 2 2 2 2 2 2 2 2];
path = 'project_data'; % In this case "datasets" folder is in this scripts folder 
b_USE_PAPER_DATA=true;
resampleFS=100; %above 100 may cause issues with iir filters
b_plot_ALL=1; % turn to 0 to stop plotting checkups
scrsz = get(groot,'ScreenSize');
addpath(genpath('utils'))


%constants
lambda = 0.0125 ;




%% 2. initialization


% we will iterate over all id's and all scenarios chosen, for ALL of the
% code from now on. it is possible to run

for indx = 1:length(IDrange)
    % Iterate all subject IDs
    ID = sprintf('GDN%04d',IDrange(indx));
    fprintf('----------- Loading %s ------------\n', ID);

    for sz = 1:length(scenarios)
        % Iterate all existing subject scenarios
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

        
%% 3. original code for reference
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
        
        time_respiration = 1/fs_z0:1/fs_z0:length(tfm_respiration)/fs_z0;
        time_ecg = 1/fs_ecg:1/fs_ecg:length(tfm_ecg)/fs_ecg;
        

%% 4. initial Radar processing- I-Q into a distance vector
% currently use thier radar_dist

%% 5. frequency domain processing
vTimeFs= 1/fs_radar:1/fs_radar:length(radar_dist)/fs_radar;

% filter & decimate signal and get a time vector
radar_dist_RsFs=decimate(radar_dist,fs_radar/resampleFS);
vTimeResample=decimate(vTimeFs,20);
% filter heart rate 
tic;
vHeartSignal=HPF_05(radar_dist_RsFs,resampleFS);
HPFtoc=toc;

tic;
vHeartSignalBand=BPF_05_3(radar_dist_RsFs,resampleFS);
BPFtoc=toc;
% filter respiration rate
tic;
vRrSignal=BPF_005_05(radar_dist_RsFs,resampleFS);
RRtoc=toc;




%% 6. time analysis to gather data
windowSeconds=5;
%1 with find peaks:
tic;
[vHrPeaks,vTimeHrPeaks] = movingWindowHR(vHeartSignalBand,resampleFS,windowSeconds,0,1);
toc
%2 with correlation:
tic
[vHrCorr,vTimeHrCorr] = movingWindowHR(vHeartSignalBand,resampleFS,windowSeconds,1,1);
toc
tic
% GT: TODO- try QRS peak finder for better GT
[vHrGT,vTimeHrGT] = movingWindowHR(tfm_ecg,fs_ecg,windowSeconds,1,1);
toc

%with QRS on reference - 
tic
[qrs_amp_raw_ref,qrs_i_raw_ref,delay_ref] = pan_tompkin(tfm_ecg,fs_ecg,0); 
HR_pan_tompkin_reference = (fs_ecg/mean(diff(qrs_i_raw_ref))) * 60
toc



%TODO: consider kalman filter
%ai: remain casual!


%% 7. print our results

% optional- plot results:
if(b_plot_ALL)
     figure('Position',[1 1 scrsz(3) scrsz(4)-80]);
            ax2(1) = subplot(1,1,1);
            hold on;

            plot(vTimeResample, vRrSignal.*1000, 'g-', 'DisplayName', 'Radar Respiration');
            plot(time_respiration, tfm_respiration, 'r-', 'DisplayName', 'TFM Respiration');
            
            hold off;
            title('Compare respiration');
            ylabel('Rel. Distance(mm)');
            xlabel('Time(s)');
            legend('show');
            grid on;


             %Heartrate 
            figure('Position',[1 1 scrsz(3) scrsz(4)-80]);
            ax2(1) = subplot(1,1,1);
            hold on;
    
            plot(vTimeResample, vHeartSignal.*1000, 'g-', 'DisplayName', 'Radar after HPF for HR');
            plot(vTimeResample, vHeartSignalBand.*1000, 'b-', 'DisplayName', 'Radar after BPF for HR');
            plot(time_ecg, tfm_ecg, 'r-', 'DisplayName', 'TFM ecg');
            hold off;
            title('Compare HR');
            ylabel('Rel. Distance(mm)');
            xlabel('Time(s)');
            legend('show');
            grid on;

end



figure(5);
hold on;
plot(vTimeHrPeaks,vHrPeaks, 'g-', 'DisplayName', 'HR every for BPF with find peaks');
plot(vTimeHrCorr,vHrCorr, 'b-', 'DisplayName', 'HR every for BPF with find peaks and correlation');
plot(vTimeHrGT, vHrGT, 'r-', 'DisplayName', 'TFM ecg HR with find peaks and correlation');
hold off;
title('Compare HR');
ylabel('HR(Bpm)');
xlabel('Time(s)');
legend('show');
grid on;

meanCorrVsGt= mean([vHrGT , vHrCorr],2);
diffCorrVsGt= diff([vHrGT , vHrCorr],1,2);
%            BlandAltman(vHrGT,vHrCorr,2,0)
%scatter(meanCorrVsGt,diffCorrVsGt,'.');
vRMseCorrVsGt=sqrt(mean((vHrCorr- vHrGT).^2))

%% 8. save results in files (png and mat)



%% end of code
    end
end