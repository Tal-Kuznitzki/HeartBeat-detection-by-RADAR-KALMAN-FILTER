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

b_CLEAR_OLD = true;
b_plot_ALL = true;

%%% TODO: 
% MAKE SURE EVERYTHING IS PRINTED AS NEEDED 
% MAKE -1 DISTINCT IN THE BA
% INTEGRATE CORRELATE_HR TO FLOW - IT SHOULD OUTPUT MEDIAN FILTER BASELINE 
% CHECK FOR APNEA/VALSALVA TIMES IN THE PAPER
% 
%STATISTICS !



IDrange = 1 ; %11:12;   
scenarios ={'Resting'};% {'Resting','Valsalva','Apnea','Tiltdown','Tiltup'}; %["Resting","Valsalva","Apnea","Tilt-down","Tilt-up"]
ECG_CHANNEL = [2 2 2 2 2 1 2 2 2 2 2 2 2 2 1 2 2 2 2 2 1 1 2 2 2 2 2 2 2 2];
path = 'project_data'; 
b_USE_PAPER_DATA=1;
resampleFS=200; 
 
scrsz = get(groot,'ScreenSize');
addpath(genpath('utils'))
windowSeconds=15; 
windowStep=1; 
saveBaseDir = 'SavedAnalysisFigures'; 
lambda = 0.0125 ;
%initiate full table so the indexes will stay the same
statisticsAPMed = statisticsClass(max(IDrange), 5); % after median, without corr to GT
statisticsPMed = statisticsClass(max(IDrange), 5); % after median and corr to GT
statisticsAPKal = statisticsClass(max(IDrange), 5); % after Kalman, without corr to GT
statisticsPKal = statisticsClass(max(IDrange), 5); % after Kalman and corr to GT
if b_CLEAR_OLD && exist(saveBaseDir,'dir')
    rmdir(saveBaseDir,'s');
end

%% 2. initialization - Loop 
% create filters
if(~exist("lpf_3"))
    [lpf_3,hpf_05]=HRfir(resampleFS); %HR filters
    [hpf_005,lpf_05]=LPF_05(resampleFS); %RR filter
end
% create a matrix for all of our data, divided by ID and scenario
dataFull=cell(length(IDrange), numel(scenarios) ); %a cell for each struct

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
        [radar_i_compensated,radar_q_compensated,phase_compensated,radar_dist] = elreko(radar_i,radar_q,measurement_info{1},0);
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
 %%% TODO: get our own radar_dist

      dataFull{indx,sz} = radarClass(ID,scenario,fs_radar,tfm_ecg,-1*radar_dist,0,tfm_respiration);
      
%% 5. frequency domain processing
        tic
        dataFull{indx,sz}.DownSampleRadar(resampleFS);
        dataFull{indx,sz}.HrFilter(lpf_3,hpf_05);
        dataFull{indx,sz}.RrFilter(lpf_05,hpf_005);
        filteringTime = toc;         
     %% 6. time analysis
       
        dataFull{indx,sz}.FindPeaks(); 
        % generates peaks: Hr, Rr , ecg(gt) ,Rr_gt
        % based solely on findPeaks() and pan_tompkin 
        % used HrSignal,RrSignal ecg_gt resp_gt

        dataFull{indx,sz}.FindRates(0); 
        % based on peaks: Hr, Rr , ecg(gt) ,Rr_gt and peaksFinal ,
        % generates rates: HrEst, HrGtEst, RrEst, RrGtEst 
        % if use_reference is true(default yes) then also generate:
        % HrEstFinal  from HrPeaksFinal
        % Hr_correlated_HrPeaks  from correlated_HrPeaks
        
        dataFull{indx,sz}.FindHrSpikes(2,1);  
        % based on  HrEst, finds missing/excess beats INDEPENDANT of ecg
        % arg1 (default is 2) is power or normalization for jitter
        % now is 2/3and 1/3 
        % arg2 (bool, default true) - whether to update HrEstFinal
        % based on the newly calculated Hr
        % generates:
        % excessLocsFromHr, missingLocsFromHr, HrEstSpikes, HrPeaksFinal
        % excess Beats      Missing Beats      corrected Hr corrected peaks

        dataFull{indx,sz}.MedianHr();
        % based on  HrEstSpikes, HrGtEst 
        % updates HrEstAfterMedian  and HrGtEstAfterMedian after median filter on each.


        % up to here we have
        % HrEstAfterMedian  - hr from HrEstSpikes after median
        % HrPeaksFinal - peaks after cleaning the spikes 
        %   
        %

        %up to here we have not used the reference. good place for KALMAN
        [~,~,v_correlate_result] = dataFull{indx,sz}.CorrelatePeaks(dataFull{indx,sz}.HrPeaksFinal);
        % get misses and excess beats with reference to the ecg
        % based on HrPeaksFinal and ecgPeaks
        % generates correlated_HrPeaks,
        % first col is reference, second is HrPeaksFinal  after taking only
        % the correlated ones. - for the incorrect ones we interpolate! 


        dataFull{indx,sz}.FindRates(1);
        % based on peaks: Hr, Rr , ecg(gt) ,Rr_gt and peaksFinal ,
        % generates rates: HrEst, HrGtEst, RrEst, RrGtEst 
        % if use_reference is true(default yes) then also generate:
        % HrEstFinal  from HrPeaksFinal
        % Hr_correlated_HrPeaks  from correlated_HrPeaks


        dataFull{indx,sz}.CalcError();
        dataFull{indx,sz}.plotCovDiag(); 
        dataFull{indx,sz}.PlotAll(true, saveBaseDir, ...
            'plot_RrSignals',false, ...
            'plot_RrRates',false);
            
        
       statisticsAPMed.updateTable(dataFull{indx,sz}.HrEstFinal,dataFull{indx,sz}.HrGtEst,indx,sz); 
    end
end



%% 8. Save Results
statisticsAPMed.exportExcel();

if ~exist(saveBaseDir, 'dir')
    mkdir(saveBaseDir);
end

matFileName = fullfile(saveBaseDir, 'Processed_Data_Full.mat');

fprintf('------------------------------------------------\n');
fprintf('Processing complete. Converting objects to structs and saving...\n');

% Convert objects to structs

save(matFileName, 'dataFull', 'IDrange', 'scenarios', '-v7.3');

% dataStructs = cell(size(dataFull));
% for i = 1:numel(dataFull)
%     if ~isempty(dataFull{i})
%         dataStructs{i} = struct(dataFull{i});
%     end
% end
% 
% % Save the data AND the mapping keys (IDrange and scenarios)
% % This allows you to know that dataStructs{i, j} corresponds to:
% % Patient ID = IDrange(i)
% % Scenario   = scenarios{j}
% 
% save(matFileName, 'dataStructs', 'IDrange', 'scenarios', '-v7.3');

fprintf('Successfully saved data to:\n %s\n', matFileName);
fprintf('------------------------------------------------\n');