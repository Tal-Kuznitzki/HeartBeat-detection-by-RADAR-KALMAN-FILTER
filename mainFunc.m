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
b_reset_filter = false;
if b_CLEAN_START
    clc; 
    close all; 
end

if b_reset_filter
    clear all; 
end


b_CLEAR_OLD = true;
b_plot_ALL = true;




IDrange = 1 ; %11:12;   
scenarios ={'Resting'};% {'Resting','Valsalva','Apnea','Tiltdown','Tiltup'}; %["Resting","Valsalva","Apnea","Tilt-down","Tilt-up"]
ECG_CHANNEL = [2 2 2 2 2 1 2 2 2 2 2 2 2 2 1 2 2 2 2 2 1 1 2 2 2 2 2 2 2 2];
path = 'project_data'; 
b_USE_PAPER_DATA=1;
resampleFS=100; 
qGridSize = length( 0.5:0.25:15);
mMseGrid = inf(qGridSize,qGridSize,length(IDrange));

scrsz = get(groot,'ScreenSize');
addpath(genpath('utils'))
windowSeconds=15; 
windowStep=1; 
saveBaseDir = 'SavedAnalysisFigures'; 
statsDirName = 'Statistics' ;
lambda = 0.0125 ;
%initiate full table so the indexes will stay the same
statisticsAPMed = statisticsClass(max(IDrange), 5,statsDirName); % after median, without corr to GT
statisticsPMed = statisticsClass(max(IDrange), 5,statsDirName); % after median and corr to GT
statisticsAPKal = statisticsClass(max(IDrange), 5,statsDirName); % after Kalman, without corr to GT
statisticsPKal = statisticsClass(max(IDrange), 5,statsDirName); % after Kalman and corr to GT

if b_CLEAR_OLD && exist(saveBaseDir,'dir')
    rmdir(saveBaseDir,'s');
end
if b_CLEAR_OLD && exist(statsDirName,'dir')
    rmdir(statsDirName,'s');
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
        %TODO: here it is still -1 
      dataFull{indx,sz} = radarClass(ID,scenario,fs_radar,tfm_ecg,-1*radar_dist,0,tfm_respiration);
      
%% 5. frequency domain processing
        tic
        dataFull{indx,sz}.DownSampleRadar(resampleFS);
        dataFull{indx,sz}.HrFilter(lpf_3,hpf_05);
        dataFull{indx,sz}.RrFilter(lpf_05,hpf_005);
        filteringTime = toc;         
     %% 6. time analysis
       
        dataFull{indx,sz}.FindPeaks(); 
        % generates peaks: HrPeaks, RrPeaks , ecgPeaks ,Rrpeaks_gt
        % based solely on findPeaks() and pan_tompkin 
        % used HrSignal,RrSignal ecg_gt resp_gt

        dataFull{indx,sz}.FindRates(); 
        % based on peaks: Hr, Rr , ecg(gt) ,Rr_gt and peaksFinal ,
        % generates rates: HrEst, HrGtEst, RrEst, RrGtEst 

        
        % based on  HrEst, HrGtEst 
        % generates HrEstAfterMedian and HrGtEstAfterMedian
        % after median filter on each.
        %dataFull{indx,sz}.FindHrSpikes(1);
        dataFull{indx,sz}.MedianHr();
%%
        %dataFull{indx,sz}.KalmanHr();
        dataFull{indx,sz}.KalmanFilterBeats();
        dataFull{indx,sz}.timeFitting(); %THIS RETURNS CORRELATED HR
        [medDelay, kalDelay] = dataFull{indx,sz}.FindMechanicalDelay();
        dataFull{indx,sz}.plot_examples();
        %%
        %TODO: xcorr to find mechanical delay
        % show all results with CorrGt and CorrKalmanHr
        
        dataFull{indx,sz}.CalcError();
        dataFull{indx,sz}.PlotAll(true, saveBaseDir, ...
           'HrEstAfterKalman',name ,...
            dataFull{indx,sz}.HrEstAfterKalman,...
            'plot_RrSignals',false, ...
            'plot_RrRates',false);

        %TODO: CHANGE AFTER WE IMPLEMENT KALMAN ! 
       statisticsAPMed.updateTable(dataFull{indx,sz}.HrEst,dataFull{indx,sz}.HrGtEst,indx,sz); 
       % for q= 0.5:0.25:15
       %   for p = 0.5:0.25:15
       %      dataFull{indx,sz}.KalmanFilterBeats(q,p);
       %      kalman=dataFull{indx,sz}.HrEstAfterKalman(:);
       %      gt= dataFull{indx,sz}.HrGtEst(:);
       %      maxlen= min(length(kalman),length(gt));
       %      mMseGrid(q*4-1, p*4-1,indx) = rmse(kalman(1:maxlen),gt(1:maxlen));
       % 
       %   end
       % end
    end
end
   
% %% CAF on different values
% [N,M,I] = size(mMseGrid);
% A2 = reshape(mMseGrid, N*M, I);   % each column = one (N,M) page
% [minVal, linIdx] = min(A2, [], 1);
% [rowIdx, colIdx] = ind2sub([N, M], linIdx);



%%


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
fprintf('Successfully saved data to:\n %s\n', matFileName);
fprintf('------------------------------------------------\n');