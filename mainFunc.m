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

folders = ["OLD_files","Functions\","original_code\","plots\",...
    "project_data\","SavedAnalysisFigures\","Statistics\","utils"];
for folderName = folders
    addpath(genpath(folderName));
end

b_CLEAN_START = false;
b_reset_filter = false;
b_CLEAR_OLD_mat=true;

if b_CLEAN_START
    clc; 
    close all; 
end

if b_reset_filter
    clear all; 
end


b_CLEAR_OLD = false;
b_plot_ALL = false;

IDrange = [52] ; %11:12;  

scenarios = ["Apnea"]; %["Resting","Valsalva","Apnea","TiltDown","TiltUp"]

ECG_CHANNEL = [2 2 2 2 2 1 2 2 2 2 2 2 2 2 1 2 2 2 2 2 1 1 2 2 2 2 2 2 2 2];
path = 'project_data'; 
resampleFS=100; 

scrsz = get(groot,'ScreenSize');
addpath(genpath('utils'))
saveBaseDir = 'SavedAnalysisFigures'; 
statsDirName = 'Statistics' ;
lambda = 0.0125 ;
%initiate full table so the indexes will stay the same
statisticsAPMed = statisticsClass(max(IDrange), 5,statsDirName); % after median, without corr to GT

if b_CLEAR_OLD && exist(saveBaseDir,'dir')
    rmdir(saveBaseDir,'s');
end
if b_CLEAR_OLD && exist(statsDirName,'dir')
    rmdir(statsDirName,'s');
end

%% 2. initialization - Loop 
% create filters
if(~exist("lpf_3"))
    [lpf_3,hpf_05,lpf_5]=HRfir(resampleFS); %HR filters
    [hpf_005,lpf_05]=LPF_05(resampleFS); %RR filter
end
% create a matrix for all of our data, divided by ID and scenario
dataFull=cell(length(IDrange), numel(scenarios)); %a cell for each struct

for indx = 1:length(IDrange)
    numericID = IDrange(indx);
    b_lab = (numericID > 40); 

    ID = sprintf('GDN%04d', numericID); % e.g., GDN0041
    path_id = fullfile(path, ID);       % e.g., project_data/GDN0041
    
    fprintf('----------- Loading %s ------------\n', ID);
    
    % Check if this subject requires S2P parsing
   

    for sz = 1:length(scenarios)
        scenario = scenarios{sz};
        fprintf('---- Scenario %s\n', scenario);
        
        current_figures = gobjects(0); 
        
 if b_lab
    %% --- NEW S2P / PPG LOGIC ---

    numericID = IDrange(indx);

    s2pFileName = fullfile(path_id, sprintf('GD%04d_%s.s2p', numericID, scenario));
    matFileName = fullfile(path_id, sprintf('GDN%04d_%s.mat', numericID, scenario));
    vidFileName = fullfile(path_id, sprintf('GD%04d_%s.mp4', numericID, scenario));

    % Cleanup for old .mat files
    if (numericID > 40) && b_CLEAR_OLD_mat
        ID = sprintf('GDN%04d', numericID);
        path_id = fullfile(path, ID);

        if exist(path_id, 'dir')
            target_files = fullfile(path_id, '*.mat');
            delete(target_files);
            fprintf('Cleared old .mat files from: %s\n', path_id);
        else
            fprintf('Folder not found, skipping: %s\n', path_id);
        end
    end

    %% --- S2P handling ---
    if exist(s2pFileName, 'file') && ~exist(matFileName, 'file')
        fprintf('Found s2pFile %s, converting to matFile %s \n', s2pFileName, matFileName);
        [~, radar_i, radar_q] = convertS2PtoMAT(s2pFileName, matFileName);

    elseif exist(matFileName, 'file')
        warning('old matFile %s found.', matFileName);

    elseif ~exist(s2pFileName, 'file')
        warning('s2pFileName file %s not found. Skipping.', s2pFileName);
        continue;
    end

    %% --- PPG INPUT ---
    if numericID > 50
        % New format: CSV PPG file
        % File only needs to contain scenario name
        % Example:
        % GDN0051_Apnea_XXXX.csv
        % GDN0051_Resting_XXXX.csv
        % GDN0051_Number_XXXX.csv

        csvPattern = fullfile(path_id, sprintf('*%s*.csv', scenario));
        csvFiles = dir(csvPattern);

        if isempty(csvFiles)
            warning('PPG csv file not found for ID %d, scenario %s. Skipping.', ...
                numericID, scenario);
            continue;
        end

        csvPath = fullfile(path_id, csvFiles(1).name);

        % Read CSV (skip 6 header lines)
        ppgMat = readmatrix(csvPath, ...
            'FileType', 'text', ...
            'NumHeaderLines', 6);

        % Column mapping from your files:
        % 1 = RED
        % 2 = IR
        % 3 = RED without ambient
        % 4 = IR without ambient

        tfm_ecg = ppgMat(:,4);   % IR without ambient channel

    else
        % Old lab format: MP4 PPG video
        mVideoPPG = VideoReader(vidFileName);
        tfm_ecg = mVideoPPG;
    end

    %% --- LOAD RADAR MAT ---
    radar = load(matFileName);

    radar_i   = radar.radar_i;
    radar_q   = radar.radar_q;
    fs_radar  = radar.fs_radar;

    dataFull{indx,sz} = radarClass( ...
        ID, ...
        scenario, ...
        fs_radar, ...
        tfm_ecg, ...
        radar_i, ...
        radar_q, ...
        0, ...
        b_lab);

    b_comp = 1;
    b_mode = 1;
    b_plot = 0;
    if b_comp
        dataFull{indx,sz}.IQcompensation(b_plot, b_mode);
    end

    dataFull{indx,sz}.calculateRadarDistFromIQ();

else
    %% --- LEGACY .MAT LOGIC ---

    files_synced_mat = dir(fullfile(path_id, sprintf('*%s*.mat', scenario)));

    if isempty(files_synced_mat)
        fprintf('---- skipped\n');
        continue;
    end

    load(fullfile(path_id, files_synced_mat(1).name));

    output.(ID).(scenario) = struct;

    [radar_i_compensated, radar_q_compensated, phase_compensated, radar_dist] = ...
        elreko(radar_i, radar_q, measurement_info{1}, 0);

    [radar_respiration, radar_pulse, radar_heartsound, tfm_respiration] = ...
        getVitalSigns(radar_dist, fs_radar, tfm_z0, fs_z0);

    if ECG_CHANNEL(numericID) == 1
        tfm_ecg = fillmissing(tfm_ecg1, 'constant', 0);
    else
        tfm_ecg = fillmissing(tfm_ecg2, 'constant', 0);
    end

    tfm_ecg = filtButter(tfm_ecg, fs_ecg, 4, [1 20], 'bandpass');

    if numericID == 10 && scenario == "TiltDown"
        tfm_ecg   = tfm_ecg(2000:end);
        radar_dist = radar_dist(2000:end);
    end

    dataFull{indx,sz} = radarClass( ...
        ID, ...
        scenario, ...
        fs_radar, ...
        tfm_ecg, ...
        radar_dist, ...
        0, ...
        tfm_respiration, ...
        b_lab);
end
        %% 5. frequency domain processing
        tic
        dataFull{indx,sz}.DownSampleRadar(resampleFS)
        dataFull{indx,sz}.HrFilter(lpf_3,hpf_05,lpf_5);
        dataFull{indx,sz}.RespFilter(lpf_05,hpf_005);

        dataFull{indx,sz}.NormalizeHrSignal(1.0); 
        if ~b_lab
        dataFull{indx,sz}.RespFilter(lpf_05,hpf_005);
 %      dataFull{indx,sz}.HrSignal = dataFull{indx,sz}.KF_HrSignal;
        end
      filtering_time=toc;
     %% 6. time analysis
       
        dataFull{indx,sz}.FindPeaks(); 
        % generates peaks: HrPeaks, RrPeaks , ecgPeaks ,Rrpeaks_gt
        % based solely on findPeaks() and pan_tompkin 
        % used HrSignal,RrSignal ecg_gt resp_gt
        

        [peaksDelay,sign] = dataFull{indx,sz}.FindMechDelay();
        
        %finds the delay between the  peaks and ECG
        dataFull{indx,sz}.radar_dist = sign.* dataFull{indx,sz}.radar_dist ;
        if(sign==-1 || peaksDelay<0.25 || peaksDelay > 0.45)
             tic
            dataFull{indx,sz}.radar_dist = -1.* dataFull{indx,sz}.radar_dist ;
        
            dataFull{indx,sz}.DownSampleRadar(resampleFS);
            dataFull{indx,sz}.HrFilter(lpf_3,hpf_05,lpf_5);
            

            filteringTime = toc;         
            dataFull{indx,sz}.NormalizeHrSignal(1.0);

            %dataFull{indx,sz}.kalmanSmoothRadarDist();
            
           
    
        end
        dataFull{indx,sz}.FindPeaks(); 
        % generates peaks: HrPeaks, RrPeaks , ecgPeaks ,Rrpeaks_gt
        % based solely on findPeaks() and pan_tompkin 
        % used HrSignal,RrSignal ecg_gt resp_gt

        dataFull{indx,sz}.FindRates(); 
        % based on peaks: Hr, Rr , ecg(gt) ,Rr_gt and peaksFinal ,
        % generates rates: HrEst, HrGtEst, RrEst, RrGtEst 
        dataFull{indx,sz}.SmoothSpikesHr(1.4);
        dataFull{indx,sz}.ComputePreFilterStats();




        
        % based on  HrEst, HrGtEst 
        % generates HrEstAfterMedian and HrGtEstAfterMedian
        % after median filter on each.


        % % [Q,R] = dataFull{indx,sz}.ProduceKalmanCoeff(); 
        % Q
        % R
        % [q_auto, r_auto] = dataFull{indx,sz}.EstimateKalmanCoeffs('BiState');
        % % Apply them
        % q_auto
        % r_auto
        
        %dataFull{indx,sz}.kalmanFilterBeats_nH(q_auto,r_auto) 
        %dataFull{indx,sz}.KalmanFilterHrGrid(0); %1 to draw CAF NEW
        
        %dataFull{indx,sz}.OptimizeKalman_Innovation(850,0);
        dataFull{indx,sz}.OptimizeKalman_NSubSignals(50,true,2);
        %dataFull{indx,sz}.OptimizeKalman_Innovation_AdaptiveR(850,0);
        dataFull{indx,sz}.MedianHr(); 
        %dataFull{indx,sz}.KalmanSmooth_BiDir();
        % generates HrPeaksAfterKalman and HrEstAfterKalman

        % or use Bistate - 
        %obj.kalmanFilterBistate(q_auto, r_auto);



        %[medDelay, kalDelay] = dataFull{indx,sz}.FindMechanicalDelay();

        dataFull{indx,sz}.timeFitting(); %generates CORRELATED HR
      
         dataFull{indx,sz}.plot_examples();
        dataFull{indx,sz}.plotRespRates();
        dataFull{indx,sz}.plotRespSignals();
        %%
        % show all results with CorrGt and CorrKalmanHr
        dataFull{indx,sz}.CalcError(dataFull{indx,sz}.CorrKalmanHr_on_gt_time);
        dataFull{indx,sz}.PlotHrCovAndBA();
        % % dataFull{indx,sz}.PlotAll(true, saveBaseDir, ...
        % %    'HR after Kalman & time fit',...
        % %     dataFull{indx,sz}.CorrKalmanHr_on_gt_time,...
        % %     dataFull{indx,sz}.HrPeaksAfterKalman,...
        % %     dataFull{indx,sz}.corrtime,... %time vector after fitting
        % %     'plot_RrSignals',false, ...
        % %     'plot_RrRates',false);



   
       statisticsAPMed.updateTable...
       (dataFull{indx,sz}.CorrKalmanHr_on_gt_time,dataFull{indx,sz}.CorrGt,indx,sz); 
       % for q= 0.5:0.25:15
       %   for p = 0.5:0.25:15
       %      dataFull{indx,sz}.KalmanFilterBeats(q,p);
       %      kalman=dataFull{indx,sz}.HrEstAfterKalman(:);
       %      gt= dataFull{indx,sz}.HrGtEst(:);
       %      maxlen= min(length(kalman),length(gt));
       %      mMseGrid(q*4-1, p*4-1,indx) = rmse(kalman(1:maxlen),gt(1:maxlen));
       % 
       %   end
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
fprintf('Preparing to save .mat file...\n');

% --- OPTIMIZATION: Clear heavy raw data to fix the "Hang" ---
fprintf('Clearing raw radar signals from memory to speed up saving...\n');
for i = 1:numel(dataFull)
    if ~isempty(dataFull{i})
        % Remove the massive raw vectors. 
        % We only need the processed Heart Rate estimates and stats.
        dataFull{i}.radar_i = [];
        dataFull{i}.radar_q = [];
        dataFull{i}.radar_dist = []; 
        
        % If you have other huge intermediate vectors, clear them too:
        % dataFull{i}.radar_decimated = []; 
    end
end
fprintf('Raw data cleared.\n');

fprintf('Saving variable "dataFull" to %s...\n', matFileName);
fprintf('This might still take a moment, but should not hang.\n');

% Save
%save(matFileName, 'dataFull', 'IDrange', 'scenarios', '-v7.3');

fprintf('Successfully saved data to:\n %s\n', matFileName);
fprintf('------------------------------------------------\n');