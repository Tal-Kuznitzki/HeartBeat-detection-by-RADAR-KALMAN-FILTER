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





IDrange = 1 ; %11:12;   
scenarios ={'Resting','Valsalva'};% {'Resting','Valsalva','Apnea','Tilt-down','Tilt-up'}; %["Resting","Valsalva","Apnea","Tilt-down","Tilt-up"]
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

if b_CLEAR_OLD && exist(saveBaseDir,'dir')
    rmdir(saveBaseDir,'s');
end

%% 2. initialization - Loop 
% create filters
[lpf_3,hpf_05]=HRfir(resampleFS); %HR filters
[hpf_005,lpf_05]=LPF_05(resampleFS); %RR filter
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
      dataFull{indx,sz} = radarClass(ID,scenario,fs_radar,tfm_ecg,radar_dist,0,tfm_respiration);
      
%% 5. frequency domain processing
        tic
        dataFull{indx,sz}.DownSampleRadar(resampleFS);
        dataFull{indx,sz}.HrFilter(lpf_3,hpf_05);
        dataFull{indx,sz}.RrFilter(lpf_05,hpf_005);
        filteringTime = toc;
        
        

        
   
           
        
%% 6. time analysis
       
        dataFull{indx,sz}.FindPeaks();
        dataFull{indx,sz}.FindRates();
        dataFull{indx,sz}.HrEst();
        dataFull{indx,sz}.RrEst();
        
        %TODO: update clearFalsePos to make a new vector, use findRates and
        %HrEst on the new vector, so we keep the data
        dataFull{indx,sz}.FindHrSpikes(2); % power or normalization for jitter
        dataFull{indx,sz}.clearFalsePos();
        dataFull{indx,sz}.CorrelatePeaks();
        dataFull{indx,sz}.FindRates();
        dataFull{indx,sz}.HrEst();
%
        dataFull{indx,sz}.HrEstFinal=dataFull{indx,sz}.HrEstSpikes;
%
        dataFull{indx,sz}.CalcError();
        dataFull{indx,sz}.Plot_all(true, saveBaseDir);
        
        
    end
end



%% 8. Save Results
if ~exist(saveBaseDir, 'dir')
    mkdir(saveBaseDir);
end

matFileName = fullfile(saveBaseDir, 'Processed_Data_Full.mat');

fprintf('------------------------------------------------\n');
fprintf('Processing complete. Converting objects to structs and saving...\n');

% Convert objects to structs
dataStructs = cell(size(dataFull));
for i = 1:numel(dataFull)
    if ~isempty(dataFull{i})
        dataStructs{i} = struct(dataFull{i});
    end
end

% Save the data AND the mapping keys (IDrange and scenarios)
% This allows you to know that dataStructs{i, j} corresponds to:
% Patient ID = IDrange(i)
% Scenario   = scenarios{j}
save(matFileName, 'dataStructs', 'IDrange', 'scenarios', '-v7.3');

fprintf('Successfully saved data to:\n %s\n', matFileName);
fprintf('------------------------------------------------\n');