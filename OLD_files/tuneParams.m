%% TuneParameters_FineGrid.m
% GOAL: Fine-grained optimization of Q and R using internal object methods.
% OUTPUT: Sorted Excel/CSV file with Best Parameters and Signal Statistics.

% --- PROTECTION: Check Data ---
if ~exist('dataFull', 'var')
    error('Error: Variable "dataFull" is missing. Please run mainFunc.m first.');
end

clc;
fprintf('Starting Fine-Grained Parameter Tuning...\n');
fprintf('Method: Internal Object Calls (Filter + TimeFitting)\n');
fprintf('Grid Size: 40x40 (1600 iterations per filter type per file)\n');
fprintf('This may take time. Progress will be saved incrementally.\n\n');

% Define temp path
tempPath = fullfile(tempdir, 'kalman_fine_progress.mat');

% --- 1. CONFIGURATION ---
[nRows, nCols] = size(dataFull);
outputCsv = 'Fine_Kalman_Tuning.csv';

% --- 2. DEFINE FINE SEARCH GRID ---
% Logarithmic spacing from 10^-5 to 10^5 to catch everything
% 40 points per parameter = 1600 combinations
Q_list = logspace(-5, 5, 40); 
R_list = logspace(-3, 5, 40);

% Prepare Results Container
final_results = {}; 
stats_counter = struct('Processed', 0, 'Skipped', 0);

startTotal = tic;
wb = waitbar(0, 'Initializing...');

% --- 3. MAIN LOOP ---
for i = 1:nRows
    for j = 1:nCols
        
        % Update Waitbar
        progress = ((i-1)*nCols + j) / (nRows*nCols);
        waitbar(progress, wb, sprintf('Patient %d/%d', i, nRows));
        
        % --- Validity Checks ---
        if isempty(dataFull{i,j})
            stats_counter.Skipped = stats_counter.Skipped + 1;
            continue; 
        end
        obj = dataFull{i,j};
        
        % ID & Scenario Extraction
        try
            ID_str = string(obj.ID);
            if isprop(obj, 'scenario'), Scen_str = string(obj.scenario); 
            elseif isprop(obj, 'sceneario'), Scen_str = string(obj.sceneario);
            else, Scen_str = "Unknown"; end
        catch
            ID_str = "Unknown"; Scen_str = "Unknown";
        end

        % --- A. CHECK INPUT DATA ---
        % We need the Median HR Estimate to run the filter
        if isempty(obj.HrEstAfterMedian) || all(isnan(obj.HrEstAfterMedian))
             if isprop(obj, 'medianHR') && ~isempty(obj.medianHR)
                 InputSig = obj.medianHR;
             else
                 % No data to process
                 stats_counter.Skipped = stats_counter.Skipped + 1;
                 continue; 
             end
        else
            InputSig = obj.HrEstAfterMedian;
        end
        
        % Ensure we have Ground Truth for correlation checking
        if isempty(obj.CorrGt)
            try, obj.timeFitting(); catch, end
        end
        if isempty(obj.CorrGt) || sum(~isnan(obj.CorrGt)) < 10
             stats_counter.Skipped = stats_counter.Skipped + 1;
             fprintf('Skipping %s %s: No GT available.\n', ID_str, Scen_str);
             continue;
        end

        stats_counter.Processed = stats_counter.Processed + 1;
        fprintf('Tuning: %s - %s ... ', ID_str, Scen_str);

        % --- B. CALCULATE STATISTICS (On Input Signal) ---
        valid_input = InputSig(~isnan(InputSig));
        
        feat_Var = var(valid_input);
        feat_Roughness = mean(abs(diff(valid_input)));
        feat_PkPk = max(valid_input) - min(valid_input);
        feat_Kurtosis = kurtosis(valid_input);
        feat_Skew = skewness(valid_input);
        
        % Initialize Best Trackers
        best_std = struct('Q', NaN, 'R', NaN, 'RMSE', Inf, 'Corr', -2);
        best_bi  = struct('Q', NaN, 'R', NaN, 'RMSE', Inf, 'Corr', -2);
        
        % --- 4. GRID SEARCH LOOP ---
        % Optimization: Loop R inside Q to minimize Q-matrix re-allocations (minor help)
        
        for iQ = 1:length(Q_list)
            q_val = Q_list(iQ);
            
            for iR = 1:length(R_list)
                r_val = R_list(iR);
                
                % --- FILTER 1: STANDARD (1-State) ---
                try
                    obj.kalmanFilterBeats(q_val, r_val);
                    obj.timeFitting(); % Align to GT
                    
                    c_val = obj.kalmanCorrValue;
                    if c_val > best_std.Corr
                        % Calc RMSE only if Corr is good (optimization)
                        diff_v = obj.CorrKalmanHr_on_gt_time - obj.CorrGt;
                        r_val_rmse = sqrt(nanmean(diff_v.^2));
                        
                        best_std.Corr = c_val;
                        best_std.RMSE = r_val_rmse;
                        best_std.Q = q_val;
                        best_std.R = r_val;
                    end
                catch
                end
                
                % --- FILTER 2: BI-STATE (2-State) ---
                try
                    obj.kalmanFilterBistate(q_val, r_val);
                    obj.timeFitting();
                    
                    c_val = obj.kalmanCorrValue;
                    if c_val > best_bi.Corr
                        diff_v = obj.CorrKalmanHr_on_gt_time - obj.CorrGt;
                        r_val_rmse = sqrt(nanmean(diff_v.^2));
                        
                        best_bi.Corr = c_val;
                        best_bi.RMSE = r_val_rmse;
                        best_bi.Q = q_val;
                        best_bi.R = r_val;
                    end
                catch
                end
                
            end 
        end
        
        fprintf('Done. (Best Std Corr: %.2f | Bi Corr: %.2f)\n', best_std.Corr, best_bi.Corr);
        
        % --- 5. STORE RESULTS ---
        new_row = { ...
            ID_str, Scen_str, ...
            feat_Var, feat_Roughness, feat_PkPk, feat_Kurtosis, feat_Skew, ...
            best_std.Q, best_std.R, best_std.RMSE, best_std.Corr, ...
            best_bi.Q,  best_bi.R,  best_bi.RMSE,  best_bi.Corr ...
        };
       
        final_results = [final_results; new_row];
        
        % Periodic Save
        if mod(length(final_results), 1) == 0
            save(tempPath, 'final_results', 'stats_counter');
        end
    end
end
close(wb);

totalTime = toc(startTotal);
fprintf('\n=================================================\n');
fprintf('Processing Finished in %.2f minutes.\n', totalTime/60);
fprintf('Processed: %d Files | Skipped: %d Files\n', stats_counter.Processed, stats_counter.Skipped);
fprintf('=================================================\n');

% --- 6. EXPORT AND SORT ---
if ~isempty(final_results)
    varNames = {'ID', 'Scenario', ...
                'In_Variance', 'In_Roughness', 'In_PkPk', 'In_Kurtosis', 'In_Skew', ...
                'Std_BestQ', 'Std_BestR', 'Std_MinRMSE', 'Std_MaxCorr', ...
                'Bi_BestQ',  'Bi_BestR',  'Bi_MinRMSE',  'Bi_MaxCorr'};

    ResultTable = cell2table(final_results, 'VariableNames', varNames);
    
    % --- SORTING LOGIC ---
    % 1. Sort by ID (GDN0001, GDN0002...)
    % 2. Then by Scenario
    try
        ResultTable = sortrows(ResultTable, {'ID', 'Scenario'});
    catch
        fprintf('Warning: Sorting failed (check ID format). Saving unsorted.\n');
    end

    writetable(ResultTable, outputCsv);
    fprintf('Results saved to %s\n', outputCsv);
else
    fprintf('No results generated. Check if dataFull has valid data.\n');
end