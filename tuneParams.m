%% Optimize Kalman Parameters (Overnight Mode - Safe Version)
% GOAL: Find best parameters for EACH (ID, Scenario) pair.
% FEATURES: 
%   - Tests Standard vs Bi-Directional filters.
%   - Ultra-fine grid search (~5000 combos per file).
%   - Robust protections against missing data or empty files.
%   - Saves progress automatically to 'temp_optimization_results.mat'.

% --- PROTECTION 1: Check if Data Exists ---
if ~exist('dataFull', 'var')
    error('Error: Variable "dataFull" is missing. Please run mainFunc.m to load your data first.');
end

clc;
fprintf('Starting Robust Overnight Parameter Optimization...\n');
fprintf('Results will be saved progressively to "temp_optimization_results.mat".\n\n');

% --- 1. DEFINITIONS ---
[nRows, nCols] = size(dataFull);
P_fixed = 5.0; % Fixed Covariance

% --- 2. ULTRA-FINE GRID ---
% Q: Process Noise (0.001 to 3000)
q_low = logspace(-3, -1, 15);       % 0.001 -> 0.1
q_mid = logspace(log10(0.15), log10(200), 50); % 0.15 -> 200 (Sweet spot)
q_high = linspace(250, 3000, 20);   % 250 -> 3000
Q_list = unique([q_low, q_mid, q_high]);

% R: Measurement Noise (1 to 10000)
r_low = logspace(0, 1.5, 15);       % 1 -> 30
r_mid = logspace(log10(35), log10(1500), 50); % 35 -> 1500 (Sweet spot)
r_high = linspace(1600, 10000, 20); % 1600 -> 10000
R_list = unique([r_low, r_mid, r_high]);

filterTypes = {'Standard', 'BiDirectional'};

% Prepare Results Container
% Columns: ID, Scenario, BestFilter, BestQ, BestR, MaxCorr, MinRMSE
final_results = {}; 

startTotal = tic;

% --- 3. MAIN LOOP (Per File) ---
for i = 1:nRows
    for j = 1:nCols
        
        % --- PROTECTION 2: Check Object Validity ---
        if isempty(dataFull{i,j})
            continue; % Skip empty cells
        end
        
        obj = dataFull{i,j};
        ID_str = string(obj.ID);
        Scen_str = string(obj.sceneario);
        
        % --- PROTECTION 3: Check Data Content ---
        % If we have no Radar peaks or no ECG ground truth, we can't optimize.
        if isempty(obj.HrPeaks) || length(obj.HrPeaks) < 5
            fprintf('Skipping %s - %s: Not enough Radar peaks.\n', ID_str, Scen_str);
            continue;
        end
        if isempty(obj.ecgPeaks) || length(obj.ecgPeaks) < 5
            fprintf('Skipping %s - %s: Not enough ECG (GT) peaks.\n', ID_str, Scen_str);
            continue;
        end
        
        fprintf('Processing: %s - %s ...\n', ID_str, Scen_str);
        
        % Initialize Best-for-File trackers
        best_for_file.Corr = -2; % Correlation can be -1, so start lower
        best_for_file.Q = NaN;
        best_for_file.R = NaN;
        best_for_file.Type = "None";
        best_for_file.RMSE = Inf;
        
        % --- 4. FILTER TYPE LOOP ---
        for fType = 1:length(filterTypes)
            currentType = filterTypes{fType};
            
            hWait = waitbar(0, sprintf('%s-%s (%s)', ID_str, Scen_str, currentType));
            
            local_best_corr = -2;
            
            % --- 5. GRID SEARCH LOOP ---
            for q_idx = 1:length(Q_list)
                q = Q_list(q_idx);
                
                for r = R_list
                    % --- PROTECTION 4: Mathematical Safety ---
                    try
                        % 1. Run Filter
                        if strcmp(currentType, 'Standard')
                            obj.KalmanFilterBeats(q, r);
                        else
                            obj.KalmanSmooth_BiDir(q, r, P_fixed);
                        end
                        
                        % 2. Run Time Fitting
                        obj.timeFitting();
                        
                        % 3. Check Metrics
                        est = obj.CorrKalmanHr_on_gt_time;
                        gt  = obj.CorrGt;
                        
                        % Ensure vectors are valid and match
                        valid = ~isnan(est) & ~isnan(gt);
                        
                        if sum(valid) > 15
                            curr_corr = obj.kalmanCorrValue;
                            
                            % Check for Local Best (Optimization Step)
                            if curr_corr > local_best_corr
                                local_best_corr = curr_corr;
                                
                                err = est(valid) - gt(valid);
                                curr_rmse = sqrt(mean(err.^2));
                                
                                % Check for Global File Best
                                if curr_corr > best_for_file.Corr
                                    best_for_file.Corr = curr_corr;
                                    best_for_file.Q = q;
                                    best_for_file.R = r;
                                    best_for_file.Type = string(currentType);
                                    best_for_file.RMSE = curr_rmse;
                                end
                            end
                        end
                        
                    catch ME
                        % If Kalman crashes (e.g. singular matrix), we log nothing and continue.
                        % Uncomment below to debug specific errors:
                        % fprintf('Error on Q=%.2f R=%.2f: %s\n', q, r, ME.message);
                    end
                end 
                % Update progress bar per Q-row
                waitbar(q_idx/length(Q_list), hWait);
            end
            close(hWait);
        end
        
        % --- 6. SAVE RESULTS FOR THIS FILE ---
        if best_for_file.Corr > -2
            fprintf('  -> WINNER: %s (Q=%.2f, R=%.0f) | Corr: %.4f\n', ...
                best_for_file.Type, best_for_file.Q, best_for_file.R, best_for_file.Corr);
            
            % Permanently apply best parameters to the object in memory
            if strcmp(best_for_file.Type, 'Standard')
                obj.KalmanFilterBeats(best_for_file.Q, best_for_file.R);
            else
                obj.KalmanSmooth_BiDir(best_for_file.Q, best_for_file.R, P_fixed);
            end
            obj.timeFitting(); 
            
            % Store in list
            new_row = {ID_str, Scen_str, best_for_file.Type, ...
                       best_for_file.Q, best_for_file.R, ...
                       best_for_file.Corr, best_for_file.RMSE};
            final_results = [final_results; new_row];
        else
            fprintf('  -> NO VALID RESULTS FOUND for %s - %s\n', ID_str, Scen_str);
        end
        
        % PROTECTION 5: Periodic Auto-Save
        save('temp_optimization_results.mat', 'final_results', 'dataFull');
    end
end

totalTime = toc(startTotal);
fprintf('\n======================================================\n');
fprintf('OPTIMIZATION COMPLETE in %.2f hours.\n', totalTime/3600);
fprintf('======================================================\n');

% --- 7. FINAL TABLE GENERATION ---
if ~isempty(final_results)
    ResultTable = cell2table(final_results, ...
        'VariableNames', {'ID', 'Scenario', 'BestFilter', 'BestQ', 'BestR', 'MaxCorr', 'MinRMSE'});

    % Sort
    ResultTable = sortrows(ResultTable, {'ID', 'Scenario'});
    
    disp(ResultTable);

    % Save to CSV
    writetable(ResultTable, 'Best_Kalman_Parameters.csv');
    fprintf('Results saved to Best_Kalman_Parameters.csv\n');
else
    fprintf('No results were generated. Please check your data quality.\n');
end