%% Optimize Kalman Parameters (Lightweight & OneDrive Safe)
% GOAL: Find best parameters for EACH (ID, Scenario) pair.
% OUTPUT: A CSV file with the best Q, R, and Filter Type for each file.
% DOES NOT save the heavy 'dataFull' object.

% --- PROTECTION 1: Check if Data Exists ---
if ~exist('dataFull', 'var')
    error('Error: Variable "dataFull" is missing. Please run mainFunc.m to load your data first.');
end

clc;
fprintf('Starting Robust Overnight Parameter Optimization...\n');

% Define a safe temporary file path (bypassing OneDrive)
tempPath = fullfile(tempdir, 'kalman_optimization_progress.mat');
fprintf('Progress will be saved to local temp: %s\n', tempPath);
fprintf('If MATLAB crashes, load that file to recover "final_results".\n\n');

% --- 1. DEFINITIONS ---
[nRows, nCols] = size(dataFull);
P_fixed = 5.0; % Fixed Covariance

% --- 2. ULTRA-FINE GRID ---
% Q: Process Noise (0.001 to 3000)
q_low = logspace(-3, -1, 15);       % 0.001 -> 0.1
q_mid = logspace(log10(0.15), log10(200), 50); % Sweet spot
q_high = linspace(250, 3000, 20);   
Q_list = unique([q_low, q_mid, q_high]);

% R: Measurement Noise (1 to 10000)
r_low = logspace(0, 1.5, 15);       
r_mid = logspace(log10(35), log10(1500), 50); % Sweet spot
r_high = linspace(1600, 10000, 20); 
R_list = unique([r_low, r_mid, r_high]);

filterTypes = {'Standard', 'BiDirectional'};

% Prepare Results Container
final_results = {}; 

startTotal = tic;

% --- 3. MAIN LOOP (Per File) ---
for i = 1:nRows
    for j = 1:nCols
        
        % --- PROTECTION 2: Check Object Validity ---
        if isempty(dataFull{i,j}), continue; end
        
        obj = dataFull{i,j};
        ID_str = string(obj.ID);
        Scen_str = string(obj.sceneario);
        
        % --- PROTECTION 3: Check Data Content ---
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
        best_for_file.Corr = -2; 
        best_for_file.Q = NaN;
        best_for_file.R = NaN;
        best_for_file.Type = "None";
        best_for_file.RMSE = Inf;
        
        % --- 4. FILTER TYPE LOOP ---
        for fType = 1:length(filterTypes)
            currentType = filterTypes{fType};
            
            % Reset local best for this specific filter type
            local_best_corr = -2;
            
            % Progress bar for this file & filter type
            hWait = waitbar(0, sprintf('%s-%s (%s)', ID_str, Scen_str, currentType));
            
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
                        valid = ~isnan(est) & ~isnan(gt);
                        
                        if sum(valid) > 15
                            curr_corr = obj.kalmanCorrValue;
                            
                            % Optimization Check
                            if curr_corr > local_best_corr
                                local_best_corr = curr_corr;
                                
                                % Check Global Best for File
                                if curr_corr > best_for_file.Corr
                                    err = est(valid) - gt(valid);
                                    curr_rmse = sqrt(mean(err.^2));
                                    
                                    best_for_file.Corr = curr_corr;
                                    best_for_file.Q = q;
                                    best_for_file.R = r;
                                    best_for_file.Type = string(currentType);
                                    best_for_file.RMSE = curr_rmse;
                                end
                            end
                        end
                    catch
                        % Ignore math errors
                    end
                end 
                % Update progress bar
                waitbar(q_idx/length(Q_list), hWait);
            end
            close(hWait);
        end
        
        % --- 6. UPDATE OBJECT & SAVE RESULTS ---
        if best_for_file.Corr > -2
            fprintf('  -> WINNER: %s (Q=%.2f, R=%.0f) | Corr: %.4f\n', ...
                best_for_file.Type, best_for_file.Q, best_for_file.R, best_for_file.Corr);
            
            % Update object in memory (so you can plot it later if you want)
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
            fprintf('  -> NO VALID RESULTS for %s - %s\n', ID_str, Scen_str);
        end
        
        % --- PROTECTION 5: Lightweight Incremental Save ---
        % Only saves the parameter list to TEMP folder.
        try
            save(tempPath, 'final_results');
        catch ME
            fprintf(2, 'Warning: Could not save temp file. Continuing... (%s)\n', ME.message);
        end
    end
end

totalTime = toc(startTotal);
fprintf('\n======================================================\n');
fprintf('OPTIMIZATION COMPLETE in %.2f hours.\n', totalTime/3600);
fprintf('======================================================\n');

% --- 7. FINAL EXPORT (Small Files Only) ---
if ~isempty(final_results)
    ResultTable = cell2table(final_results, ...
        'VariableNames', {'ID', 'Scenario', 'BestFilter', 'BestQ', 'BestR', 'MaxCorr', 'MinRMSE'});
    ResultTable = sortrows(ResultTable, {'ID', 'Scenario'});
    
    disp(ResultTable);

    % Save Table to CSV
    csvName = 'Best_Kalman_Parameters.csv';
    writetable(ResultTable, csvName);
    
    % Save Table to MAT (Just the table, not the data)
    matName = 'Optimization_Results_Table.mat';
    save(matName, 'ResultTable'); 
    
    fprintf('Success! Saved parameters to:\n  1. %s\n  2. %s\n', csvName, matName);
else
    fprintf('No results were generated.\n');
end