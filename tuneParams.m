%% Optimize Kalman Parameters (Massive Grid Search)
% ASSUMPTION: 'dataFull' variable exists in the workspace.
% GOAL: Find global maximum for Correlation (primary) and minimum RMSE (secondary).

if ~exist('dataFull', 'var')
    error('Variable "dataFull" not found. Please run mainFunc.m once to load data.');
end

clc;
fprintf('Starting Massive Parameter Optimization on "dataFull"...\n');

% --- 1. EXPANDED PARAMETER GRID ---
% Q: Process Noise (Lower = smoother/stiffer, Higher = more reactive/jumpy)
Q_list = [0.01, 0.1, 0.5, 1, 2.5, 5, 10, 15, 25, 40, 55, 75, 100, 150, 250, 500];

% R: Measurement Noise (Lower = trust Radar peaks, Higher = trust Kalman model)
R_list = [1, 5, 10, 25, 50, 75, 100, 125, 150, 200, 300, 450, 600, 800, 1000, 1500, 2000];

% P: Initial Error Covariance (Convergence speed)
% Storage for results: [Q, R, P, Mean_RMSE, Mean_Corr]
results = []; 

% Get size of data
[nRows, nCols] = size(dataFull);

% Total iterations
total_iters = length(Q_list) * length(R_list) * length(P_list);
fprintf('Total combinations to test: %d\n', total_iters);

iter_count = 0;
hWait = waitbar(0, 'Initializing Grid Search...');
startTime = tic;

% --- 2. OPTIMIZATION LOOP ---
for q = Q_list
    for r = R_list
            iter_count = iter_count + 1;
            
            % Accumulators
            sum_rmse = 0;
            sum_corr = 0;
            count = 0;
            
            % Loop through every patient/scenario
            for i = 1:nRows
                for j = 1:nCols
                    obj = dataFull{i,j};
                    if isempty(obj), continue; end
                    
                    % --- A. RUN KALMAN ---
                    %obj.KalmanSmooth_BiDir(q, r, p);
                    obj.KalmanFilterBeats(q,r);
                    % --- B. RUN TIME FITTING ---
                    obj.timeFitting();
                    
                    % --- C. CALCULATE METRICS ---
                    est = obj.CorrKalmanHr_on_gt_time;
                    gt  = obj.CorrGt;
                    
                    % Filter NaNs
                    valid = ~isnan(est) & ~isnan(gt);
                    
                    if sum(valid) > 10
                        % RMSE
                        err = est(valid) - gt(valid);
                        rmse = sqrt(mean(err.^2));
                        
                        % Accumulate
                        sum_rmse = sum_rmse + rmse;
                        sum_corr = sum_corr + obj.kalmanCorrValue;
                        count = count + 1;
                  end
            end
            
            % Save averages
            if count > 0
                avg_rmse = sum_rmse / count;
                avg_corr = sum_corr / count;
                results = [results; q, r, p, avg_rmse, avg_corr];
            end
            
            % Update progress bar (every 20 iterations to save UI time)
            if mod(iter_count, 20) == 0
                elapsed = toc(startTime);
                avg_time_per_iter = elapsed / iter_count;
                remain_time = (total_iters - iter_count) * avg_time_per_iter;
                
                % Find best correlation so far for display
                current_best_corr = 0;
                if ~isempty(results)
                    current_best_corr = max(results(:,5));
                end

                waitbar(iter_count/total_iters, hWait, ...
                    sprintf('Progress: %.1f%% | Best Corr: %.4f | ETA: %.0fs', ...
                    (iter_count/total_iters)*100, current_best_corr, remain_time));
            end
        end
    end
end
close(hWait);
totalTime = toc(startTime);
fprintf('Optimization finished in %.1f seconds.\n', totalTime);

% --- 3. ANALYZE RESULTS ---
T = array2table(results, 'VariableNames', {'Q', 'R', 'P', 'RMSE', 'Correlation'});

% --- SORTING: Correlation DESCENDING, then RMSE ASCENDING ---
T_sorted = sortrows(T, {'Correlation', 'RMSE'}, {'descend', 'ascend'});

fprintf('\n=========================================\n');
fprintf('       TOP 20 CONFIGURATIONS             \n');
fprintf('       (Sorted by Max Correlation)       \n');
fprintf('=========================================\n');
disp(T_sorted(1:min(20, height(T_sorted)), :));

% Extract Winner
best_Q = T_sorted.Q(1);
best_R = T_sorted.R(1);
best_P = T_sorted.P(1);
best_Corr = T_sorted.Correlation(1);
best_RMSE = T_sorted.RMSE(1);

fprintf('\nWINNER: Q=%.2f, R=%.2f, P=%.2f\n', best_Q, best_R, best_P);
fprintf('Metrics: Correlation = %.4f | RMSE = %.4f\n', best_Corr, best_RMSE);
fprintf('Applying best parameters to all objects in dataFull...\n');

% --- 4. APPLY BEST PARAMS ---
for i = 1:nRows
    for j = 1:nCols
        if ~isempty(dataFull{i,j})
            dataFull{i,j}.KalmanSmooth_BiDir(best_Q, best_R, best_P);
            dataFull{i,j}.timeFitting();
        end
    end
end
fprintf('Done. dataFull updated.\n');

% --- 5. VISUALIZE SEARCH SPACE (Heatmap) ---
% We visualize R vs Q for the specific optimal P value found.
subset = results(results(:,3) == best_P, :);

if ~isempty(subset)
    uQ = unique(subset(:,1));
    uR = unique(subset(:,2));
    corr_grid = zeros(length(uR), length(uQ));
    
    for i = 1:length(uQ)
        for j = 1:length(uR)
            % Find the row matching this Q, R (and best P)
            row = subset(subset(:,1)==uQ(i) & subset(:,2)==uR(j), :);
            if ~isempty(row)
                corr_grid(j,i) = row(5); % Correlation
            else
                corr_grid(j,i) = NaN;
            end
        end
    end
    
    figure('Name', 'Correlation Heatmap (Expanded)', 'Color', 'w');
    % Use imagesc but force axes to be categorical strings for readability
    imagesc(corr_grid);
    colorbar;
    
    % Beautify axes
    set(gca, 'XTick', 1:length(uQ), 'XTickLabel', string(uQ));
    set(gca, 'YTick', 1:length(uR), 'YTickLabel', string(uR));
    
    xtickangle(45);
    xlabel('Process Noise (Q)');
    ylabel('Measurement Noise (R)');
    title(sprintf('Correlation Landscape (Fixed P=%.1f)\nBrighter is Better', best_P));
    axis xy;
end