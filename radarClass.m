% this is a class for the signal and it's relevant information:
classdef radarClass < handle
    properties
        %instance details
        ID
        sceneario
        %basic signals
        radar_i
        radar_q
        radar_dist
        ecg_gt
        
        %proccessed signals
        radar_decimated
        HrSignal
        RrSignal
        HrPeaks
        HrPeaksFinal
        RrPeaks
        ecgPeaks
        correlated_HrPeaks %saving the corresponding peaks of gt and Hr
        % correlated peak(i,2)=-1 if the peak in (i,1) had no match

        %results
        HrEst
        HrEstFinal
        RrEst
        GtEst
        HrGtMean
        HrGtDiff
        mseRaw
        mseFitted
        maeRaw
        maeFitted
        mse2HrRaw
        mse2HrFitted

        %additional information
        vTimeOriginal
        vTimeNew
        
        fs_radar
        fs_new
        fs_ecg
        % filters?
        
    end
    properties(Constant)
        MaxDiffFromGt = 0.3; %maximum difference between beats of radar and GT allowed 

    end
    methods
        % constructor: call when making a new object
        function obj = radarClass(ID,scenario,fs,ecg,radar_i, radar_q)
            arguments
                ID {mustBeNonnegative}                 
                scenario (1,1) string {mustBeMember(scenario, ["resting","valsalva","apnea","tilt-down","tilt-up"])} 
                fs {mustBeNonnegative}
                ecg {mustBeColumn}
                radar_i {mustBeColumn} 
                radar_q {mustBeColumn} = 0 %optional, in case we get radar_dist
                
            end

            obj.ID = ID;
            obj.fs_radar = fs;
            obj.fs_ecg = fs;
            obj.sceneario = scenario;
            obj.ecg_gt = ecg;
            if(radar_q==0)
                obj.radar_dist=radar_i;
                sprintf('Only one signal detected, saved as radar_dist')
            else
            obj.radar_i = radar_i;
            obj.radar_q = radar_q;
            end
            obj.vTimeOriginal= 1/fs:1/fs:(length(ecg))/fs; % len-1?

        end
        
        %I-Q compentasion- create radar_dist from radar_i and radar_Q
            %TODO

        % downSample radar dist- accept new fs or use default. change
        % fs_radar accordingly
        function DS = DownSampleRadar(obj,fs)
            if(nargin<2 || isempty(fs))
                fs=100;

            end
            DS=obj.fs_radar/fs;
            sprintf('Set new radar fs to %d',fs)
            obj.radar_decimated = decimate(obj.radar_dist, DS);
            obj.fs_new = fs;
            obj.vTimeNew = 1/fs:1/fs:length(obj.radar_decimated)/fs; %len-1?
        end
    %end 
    %methods(Static)
        % apply HR filters and make %TODO use outside filter so we wont
        % designfilt each iteration
        function HrFilter(obj, LPF, HPF)
            if nargin < 3
                 HPF = [];
            end
            if nargin < 2
                LPF = [];
            end
            if((obj.fs_new)==0)
                obj.fs_new=obj.fs_radar;
                sprintf('warning: signal not yet decimated')
            end
            if(isempty(HPF) || nargin<3)
                firH = designfilt('highpassfir','StopbandFrequency',0.4,...
                'PassbandFrequency',0.7,'StopbandAttenuation',60, ...
                'SampleRate',obj.fs_new);
            else
                firH=HPF;
            end
            if(isempty(LPF) || nargin<2)
                firL = designfilt('lowpassfir','PassbandFrequency',3,...
                'StopbandFrequency',4,'StopbandAttenuation',80, ...
                'SampleRate',obj.fs_new);
            else
                firL=LPF;
            end
            % Apply the filters to the decimated radar signal
            obj.HrSignal = filtfilt(firL, obj.radar_decimated);
            obj.HrSignal = filtfilt(firH, obj.HrSignal);
        end

        %
        function RrFilter(obj, LPF) %currently without highpass. maybe pass through median filter. %TODO use outside filter so we wont
        % designfilt each iteration
          if nargin < 2
                LPF = [];
          end
          if(isnan(obj.fs_new))
              obj.fs_new=obj.fs_radar;
              sprintf('warning: signal not yet decimated')
          end
          
            % firH = designfilt('highpassfir','StopbandFrequency',0.02,...
            % 'PassbandFrequency',0.1,'StopbandAttenuation',40, ...
            % 'SampleRate',obj.fs_new);
            if(isempty(LPF) || nargin<2)
                firL = designfilt('lowpassfir','PassbandFrequency',0.4,...
                'StopbandFrequency',0.8,'StopbandAttenuation',40, ...
                'SampleRate',obj.fs_new);
            else
            firL=LPF;
            end
            % Apply the filters to the decimated radar signal
            obj.RrSignal = filtfilt(firL, obj.radar_decimated);
            % obj.RrSignal = filter(firH, obj.RrSignal);
        end
        %finding the peaks and normalize to seconds
        function FindPeaks(obj)
            thresholdHr= mean(abs((obj.HrSignal)))*0.05;
            thresholdRr= mean(abs((obj.RrSignal)))*0.05;
            [~,obj.HrPeaks, ~,~] = findpeaks(obj.HrSignal, "MinPeakHeight",...
        thresholdHr,'MinPeakDistance',0.33*obj.fs_new);
            [~,obj.RrPeaks, ~,~] = findpeaks(obj.RrSignal, "MinPeakHeight",...
        thresholdRr,'MinPeakDistance',0.33*obj.fs_new);
            [~,obj.ecgPeaks,~] = pan_tompkin(obj.ecg_gt,obj.fs_ecg,0); 
            
            obj.HrPeaks = obj.HrPeaks / obj.fs_new;
            obj.RrPeaks = obj.RrPeaks / obj.fs_new;
            obj.ecgPeaks = obj.ecgPeaks / obj.fs_ecg;
            % TODO: add peaks for GT Rr
        end

        % finding the rates
        function FindRates(obj)
            obj.HrEst = 60 ./  diff(obj.HrPeaks);
            obj.GtEst = 60 ./  diff(obj.ecgPeaks); 
            obj.RrEst = 60 ./  diff(obj.RrPeaks);
             % TODO: add rates for GT Rr
            if(~isempty(obj.HrPeaksFinal))
                obj.HrEstFinal = 60 ./  diff(obj.HrPeaksFinal);
            end
        end
        % align HR detected peaks to gt peaks- validation
        % also returns a vector with missed beats and their offset,
        % and a vector of false positives
        function [missed,excess] = CorrelatePeaks(obj)
            if(isempty(obj.HrPeaksFinal))
                obj.HrPeaksFinal=obj.HrPeaks;
            end
            final=zeros(length(obj.ecgPeaks),2);
            localHr=obj.HrPeaksFinal;
            missed=zeros(length(obj.ecgPeaks),2);
            excess = ones(length(obj.HrPeaksFinal),1);
            for i=1:length(obj.ecgPeaks)
                diffI = abs(obj.ecgPeaks(i)-localHr);
                [val,loc]=min(diffI);
                if(val>obj.MaxDiffFromGt) %arbitrary value for max error
                    missed(i,1) = obj.ecgPeaks(i);
                    missed(i,2) = val;
                    final(i,1) = obj.ecgPeaks(i);
                    final(i,2) = -1;
                else 
                    excess(loc) = 0;
                    final(i,1) = obj.ecgPeaks(i);
                    final(i,2) = obj.HrPeaksFinal(loc);
                    localHr(loc) = inf; %making sure this value wont apply to 2 different beats.
                end
            end

            obj.correlated_HrPeaks = final;
            obj.HrGtMean = mean(final,2);
            obj.HrGtDiff = final(:,1) - final(:,2);
        end
        % function to find missing beats (positives)
        function [additions] = findMissingBeats(obj)
            if(isempty(obj.HrPeaksFinal))
                obj.HrPeaksFinal=obj.HrPeaks;
            end
            hrvec = obj.HrPeaksFinal(:);
            diffs= diff(obj.HrPeaksFinal);
            additions = [];
            sumAdded=0;
            hrLen=length(diffs)+1;
            % go over diffs, if we find a diff that is +-10% x2 than the
            % ones before and after it, wwe missed a bit, make sure the
            % others make sense and are in the HR range.
            %{
                [1 2 4 5] %missing between i=2 and i=3
                [1 2 1] i=2 is sus
                avg = (hrvec(i+sumAdded)+hrvec(i+1+sumAdded))/2
                hrvec= [hrvec(1:i+sumAdded) (avg) hrvec(i+1+sumAdded:len)]
                len += 1;
                sumAdded +=1;
                additions = [additions; avg i]
                i+=1; 
            %}
            for i=2:length(diffs)-1
                if(diffs(i)> 0.9*(diffs(i-1) +diffs(i+1)))... %find anomalies
                     &&(diffs(i)>1.5*(diffs(i-1)&&diffs(i)>1.5*diffs(i+1)))...   
                    && (diffs(i-1)>0.33&&diffs(i-1)<2) &&...
                        (diffs(i+1)>0.33&&diffs(i+1)<2) &&... %make sure peaks are valid
                        (diffs(i)>0.6 &&diffs(i)<4 )
                        % add a peak equal to the average of the last two
                        %peaks, and rearragne the vectors, choose new i 
                        avg = (hrvec(i+sumAdded)+hrvec(i+1+sumAdded))/2;
                        hrvec= [hrvec(1:i+sumAdded) ; (avg) ; hrvec(i+1+sumAdded:hrLen)];
                        additions = [additions; i+1+sumAdded avg];
                        hrLen = hrLen + 1;
                        sumAdded = sumAdded+1;
                        
                    
                end
            end
            obj.HrPeaksFinal = hrvec;

        end
        % function to find false positives independently for radar
        function [removals] = clearFalsePos(obj)
             if(isempty(obj.HrPeaksFinal))
                obj.HrPeaksFinal=obj.HrPeaks;
            end
            % go over obj.HrPeaks, find beats that their diff is suspicous
            diffs= diff(obj.HrPeaksFinal);
            removals= zeros(size(obj.HrPeaksFinal));
            for i=2:length(diffs)-2
                %scenario 1: single phase shift
                %scenario 2: false positive
                if(diffs(i)<0.6*diffs(i-1) || diffs(i+1)<0.6*diffs(i+2))
                   %suspected false positive in index i+1. 
                   % [1 2 @2.5 3 4 5] error in i=3
                   % [1 0.5 0.5 1 1 ] small diff found in i=2
                   % now check if diff(i) + diff(i+1) is around a normal
                   % diff
                   potDiff=diffs(i)+diffs(i+1);
                   avgDiff = (diffs(i-1)+diffs(i+2))/2;
                   if(potDiff<1.5*(avgDiff) && avgDiff<2 && avgDiff>0.33) %decide to remove extra beat
                       removals(i+1) = obj.HrPeaksFinal(i+1);
                      
                       obj.HrPeaksFinal(i+1) = 0;  
                       
                   end
                end
                
            end %end of loop
            obj.HrPeaksFinal= obj.HrPeaksFinal(obj.HrPeaksFinal>0); %remove removed beats 
            for i=2:length(diffs)-2
                %scenario 1: single phase shift:
                 minI=max(2,i-10);
                 meanWin= mean(diffs(minI:i));
                 if ((diffs(i)>1.2*meanWin && diffs(i+1)<1.2*meanWin)...
                  || (diffs(i)<1.2*diffs(i-1) && diffs(i+1)>1.2*diffs(i+2)))
                     %beat i came late
                     if(abs((diffs(i)+diffs(i+1)/2)-meanWin)<meanWin*0.5)
                        obj.HrPeaksFinal(i+1) =...
                       (obj.HrPeaksFinal(i+2)+obj.HrPeaksFinal(i))/2;
                     end
                 end
             end
            %after this function, you should run FindRates again to save he
            %new Heart rate estimation.
           
            
        end 
        
        function [] = CalcError(obj)
            
            v1=obj.HrEstFinal;
            v2=obj.GtEst;
            mseraw=v1.^2-v2.^2;
            obj.mseRaw = sum(mseraw);
            obj.maeRaw = sum(abs(v1-v2));
            obj.mse2HrRaw = sortrows([v2, mseraw]); % N,2, acending error per beat rate
            
            v1=obj.correlated_HrPeaks(:,1);
            v2=obj.correlated_HrPeaks(:,2);
            mseraw=v1.^2-v2.^2;
            obj.mseFitted = sum(mseraw);
            obj.maeFitted = sum(abs(v1-v2));
            obj.mse2HrFitted = sortrows([v2, mseraw]); % N,2, acending error per beat rate
            
            %calculate the MSE, MAE, MRE between v1 and v2:
            


            
        end

       % plotting functions       
        %plot the Hr peaks and the signals
       % function [] = plotHrpeaks(obj) % plot only the HRpeaks
       % function [] = plotBA(obj) %plot BlandAltman plot
       % function [] = plotRR(obj) % plots resipration rate plot 
       % function [] plotDashBoard(obj) 
            %plot the following 
            % 1. radar heart signal after filter(with peaks), 
            % 2. ecg refernce signal with peaks 
            % 3.Heart rate comparison between the GT and the calculated IBI. (all three are linked axis.
            % 4. bland altman plot 
        
% ---------------------------------------------------------
        % Plotting Functions
        % ---------------------------------------------------------

        %% Plot only the HR peaks and the signals (Radar vs ECG)
        function h = plotHrpeaks(obj)
            h = figure('Name', 'Radar_ECG_Peaks', 'Color', 'w');
            
            % Subplot 1: Radar Signal
            ax(1) = subplot(2,1,1);
            hold on;
            title(sprintf('Radar HR Filters & Peak Finder - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
            
            plot(obj.vTimeNew, obj.HrSignal*1e4, 'b', 'DisplayName', 'Filtered Radar');
            
            % Interpolate to find amplitude at peak times
            if ~isempty(obj.HrPeaks)
                peakAmps = interp1(obj.vTimeNew, obj.HrSignal, obj.HrPeaks); 
                plot(obj.HrPeaks, peakAmps*1e4, 'r*', 'MarkerSize', 8, 'DisplayName', 'Radar Peaks');
            end
            
            ylabel('Amp (scaled)');
            legend('show'); grid on; hold off;

            % Subplot 2: ECG Reference
            ax(2) = subplot(2,1,2);
            hold on;
            title(sprintf('ECG Reference Signal - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
            
            plot(obj.vTimeOriginal, obj.ecg_gt, 'k', 'DisplayName', 'ECG Signal');
            
            if ~isempty(obj.ecgPeaks)
                peakAmpsEcg = interp1(obj.vTimeOriginal, obj.ecg_gt, obj.ecgPeaks);
                plot(obj.ecgPeaks, peakAmpsEcg, 'r*', 'MarkerSize', 8, 'DisplayName', 'ECG Peaks');
            end
            
            xlabel('Time(s)'); ylabel('Amp');
            legend('show'); grid on; hold off;
            
            linkaxes(ax, 'x');
        end

        %% Plot Bland-Altman plot
        function h = plotBA(obj)
            h = []; % Default empty if no data
            % Check if estimates exist
            if isempty(obj.HrEst) || isempty(obj.GtEst)
                warning('HR Estimates are empty. Run FindRates() first.');
                return;
            end

            h = figure('Name', 'Bland_Altman_Analysis', 'Color', 'w');
            
            % Sync lengths
            min_len = min(length(obj.GtEst), length(obj.HrEst));
            vec_ecg = obj.GtEst(1:min_len);
            vec_radar = obj.HrEst(1:min_len);

            if exist('BlandAltman', 'file')
                BlandAltman(vec_ecg, vec_radar, 2, 0);
            else
                % Fallback
                diffs = vec_ecg - vec_radar;
                means = (vec_ecg + vec_radar) / 2;
                plot(means, diffs, 'o');
                yline(mean(diffs), '-r', 'Mean Diff');
                yline(mean(diffs) + 1.96*std(diffs), '--r');
                yline(mean(diffs) - 1.96*std(diffs), '--r');
                xlabel('Mean of methods'); ylabel('Difference');
            end
            
            title(sprintf('Bland-Altman - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
        end

        %% Plots respiration rate signal
        function h = plotRR(obj) 
            h = figure('Name', 'Respiration_Signal', 'Color', 'w');
            hold on;
            
            plot(obj.vTimeNew, obj.RrSignal*1000, 'r-', 'DisplayName', 'Radar Respiration');
            if ~isempty(obj.resp_gt)
                 % Assuming resp_gt matches vTimeNew or needs interpolation. 
                 % If lengths differ, plotting might fail, so usually we need a time vector for GT or assume alignment.
                 % Here assuming alignment based on previous code:
                 plot(obj.vTimeNew, obj.resp_gt*1000, 'b-', 'DisplayName', 'TFM Respiration');           
            end
            
            title(sprintf('Respiration Signal - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
            ylabel('Rel. Distance(mm)');
            xlabel('Time(s)');
            legend('show'); grid on; hold off;
        end

        %% Plot the full dashboard (Summary)
        function h = plotDashBoard(obj) 
            h = figure('Name', 'Summary_Dashboard', 'Units', 'normalized', 'Position', [0.1 0.1 0.6 0.8], 'Color', 'w');
            
            ax_link = [];

            % 1. Radar Heart signal WITH PEAKS
            ax_link(1) = subplot(4,1,1);
            hold on;
            plot(obj.vTimeNew, obj.HrSignal, 'b', 'DisplayName', 'Radar Signal');
            if ~isempty(obj.HrPeaks)
                peakAmps = interp1(obj.vTimeNew, obj.HrSignal, obj.HrPeaks);
                plot(obj.HrPeaks, peakAmps, 'r*', 'MarkerSize', 8, 'DisplayName', 'Radar Peaks');
            end
            title(sprintf('Radar Signal - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
            ylabel('Amp'); legend('show'); grid on; axis tight; hold off;

            % 2. ECG reference WITH PEAKS
            ax_link(2) = subplot(4,1,2);
            hold on;
            plot(obj.vTimeOriginal, obj.ecg_gt, 'k', 'DisplayName', 'ECG Signal');
            if ~isempty(obj.ecgPeaks)
                peakAmpsEcg = interp1(obj.vTimeOriginal, obj.ecg_gt, obj.ecgPeaks);
                plot(obj.ecgPeaks, peakAmpsEcg, 'r*', 'MarkerSize', 8, 'DisplayName', 'QRS Peaks');
            end
            title(sprintf('ECG Reference - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
            ylabel('Amp'); legend('show'); grid on; axis tight; hold off;

            % 3. Heart Rate Comparison (BPM vs Time)
            ax_link(3) = subplot(4,1,3);
            hold on;
            if ~isempty(obj.GtEst)
                time_gt_bpm = obj.ecgPeaks(2:end); 
                plot(time_gt_bpm, obj.GtEst, 'r.-', 'LineWidth', 1.5, 'DisplayName', 'ECG GT');
            end
            if ~isempty(obj.HrEst)
                time_radar_bpm = obj.HrPeaks(2:end);
                plot(time_radar_bpm, obj.HrEst, 'b.--', 'LineWidth', 1.2, 'DisplayName', 'Radar Est');
            end
            title(sprintf('Heart Rate (BPM) - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
            xlabel('Time (s)'); ylabel('BPM'); legend('Location', 'best'); grid on; axis tight; hold off;

            % 4. Bland-Altman (NOT LINKED)
            subplot(4,1,4);
            if ~isempty(obj.GtEst) && ~isempty(obj.HrEst)
                min_len = min(length(obj.GtEst), length(obj.HrEst));
                vec_ecg = obj.GtEst(1:min_len);
                vec_radar = obj.HrEst(1:min_len);
                if exist('BlandAltman', 'file')
                    BlandAltman(vec_ecg, vec_radar, 2, 0);
                else
                    plot(vec_ecg, vec_radar, 'o'); xlabel('ECG'); ylabel('Radar');
                end
                title(sprintf('Bland-Altman - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
            end

            linkaxes(ax_link, 'x');            
        end

        %% Save Figures
        function [] = saveFigures(obj, figHandles, saveDir)
             % Set default directory if not provided
            if nargin < 3
                saveDir = 'SavedAnalysisFigures';
            end

            % Filter out invalid handles
            if isempty(figHandles)
                fprintf('No figure handles provided for ID %s.\n', string(obj.ID));
                return;
            end

            % Create directory if it doesn't exist
            if ~exist(saveDir, 'dir')
                mkdir(saveDir);
                fprintf('Created save directory: %s\n', saveDir);
            end

            fprintf('Saving %d figures to: %s\n', length(figHandles), saveDir);

            for i = 1:length(figHandles)
                hFig = figHandles(i);
                
                % Use isgraphics to ensure it's a valid figure
                if ~isgraphics(hFig)
                    continue; 
                end

                % Get Name
                figName = get(hFig, 'Name');
                if isempty(figName)
                    figName = sprintf('Figure%d', get(hFig, 'Number'));
                end

                % Clean filename
                safeFigName = strrep(figName, ' ', '_');
                safeFigName = regexprep(safeFigName, '[^a-zA-Z0-9_]', '');
                
                % Construct full filename: ID_Scenario_FigName
                baseFileName = sprintf('%s_%s_%s', string(obj.ID), obj.sceneario, safeFigName);
                
                pngFileName = fullfile(saveDir, [baseFileName, '.png']);
                figFileName = fullfile(saveDir, [baseFileName, '.fig']);

                try
                    saveas(hFig, pngFileName, 'png'); 
                    saveas(hFig, figFileName, 'fig');
                    fprintf('  -> Saved: %s\n', baseFileName);
                catch ME
                    fprintf(2, 'Warning: Could not save %s. Error: %s\n', baseFileName, ME.message);
                end
            end
        end

        %% Plot All and Optionally Save
        function [] = plot_all_(obj, bsave, saveDir) 
            % bsave - boolean: save all figures?
            % saveDir - optional: string path to save folder
            
            if nargin < 3
                saveDir = 'SavedAnalysisFigures';
            end
            if nargin < 2
                bsave = false;
            end

            figHandles = gobjects(0); % Initialize empty graphics array

            % Generate plots and collect handles
            % 1. HR Peaks
            h1 = obj.plotHrpeaks();
            if isgraphics(h1), figHandles(end+1) = h1; end
            
            % 2. Respiration
            h2 = obj.plotRR();
            if isgraphics(h2), figHandles(end+1) = h2; end

            % 3. Bland Altman (only if data exists)
            h3 = obj.plotBA();
            if isgraphics(h3), figHandles(end+1) = h3; end

            % 4. Dashboard
            h4 = obj.plotDashBoard();
            if isgraphics(h4), figHandles(end+1) = h4; end

            % Save if requested
            if bsave
                obj.saveFigures(figHandles, saveDir);
                
                % Optional: Close figures after saving to prevent memory buildup during loops
                % close(figHandles); 
            end
        end
    end
end


