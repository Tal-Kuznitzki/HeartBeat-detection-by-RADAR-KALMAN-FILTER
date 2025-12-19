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
        resp_gt
        
        %proccessed signals
        radar_decimated
        HrSignal
        RrSignal
        HrPeaks
        HrPeaksFinal
        RrPeaks
        ecgPeaks
        Rrpeaks_gt
        correlated_HrPeaks %saving the corresponding peaks of gt and Hr
        % correlated peak(i,2)=-1 if the peak in (i,1) had no match
       
        %results
        HrEst
        HrEstFinal
        HrEstSpikes
        excessLocsFromHr
        missingLocsFromHr
        RrEst
        HrGtEst
        RrGtEst
        HrGtMean
        HrGtDiff
        mseRaw
        mseFitted
        maeRaw
        maeFitted
        mse2HrRaw
        mse2HrFitted
        mae2HrFitted
        error2sSorted
        correlation_misses
        correlation_excess

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
        function obj = radarClass(ID,scenario,fs,ecg,radar_i, radar_q,gt_resp)
            arguments
                ID            
                scenario (1,1) string {mustBeMember(scenario, ["Resting","Valsalva","Apnea","Tilt-down","Tilt-up"])} 
                fs {mustBeNonnegative}
                ecg {mustBeColumn}
                radar_i {mustBeColumn} 
                radar_q {mustBeColumn} = 0 %optional, in case we get radar_dist
                gt_resp{mustBeColumn} = 0 
            end
                 %TODO : create the process for getting TFM_respiration

            %similliar to getvitalsigns() instead of getting it directly
            obj.ID = string(ID);
            obj.fs_radar = fs;
            obj.fs_ecg = fs;
            obj.sceneario = scenario;
            obj.ecg_gt = ecg;
            obj.resp_gt = gt_resp;
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
        function RrFilter(obj, LPF, HPF) %currently without highpass. maybe pass through median filter. %TODO use outside filter so we wont
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
            firH=HPF;
            end
            % Apply the filters to the decimated radar signal
            obj.RrSignal = filtfilt(firL, obj.radar_decimated);
            obj.RrSignal = filtfilt(firH, obj.RrSignal);
        end
        %finding the peaks and normalize to seconds
        function FindPeaks(obj)
            thresholdHr= mean(abs((obj.HrSignal)))*0.05;
            thresholdRr= mean(abs((obj.RrSignal)))*0.05;
            [~,obj.HrPeaks, ~,~] = findpeaks(obj.HrSignal, "MinPeakHeight",...
        thresholdHr,'MinPeakDistance',0.33*obj.fs_new);
            [~,obj.RrPeaks, ~,~] = findpeaks(obj.RrSignal, "MinPeakHeight",...
        thresholdRr,'MinPeakDistance',2*obj.fs_new);
            [~,obj.ecgPeaks,~] = pan_tompkin(obj.ecg_gt,obj.fs_ecg,0); 
            [~,obj.Rrpeaks_gt, ~,~] = findpeaks(obj.resp_gt, "MinPeakHeight",...
        thresholdRr,'MinPeakDistance',2*(obj.fs_new/2.5));

            
            obj.HrPeaks = obj.HrPeaks / obj.fs_new;
            obj.RrPeaks = obj.RrPeaks / obj.fs_new;
            obj.ecgPeaks = obj.ecgPeaks / obj.fs_ecg;
            obj.Rrpeaks_gt = obj.Rrpeaks_gt /(100);
        end

        % finding the rates
        function FindRates(obj)
            obj.HrEst = 60 ./  diff(obj.HrPeaks);
            obj.HrGtEst = 60 ./  diff(obj.ecgPeaks); 
            obj.RrEst = 60 ./  diff(obj.RrPeaks);
            obj.RrGtEst = 60 ./ diff(obj.Rrpeaks_gt);            
             % TODO: add rates for GT Rr
            if(~isempty(obj.HrPeaksFinal))
                obj.HrEstFinal = 60 ./  diff(obj.HrPeaksFinal);
            end
        end
        % align HR detected peaks to gt peaks- validation
        % also returns a vector with missed beats and their offset,
        % and a vector of false positives
        function [missed,excess] = CorrelatePeaks(obj) %% correlated_HrPeaks
                     
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
            obj.correlated_HrPeaks = final(final(:,2) ~= -1, :);
            obj.HrGtMean = mean(final,2);
            obj.HrGtDiff = final(:,1) - final(:,2);
            obj.correlation_misses=nnz( missed(:,1) ) ;
            obj.correlation_excess = sum(excess(:,1)) ;
            
        end
        % function to find missing beats (positives)
        function [additions] = findMissingBeats(obj)    %% HrPeaksFinal
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
        % function to fix hr spikes on the hrEst
        function [missingBeats, excessBeats] = FindHrSpikes(obj,norm) %%HrEstSpikes
            % this function finds anomalies:
            % excess peaks removed
            % missing peaks added
            % high jitter reduced
            % input: power of jitter reduction (N/N+1)
            % output: Nx2 of excess or missing beats, 
            % first col for old index and 2nd col for new

            hr = obj.HrEst;
            N  = length(hr);
            sizeDiff=0;
            for i=2:N-2
                % check if hr is out of bounds, if its too low, or 
                % irregularly low, we are missing a beat, 
                % double the value and concatenate a copy
                % of it.
                % 
                % if its too high, or irregularly high,
                % remove this sample. 
                
                missingBeats = [];
                excessBeats = [];
                effI= i+sizeDiff;
                %high
                if (hr(effI) > 180 || hr(effI)+hr(effI+1)>1.7*(hr(effI-1)+hr(effI+2)))&&...
                    (hr(effI) > hr(effI-1) && hr(effI) > hr(effI+2))...
                    && (hr(effI+1) > hr(effI-1) && hr(effI+1) > hr(effI+2))

                    hr = [hr(1:effI-1);1/(1/(hr(effI)) + 1/(hr(effI+1)));...
                                                hr(effI+2:N+sizeDiff)];
                    sizeDiff=sizeDiff-1;
                    excessBeats = [excessBeats ; i , effI-1]; 
                    continue;
                end
                %find high var couple to skip
                if(hr(effI)<0.65*hr(effI-1) && hr(effI+1)>1.35*hr(effI+2)) ||...
                  (hr(effI)>1.35*hr(effI-1) && hr(effI+1)<0.65*hr(effI+2))
                    
                    hr(effI)   = (norm*hr(effI-1) + hr(effI))  /(norm+1);
                    hr(effI+1) = (norm*hr(effI+2) + hr(effI+1))/(norm+1);

                    continue;
                end
                %low
                if ( hr(effI)<35 || hr(effI)<0.55*(hr(effI-1)+hr(effI+1))/2)
                    % TODO, check if irregular high caused this
                    hr(effI) = (hr(effI-1)+hr(effI+1))/2;
                    hr = [hr(1:effI);hr(effI);hr(effI+1:N+sizeDiff)];
                    sizeDiff=sizeDiff+1;
                    missingBeats = [missingBeats; i , effI+1];
                    continue;
                end
                
            end
            obj.excessLocsFromHr = excessBeats;
            obj.missingLocsFromHr = missingBeats;
            obj.HrEstSpikes = hr;
        end
        % function to find false positives independently for radar
        function [removals] = clearFalsePos(obj)                      %HrPeaksFinal        
             if(isempty(obj.HrPeaksFinal))
                obj.HrPeaksFinal=obj.HrPeaks;
            end
            % go over obj.HrPeaks, find beats that their diff is suspicous
            diffs= diff(obj.HrPeaksFinal);
            removals = zeros(size(obj.HrPeaksFinal));
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
            
            if length(obj.HrPeaksFinal) < 3
                return; % Not enough peaks to continue analysis
            end
            diffs = diff(obj.HrPeaksFinal);
            
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
        function [] = CorrelateHr(obj)

          gtHr=obj.HrGtEst(:);
          calculatedHr=obj.HrEstSpikes(:);
          minlen=min(length(gtHr),length(calculatedHr));
          gtHr=obj.HrGtEst(1:minlen);
          calculatedHr=obj.HrEstSpikes(1:minlen);
          gt_after_median_filt=medfilt1(gtHr,7);
          calculatedHr_after_median_filt=medfilt1(calculatedHr,7);


            %FOR RADAR SIGNAL
         indx = find(calculatedHr > 1.4 * calculatedHr_after_median_filt | calculatedHr < 0.6 * calculatedHr_after_median_filt); 
            
          hr_replaced_with_median_in_outliers = calculatedHr;
          hr_replaced_with_median_in_outliers(indx) =  calculatedHr_after_median_filt(indx);
          
           %FOR GT
          indx = find(gtHr > 1.4 * gt_after_median_filt | gtHr < 0.6 * gt_after_median_filt); 
            
          gt_replaced_with_median_in_outliers = gtHr;
          gt_replaced_with_median_in_outliers(indx) =  gt_after_median_filt(indx);
          
        hr_replaced_with_median_in_outliers=hr_replaced_with_median_in_outliers(:);
        gt_replaced_with_median_in_outliers=gt_replaced_with_median_in_outliers(:);
        cor=corr(hr_replaced_with_median_in_outliers,gt_replaced_with_median_in_outliers)
        R = corrcoef(hr_replaced_with_median_in_outliers,gt_replaced_with_median_in_outliers);
        correlation_value = R(1,2);
        fprintf('The correlation coefficient is: %.4f\n', correlation_value);
          





          h = figure('Name', 'diff', 'Color', 'w');
          legend('show', 'Location', 'best');
          subplot(2,1,1);
           plot(gtHr,'Color', 'b', 'DisplayName', 'GT');
           hold on;
           plot(gt_after_median_filt,'Color', 'r', 'DisplayName', 'gt filtered');
           plot(gt_replaced_with_median_in_outliers,'Color', 'g', 'DisplayName', 'after median fix');     
           legend('show', 'Location', 'best');
           subplot(2,1,2);
           plot(calculatedHr,'Color', 'b', 'DisplayName', 'Radar Signal');
          hold on;
          plot(calculatedHr_after_median_filt,'Color', 'r', 'DisplayName', 'radar filtered');
          plot(hr_replaced_with_median_in_outliers,'Color', 'g', 'DisplayName', 'after median fix');
          
          
          legend('show', 'Location', 'best');


          linkaxes







        end
        function [] = CalcError(obj)
            
            % pre validation with clearFalsePos, CorrelatePeaks ,
            % FindHrSpikes








            min_len = min(length(obj.HrGtEst), length(obj.HrEstFinal));
            v1=obj.HrEstFinal(1:min_len);
            v2=obj.HrGtEst(1:min_len);
            v2=v2(:);
            mseraw= (v1-v2).^2;
            obj.mseRaw = mean(mseraw);
            obj.maeRaw = mean(abs(v1-v2));
            obj.mse2HrRaw = sortrows([v2, mseraw]); % N,2, acending error per beat rate
            
            v1_gt = 60 ./ diff(obj.correlated_HrPeaks(:,1) );% first colum is GT in this matrix
            v2_radar= 60 ./ diff(obj.correlated_HrPeaks(:,2) );
            
            raw_err = v1_gt-v2_radar;
            mseraw=(raw_err).^2;
            abs_err = abs(raw_err);
            


            obj.mseFitted = mean(mseraw);
            obj.maeFitted = mean(abs_err);
            gt_vs_mse_raw_matrix = [v1_gt, mseraw];
            gt_vs_rawerr_raw_matrix = [v1_gt, raw_err];
            gt_vs_abserr_raw_matrix = [v1_gt, abs_err];
            obj.mse2HrFitted = sortrows(gt_vs_mse_raw_matrix); % N,2, acending error per beat rate
            obj.error2sSorted = sortrows(gt_vs_rawerr_raw_matrix ); 
            obj.mae2HrFitted = sortrows(gt_vs_abserr_raw_matrix );




            
            


            
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
        function h = plotCovDiag(obj)
            h = []; % Default empty if no data
            % Check if estimates exist
            if isempty(obj.HrEst) || isempty(obj.HrGtEst)
                warning('HR Estimates are empty. Run FindRates() first.');
                return;
            end

            h = figure('Name', 'Covariance_Diag_Analysis', 'Color', 'w');
            
            % Sync lengths
          
            Hr_after_corr_fix = obj.HrEstFinal;% 60 ./ diff( obj.correlated_HrPeaks(:,2) );
            Hr_gt_after_corr_fix = obj.HrGtEst;%60 ./ diff( obj.correlated_HrPeaks(:,1) );
            hold on;
            plot(Hr_gt_after_corr_fix,Hr_after_corr_fix,"DisplayName",'Ground truth to RADAR heart rate'...
                ,'Marker','*','LineStyle','none')
            xlabel('ground truth HR')
            ylabel('RADAR estimated HR')
            
            plot((30:1:150),(30:1:150));
            %title(sprintf('Bland-Altman - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
            legend('Measurements', 'Cov Diagonal'); grid on; hold off;
        end       
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
            legend('show', 'Location', 'best'); grid on; hold off;

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
            legend('show', 'Location', 'best'); grid on; hold off;
            
            linkaxes(ax, 'x');
        end

        %% Plot Bland-Altman plot
        function h = plotBA(obj)
            h = []; % Default empty if no data
            % Check if estimates exist
            if isempty(obj.HrEst) || isempty(obj.HrGtEst)
                warning('HR Estimates are empty. Run FindRates() first.');
                return;
            end

            h = figure('Name', 'Bland_Altman_Analysis', 'Color', 'w');
            
            % Sync lengths
            if exist('BlandAltman', 'file')


                Hr_after_corr_fix = 60 ./ diff( obj.correlated_HrPeaks(:,2) );
                Hr_gt_after_corr_fix = 60 ./ diff( obj.correlated_HrPeaks(:,1) );

                BlandAltman(Hr_gt_after_corr_fix, Hr_after_corr_fix, 2, 0);
            else
                % Fallback
                diffs = vec_ecg - vec_radar;
                means = (vec_ecg + vec_radar) / 2;
                plot(means, diffs, 'o', 'DisplayName', 'Data Points');
                hold on;
                yline(mean(diffs), '-r', 'Mean Diff', 'DisplayName', 'Mean Difference');
                yline(mean(diffs) + 1.96*std(diffs), '--r', 'DisplayName', '+1.96 SD');
                yline(mean(diffs) - 1.96*std(diffs), '--r', 'DisplayName', '-1.96 SD');
                xlabel('Mean of methods'); ylabel('Difference');
                legend('show', 'Location', 'best');
                hold off;
            end
            
            title(sprintf('Bland-Altman - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
        end

       %% Plots respiration rate signal (Split into 2 Figures)
            %% Plot Respiration Rates (Trends)
        function h = plotRrRates(obj)
            h = figure('Name', 'Respiration_Rates', 'Color', 'w');
            hold on;
            title(sprintf('Respiration Rate Comparison - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
           
            if ~isempty(obj.RrEst)
                plot(obj.RrEst, 'b.-', 'LineWidth', 1.5, 'DisplayName', 'Radar Respiration');
            end

            if ~isempty(obj.RrGtEst)
                 plot(obj.RrGtEst, 'r.-', 'LineWidth', 1.5, 'DisplayName', 'TFM Respiration');           
            end
            
            ylabel('Breaths Per Minute'); 
            xlabel('Window Index / Time');
            legend('show', 'Location', 'best'); 
            grid on; hold off;
        end

        %% Plot Respiration Signals (Time Domain)
        function h = plotRrSignals(obj)
            h = figure('Name', 'Respiration_Signals_Comparison', 'Color', 'w');
            
            ax = [];
            
            % Subplot 1: Radar Respiration Signal
            ax(1) = subplot(2,1,1);
            hold on;
            title(sprintf('Radar Respiration Signal - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
            
            plot(obj.vTimeNew, obj.RrSignal, 'b-', 'DisplayName', 'Respiration Signal');
            if ~isempty(obj.RrPeaks)
                peakAmps = interp1(obj.vTimeNew, obj.RrSignal, obj.RrPeaks);
                plot(obj.RrPeaks, peakAmps, 'k*', 'MarkerSize', 8, 'DisplayName', 'Radar Peaks');
            end
            ylabel('Rel. Distance(mm)');
            legend('show', 'Location', 'best'); grid on; hold off;

            % Subplot 2: GT Respiration Signal
            ax(2) = subplot(2,1,2);
            hold on;
            title(sprintf('GT Respiration Signal (TFM) - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));

            % Create time vector for GT if not stored
            time_respiration = 1/100:1/100:length(obj.resp_gt)/100;
            
            plot(time_respiration, obj.resp_gt, 'r-', 'DisplayName', 'GT Signal');
            if ~isempty(obj.Rrpeaks_gt)
                peakAmps = interp1(time_respiration, obj.resp_gt, obj.Rrpeaks_gt);
                plot(obj.Rrpeaks_gt, peakAmps , 'k*', 'MarkerSize', 8, 'DisplayName', 'GT Peaks');
            end
            
            ylabel('Rel. Distance(mm)');
            xlabel('Time(s)');
            legend('show', 'Location', 'best'); grid on; hold off;          
            
            linkaxes(ax, 'x');
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
            ylabel('Amp'); legend('show', 'Location', 'best'); grid on; axis tight; hold off;

            % 2. ECG reference WITH PEAKS
            ax_link(2) = subplot(4,1,2);
            hold on;
            plot(obj.vTimeOriginal, obj.ecg_gt, 'k', 'DisplayName', 'ECG Signal');
            if ~isempty(obj.ecgPeaks)
                peakAmpsEcg = interp1(obj.vTimeOriginal, obj.ecg_gt, obj.ecgPeaks);
                plot(obj.ecgPeaks, peakAmpsEcg, 'r*', 'MarkerSize', 8, 'DisplayName', 'QRS Peaks');
            end
            title(sprintf('ECG Reference - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
            ylabel('Amp'); legend('show', 'Location', 'best'); grid on; axis tight; hold off;

            % 3. Heart Rate Comparison (BPM vs Time)
            ax_link(3) = subplot(4,1,3);
            hold on;
            if ~isempty(obj.HrGtEst)
                time_gt_bpm = obj.ecgPeaks(2:end); 
                plot(time_gt_bpm, obj.HrGtEst, 'r.-', 'LineWidth', 1.5, 'DisplayName', 'ECG GT');
            end
            if ~isempty(obj.HrEst)
                time_radar_bpm = obj.HrPeaks(2:end);
                plot(time_radar_bpm, obj.HrEst, 'b.--', 'LineWidth', 1.2, 'DisplayName', 'Radar Est');
            end
            title(sprintf('Heart Rate (BPM) - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
            xlabel('Time (s)'); ylabel('BPM'); legend('show', 'Location', 'best'); grid on; axis tight; hold off;

            % 4. Bland-Altman (NOT LINKED)
            subplot(4,1,4);
            if ~isempty(obj.HrGtEst) && ~isempty(obj.HrEst)


                Hr_after_corr_fix = 60 ./ diff( obj.correlated_HrPeaks(:,2) );
                Hr_gt_after_corr_fix = 60 ./ diff( obj.correlated_HrPeaks(:,1) );
         
                if exist('BlandAltman', 'file')
                   BlandAltman(Hr_gt_after_corr_fix, Hr_after_corr_fix, 2, 0);
                else
                    plot(vec_ecg, vec_radar, 'o', 'DisplayName', 'Data Points'); 
                    xlabel('ECG'); ylabel('Radar');
                    legend('show', 'Location', 'best');
                end
                title(sprintf('Bland-Altman - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
            end

            linkaxes(ax_link, 'x');            
        end
%% Plot Error Analysis (Corrected Titles)
    function h = plotErrors(obj)
            h = figure('Name', 'HR_Error_Analysis', 'Color', 'w');
            
            % --- Subplot 1: Raw Error vs HR ---
            subplot(3,1,1);
            hold on; grid on;
            if ~isempty(obj.error2sSorted)
                plot(obj.error2sSorted(:,1), obj.error2sSorted(:,2), 'ko', 'MarkerFaceColor', 'b', 'MarkerSize', 4);
            end
            yline(0, 'r--', 'LineWidth', 1.5);
            % Updated Title with ID and Scenario
            title(sprintf('Raw Error (Bias) vs ECG HR - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
            ylabel('Error (BPM)');
            legend('Raw Error', 'Zero Line', 'Location', 'best');
            
            % --- Subplot 2: Absolute Error vs HR (The "MAE" request) ---
            subplot(3,1,2);
            hold on; grid on;
            if ~isempty(obj.mae2HrFitted)
                plot(obj.mae2HrFitted(:,1), obj.mae2HrFitted(:,2), 'ko', 'MarkerFaceColor', 'g', 'MarkerSize', 4);
            end
            % Updated Title with ID and Scenario
            title(sprintf('Absolute Error (Magnitude) vs ECG HR - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
            ylabel('|Error| (BPM)');
            yline(mean(obj.mae2HrFitted(:,2)), 'm--', 'DisplayName', 'Mean (MAE)');
            legend('Abs Error', 'MAE Level', 'Location', 'best');

            % --- Subplot 3: Squared Error vs HR (MSE) ---
            subplot(3,1,3);
            hold on; grid on;
            if ~isempty(obj.mse2HrFitted)
                plot(obj.mse2HrFitted(:,1), obj.mse2HrFitted(:,2), 'ko', 'MarkerFaceColor', 'r', 'MarkerSize', 4);
            end
            % Updated Title with ID and Scenario
            title(sprintf('Squared Error (Outliers) vs ECG HR - ID: %s, Scenario: %s', string(obj.ID), obj.sceneario));
            xlabel('ECG Heart Rate (BPM)');
            ylabel('Error^2');
            
            % Print Stats
            fprintf('------------------------------------------------\n');
            fprintf('Error Analysis for ID: %s, Scenario: %s\n', obj.ID, obj.sceneario);
            fprintf('Mean Absolute Error (MAE) -> RAW: %.2f | Fitted: %.2f\n', obj.maeRaw, obj.maeFitted);
            fprintf('Mean Squared Error (MSE)  -> RAW: %.2f | Fitted: %.2f\n', obj.mseRaw, obj.mseFitted);
            fprintf('# Missed beats: %d\n', obj.correlation_misses);
            fprintf('# Excess beats: %d\n', obj.correlation_excess);
            fprintf('------------------------------------------------\n');
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
        function [] = Plot_all(obj, bsave, saveDir) 
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
            h2 = obj.plotRrRates();
            if isgraphics(h2), figHandles(end+1) = h2; end
            h3 = obj.plotRrSignals();
            if isgraphics(h3), figHandles(end+1) = h3; end

            % 3. Bland Altman (only if data exists)
            h4 = obj.plotBA();
            if isgraphics(h4), figHandles(end+1) = h4; end

            % 4. Dashboard
            h5 = obj.plotDashBoard();
            if isgraphics(h5), figHandles(end+1) = h5; end

            % 5. Error Analysis
            h6 = obj.plotErrors();
            if isgraphics(h6), figHandles(end+1) = h6; end

            % Save if requested
            if bsave
                obj.saveFigures(figHandles, saveDir);
                
                % Optional: Close figures after saving to prevent memory buildup during loops
                % close(figHandles); 
            end
        end
    end
end