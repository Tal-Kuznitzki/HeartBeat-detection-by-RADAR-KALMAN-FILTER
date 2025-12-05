% this is a class for the signal and it's relevant information:
classdef radarClass
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
        RrPeaks
        ecgPeaks
        correlated_HrPeaks %saving the corresponding peaks of gt and Hr
        % correlated peak(i,2)=-1 if the peak in (i,1) had no match

        %results
        HrEst
        RrEst
        GtEst
        HrGtMean
        HrGtDiff

        %additional information
        vTimeRadar
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
                ID mustBeNumeric                 
                scenario (1,1) string {mustBeMember(scenario, ["resting","valsalva","apnea","tilt-down","tilt-up"])} 
                fs mustBeNonnegative
                ecg mustBeColumn
                radar_i mustBeColumn 
                radar_q mustBeColumn = 0 %optional, in case we get radar_dist
                
            end

            obj.ID = ID;
            obj.fs_radar = fs;
            obj.sceneario = scenario;
            obj.ecg_gt = ecg;
            if(radar_q==0)
                obj.radar_dist=radar_i;
                sprintf('Only one signal detected, saved as radar_dist')
            else
            obj.radar_i = radar_i;
            obj.radar_q = radar_q;
            end
        end
        
        %I-Q compentasion- create radar_dist from radar_i and radar_Q
            %TODO

        % downSample radar dist- accept new fs or use default. change
        % fs_radar accordingly
        function [DS,obj] = DownSampleRadar(obj,fs)
            if(nargin<1)
                fs=20;

            end
            DS=obj.fs_radar/fs;
            sprintf('Set new radar fs to %d',fs)
            obj.radar_decimated = decimate(obj.radar_dist, DS);
            obj.fs_new = fs;

        end
    %end 
    %methods(Static)
        % apply HR filters and make %TODO use outside filter so we wont
        % designfilt each iteration
        function [obj] = HrFilter(obj)
            if(isnan(obj.fs_new))
                obj.fs_new=obj.fs_radar;
                sprintf('warning: signal not yet decimated')
            end
            firH = designfilt('highpassfir','StopbandFrequency',0.4,...
            'PassbandFrequency',0.7,'StopbandAttenuation',60, ...
            'SampleRate',obj.fs_new);

            firL = designfilt('lowpassfir','PassbandFrequency',3,...
            'StopbandFrequency',4,'StopbandAttenuation',80, ...
            'SampleRate',obj.fs_new);

            % Apply the filters to the decimated radar signal
            obj.HrSignal = filter(firL, obj.radar_decimated);
            obj.HrSignal = filter(firH, obj.HrSignal);
        end

        %
        function [obj] = RrFilter(obj) %currently without highpass. maybe pass through median filter. %TODO use outside filter so we wont
        % designfilt each iteration
            if(isnan(obj.fs_new))
                obj.fs_new=obj.fs_radar;
                sprintf('warning: signal not yet decimated')
            end
            % firH = designfilt('highpassfir','StopbandFrequency',0.02,...
            % 'PassbandFrequency',0.1,'StopbandAttenuation',40, ...
            % 'SampleRate',obj.fs_new);

            firL = designfilt('lowpassfir','PassbandFrequency',0.5,...
            'StopbandFrequency',0.8,'StopbandAttenuation',60, ...
            'SampleRate',obj.fs_new);

            % Apply the filters to the decimated radar signal
            obj.RrSignal = filter(firL, obj.radar_decimated);
            % obj.RrSignal = filter(firH, obj.RrSignal);
        end
        %finding the peaks and normalize to seconds
        function [obj] = FindPeaks(obj)
            thresholdHr= mean(abs((obj.HrSignal)))*0.05;
            thresholdRr= mean(abs((obj.HrSignal)))*0.05;
            [~,obj.HrPeaks, ~,~] = findpeaks(obj.HrSignal, "MinPeakHeight",...
        thresholdHr,'MinPeakDistance',0.33*obj.fs_new);
            [~,obj.RrPeaks, ~,~] = findpeaks(obj.RrSignal, "MinPeakHeight",...
        thresholdRr,'MinPeakDistance',0.33*obj.fs_new);
            [~,obj.ecgPeaks,~] = pan_tompkin(obj.ecg_gt,obj.fs_ecg,0); 
            
            obj.HrPeaks = obj.HrPeaks / obj.fs_new;
            obj.RrPeaks = obj.RrPeaks / obj.fs_new;
            obj.ecgPeaks = obj.ecgPeaks / obj.fs_ecg;
        end

        % finding the rates
        function [obj] = FindRates(obj)
            obj.HrEst = 60 /  diff(obj.HrPeaks);
            obj.GtEst = 60 /  diff(obj.ecgPeaks); 
            obj.RrEst = 60 /  diff(obj.RrPeaks);
        end
        % align HR detected peaks to gt peaks- validation
        % also returns a vector with missed beats and their offset,
        % and a vector of false positives
        function [missed,excess ,obj] = CorrelatePeaks(obj)
            final=zeros(length(obj.ecgPeaks),2);
            localHr=obj.HrPeaks;
            missed=zeros(length(obj.ecgPeaks),2);
            excess = ones(length(obj.HrPeaks),1);
            for i=1:length(obj.ecgPeaks)
                diffI = abs(obj.ecgPeaks(i)-localHr);
                [val,loc]=min(diffI);
                if(val>obj.MaxDiffFromGt) %arbitrary value for max error
                    missed(i,1) = obj.ecgPeaks(i);
                    missed(i,2) = val;
                    final(i,1) = obj.ecgPeaks(i);
                    final(i,2) = -1;
                else 
                    excess(i) = 0;
                    final(i,1) = obj.ecgPeaks(i);
                    final(i,2) = obj.HrPeaks(loc);
                    localHr(loc) = inf; %making sure this value wont apply to 2 different beats.
                end
            end

            obj.correlated_HrPeaks = final;
            obj.HrGtMean
            obj.HrGtDiff
        end
        % function to find false positives independently for radar
        function [obj,removals] = clearFalsePos(obj)
            % go over obj.HrPeaks, find beats that their diff is suspicous
            diffs= diff(obj.HrPeaks);
            removals= zeros(size(obj.HrPeaks));
            for i=2:length(diffs)-2
                if(diffs(i)>1.5*diffs(i-1) || diffs(i+1)>1.5*diffs(i+2))
                    %suspected false positive in index i+1. 
                   % [1 2 @2.5 3 4 5] error in i=3
                   % [1 0.5 0.5 1 1 ] small diff found in i=2
                   %now check if diff(i) + diff(i+1) is around a normal
                   %diff
                   potDiff=diffs(i)+diffs(i+1);
                   avgDiff = (diff(i-1)+diff(i+2))/2;
                   if(potDiff<1.3*(avgDiff)) %decide to remove extra beat
                        obj.HrPeaks(i+1) = 0;  
                   end
                end
            end %end of loop
            obj.HrPeaks= obj.HrPeaks(obj.HrPeaks>0); %remove removed beats 
            %after this function, you should run FindRates again to save he
            %new Heart rate estimation.
           
            
        end 

        % plotting functions
        
        %plot the Hr peaks and the signal
        function [] = plotHrPeaks(obj)
            
        end
    end
end


