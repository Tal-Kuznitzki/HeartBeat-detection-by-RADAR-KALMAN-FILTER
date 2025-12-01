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
        
        %proccessed signals
        radar_decimated
        HrSignal
        RrSignal
        HrPeaks
        RrPeaks

        %results
        HrEst
        RrEst
        
        %additional information
        vTimeRadar
        vTimeNew
  
        fs_radar
        fs_new
        
        % filters?
        
    end
    methods
        % constructor: call when making a new object
        function obj = radarClass(ID,scenario,fs,radar_i, radar_q)
            arguments
                ID mustBeNumeric                 
                scenario (1,1) string {mustBeMember(scenario, ["resting","valsalva","apnea","tilt-down","tilt-up"])} 
                fs mustBeNonnegative
                radar_i mustBeColumn 
                radar_q mustBeColumn = 0 %optional, in case we get radar_dist
            end

            obj.ID = ID;
            obj.fs_radar = fs;
            obj.sceneario = scenario;
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
        function DS = DownSampleRadar(fs)
            if(nargin<1)
                fs=20;

            end
            DS=obj.fs_radar/fs;
            sprintf('Set new radar fs to %d',fs)
            obj.radar_decimated = decimate(obj.radar_dist, DS);
            obj.fs_new = fs;

        end
    end 
    methods(Static)
        % apply HR filters and make %TODO use outside filter so we wont
        % designfilt each iteration
        function [] = HrFilter()
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
        function [] = RrFilter() %currently without highpass. maybe pass through median filter. %TODO use outside filter so we wont
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
            obj.RrSignal = filter(firL, obj.redar_decimated);
            % obj.RrSignal = filter(firH, obj.RrSignal);
        end
        %finding the peaks
        function [] = FindPeaks()
            thresholdHr= mean(abs((obj.HrSignal)))*0.05;
            thresholdRr= mean(abs((obj.HrSignal)))*0.05;
            [~,obj.HrPeaks, ~,~] = findpeaks(obj.HrSignal, "MinPeakHeight",...
        thresholdHr,'MinPeakDistance',0.33*fs_new);
            [~,obj.RrPeaks, ~,~] = findpeaks(obj.RrSignal, "MinPeakHeight",...
        thresholdRr,'MinPeakDistance',0.33*fs_new);
        end

        % finding the rates
        function [] = FindRates()
            obj.HrEst = obj.fs_new * (60 /  diff(obj.HrPeaks));
            obj.RrEst = obj.fs_new * (60 /  diff(obj.RrPeaks));
        end
        %

        % plotting functions
        
        %plot the Hr peaks and the signal
        function [] = plotHrPeaks()
            
        end
    end
end


