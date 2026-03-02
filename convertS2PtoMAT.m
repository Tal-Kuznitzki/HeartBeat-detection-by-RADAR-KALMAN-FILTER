function [matFileName,radar_i,radar_q] = convertS2PtoMAT(s2pFilePath, matFileName)
            % convertS2PtoMAT: Parses a 2-port Touchstone .s2p file and 
            % extracts I and Q from the S22 parameter, saving to a .mat file.
            
            % Read the numeric data matrix, ignoring '!' and '#' comment lines
            opts = detectImportOptions(s2pFilePath, 'FileType', 'text');
            opts.CommentStyle = {'!', '#'};
            data = readmatrix(s2pFilePath, opts);
            
            % Check that we have the 9 columns (Freq, S11_MA, S21_MA, S12_MA, S22_MA)
            if size(data, 2) >= 9
                data(any(isnan(data), 2), :) = [];
                S22_mag = data(:, 8);
                S22_ang_deg = data(:, 9);
            else
                error('Unexpected .s2p data format. Expected 9 columns.');
            end
            
            % Convert MA (Magnitude-Angle) to I and Q components
            radar_i = S22_mag .* cosd(S22_ang_deg);
            radar_q = S22_mag .* sind(S22_ang_deg);
            
            % Mock ground truth arrays to keep the constructor happy
            signal_gt = zeros(size(obj.radar_i));
            resp_gt = zeros(size(obj.radar_i));
            
            % Default fs (you can adjust this if your s2p samples are captured at a different rate)
            fs_radar = 100;             
            % Save to the standard .mat format expected by the rest of the script
            save(matFileName, 'radar_i', 'radar_q', 'signal_gt', 'resp_gt', 'fs_radar');
            fprintf('Successfully converted %s to %s\n', s2pFilePath, matFileName);
        end