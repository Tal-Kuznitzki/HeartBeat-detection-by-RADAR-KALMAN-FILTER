function saveFigures(figHandles, ID, scenario, saveDir)
% SAVEFIGURES Saves a specific list of figure handles to the specified directory.
    
    % Filter out any empty or invalid entries immediately
    if isempty(figHandles)
        fprintf('No figure handles provided to save for ID %s, Scenario %s.\n', ID, scenario);
        return;
    end

    % 2. Verify and create the single folder
    if ~exist(saveDir, 'dir')
        [status, msg, msgID] = mkdir(saveDir);
        if ~status
            fprintf(2, 'Error creating save directory: %s\n', msg);
            return;
        end
        fprintf('Created save directory: %s\n', saveDir);
    end

    % 3. Loop through and save each figure
    fprintf('Saving %d figures to: %s\n', length(figHandles), saveDir);

    for i = 1:length(figHandles)
        hFig = figHandles(i);

        % --- FIX: Use isgraphics() ---
        % isgraphics is robust and works on both double handles and objects.
        % isvalid() crashes if hFig is accidentally a double (number).
        if ~isgraphics(hFig)
             fprintf(2, 'Warning: Skipped an invalid or non-graphic handle.\n');
             continue;
        end
        
        % Get the plot title/description from the figure's Name property
        figName = get(hFig, 'Name');
        if isempty(figName)
            figName = sprintf('Figure%d', get(hFig, 'Number'));
        end

        safeFigName = strrep(figName, ' ', '_');
        safeFigName = regexprep(safeFigName, '[^a-zA-Z0-9_]', '');
        baseFileName = sprintf('%s_%s_%s', ID, scenario, safeFigName);
        
        pngFileName = fullfile(saveDir, [baseFileName, '.png']);
        figFileName = fullfile(saveDir, [baseFileName, '.fig']);

        try
            saveas(hFig, pngFileName, 'png'); 
            saveas(hFig, figFileName, 'fig');
            fprintf('  -> Overwriting: %s.png and .fig\n', baseFileName);
        catch ME
            fprintf(2, 'Warning: Could not save figure %s. Error: %s\n', figName, ME.message);
        end
    end
    fprintf('Figure saving complete.\n');
end