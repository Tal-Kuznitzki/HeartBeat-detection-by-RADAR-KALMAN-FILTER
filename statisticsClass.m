classdef statisticsClass < handle
    % statisticsClass
    % Stores per-patient/per-scenario metrics (Cov/MSE/MAE) and exports
    % raw estimated/GT HR signals into Excel for later statistical analysis.
    %
    % Scenario index mapping:
    %   j=1 resting, j=2 valsalva, j=3 apnea, j=4 tiltUp, j=5 tiltDown

    properties
        HrCovTable   % (nPatients x nScenarios) covariance between GT and estimate
        HrMseTable   % (nPatients x nScenarios) mean squared error
        HrMaeTable   % (nPatients x nScenarios) mean absolute error

        restingCols
        valsalvaCols
        apneaCols
        tiltUpCols
        tiltDownCols
    end

    properties (Access = private)
        nPatients
        nScenarios
        dirName
    end

    methods
        function obj = statisticsClass(nPatients, nScenarios,dirName)
            % Constructor: prealloc tables and scenario signal storage
            obj.nPatients  = nPatients;
            obj.nScenarios = nScenarios;
            obj.dirName = dirName;
            obj.HrCovTable = nan(nPatients, nScenarios);
            obj.HrMseTable = nan(nPatients, nScenarios);
            obj.HrMaeTable = nan(nPatients, nScenarios);

            blank = repmat(struct('est', [], 'gt', []), nPatients, 1);

            % Always allocate the 5 scenario holders (even if nScenarios < 5)
            obj.restingCols   = blank;
            obj.valsalvaCols  = blank;
            obj.apneaCols     = blank;
            obj.tiltUpCols    = blank;
            obj.tiltDownCols  = blank;
        end

        function updateTable(obj, gtSig, radarSig, i, j)
            % updateTable(gtSig, radarSig, i, j)
            % Stores raw signals into the scenario "Cols" field and updates
            % Cov/MSE/MAE at (i,j) in the metric tables.

            if i < 1 || i > obj.nPatients
                error('updateTable:BadI', 'Patient index i out of range.');
            end
            if j < 1 || j > obj.nScenarios
                error('updateTable:BadJ', 'Scenario index j out of range.');
            end

            gtSig    = gtSig(:);
            radarSig = radarSig(:);

            % Save raw signals for export
            S = obj.getScenarioStruct(j);
            S(i).est = radarSig;
            S(i).gt  = gtSig;
            obj.setScenarioStruct(j, S);

            % Compute metrics on overlapped valid region
            L = min(numel(gtSig), numel(radarSig));
            if L < 2
                obj.HrCovTable(i,j) = nan;
                obj.HrMseTable(i,j) = nan;
                obj.HrMaeTable(i,j) = nan;
                return;
            end

            g = gtSig(1:L);
            r = radarSig(1:L);

            valid = isfinite(g) & isfinite(r);
            if nnz(valid) < 2
                obj.HrCovTable(i,j) = nan;
                obj.HrMseTable(i,j) = nan;
                obj.HrMaeTable(i,j) = nan;
                return;
            end

            g = g(valid);
            r = r(valid);

            C = corrcoef(r, g);               % 2x2
            obj.HrCovTable(i,j) = C(1,2);
            obj.HrMseTable(i,j) = mean((r - g).^2);
            obj.HrMaeTable(i,j) = mean(abs(r - g));
        end

        function filename = exportExcel(obj, filename)
            % exportExcel()
            % exportExcel(filename)
            %
            % Exports:
            %  - HrCovTable, HrMseTable, HrMaeTable  (patients x scenarios)
            %  - One "Signals" sheet per scenario with 2*nPatients columns
            %
            % Returns: filename actually written.
            


            if nargin < 2 || isempty(filename)
                ts = datestr(now, 'yyyymmdd_HHMMSS');
                filename = sprintf('RadarStatistics_%s.xlsx', ts);
            end
            
            outDir = fullfile(pwd, obj.dirName);
            if ~exist(outDir, 'dir')
                mkdir(outDir);
            end
            filename = fullfile(outDir, filename);
            scenNames = obj.getScenarioNames();

            % ---- Metric sheets (patients x scenarios) with headers ----
            obj.writeMetricSheet(filename, 'HrCov', obj.HrCovTable, scenNames);
            obj.writeMetricSheet(filename, 'HrMSE', obj.HrMseTable, scenNames);
            obj.writeMetricSheet(filename, 'HrMAE', obj.HrMaeTable, scenNames);

            % ---- Raw signal sheets per scenario ----
            for j = 1:obj.nScenarios
                sheetName = [scenNames{j} '_Signals'];
                S = obj.getScenarioStruct(j);
                [M, headers] = obj.buildSignalMatrix(S);
                obj.writeSignalSheet(filename, sheetName, M, headers);
            end

            % Optional: meta sheet
            meta = {
                'Created', datestr(now)
                'nPatients', obj.nPatients
                'nScenarios', obj.nScenarios
                'Scenario mapping', '1=Resting, 2=Valsalva, 3=Apnea, 4=TiltUp, 5=TiltDown'
            };
            writecell(meta, filename, 'Sheet', 'Meta', 'Range', 'A1');

        end
    end

    methods (Access = private)

        function scenNames = getScenarioNames(obj)
            % Returns scenario names up to obj.nScenarios
            base = {'Resting','Valsalva','Apnea','Tilt-Up','Tilt-Down'};
            if obj.nScenarios <= numel(base)
                scenNames = base(1:obj.nScenarios);
            else
                % If you ever extend beyond 5, name generically
                scenNames = [base, arrayfun(@(k)sprintf('Scenario%d',k), 6:obj.nScenarios, 'UniformOutput', false)];
            end
        end

        function S = getScenarioStruct(obj, j)
            switch j
                case 1, S = obj.restingCols;
                case 2, S = obj.valsalvaCols;
                case 3, S = obj.apneaCols;
                case 4, S = obj.tiltUpCols;
                case 5, S = obj.tiltDownCols;
                otherwise
                    error('getScenarioStruct:BadJ','Scenario j must be 1..5 (or extend mapping).');
            end
        end

        function setScenarioStruct(obj, j, S)
            switch j
                case 1, obj.restingCols  = S;
                case 2, obj.valsalvaCols = S;
                case 3, obj.apneaCols    = S;
                case 4, obj.tiltUpCols   = S;
                case 5, obj.tiltDownCols = S;
                otherwise
                    error('setScenarioStruct:BadJ','Scenario j must be 1..5 (or extend mapping).');
            end
        end

        function writeMetricSheet(~, filename, sheetName, M, scenNames)
            nPatients  = size(M,1);
            nScenarios = size(M,2);

            hdr = [{'Patient'}, scenNames(1:nScenarios)];
            out = cell(nPatients+1, nScenarios+1);
            out(1,:) = hdr;

            out(2:end,1) = num2cell((1:nPatients).');
            out(2:end,2:end) = num2cell(M);

            writecell(out, filename, 'Sheet', sheetName, 'Range', 'A1');
        end

        function [M, headers] = buildSignalMatrix(obj, S)
            % Builds a numeric matrix (maxLen x 2*nPatients) where
            % col 2*i-1 = est, col 2*i = gt
            n = obj.nPatients;
            maxLen = 0;
            for i = 1:n
                maxLen = max(maxLen, numel(S(i).est));
                maxLen = max(maxLen, numel(S(i).gt));
            end

            M = nan(maxLen, 2*n);
            headers = cell(1, 2*n);

            for i = 1:n
                headers{2*i-1} = sprintf('P%02d_Est', i);
                headers{2*i}   = sprintf('P%02d_GT',  i);

                e = S(i).est(:);
                g = S(i).gt(:);

                if ~isempty(e)
                    M(1:numel(e), 2*i-1) = e;
                end
                if ~isempty(g)
                    M(1:numel(g), 2*i) = g;
                end
            end
        end

        function writeSignalSheet(~, filename, sheetName, M, headers)
            out = cell(size(M,1)+1, size(M,2));
            out(1,:) = headers;
            out(2:end,:) = num2cell(M);
            writecell(out, filename, 'Sheet', sheetName, 'Range', 'A1');
        end

    end
end
