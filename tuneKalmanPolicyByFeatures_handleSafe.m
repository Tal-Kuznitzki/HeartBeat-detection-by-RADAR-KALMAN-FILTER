function WinnerTable = tuneKalmanPolicyByFeatures_handleSafe(dataFull, outCsv, varargin)
% Handle-safe policy search: Q policies (~20) x R policies (5) per scenario.
% Restores object state after each run so handle mutation doesn't accumulate.
%
% Assumes you will add these fields to radarClass later:
%   HrEstVar, HrEstRoughNorm, IQRadCV, DistSNR_HR_dB, DistPeakToFloor
%
% Assumes these methods exist:
%   obj.KalmanFilterBeats(Q,R)  (or adjust below)
%   obj.timeFitting()
% Uses:
%   obj.kalmanCorrValue  (primary score)

p = inputParser;
p.addParameter('FilterRunner', @runStandard, @(f)isa(f,'function_handle'));
p.addParameter('QClamp', [1e-2, 3e3], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('RClamp', [1, 1e4],   @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('Verbose', true, @(x)islogical(x)&&isscalar(x));
p.parse(varargin{:});
opt = p.Results;

if nargin < 2 || isempty(outCsv), outCsv = "PolicySearch_Results.csv"; end
assert(iscell(dataFull), 'dataFull must be cell array data{i,j}.');

% --- Flatten objects and scenarios ---
objs = {};
scens = strings(0,1);

[nI,nJ] = size(dataFull);
for i=1:nI
    for j=1:nJ
        if isempty(dataFull{i,j}), continue; end
        obj = dataFull{i,j};

        % keep only ones with peaks etc. (same spirit as your tuneParams)
        if isempty(obj.HrPeaks) || numel(obj.HrPeaks) < 2, continue; end
        if isempty(obj.ecgPeaks) || numel(obj.ecgPeaks) < 2, continue; end

        objs{end+1,1} = obj; %#ok<AGROW>
        scens(end+1,1) = string(obj.sceneario);
    end
end

scenarioList = unique(scens);
if opt.Verbose
    fprintf('Found %d valid objects across %d scenarios.\n', numel(objs), numel(scenarioList));
end

% --- Build policy families (same for each scenario; selection happens by results) ---
[Qfuns, Qdesc] = makeQPolicies(opt.QClamp);
[Rfuns, Rdesc] = makeRPolicies(opt.RClamp);

% --- Run per scenario ---
AllRows = {};
for s = scenarioList.'
    idx = find(scens == s);
    if opt.Verbose
        fprintf('\nScenario %s: %d objs | %d Q x %d R = %d combos\n', s, numel(idx), numel(Qfuns), numel(Rfuns), numel(Qfuns)*numel(Rfuns));
    end

    for qi=1:numel(Qfuns)
        for ri=1:numel(Rfuns)
            corrs = nan(numel(idx),1);

            for t=1:numel(idx)
                obj = objs{idx(t)};
                feats = getPreFilterFeatures(obj); % requires fields to exist

                Q = Qfuns{qi}(feats);
                R = Rfuns{ri}(feats);

                snap = snapshotForKalman(obj); % handle-safe
                try
                    corrs(t) = opt.FilterRunner(obj, Q, R);
                catch
                    corrs(t) = NaN;
                end
                restoreSnapshot(obj, snap);
            end

            meanCorr = mean(corrs,'omitnan');
            medCorr  = median(corrs,'omitnan');
            stdCorr  = std(corrs,'omitnan');
            nOK      = sum(~isnan(corrs));

            AllRows(end+1,:) = {string(s), qi, string(Qdesc{qi}), ri, string(Rdesc{ri}), ...
                                meanCorr, medCorr, stdCorr, nOK, numel(idx)}; %#ok<AGROW>
        end
    end
end

T = cell2table(AllRows, 'VariableNames', ...
    {'Scenario','QPolicyID','QPolicyDesc','RPolicyID','RPolicyDesc','MeanCorr','MedianCorr','StdCorr','NValid','NTotal'});
writetable(T, outCsv);

% --- Winners per scenario (max MeanCorr) ---
W = [];
for s = scenarioList.'
    rows = T(T.Scenario==s,:);
    if isempty(rows), continue; end
    [~,ii] = max(rows.MeanCorr);
    W = [W; rows(ii,:)]; %#ok<AGROW>
end
WinnerTable = W;
[pth,nm,~] = fileparts(outCsv);
writetable(WinnerTable, fullfile(pth, nm + "_WINNERS.csv"));

if opt.Verbose
    disp(WinnerTable(:,{'Scenario','QPolicyID','RPolicyID','MeanCorr','NValid'}));
end
end

% ==================== Policy definitions ====================

function [Qfuns, Qdesc] = makeQPolicies(QClamp)
% ~20 options using pairs among: HrEstVar, HrEstRoughNorm, IQRadCV
% feats = [Var, RoughNorm, CV, SNR, P2F]
clamp = @(x) min(QClamp(2), max(QClamp(1), x));
V=1; RN=2; CV=3;

Q0 = 30; % base anchor; you can change globally later

% Utility: safe log mapping
lz = @(x) log10(max(x, eps));

Qfuns = {};
Qdesc = {};

pairs = {[V,RN],[V,CV],[RN,CV]};
weights = [ ...
    1.0 0.0;
    0.0 1.0;
    0.7 0.3;
    0.3 0.7;
    0.5 0.5;
    0.8 0.2;
    0.2 0.8];

% Build 18 from pair+weight with log mapping
for p=1:numel(pairs)
    a = pairs{p}(1); b = pairs{p}(2);
    for w=1:size(weights,1)
        wa = weights(w,1); wb = weights(w,2);
        Qfuns{end+1} = @(f) clamp( Q0 * 10^( wa*lz(f(a)) + wb*lz(f(b)) ) ); %#ok<AGROW>
        Qdesc{end+1} = sprintf('Q=Q0*10^(%.1f*logVar(Feat%d)+%.1f*logVar(Feat%d))', wa,a,wb,b); %#ok<AGROW>
        if numel(Qfuns) >= 18, break; end
    end
    if numel(Qfuns) >= 18, break; end
end

% Add 2 "max/robust" variants (to reach 20)
Qfuns{end+1} = @(f) clamp( Q0 * 10^( max(lz(f(V)), lz(f(RN))) ) );
Qdesc{end+1} = 'Q=Q0*10^(max(logVar,logRoughNorm))';

Qfuns{end+1} = @(f) clamp( Q0 * 10^( 0.7*lz(f(V)) + 0.3*abs(lz(f(CV))) ) );
Qdesc{end+1} = 'Q=Q0*10^(0.7*logVar+0.3*abs(logCV))';

Qfuns = Qfuns(1:20);
Qdesc = Qdesc(1:20);
end

function [Rfuns, Rdesc] = makeRPolicies(RClamp)
% 5 options using pairs among: HrEstVar, DistSNR_HR_dB, DistPeakToFloor
% feats = [Var, RoughNorm, CV, SNR, P2F]
clamp = @(x) min(RClamp(2), max(RClamp(1), x));
V=1; SNR=4; P2F=5;

R0 = 200;

Rfuns = cell(5,1);
Rdesc = cell(5,1);

% Higher SNR/P2F => lower R
Rfuns{1} = @(f) clamp( R0 * 10^(-0.8*(f(SNR)-5.8)) );   Rdesc{1}='R=R0*10^(-0.8*(SNR-const))';
Rfuns{2} = @(f) clamp( R0 * 10^(-0.8*(f(P2F)-3.8)) );   Rdesc{2}='R=R0*10^(-0.8*(P2F-const))';
% Var-driven
Rfuns{3} = @(f) clamp( R0 * 10^(+0.6*log10(max(f(V),eps))) ); Rdesc{3}='R=R0*10^(0.6*logVar)';
% Combined
Rfuns{4} = @(f) clamp( R0 * 10^(-0.5*(f(SNR)-5.8) -0.5*(f(P2F)-3.8)) ); Rdesc{4}='R=R0*10^(-0.5*(SNR-const)-0.5*(P2F-const))';
Rfuns{5} = @(f) clamp( R0 * 10^(-0.6*(f(SNR)-5.8) +0.2*log10(max(f(V),eps))) ); Rdesc{5}='R=R0*10^(-0.6*(SNR-const)+0.2*logVar)';
end

% ==================== Runner + snapshot/restore ====================

function corrVal = runStandard(obj, Q, R)
% Default runner: adjust here if your tuneParams uses extra steps
obj.param_Q = Q;
obj.param_R = R;

obj.KalmanFilterBeats(Q, R);
obj.timeFitting();

corrVal = obj.kalmanCorrValue;
end

function feats = getPreFilterFeatures(obj)
% Requires these fields to exist on obj (you said you'll add them)
req = ["HrEstVar","HrEstRoughNorm","IQRadCV","DistSNR_HR_dB","DistPeakToFloor"];
for k=1:numel(req)
    if ~isprop(obj, req(k))
        error('Missing required pre-filter feature field "%s" on radarClass.', req(k));
    end
end
feats = [obj.HrEstVar, obj.HrEstRoughNorm, obj.IQRadCV, obj.DistSNR_HR_dB, obj.DistPeakToFloor];
end

function snap = snapshotForKalman(obj)
% Minimal snapshot of fields commonly overwritten.
names = { ...
    'param_Q','param_R','param_P', ...
    'KF_HrSignal', ...
    'HrPeaksAfterKalman','HrEstAfterKalman', ...
    'CorrKalmanHr_on_gt_time','kalmanCorrValue', ...
    'mechanicalDelayKalman','corrtime' ...
    };

snap = struct();
for i=1:numel(names)
    nm = names{i};
    if isprop(obj, nm)
        snap.(nm) = obj.(nm);
    end
end
end

function restoreSnapshot(obj, snap)
f = fieldnames(snap);
for i=1:numel(f)
    obj.(f{i}) = snap.(f{i});
end
end
