function out = kalmanSmoothRadarDist(obj, config)
% kalmanSmoothRadarDist
% Kalman smoothing for raw radar_dist at high fs (e.g. 2000 Hz).
% Supports 2-state (x, xdot) and 3-state (x, xdot, xddot).
%
% INPUTS
%   radar_dist : Nx1 (or 1xN) double, raw distance (meters or arbitrary units)
%   fs         : sampling rate in Hz (e.g. 2000)
%   config     : struct with optional fields:
%       .order            (2 or 3) default 3
%       .useRTS           (true/false) default true (Rauch-Tung-Striebel smoother)
%       .initPScale       default 10
%       .sigmaMeas        measurement noise std (same units as radar_dist)
%       .sigmaAcc         process accel noise std (units/s^2) for order=3
%       .sigmaJerk        process jerk noise std (units/s^3) for order=2
%       .robustDetrend    (true/false) default true
%       .hpCutHz          detrend highpass cutoff (Hz) default 0.05
%       .nanPolicy        'interp' or 'omit' default 'interp'
%
% OUTPUT (struct)
%   out.x, out.xdot, out.xddot (if order=3)
%   out.x_filt (filtered forward pass), out.x_smooth (RTS)
%   out.innov, out.S (innovation & its variance)
%   out.params (A,H,Q,R,dt,order, estimated sigmas)

arguments
    obj {mustBeNonempty}
    config.order (1,1) double {mustBeMember(config.order,[2 3])} = 3
    config.useRTS (1,1) logical = true
    config.initPScale (1,1) double {mustBePositive} = 10
    config.sigmaMeas (1,1) double = NaN
    config.sigmaAcc (1,1) double = NaN
    config.sigmaJerk (1,1) double = NaN
    config.robustDetrend (1,1) logical = true
    config.hpCutHz (1,1) double {mustBeNonnegative} = 0.05
    config.nanPolicy (1,:) char {mustBeMember(config.nanPolicy,{'interp','omit'})} = 'interp'
end


x = obj.HrSignal(:);
fs = obj.fs_radar;
N = numel(x);
dt = 1/fs;

% --- NaN handling ---
nanMask = isnan(x);
if any(nanMask)
    switch config.nanPolicy
        case 'interp'
            t = (0:N-1)'*dt;
            x(nanMask) = interp1(t(~nanMask), x(~nanMask), t(nanMask), 'linear', 'extrap');
        case 'omit'
            % For 'omit', we keep NaNs and skip measurement update when NaN
            % (implemented below)
    end
end

% --- Robust detrend / drift removal (VERY important at 2000 Hz) ---
x0 = x;
if config.robustDetrend
    % remove slow drift with a gentle highpass implemented via moving median
    % window length ~ 1/config.hpCutHz seconds (clamped)
    wSec = max(5, 1/max(config.hpCutHz,1e-3)); % seconds
    w = max(201, 2*floor((wSec*fs)/2)+1);      % odd, >=201
    trend = movmedian(x, w, 'omitnan');
    x = x - trend;
end

% --- Estimate measurement noise if not provided ---
% Use robust estimate on high-frequency component: diff(x) is dominated by noise at high fs.
if isnan(config.sigmaMeas)
    dx = diff(x);
    sigmaMeas = mad(dx,1) / sqrt(2);  % robust, approx for white noise
    if sigmaMeas == 0 || ~isfinite(sigmaMeas)
        sigmaMeas = std(dx)/sqrt(2) + eps;
    end
else
    sigmaMeas = config.sigmaMeas;
end

% --- Build state-space model ---
order = config.order;
H = zeros(1, order); H(1) = 1; % measure position only

if order == 2
    % State: [x; xdot]
    A = [1 dt;
         0 1];

    % Process noise driven by "jerk" (random accel changes) on xdot
    % Continuous-time white noise on acceleration integrated gives Q ~ sigmaJerk^2 * [...]
    if isnan(config.sigmaJerk)
        % crude default: set jerk so that model is moderately flexible
        sigmaJerk = 10*sigmaMeas / (dt^2); % heuristic
    else
        sigmaJerk = config.sigmaJerk;
    end
    q = sigmaJerk^2;
    Q = q * [dt^3/3, dt^2/2;
             dt^2/2, dt];

elseif order == 3
    % State: [x; xdot; xddot]
    A = [1 dt 0.5*dt^2;
         0 1  dt;
         0 0  1];

    % Process noise driven by "acceleration random walk" (white noise on jerk)
    if isnan(config.sigmaAcc)
        % default: allow moderate accel changes; scale from sigmaMeas
        sigmaAcc = 5*sigmaMeas / (dt^2); % heuristic (units/s^2)
    else
        sigmaAcc = config.sigmaAcc;
    end
    q = sigmaAcc^2;
    % Standard discrete Q for constant-acceleration model with white noise on acceleration
    Q = q * [dt^5/20, dt^4/8,  dt^3/6;
             dt^4/8,  dt^3/3,  dt^2/2;
             dt^3/6,  dt^2/2,  dt];
end

R = sigmaMeas^2;

% --- Initialize ---
xhat = zeros(order, N);
P = eye(order) * config.initPScale;

% Initialize state from first samples
xhat(1,1) = x(1);
if N >= 2
    xhat(2,1) = (x(2)-x(1))/dt; % crude derivative init
end
if order == 3 && N >= 3
    xhat(3,1) = ((x(3)-x(2)) - (x(2)-x(1))) / (dt^2);
end

innov = nan(1,N);
S = nan(1,N);

% Store forward pass for RTS
xpred_store = zeros(order,N);
Ppred_store = zeros(order,order,N);
P_store     = zeros(order,order,N);

% --- Forward Kalman filter ---
for k = 2:N
    % Predict
    xpred = A * xhat(:,k-1);
    Ppred = A * P * A' + Q;

    xpred_store(:,k) = xpred;
    Ppred_store(:,:,k) = Ppred;

    zk = x0(k);
    hasMeas = ~(isnan(zk) && strcmp(config.nanPolicy,'omit'));

    if hasMeas
        % Use detrended measurement if detrending enabled, else raw
        z = x(k);

        % Update
        y = z - H*xpred;             % innovation
        Sk = H*Ppred*H' + R;         % innovation variance (scalar)
        K = (Ppred*H') / Sk;         % Kalman gain

        xnew = xpred + K*y;
        Pnew = (eye(order) - K*H) * Ppred;

        innov(k) = y;
        S(k) = Sk;
    else
        xnew = xpred;
        Pnew = Ppred;
    end

    xhat(:,k) = xnew;
    P = Pnew;
    P_store(:,:,k) = P;
end

% --- RTS smoother (optional) ---
xs = xhat;
Ps = P_store;
if config.useRTS
    for k = N-1:-1:1
        Pk = P_store(:,:,k);
        Ppred_next = Ppred_store(:,:,k+1);
        if all(Ppred_next(:)==0)
            continue
        end
        G = (Pk*A') / Ppred_next; % smoother gain
        xs(:,k) = xhat(:,k) + G*(xs(:,k+1) - xpred_store(:,k+1));
        Ps(:,:,k) = Pk + G*(Ps(:,:,k+1) - Ppred_next)*G';
    end
end

% --- Outputs ---
out = struct();
out.params = struct('A',A,'H',H,'Q',Q,'R',R,'dt',dt,'fs',fs,'order',order, ...
                    'sigmaMeas',sigmaMeas);

out.x_filt = xhat(1,:).';
out.x_smooth = xs(1,:).';
out.innov = innov(:);
out.S = S(:);

if order >= 2
    out.xdot_filt = xhat(2,:).';
    out.xdot_smooth = xs(2,:).';
end
if order == 3
    out.xddot_filt = xhat(3,:).';
    out.xddot_smooth = xs(3,:).';
end
obj.KF_HrSignal = out.x_filt(:);
end
