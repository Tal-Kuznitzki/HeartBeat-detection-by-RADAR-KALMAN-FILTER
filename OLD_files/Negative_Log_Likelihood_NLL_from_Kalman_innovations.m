%% =========================================================
% Offline RR Estimation
% 2-State Kalman (RR + dRR) + RTS Smoother
% Parameter tuning using NLL (no ground truth)
% Robustness to unknown measurement noise R
%% =========================================================
%clear; close all; clc;

set(groot,'defaultAxesFontSize',14);
set(groot,'defaultTextFontSize',16);
set(groot,'defaultLegendFontSize',14);

rng(1);

%% =========================================================
% PARAMETERS
%% =========================================================
N = 240;
RR0 = 1.0;

% --- TRUE vs ASSUMED measurement noise ---
R_rr_true    = (0.0016)^2;   % Used ONLY to generate data (unknown in practice)
R_rr_assumed = (1e-6)^1;   % Used by Kalman filter (imperfect assumption)
R_info_str = sprintf('R_{rr}^{true}=%.3f  \\neq  R_{rr}^{assumed}=%.3f', ...
                     R_rr_true, R_rr_assumed);



MED_WIN = 9;               % Median window (odd)

% Parameter sweep for Kalman acceleration noise
q_acc_vec = logspace(-8, -1, 25);

fprintf('True R_rr = %.4f | Assumed R_rr = %.4f\n', ...
        R_rr_true, R_rr_assumed);

%% =========================================================
% TRUE RR WITH REALISTIC DYNAMICS (FOR EVALUATION ONLY)
%% =========================================================
RR_true = zeros(N,1);
RR_true(1) = RR0;

drift = zeros(N,1);
drift(1:80)    = linspace(0, -0.10, 80);
drift(81:140)  = drift(80);
drift(141:240) = drift(80) + 0.14*((1:100)'/100).^1.5;

Q_true_jitter = (0.01)^2;

for k = 2:N
    RR_det = RR0 + drift(k);
    RR_true(k) = RR_true(k-1) + ...
        (RR_det - (RR0 + drift(k-1))) + sqrt(Q_true_jitter)*randn;
end

RR_true = min(max(RR_true,0.55),1.60);
%%%%
RR_true = dataFull{indx,sz}.RrGtEst;
RR_true = RR_true./vecnorm(RR_true);
%%%%
%% =========================================================
% NOISY MEASUREMENTS (TRUE R_rr USED HERE ONLY)
%% =========================================================
%z = RR_true + sqrt(R_rr_true)*randn(N,1);
%%%%
z = dataFull{indx,sz}.RrEst;
z = z./vecnorm(z);
%%%%

%% =========================================================
% OFFLINE MEDIAN (NON-CAUSAL BASELINE)
%% =========================================================
RR_med = medfilt1(z, MED_WIN, 'truncate');

%% =========================================================
% PARAMETER TUNING LOOP (NLL — NO GROUND TRUTH)
%% =========================================================
nll_vec   = zeros(size(q_acc_vec));
RR_kf_all = cell(length(q_acc_vec),1);

for i = 1:length(q_acc_vec)

    q_acc = q_acc_vec(i);

    % Forward Kalman (uses ASSUMED R_rr)
    [xf, Pf, innov, S] = kalman_forward_2state_nll(z, R_rr_assumed, q_acc);

    % Negative log-likelihood from innovations
    nll_vec(i) = 0.5 * sum( log(S) + (innov.^2)./S );

    % RTS smoother (offline)
    [xs, ~] = rts_smoother_2state(xf, Pf, q_acc);
    RR_kf_all{i} = xs(1,:)';
end

% Optimal parameter (maximum likelihood)
[nll_best, idx_best] = min(nll_vec);
q_acc_best = q_acc_vec(idx_best);
RR_kf_best = RR_kf_all{idx_best};

fprintf('Best q_acc = %.2e | NLL = %.2f\n', q_acc_best, nll_best);

%% =========================================================
% FIGURE 0 — PARAMETER TUNING CURVE (NLL)
%% =========================================================
figure('Color','k','Position',[100 100 700 350]);
hold on; grid on;

ax = gca;
ax.Color = 'k';
ax.XColor = 'w';
ax.YColor = 'w';
ax.GridColor = [0.4 0.4 0.4];
ax.GridAlpha = 0.35;

semilogx(q_acc_vec, nll_vec, 'w-o','LineWidth',1.6);
plot(q_acc_best, nll_best,'ro','MarkerSize',10,'LineWidth',2);

xlabel('q_{acc}','Color','w');
ylabel('Negative Log-Likelihood','Color','w');
title({ ...
    'Kalman Parameter Tuning via Innovation Likelihood', ...
    R_info_str}, 'Color','w');

%% =========================================================
% FIGURE 1 — RR TIME SERIES (FINAL COMPARISON)
%% =========================================================
figure('Color','k','Position',[100 480 1050 360]);
hold on; grid on;

ax = gca;
ax.Color = 'k';
ax.XColor = 'w';
ax.YColor = 'w';
ax.GridColor = [0.4 0.4 0.4];
ax.GridAlpha = 0.35;

plot(RR_true,'g','LineWidth',2.0);
plot(z,'.','Color',[0.6 0.6 0.6]);
plot(RR_kf_best,'r','LineWidth',2.2);
plot(RR_med,'b','LineWidth',1.6);

xlabel('Beat Index','Color','w');
ylabel('RR Interval [s]','Color','w');
title({ ...
    'Offline RR Estimation (RR Domain)', ...
    R_info_str}, 'Color','w');

legend({ ...
    'RR True', ...
    'RR Measured', ...
    sprintf('RR Kalman (q_{acc}=%.1e)', q_acc_best), ...
    'RR Median (non-causal)'}, ...
    'TextColor','w','Location','best','Box','off');



%% =========================================================
% BLAND–ALTMAN (RR DOMAIN)
%% =========================================================
diff_kf = RR_kf_best - RR_true;
mean_kf = (RR_kf_best + RR_true)/2;
bias_kf = mean(diff_kf);
loa_kf  = 1.96 * std(diff_kf);

diff_med = RR_med - RR_true;
mean_med = (RR_med + RR_true)/2;
bias_med = mean(diff_med);
loa_med  = 1.96 * std(diff_med);

figure('Color','k','Position',[100 860 950 360]);
hold on; grid on;

ax = gca;
ax.Color = 'k';
ax.XColor = 'w';
ax.YColor = 'w';
ax.GridColor = [0.4 0.4 0.4];
ax.GridAlpha = 0.35;

% Kalman (red)
h_kf_pts  = scatter(mean_kf, diff_kf, 24, 'filled', ...
    'MarkerFaceColor',[1 0.2 0.2], 'MarkerEdgeColor','none');
h_kf_bias = yline(bias_kf,'-','Color',[1 0.4 0.4],'LineWidth',1.6);
h_kf_loa  = yline(bias_kf+loa_kf,'--','Color',[1 0.6 0.6],'LineWidth',1.2);
yline(bias_kf-loa_kf,'--','Color',[1 0.6 0.6],'LineWidth',1.2);

% Median (light blue)
median_color = [0.3 0.75 1.0];
h_med_pts  = scatter(mean_med, diff_med, 18, 'filled', ...
    'MarkerFaceColor',median_color, 'MarkerEdgeColor','none');
h_med_bias = yline(bias_med,'-','Color',median_color,'LineWidth',1.4);
h_med_loa  = yline(bias_med+loa_med,'--','Color',median_color,'LineWidth',1.1);
yline(bias_med-loa_med,'--','Color',median_color,'LineWidth',1.1);

xlabel('Mean RR [s]','Color','w');
ylabel('RR Difference (Estimate − True) [s]','Color','w');
title({ ...
    'Bland–Altman Analysis (RR Domain)', ...
    R_info_str}, 'Color','w');

legend([h_kf_pts, h_kf_bias, h_kf_loa, h_med_pts, h_med_bias, h_med_loa], ...
    { ...
    'Kalman estimates', ...
    sprintf('Kalman bias = %.4f s', bias_kf), ...
    sprintf('Kalman LoA = [%.4f, %.4f] s', bias_kf-loa_kf, bias_kf+loa_kf), ...
    'Median estimates', ...
    sprintf('Median bias = %.4f s', bias_med), ...
    sprintf('Median LoA = [%.4f, %.4f] s', bias_med-loa_med, bias_med+loa_med) ...
    }, ...
    'TextColor','w','Location','northeast','Box','off');

%% =========================================================
% LOCAL FUNCTIONS
%% =========================================================
function [xf, Pf, innov, S] = kalman_forward_2state_nll(z, R, q_acc)
    N = length(z);
    F = [1 1; 0 1];
    H = [1 0];
    Q = q_acc * [0.25 0.5; 0.5 1];

    xf = zeros(2,N);
    Pf = zeros(2,2,N);
    innov = zeros(N,1);
    S = zeros(N,1);

    xf(:,1) = [z(1); 0];
    Pf(:,:,1) = diag([0.1 0.01]);

    for k = 2:N
        xp = F*xf(:,k-1);
        Pp = F*Pf(:,:,k-1)*F' + Q;

        vk = z(k) - H*xp;
        Sk = H*Pp*H' + R;

        K = (Pp*H')/Sk;
        xf(:,k) = xp + K*vk;
        Pf(:,:,k) = (eye(2) - K*H)*Pp;

        innov(k) = vk;
        S(k) = Sk;
    end

    innov(1) = 0;
    S(1) = H*Pf(:,:,1)*H' + R;
end

function [xs, Ps] = rts_smoother_2state(xf, Pf, q_acc)
    [~,N] = size(xf);
    F = [1 1; 0 1];
    Q = q_acc * [0.25 0.5; 0.5 1];

    xs = xf;
    Ps = Pf;

    for k = N-1:-1:1
        Pp = F*Pf(:,:,k)*F' + Q;
        G  = Pf(:,:,k)*F'/Pp;
        xs(:,k) = xf(:,k) + G*(xs(:,k+1) - F*xf(:,k));
        Ps(:,:,k) = Pf(:,:,k) + G*(Ps(:,:,k+1) - Pp)*G';
    end
end
