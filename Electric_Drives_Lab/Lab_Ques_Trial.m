%% Your system in the 'design skeleton' style (margins + step metrics)
clear; clc; close all;
s = tf('s');

% ---- 1) Plant blocks ----

gs = 0.76012*(2.4)/(1+1.054*(s/46.35)+(s/46.35)*(s/46.35));
hs = 5481.590626/(s/0.004008+1);

Gp = gs * hs;                   % treat this as open-loop plant L(s) if C(s)=1

% ---- 2) Controller (identity here) ----
C = 1;                          % no extra compensation; keep exactly your loop
% ---- 3) Loop and closed-loop ----
L  = C * Gp;                    % open-loop
T  = feedback(L, 1);            % unity negative feedback closed-loop

% ---- 4) zpk display in frequency format ----
L_zpk  = zpk(L);  L_zpk.DisplayFormat  = 'frequency';
T_zpk  = zpk(T);  T_zpk.DisplayFormat  = 'frequency';
disp('Closed-loop (T) as ZPK:');
disp(T_zpk);

% ---- 5) Bode with margins for L ----
figure('Name','Open-Loop Bode with Margins');
margin(L); grid on;
try
    rp = gcr; opt = getoptions(rp);
    opt.MagUnits  = 'dB';
    opt.FreqUnits = 'rad/s';
    opt.Grid      = 'on';
    setoptions(rp,opt);
catch, end
title('Open-Loop L(s) = gs \times hs (Bode with GM/PM markers)');

% ---- 6) Also show Bode of closed-loop T ----
figure('Name','Closed-Loop Bode');
bode(T); grid on;
try
    rp = gcr; opt = getoptions(rp);
    opt.MagUnits  = 'dB';
    opt.FreqUnits = 'rad/s';
    opt.Grid      = 'on';
    setoptions(rp,opt);
catch, end
title('Closed-Loop T(s) = L/(1+L) Bode');

% ---- 7) Print ALL margins ----
AM = allmargin(L);
fprintf('\n===== ALL Stability Margins (Open-Loop L) =====\n');
if ~isempty(AM.GainMargin)
    for k=1:numel(AM.GainMargin)
        gm = AM.GainMargin(k); w = AM.GMFrequency(k);
        if isinf(gm), gmdb = Inf; else gmdb = 20*log10(gm); end
        fprintf('Gain Margin #%d: %s dB at ω = %.4g rad/s\n', ...
            k, num2str(gmdb), w);
    end
else
    fprintf('No finite Gain Margins.\n');
end
if ~isempty(AM.PhaseMargin)
    for k=1:numel(AM.PhaseMargin)
        pm = AM.PhaseMargin(k); w = AM.PMFrequency(k);
        fprintf('Phase Margin #%d: %.2f° at ω = %.4g rad/s\n', k, pm, w);
    end
else
    fprintf('No Phase Margins.\n');
end
if isfield(AM,'DelayMargin') && ~isempty(AM.DelayMargin)
    fprintf('Delay Margin: %.4g s at ω = %.4g rad/s\n', AM.DelayMargin, AM.DMFrequency);
end

% ---- 8) Closed-loop step + metrics ----
figure('Name','Closed-Loop Step');
[y,t] = step(T); plot(t,y,'LineWidth',1.6); grid on;
xlabel('Time (s)'); ylabel('Amplitude'); title('Step: Closed-Loop T(s)');
info = stepinfo(T);
ess  = abs(1 - dcgain(T));

fprintf('\n===== Closed-Loop Step Metrics =====\n');
fprintf('Rise Time     : %.6g s\n', info.RiseTime);
fprintf('Settling Time : %.6g s\n', info.SettlingTime);
fprintf('Overshoot     : %.2f %%\n', info.Overshoot);
fprintf('Peak          : %.6g at t = %.6g s\n', info.Peak, info.PeakTime);
fprintf('Steady-State Error to 1-step: %.6g\n', ess);


figure('Name','Overview (L and T)');
subplot(2,2,1); margin(L); grid on; title('Open-Loop L: Bode+Margins');
subplot(2,2,2); bode(T);   grid on; title('Closed-Loop T: Bode');
subplot(2,2,[3 4]); step(T); grid on; title('Closed-Loop T: Step');

