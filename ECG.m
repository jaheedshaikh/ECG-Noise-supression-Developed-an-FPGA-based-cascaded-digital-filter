clc; clear; close all;

%% ========================================================================
%  STEP 1: Generate Synthetic ECG (Realistic Morphology)
%  ========================================================================
Fs = 360;                       
t = 0:1/Fs:5;                   

% Synthetic ECG components
hr_freq = 1.2; % ~72 BPM
r_peaks = 1.5 * exp(-((mod(t, 1/hr_freq) - 0.41/hr_freq)/(0.015)).^2);
p_wave  = 0.2 * exp(-((mod(t, 1/hr_freq) - 0.35/hr_freq)/(0.03)).^2);
t_wave  = 0.4 * exp(-((mod(t, 1/hr_freq) - 0.60/hr_freq)/(0.05)).^2);
ecg_clean = r_peaks + p_wave + t_wave;

%% ========================================================================
%  STEP 2: Add Realistic Noise
%  ========================================================================
bw_noise = 0.6 * sin(2*pi*0.3*t);        % Baseline Wander
pli_noise = 0.35 * sin(2*pi*50*t);       % 50 Hz Powerline
emg_noise = 0.12 * randn(size(t));       % High-Freq Muscle Noise

ecg_noisy = ecg_clean + bw_noise + pli_noise + emg_noise;

%% ========================================================================
%  STEP 3: Multi-Stage Denoising (Manual Implementation)
%  ========================================================================

% --- Stage A: Baseline Removal (Simple Moving Average High-Pass) ---
% Using a larger window (L) to capture the slow drift
L = round(0.5 * Fs); 
base_est = conv(ecg_noisy, ones(1,L)/L, 'same');
ecg_detrended = ecg_noisy - base_est;

% --- Stage B: 50 Hz Notch Filter (Manual Design) ---
f0 = 50;                % Notch frequency
w0 = 2*pi*f0/Fs;        % Angular frequency
r = 0.96;               % Pole radius (closer to 1 = narrower notch)

% Transfer Function: H(z) = (1 - 2cos(w0)z^-1 + z^-2) / (1 - 2rcos(w0)z^-1 + r^2z^-2)
b_notch = [1, -2*cos(w0), 1];
a_notch = [1, -2*r*cos(w0), r^2];
ecg_notched = filter(b_notch, a_notch, ecg_detrended);

% --- Stage C: Low-Pass Filter (Simple Smoothing) ---
% A 5-point moving average for final high-freq cleanup
b_lp = ones(1,5)/5;
ecg_final = filter(b_lp, 1, ecg_notched);

%% ========================================================================
%  STEP 4: Visualization
%  ========================================================================
figure('Color', 'w', 'Position', [100, 100, 900, 600]);

subplot(3,1,1);
plot(t, ecg_noisy, 'r'); hold on;
plot(t, ecg_clean, 'k', 'LineWidth', 1.2);
title('Input: Noisy ECG vs. Ground Truth');
legend('Noisy', 'Original'); grid on;

subplot(3,1,2);
plot(t, ecg_detrended, 'Color', [0.4 0.4 0.8]);
title('After Baseline Removal & Detrending');
grid on;

subplot(3,1,3);
plot(t, ecg_final, 'b', 'LineWidth', 1.5); hold on;
plot(t, ecg_clean, 'k--', 'LineWidth', 1);
title('Final Denoised Signal');
xlabel('Time (s)'); legend('Denoised', 'Original Target');
grid on;

% Performance check
snr_val = 10 * log10(var(ecg_clean) / var(ecg_clean - ecg_final));
fprintf('Final SNR Improvement: %.2f dB\n', snr_val);