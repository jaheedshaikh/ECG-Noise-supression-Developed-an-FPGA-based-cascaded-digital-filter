%% Original signal

clc
clear
close all

fs = 360;                         % Sampling frequency

ecg = read_mitbih('100.dat');     % Read MIT-BIH ECG record

ecg = ecg(1:3600);                % First 10 seconds of signal

t = (0:length(ecg)-1)/fs;         % Time axis

figure
subplot(1,1,1)
plot(t, ecg, 'LineWidth', 1.25)    % Line width applied here
title('Original MIT-BIH ECG')
xlabel('Time (seconds)')
ylabel('Amplitude')
legend('ECG', 'Location', 'northeast')
grid on

%% Figure 2: ECG signal with added baseline wander and 50 Hz interference and emg noise for testing 

% Normalize ECG
ecg = ecg / max(abs(ecg));

%power_noise   = 0.02*sin(2*pi*50*t)';

power_noise = 0.03*sin(2*pi*50*t)';
baseline_noise = 0.03*sin(2*pi*0.3*t)';
emg_noise      = 0.01*randn(length(ecg),1);

total_noise = power_noise + baseline_noise + emg_noise;

ecg_noisy = ecg + total_noise;


figure

subplot(2,1,1)
plot(t,ecg,'LineWidth',1.25)
title('Normalized ECG')
xlabel('Time (seconds)')
ylabel('Amplitude')
grid on

subplot(2,1,2)
plot(t,ecg_noisy)
title('ECG with Noise')
xlabel('Time (seconds)')
ylabel('Amplitude')
grid on

%% 3.Baseline Wander Removal, Moving average high-pass.

fs = 360;   % MIT-BIH sampling frequency

% Window size (about 0.6 sec recommended)
window_size = round(0.6*fs);

% Moving average filter (baseline estimation)
baseline_est = movmean(ecg_noisy, window_size);

% High-pass output
ecg_baseline_removed = ecg_noisy - baseline_est;

figure

subplot(3,1,1)
plot(ecg_noisy)
title('Noisy ECG')
xlabel('Samples')
ylabel('Amplitude')
grid on

subplot(3,1,2)
plot(baseline_est)
title('Estimated Baseline (Moving Average)')
xlabel('Samples')
ylabel('Amplitude')
grid on

subplot(3,1,3)
plot(ecg_baseline_removed)
title('ECG after Baseline Wander Removal')
xlabel('Samples')
ylabel('Amplitude')
grid on


%% 4. Powerline Noise Removal (50 Hz Notch Filter)

f0 = 50;        % Powerline frequency
bw = 1.1;         % Bandwidth around 50 Hz

low  = (f0-bw)/(fs/2);   % Lower normalized frequency
high = (f0+bw)/(fs/2);   % Upper normalized frequency

[b,a] = butter(2,[low high],'stop');   % Band-stop filter

ecg_notch = filtfilt(b,a,ecg_baseline_removed);

figure

subplot(2,1,1)
plot(t, ecg_baseline_removed)
title('ECG after Baseline Removal')
xlabel('Time (seconds)')
ylabel('Amplitude')
grid on

subplot(2,1,2)
plot(t, ecg_notch)
title('ECG after 50 Hz Powerline Removal')
xlabel('Time (seconds)')
ylabel('Amplitude')
grid on

%% Frequency Spectrum Comparison (Before and After Notch Filter)

% In the top spectrum plot: A clear peak near 50 Hz (powerline interference)

% In the bottom spectrum plot:The 50 Hz peak should be reduced or disappear

figure
N = length(ecg_baseline_removed);

f = (0:N-1)*(fs/N);

Y1 = abs(fft(ecg_baseline_removed));
Y2 = abs(fft(ecg_notch));

subplot(2,1,1)
plot(f,Y1)
xlim([0 100])
title('Spectrum Before Notch Filter')

subplot(2,1,2)
plot(f,Y2)
xlim([0 100])
title('Spectrum After Notch Filter')


%% 5. FIR Low Pass Filter to Remove EMG Noise

cutoff = 45;                     % ECG useful frequency range (~0–40 Hz)
order = 50;                     % Filter order

% Normalized cutoff frequency
fc = cutoff/(fs/2);

% Design FIR filter
b = fir1(order, fc, 'low');

% Apply filter
ecg_clean = filtfilt(b,1,ecg_notch);

figure

subplot(3,1,1)
plot(t, ecg,'LineWidth',1.2)
title('Original ECG for Comparison')
xlabel('Time (seconds)')
ylabel('Amplitude')
grid on

subplot(3,1,2)
plot(t, ecg_notch)
title('ECG after 50 Hz Powerline Removal')
xlabel('Time (seconds)')
ylabel('Amplitude')
grid on

subplot(3,1,3)
plot(t, ecg_clean)
title('ECG after FIR Low Pass Filter (EMG Removed)')
xlabel('Time (seconds)')
ylabel('Amplitude')
grid on


%% 6. Savitzky-Golay Filter, Smooth ECG

sgolay_order = 3;      % Polynomial order
frame_len = 11;        % Frame length (must be odd)

ecg_smooth = sgolayfilt(ecg_clean, sgolay_order, frame_len);

figure

subplot(3,1,1)
plot(t, ecg_clean,'LineWidth',1.25)
title('ECG after FIR Low Pass Filter')
xlabel('Time (seconds)')
ylabel('Amplitude')
grid on

subplot(3,1,2)
plot(t, ecg_smooth,'LineWidth',1.25)
title('ECG after Savitzky-Golay Smoothing')
xlabel('Time (seconds)')
ylabel('Amplitude')
grid on

subplot(3,1,3)
plot(t, ecg,'LineWidth',1.2)
title('Original ECG')
xlabel('Time (seconds)')
ylabel('Amplitude')
grid on

%% Final ECG Alignment (Shift signal to original level)

ecg_final = ecg_smooth + mean(ecg);

%% Comparison: Original vs Noisy vs Final ECG

figure

subplot(3,1,1)
plot(t, ecg,'LineWidth',1.2)
title('Original ECG')
xlabel('Time (seconds)')
ylabel('Amplitude')
grid on

subplot(3,1,2)
plot(t, ecg_noisy)
title('Noisy ECG (Baseline + Powerline + EMG)')
xlabel('Time (seconds)')
ylabel('Amplitude')
grid on

subplot(3,1,3)
plot(t, ecg_final,'LineWidth',1.2)
title('Final Denoised ECG')
xlabel('Time (seconds)')
ylabel('Amplitude')
grid on

%% COMPUTATION 

%% ===============================
% Performance Metrics Computation
% ===============================

% Noise before filtering
noise_before = ecg_noisy - ecg;

% Noise after filtering
noise_after = ecg_final - ecg;

% Signal Power
signal_power = sum(ecg.^2);

% Noise Power
noise_power_before = sum(noise_before.^2);
noise_power_after  = sum(noise_after.^2);

% SNR before filtering
SNR_before = 10*log10(signal_power / noise_power_before);

% SNR after filtering
SNR_after = 10*log10(signal_power / noise_power_after);

% SNR Improvement
SNR_improvement = SNR_after - SNR_before;

% Mean Square Error
MSE = mean((ecg - ecg_final).^2);

% Root Mean Square Error
RMSE = sqrt(MSE);

% Correlation Coefficient
CC = corr(ecg - mean(ecg), ecg_final - mean(ecg_final));
% Display results
fprintf('\n===== Performance Metrics =====\n');
fprintf('SNR Before Filtering : %.2f dB\n', SNR_before);
fprintf('SNR After Filtering  : %.2f dB\n', SNR_after);
fprintf('SNR Improvement      : %.2f dB\n', SNR_improvement);
fprintf('MSE                  : %.6f\n', MSE);
fprintf('RMSE                 : %.6f\n', RMSE);
fprintf('Correlation Coefficient (CC) : %.4f\n', CC);


%% ===============================
% Export for FPGA (Verilog)
% ===============================

scale = 2^14;

% ----- FIR coefficients -----
b_fir = fir1(order, fc, 'low');   % ensure correct variable
b_fixed = round(b_fir * scale);
writematrix(b_fixed, 'fir_coefficients.txt');

% ----- Notch filter coefficients -----
[b_notch,a_notch] = butter(2,[low high],'stop');

b_notch_fixed = round(b_notch * scale);
a_notch_fixed = round(a_notch * scale);

writematrix(b_notch_fixed, 'notch_b.txt');
writematrix(a_notch_fixed, 'notch_a.txt');

% ----- Input ECG (noisy signal) -----
ecg_fixed = round(ecg_noisy * scale);
writematrix(ecg_fixed, 'ecg_input.txt');

% ----- Expected output (filtered ECG) -----
ecg_out_fixed = round(ecg_final * scale);
writematrix(ecg_out_fixed, 'ecg_expected.txt');

disp('Export completed for FPGA implementation');