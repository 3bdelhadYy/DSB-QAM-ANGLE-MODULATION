clear; clc;

samplingFreq = 1e6;
Time = 4e-3; % snap shot of the signal

% sampling the Time required into segments based on the sampling frequency
% for it to be discrete
timeVector = 0:1/samplingFreq:Time;

%freq of the first message
TP1 = 1e-3; % Time period for one cycle
f1 = 1/TP1;

%=====================(1)=====================

% sawtooth function:
% - inverted
% - phase shiftied forward about 180 degrees or pi rads
m1 = - sawtooth((2*pi*f1*timeVector) - pi);

figure;
plot(timeVector * 1e3, m1, 'Color', 'r', 'LineWidth', 2);
xlabel('Time (ms)');
ylabel('M1(t)');
title('Sawtooth Message Signal M1');
grid on;
%savefig('m_1(t).fig');

% time stamps where the value of M2 change
stampsM2 = [0, 1e-6, 0.5, 1, 1.5, 2];
% values of M2 at each time stamp
m2 = [0, 1, 0.5, -0.5, -1, 0];

% New time vector with 4001 points
t_expanded = linspace(stampsM2(1), stampsM2(end), 4001);

% Interpolate using step-wise hold (previous value)
M2 = interp1(stampsM2, m2, t_expanded, 'previous');

% Plot to visualize
figure;
stairs(timeVector * 1e3, M2, 'Color', 'b', 'LineWidth', 2);
xlabel('Time (ms)');
ylabel('M2(t)');
title('Step-like Message Signal M2');
grid on;
%savefig('m_2(t).fig');

%=====================(2)=====================

%Carrier signals:
amp = 5;
carrierFreq = 5e3;
% define
inPhase = amp * cos(2 * pi * carrierFreq * timeVector);
quadcarrier = amp * sin(2 * pi * carrierFreq * timeVector);

%plot signals
figure;
plot(timeVector * 1e3, quadcarrier, 'g', 'LineWidth', 1.5); hold on;
plot(timeVector * 1e3, inPhase, 'b', 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Amplitude');
title('Carrier Signals');
legend('Quadrature (sin)', 'In-phase (cos)');
xlim([0 1]);
grid on;
%savefig('Carrier_signals.fig');

modulated = (inPhase .* m1) + (quadcarrier .* M2);

%plotting the modulated signal
figure;
plot(timeVector * 1e3, modulated, 'm', 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Amplitude');
title('QAM Modulated Signal');
grid on;
%savefig('QAM_Modulated_Signal.fig');

%=====================(3)=====================

%Demodulating:
% Multiply with carriers
demod_I = modulated .* inPhase;
demod_Q = modulated .* quadcarrier;

% Low-pass filtering
cutoff_freq = 4e3; % 10kHz cutoff (2x message bandwidth)
filter_order = 6;
[b, a] = butter(filter_order, cutoff_freq/(samplingFreq/2), 'low');

% Zero-phase filtering
filtered_I = filtfilt(b, a, demod_I);
filtered_Q = filtfilt(b, a, demod_Q);

% Amplitude compensation
scale_factor = (amp^2)/2; % From: (A·cosθ × A·cosθ) -> A²/2
recovered_m1 = (2/amp^2) * filtered_I;
recovered_m2 = (2/amp^2) * filtered_Q;

%%Original vs Recovered
figure;
subplot(2,1,1);
plot(timeVector * 1e3, m1, 'r', 'LineWidth', 1.5); hold on;
plot(timeVector * 1e3, recovered_m1, 'w', 'LineWidth', 1);
title('M1 Recovery');
xlabel('Time (ms)');
ylabel('Amplitude');
legend('Original M1', 'Recovered M1');
grid on;

subplot(2,1,2);
stairs(timeVector * 1e3, M2, 'b', 'LineWidth', 1.5); hold on;
plot(timeVector * 1e3, recovered_m2, 'w', 'LineWidth', 1);
title('M2 Recovery');
xlabel('Time (ms)');
ylabel('Amplitude');
legend('Original M2', 'Recovered M2');
grid on;
%savefig('Original_vs_Recovered_signal.fig');

%=====================(4)=====================

%Shifted Phase at receiver -> cos(2πfct + π/3)
shifted_Inphase = amp * cos(2*pi*carrierFreq*timeVector + pi/3);
shifted_Quad = amp * sin(2*pi*carrierFreq*timeVector + pi/3);

%phase error modulation
demod_I_shifted = modulated .* shifted_Inphase;
demod_Q_shifted = modulated .* shifted_Quad;

%Filtering
filtered_I_shifted = filtfilt(b, a, demod_I_shifted);
filtered_Q_shifted = filtfilt(b, a, demod_Q_shifted);

%recoverd signals
recovered_m1_shifted = (2 / amp^2) * filtered_I_shifted;
recovered_m2_shifted = (2 / amp^2) * filtered_Q_shifted;

figure;
subplot(2,1,1);
plot(timeVector * 1e3, m1, 'r', 'LineWidth', 1.5); hold on;
plot(timeVector * 1e3, recovered_m1, 'w'); % Ideal
plot(timeVector * 1e3, recovered_m1_shifted, 'b'); % π/3 phase shift
title('M1 Recovery: Ideal vs Shifted vs Offset');
legend('Original M1', 'Ideal', 'Phase Shift π/3');
grid on;

subplot(2,1,2);
stairs(timeVector * 1e3, M2, 'b', 'LineWidth', 1.5); hold on;
plot(timeVector * 1e3, recovered_m2, 'w'); % Ideal
plot(timeVector * 1e3, recovered_m2_shifted, 'r'); % π/3 phase shift
title('M2 Recovery: Ideal vs Shifted vs Offset');
legend('Original M2', 'Ideal', 'Phase Shift π/3');
grid on;
%savefig('Shifted_Phase_at_receiver.fig');

%=====================(5)=====================

%Frequency offset -> cos(2.02πfct)
freqOffset_I = amp * cos(2.02*pi*carrierFreq*timeVector);
freqOffset_Q = amp * sin(2.02*pi*carrierFreq*timeVector);

%Demodulation
demod_I_freqOffset = modulated .* freqOffset_I;
demod_Q_freqOffset = modulated .* freqOffset_Q;

%Filter
filtered_I_freqOffset = filtfilt(b, a, demod_I_freqOffset);
filtered_Q_freqOffset = filtfilt(b, a, demod_Q_freqOffset);

%Recover
recovered_m1_freqOffset = (2 / amp^2) * filtered_I_freqOffset;
recovered_m2_freqOffset = (2 / amp^2) * filtered_Q_freqOffset;

%Plot
figure;
subplot(2,1,1);
plot(timeVector * 1e3, m1, 'r', 'LineWidth', 1.5); hold on;
plot(timeVector * 1e3, recovered_m1, 'w'); % Ideal
plot(timeVector * 1e3, recovered_m1_freqOffset, 'g'); % freq offset
title('M1 Recovery: Ideal vs Shifted vs Offset');
legend('Original M1', 'Ideal', 'Freq Offset');
grid on;

subplot(2,1,2);
stairs(timeVector * 1e3, M2, 'b', 'LineWidth', 1.5); hold on;
plot(timeVector * 1e3, recovered_m2, 'w'); % Ideal
plot(timeVector * 1e3, recovered_m2_freqOffset, 'g'); % freq offset
title('M2 Recovery: Ideal vs Shifted vs Offset');
legend('Original M2', 'Ideal', 'Freq Offset');
grid on;
%savefig('Frequency_offset.fig');
