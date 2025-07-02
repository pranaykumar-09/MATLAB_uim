clc;
clear all;
close all;
pkg load signal;

%% Load real-only signal data from .dat file
raw = load('custom_real_signal.dat');  % Replace with your actual .dat file

% Ensure it is a column vector
x = raw';  % real-valued samples
m=mean(x);

%% Convert to complex analytic signal using Hilbert transform
 % I + jQ approximation

%% Parameters
fs = 1000;      % Sampling frequency (adjust as needed)
N = length(x);
for i=1:1:N
  x(i)=x(i)-m;
end
analytic_signal = hilbert(x);
Nx = 100;       % Window size
step = 1;
NT = N - Nx + 1; % Number of time segments
p_order = 1;    % Polynomial order for curve fitting

IP_est = zeros(1, NT);  % Instantaneous phase per window

%% Sliding window + modified Kay method + MLE phase estimation
for idx = 1:step:NT
    seg = analytic_signal(idx : idx + Nx - 1);

    % Step 1: Unwrap instantaneous phase
    phase_unwrapped = unwrap(angle(seg));

    % Step 2: First-order phase difference
    D = diff(phase_unwrapped).';

    % Step 3: Build covariance matrix C and optimal weight vector
    H = ones(Nx - 1, 1);
    C = 2*eye(Nx - 1) - diag(ones(Nx - 2, 1), 1) - diag(ones(Nx - 2, 1), -1);
    C_inv = inv(C);
    w_vector = C_inv * H;  % Optimal weights

    % Step 4: Estimate digital angular frequency (omega)
    omega_est = sum(w_vector .* D) / sum(w_vector);
    f_est = omega_est / (2*pi) * fs;

    % Step 5: MLE estimation of instantaneous phase
    n = 0:Nx-1;
    re_part = real(seg);
    num_phase = sum(re_part .* sin(2*pi*f_est/fs * n));
    den_phase = sum(re_part .* cos(2*pi*f_est/fs * n));
    IP_est(idx) = atan2(num_phase, den_phase);
end

%% Fit polynomial trend and extract UIMIP
t_IP = (0:NT-1) / fs;
IP_unwrapped = unwrap(IP_est);

% Polynomial fit of intentional modulation trend
coeffs = polyfit(t_IP, IP_unwrapped, p_order);
IP_fit = polyval(coeffs, t_IP);

% Residual = UIMIP
UIMIP = IP_unwrapped - IP_fit;

%% Plot results
figure;
subplot(3,1,1);
plot(t_IP, IP_unwrapped, 'b');
title('Estimated Instantaneous Phase (IP)');
xlabel('Time (s)');
ylabel('Phase (rad)');
grid on;

subplot(3,1,2);
plot(t_IP, IP_fit, 'r');
title('Fitted Trend (Intentional Phase)');
xlabel('Time (s)');
ylabel('Phase (rad)');
grid on;

subplot(3,1,3);
plot(t_IP, UIMIP, 'k');
title('Extracted UIMIP (Residual)');
xlabel('Time (s)');
ylabel('Phase (rad)');
grid on;

