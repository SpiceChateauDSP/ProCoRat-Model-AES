% RatInput_Simple_DK.m
clear; clc; close all;

% Utility
Fs = 48000;
Ts = 1/Fs;

% Input/Output
Vi = [1;zeros(Fs/2 - 1, 1)];
N = length(Vi);
Vo = zeros(N,1);
Ps = 0; % Power Supply

% Components
C1 = 22e-9;
C2 = 1e-9;
R1 = Ts / (2 * C1);
R2 = Ts / (2 * C2);
R3 = 1e6;
R4 = 1e3;

% Calculating substitutions
G1 = 1/R1 + 1/R3 + 1/R4;
G2 = R4/R2 + 1;
G3 = G2 - 1/(R4 * G1);

% States
x1 = 0;
x2 = 0;

% Coefficients
b0 = 1/(R1 * G1 * G3);
b1 = -1/(G1 * G3);
b2 = R4/G3;

% DSP
for n = 1 : N
    % Transfer Function
    Vo(n, 1) = Vi(n, 1) * b0 + x1 * b1 + x2 * b2;
    
    % State Updates
    Va = Vo(n, 1) * G2 - x2 * R4;
    x1 = (2 / R1) * (Vi(n, 1) - Va) - x1;
    x2 = (2 / R2) * Vo(n, 1) - x2;
end

% Graph
[H, W] = freqz(Vo, 1, N, Fs);
semilogx(W, 20*log10(abs(H)));
axis([0 24000 -6 6]);