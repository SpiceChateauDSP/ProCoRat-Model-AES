% RatGain_DK.m
clear; clc; close all;

% Utility
Fs = 48000;
Ts = 1/Fs;

% Input/Output
Vi = [1; zeros(Fs/2 - 1, 1)];
N = length(Vi);
Vo = zeros(N, 1);

for distortion = 0:.1:1
    if (distortion == 0)
        distortion = 0.00001;
    elseif (distortion == 1)
        distortion = 0.99999;
    end
    
    % Components
    C1 = 100e-12;
    R1 = Ts / (2 * C1);
    C2 = 2.2e-6;
    R2 = Ts / (2 * C2);
    C3 = 4.7e-6;
    R3 = Ts  / (2 * C3);
    R4 = 47;
    R5 = 560;
    P1 = distortion * 100e3;

    %%% Gain Bandwidth Product
    fc = ((1-distortion) * 23e3) + 500;
    WnGBP = fc/24000;
    % Find Filter Coefficients
    [b,a] = butter(1,WnGBP);

    % G Substitutions
    G1 = 1/P1 + 1/R1;
    G2 = G1 + 1/R4 + 1/R5;
    G3 = 1/R2 + 1/R4;
    G4 = 1/R3 + 1/R5;

    % States
    x1 = 0;
    x2 = 0;
    x3 = 0;
    
    % DSP
    for n = 1 : N
        % Node Calculations
        Va = Vi(n, 1)/(R4 * G3) + x2/G3;
        Vb = Vi(n, 1)/(R5 * G4) + x3/G4;

        % Transfer Function
        Vo(n, 1) = (Vi(n,1 ) * G2)/G1 - Va/(R4 * G1) - Vb/(R5 * G1) - x1/G1;

        % State Update
        x1 = (2/R1) * (Vi(n, 1) - Vo(n, 1)) - x1;
        x2 = (2/R2) * Va - x2;
        x3 = (2/R3) * Vb - x3;
    end

    % Graph
    [H, W] = freqz(Vo, 1, N, Fs);
    semilogx(W, 20*log10(abs(H)));
    axis([20 26000 -3 70]);
    hold on;
end
hold off;