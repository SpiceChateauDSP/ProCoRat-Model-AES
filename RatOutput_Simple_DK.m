% RatOutput_Simple_DK.m
clear; clc; close all;

% Utility
Fs = 48000;
Ts = 1 / Fs;

% Input/Output
Vi = [1; zeros(Fs/2 - 1, 1)];
N = length(Vi);
Vo = zeros(N, 1);

for volume = 0:.1:1
    volume = volume*98e3 + 1e3;

    % Components
    C1 = 1e-6;
    R1 = Ts / (2 * C1);
    R2 = volume;
    R3 = (99e3 - volume);

    % G Substitutions
    G1 = 1/R1 + 1/R2;
    G2 = 1 + R2/R3;
    G3 = G2 - 1/(G1*R2);

    % States
    x1 = 0;

    % Coefficients
    b0 = 1/(G1*G3*R1);
    b1 = -1/(G1*G3);

    % DSP
    for n = 1 : N
        % Transfer Function
        Vo(n, 1) = Vi(n, 1) * b0 + x1 * b1;

        % State Updates
        Va = Vo(n,1)*G2;
        x1 = (2/R1) * (Vi(n, 1) - Va) - x1;
    end

    % Graph
    [H, W] = freqz(Vo, 1, N, Fs);
    semilogx(W, 20*log10(abs(H)));
    axis([0 24000 -60 6]);
    hold on;
end
hold off;