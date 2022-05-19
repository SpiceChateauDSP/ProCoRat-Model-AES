% RatTone_DK.m
clear; clc; close all;

% Utility
Fs = 4*48000;
Ts = 1 / Fs;

% Input/Output
Vi = [1; zeros(Fs/2 - 1, 1)];
N = length(Vi);
Vo = zeros(N, 1);

for tone = 0:.1:1
    % Components
    C1 = 3.3e-9;
    C2 = 22e-9;
    R1 = Ts / (2 * C1);
    R2 = Ts / (2 * C2);
    R3 = 1e6;
    P = tone*1e5+1.5e3; % Tone Knob

    % G Substitutions
    G1 = 1/P + 1/R1 + 1/R2;
    G2 = 1 + R2/R3;
    G3 = G2 - 1/(R2*G1);
    G4 = 1/G1 - R2;

    % States
    x1 = 0;
    x2 = 0;

    % Coefficients
    b0 = 1/(P*G1*G3);
    b1 = 1/(G1*G3);
    b2 = G4/G3;

    % DSP
    for n = 1 : N
        % Transfer Function
        Vo(n,1) = Vi(n,1)*b0 + x1*b1 + x2*b2;
        % Nodes
        Va = Vo(n,1)*G2 + x2*R2;
        % State Updates
        x1 = (2/R1)*Va - x1;
        x2 = (2/R2)*(Va-Vo(n,1)) - x2;
    end

    % Graph
    [H, W] = freqz(Vo, 1, N, Fs);
    semilogx(W, 20*log10(abs(H)));
    axis([0 24000 -60 6]);
    hold on;
end
hold off;