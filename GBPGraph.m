clear; clc; close all;

% Utility
Fs = 48000;
Ts = 1/Fs;

% Input/Output
Vi = [1; zeros(Fs/2 - 1, 1)];
N = length(Vi);
Vo = zeros(N, 1);
x1 = 0;

for distortion = 0.01:.1:1.1
    % Components
    C1 = 1.6e-7;
    R1 = Ts/(2*C1);
    R2 = 3.1623;
    P1 = distortion * 100e3;

    G1 = 1/P1 + 1/R1;
    
    % DSP
    for n = 1 : N
        gainX = Vi(n,1) * (P1/R2);

        Vo(n,1) = gainX/(P1*G1) + x1/G1;

        x1 = (2/R1)*Vo(n,1) - x1;
    end

    
    % Graph
    [H, W] = freqz(Vo, 1, N, Fs);
    semilogx(W, 20*log10(abs(H)));
    axis([0 20000 30 95]);
    title('Open Loop Frequency Response - LM308 Model');
    ylabel('Gain (dB)'); xlabel('Frequency (Hz)');
    hold on;
end
hold off;