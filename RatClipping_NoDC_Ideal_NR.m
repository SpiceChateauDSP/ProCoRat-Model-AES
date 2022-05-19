% RatClipping_NoDC_Ideal_NR.m
clear; clc; close all;

% Input Parameters
Fs = 48000;

Ts = 1/Fs;
t = [0:Ts:.1-Ts].';
f = 100;

% Generate Input
x = 4.5 * sin(2 * pi * f * t);
N = length(x);

% Diode Characteristics
Is = 10e-15;
eta = 1.000;
Vt = 0.026;

% Other Components
R1 = 1e3;

% Prepare Loop
y = zeros(N, 1);
Vo = 0;

for n = 1 : N
    Vi = x(n, 1);
    
    % Transfer Function
    p = -Vi/R1;
    num = 2 * Is * sinh(Vo/(eta * Vt)) + Vo/R1 + p;
    
    b = 1;
    count = 1;
    if (abs(num) > .0000000001 && count < 10)
        % Derivate of Transfer Function for Vo
        den = 2 * Is / (eta * Vt) * cosh(Vo / (eta * Vt)) + 1/R1;

        Vnew = Vo - b * num/den;
        num2 = 2 * Is * sinh(Vnew/(eta * Vt)) + Vnew/R1 + p;
        
        % Damp the Newton Raphson to eliminate divergence
        if (abs(num2) > abs(num))
            b = b/2;
        else
            Vo = Vnew;
            b = 1;
        end
        
        % Recalculate numerator
        num = 2 * Is * sinh(Vo/(eta * Vt)) + Vo/R1 + p;
        
        count = count + 1;
    end
    
    y(n, 1) = Vo;
end

plot(x); hold on;  plot(y); hold off;