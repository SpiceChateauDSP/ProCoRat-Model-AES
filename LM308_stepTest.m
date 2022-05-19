% LM308_stepTest.m
clear; clc; close all;

% Set Fs
Fs = 48000;
Ts = 1/Fs;
t = [0:Ts:2e-4 - Ts].';
f = 5000;

% Prepare LM308
slewRate = 2e5;
slope = slewRate/Fs;
m = 0;

% Create Step
x = 4.5 * square(2*pi*f*t);
N = length(x);
y = zeros(N,1);

% Process Step
for n = 1:N
    delta = x(n,1) - m;
    if delta > slope
        delta = slope;
    elseif delta < -slope
        delta = -slope;
    end
    y(n,1) = m + delta;
    % Update Memory
    m = y(n,1);
end

plot(t, x); hold on; plot(t, y); hold off;
axis([0 2e-4 -5 5]);
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', xt*1e6);
title('Voltage Follower Pulse Response - LM308 Model');
legend('Input', 'Output'); legend('boxoff');
ylabel('Voltage (V)'); xlabel('Time (us)');