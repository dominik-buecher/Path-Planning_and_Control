clear; clc;
syms kr Ti % Define symbolic variables

% Given values
ku = 2.51;
T = 316e-3;
wd = pi;
T_t = 0.1;
phi_res = deg2rad(65); % Convert degrees to radians

% Define transfer functions
GR_cal = kr * (1 + Ti*1i*wd)/(Ti*1i*wd);
GS_cal = (ku /((T*1i*wd)+1))*exp(-1i*wd*T_t);

% Open-loop system
G0_cal = GR_cal * GS_cal;

% Formulate equations
eq1 = abs(G0_cal) == 1;
eq2 = angle(G0_cal) == -pi + phi_res;

% Solve equations
[sol_kr, sol_Ti] = vpasolve([eq1, eq2], [kr, Ti]);

% Extract numerical solutions
sol_kr_num = double(sol_kr);
sol_Ti_num = double(sol_Ti);

s = tf('s');
GS = (ku * (1) / (1 + s*T)) * exp(-s*T_t);
GR = sol_kr_num * (1 + s*sol_Ti_num) / (s*sol_Ti_num);

% Open-loop transfer function
G0 = GR * GS;
GW = G0 / (1 + G0);

% Bode plot and margins
[mag, phase, wout] = bode(G0);
[GM, PM, Wcg, Wcp] = margin(G0);

% Display phase and gain margins
fprintf('Phase Margin: %f degrees\n', PM);
fprintf('Gain Margin: %f dB\n', 20*log10(GM));

% Plot Bode diagram
figure;
margin(G0);
grid on;

% Plot the step response
figure;
step(GW);
title('Step Response of the Closed-Loop System (Gw)');

% Evaluate and display step response characteristics
info = stepinfo(GW);
fprintf('Rise Time: %f seconds\n', info.RiseTime);
fprintf('Settling Time: %f seconds\n', info.SettlingTime);
fprintf('Overshoot: %f%%\n', info.Overshoot);
fprintf('Peak: %f\n', info.Peak);
fprintf('Peak Time: %f seconds\n', info.PeakTime);
