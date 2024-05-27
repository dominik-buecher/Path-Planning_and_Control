% Clear workspace and command window
clc;
clear;

% Given values
e_y = 0.1;       % Set the desired steady-state error here
v_star = 0.1;    % [m/s^2] Set the value for v* here
T_i = 0.2468;    % Set the value for T_i here
T = 316e-3;      % Set the value for T here
T_t = 0.1;       % Set the value for T_t here
k_u = 2.51;      % Set the value for k_u here
k_r = 0.3440;

% Define the transfer function
syms k_p s
G_wp = ((T_i*s) + 1) / ((((T_i*T*T_t)/(k_r*k_u*k_p))*s^4) + (((T_i*(T+T_t))/(k_r*k_u*k_p))*s^3) + ((T_i/k_p)*(1 + (1/(k_r*k_u)))*s^2) + (((1/k_p)+T_i)*s) + 1);

% Laplace transform of the ramp signal
W_p_s = v_star / s^2;

% Expression for E_p
E_p_s = (1 - G_wp) * W_p_s;
E_p_s_s = s * E_p_s;

% Calculate the steady-state error g_y
g_y = limit(E_p_s_s, s, 0);

% Calculate k_p for the desired steady-state error
kp_eqn = e_y == g_y; % Consider the sign
k_p_solution = solve(kp_eqn, k_p);

% Display the calculated value of k_p
disp('Calculated value of k_p:');
disp(k_p_solution);

% Import the Control System Toolbox
s = tf('s');
k_p_num = double(k_p_solution);
G_wp_num = ((T_i*s) + 1) / ((((T_i*T*T_t)/(k_r*k_u*k_p_num))*s^4) + (((T_i*(T+T_t))/(k_r*k_u*k_p_num))*s^3) + ((T_i/k_p_num)*(1 + (1/(k_r*k_u)))*s^2) + (((1/k_p_num)+T_i)*s) + 1);

% Simulation time
t = 0:0.1:10; % 10 seconds with 0.1 second steps

% Generate the ramp signal w_p(t)
w_p = v_star * t;

% Simulate the system response
y_p = lsim(G_wp_num, w_p, t);

% Plot the signal-time diagrams
figure;
plot(t, w_p, 'b-', 'LineWidth', 2);
hold on;
plot(t, y_p, 'r--', 'LineWidth', 2);
grid on;
title('Ramp response y_p(t)');
xlabel('Time (s)');
ylabel('Signal Value');
legend('w_p(t)', 'y_p(t)');

% Plot the Root Locus
figure;
rlocus(G_wp_num);
grid on;

%% 7.1 e)
s = tf('s');
vmax = 0.5;
x0 = 0;
xs = 1;
T_a = 20e-3;
k = 2.51;
k_p = 1;

% Compute coefficients for the reference polynomial with maximum velocity
[c, te] = cd_refpoly_vmax(vmax, x0, xs);

% Compute feedforward coefficients
cff = cd_refpoly_ff(c, k, T, T_t, k_r, T_i);

t = 0:T_a:te;

% Calculate u_vp1
u_vp1_t = polyval(cff, t);

% Reference trajectory: w_p(t) = c5*(t^5) + c4*(t^4) + c3*(t^3) + c2*(t^2) + c1*(t^1) + c0
w_p = polyval(c, t);

% Compute derivatives of the reference trajectory
c_punkt = polyder(c);
c_punkt_punkt = polyder(c_punkt);
w_p_punkt = polyval(c_punkt, t);
w_p_punkt_punkt = polyval(c_punkt_punkt, t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms kr Ti % Define symbolic variables

% Given values
ku = 2.51;
wd = pi;
phi_res = deg2rad(65); % Convert degrees to radians

% Define transfer functions
GR_cal = kr * (1 + Ti*1i*wd)/(Ti*1i*wd);
GS_cal = (ku /((T*1i*wd)+1))*exp(-1i*wd*T_t);

% Open-loop system
G0_cal = GR_cal*GS_cal;

% Formulate equations
eq1 = abs(G0_cal) == 1;
eq2 = angle(G0_cal) == -pi + phi_res;

% Solve equations
[sol_kr, sol_Ti] = vpasolve([eq1, eq2], [kr, Ti]);

sol_kr_num = double(sol_kr);
sol_Ti_num = double(sol_Ti);

s = tf('s');
% GS = (ku * (1) / (1 + s*T)) * exp(-s*T_t);
GS = (ku * (1)) / ((1 + s*T) * (1+s*T_t));
GR = sol_kr_num * (1 + s*sol_Ti_num) / (s*sol_Ti_num);

% Open-loop transfer function
G0 = GR * GS;
GW = G0 / (1 + G0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = 0:T_a:te;
G_vp1_s = s/((T_i * s) + 1);
y_p = lsim(minreal(GW * G_vp1_s /s), u_vp1_t, t);

% Plot system response w_p
figure;
plot(t, w_p);
grid on;
title('System response w_p(t)');
xlabel('Time (t)');
ylabel('w_p(t)');

% Plot system response w_p_punkt
figure;
plot(t, w_p_punkt);
grid on;
title('System response w_p dot(t)');
xlabel('Time (t)');
ylabel('w_p dot(t)');

% Plot system response w_p_punkt_punkt
figure;
plot(t, w_p_punkt_punkt);
grid on;
title('System response w_p dot dot(t)');
xlabel('Time (t)');
ylabel('w_p dot dot(t)');

% Plot system response u_vp1
figure;
plot(t, u_vp1_t);
grid on;
title('System response u_{vp1}(t)');
xlabel('Time (t)');
ylabel('u_vp1(t)');

% Plot system response y_p
figure;
hold on;
plot(t, y_p);
grid on;
title('System response y_p(t)');
xlabel('Time (t)');
ylabel('y_p(t)');
plot(t, w_p);
title('System response w_p(t)');
xlabel('Time(t)');
ylabel('w_p(t)');
grid on;
hold off;

%% 7.1 f)
G_vp1_z = c2d(G_vp1_s, T_a, 'tustin');
u_vpk = lsim(G_vp1_z, u_vp1_t, t, 0);

figure;
plot(t, u_vpk);
grid on;
title('u_{vpk}');
xlabel('Time(t)');
ylabel('u_vpk');
