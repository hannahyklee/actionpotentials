%--------------------------------------------------------------------------
% Hannah Lee
% AMATH 383 Spring 2020
% University of Washington
% hh_model numerically integrates and plots solutions of V(t) for the 
% Hodgkin-Huxley model for neuron action potentials using the built-in
% MATLAB function ODE45. The gating variables for the sodium and potassium 
% conductances are also plotted, as well as the externally applied current.
%--------------------------------------------------------------------------

clear all; close all; clc;

%-------------------
% Defining Constants
%-------------------

% NOTE: all constants defined in 1952 original Hodgkin-Huxley model; Table
% 3

% Equilibrium potentials  
E_K = -12;     
E_Na = 115;     
E_leak = 10.6; 

% Maximum conductances 
g_K = 36;
g_Na = 120;
g_leak = 0.3;

c_m = 1;

% Input current (stimulus) functions

% constant current of 10mV
I_ext = @(t)(10.*(t>0));

% step current
%I_ext = @(t)(2.*(t>20 & t <= 30) + 10.*(t>50 & t<100)...
 %             + 20.*(t>150));

%-----------------------------------------------------
% Solving HH model using MATLAB builtin function ODE45
%-----------------------------------------------------

hh_eqns = @(t, vars)[
    (1/c_m).*(g_Na .* vars(2)^3 .* vars(3) .* (E_Na - vars(1))...
             + g_K .* vars(4)^4 .* (E_K - vars(1))...
             + g_leak .* (E_leak - vars(1)) + I_ext(t));
     alpha_m(vars(1)) .* (1 - vars(2)) - beta_m(vars(1)) .* vars(2);
     alpha_h(vars(1)) .* (1 - vars(3)) - beta_h(vars(1)) .* vars(3);
     alpha_n(vars(1)) .* (1 - vars(4)) - beta_n(vars(1)) .* vars(4)];
 
tspan = [0, 200];
vars_0 = [0; .1; 0.6; 0.3]; % initial values for V, m, h, n

[t, vals] = ode45(hh_eqns, tspan, vars_0);

V = vals(:, 1);
m = vals(:, 2);
h = vals(:, 3);
n = vals(:, 4);

%-----------------------------------------------------------------
% Plotting membrane potential, gating variables, and input current
%-----------------------------------------------------------------

% Membrane Potential alone
figure
plot(t, V)
xlabel('Time (ms)');
ylabel('Membrane Potential');
title('Voltage vs. Time');
savefig('hh_membrane')

% Three gating variables
figure
plot(t,m, 'b');
hold 'on'
plot(t, h, 'r');
hold 'on'
plot(t, n, 'g');
legend('m','h','n');
xlabel('Time (ms)');
title('Gating Variables vs. Time');


% Input current alone
figure
plot(t, I_ext(t));
xlabel('Time (ms)');
ylabel('Input Current (mV)');
title('Input Current vs. Time');

% Subplots of 3 gating variables and input current
figure
xlabel('Time (ms)');
subplot(4,1,1)
plot(t,m);
title('Gating Variable: m');

subplot(4,1,2)
plot(t, n);
title('Gating Variable: n');

subplot(4,1,3)
plot(t, h);
title('Gating Variable: h');

subplot(4,1,4)
plot(t, I_ext(t));
title('Input Current');

% subplots of membrane potential, 3 gating variables, and input current
figure
subplot(5,1,1)
plot(t,V);
title('Membrane Potential');

subplot(5,1,2)
plot(t,m);
title('Gating Variable: m');

subplot(5,1,3)
plot(t, n);
title('Gating Variable: n');

subplot(5,1,4)
plot(t, h);
title('Gating Variable: h');

subplot(5,1,5)
plot(t, I_ext(t));
title('Input Current');
