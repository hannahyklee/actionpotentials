%--------------------------------------------------------------------------
% Hannah Lee
% AMATH 383 Spring 2020
% University of Washington
% hh_reduced plots, solves, and analyzes a simplified model with a 
% reduction of dimensions of the Hodgkin-Huxley model for neuron action 
% potentials.
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
V_init = 0;

%-----------------------------------
% Input current (stimulus) functions
%-----------------------------------

% constant current of 10mV
I_ext = @(t)(10.*(t>0));

% step current
%I_ext = @(t)(2.*(t>20 & t <= 30) + 10.*(t>50 & t<100)...
 %             + 20.*(t>150));

%-------------------------------------------------------------
% Solving HH reduced model using MATLAB builtin function ODE45
%-------------------------------------------------------------

hh_eqns = @(t, vars)[
    (1/c_m).*(g_Na .* m_inf(vars(1))^3.* (.8 - vars(2)).* (E_Na - vars(1))...
             + g_K .* vars(2)^4 .* (E_K - vars(1))...
             + g_leak .* (E_leak - vars(1)) + I_ext(t));
     alpha_n(vars(1)) .* (1 - vars(2)) - beta_n(vars(1)) .* vars(2)];
 
tspan = [0, 200];
vars_0 = [V_init; 0.3]; % initial values for V, n

tic
[t, vals] = ode45(hh_eqns, tspan, vars_0);
toc 

V = vals(:, 1);
n = vals(:, 2);

%-----------------------------------------------------------------
% Plotting membrane potential, gating variables, and input current
%-----------------------------------------------------------------

% Membrane Potential alone
figure
plot(t, V, 'r')
xlabel('Time (ms)');
ylabel('Membrane Potential');
title('Voltage vs. Time');
savefig('reduced_membrane')

% Input current alone
figure
plot(t, I_ext(t));
xlabel('Time (ms)');
ylabel('Input Current (mV)');
title('Input Current vs. Time');

% subplots of membrane potential, gating variables combined into variable
% n, and input current
figure
subplot(3,1,1)
plot(t,V);
title('Membrane Potential');

subplot(3,1,2)
plot(t, n);
xlabel('Time (ms)');
title('Gating variables represented by n vs. time');

subplot(3,1,3)
plot(t, I_ext(t));
title('Input Current');

%function to calculate the approximation for m
function y = m_inf(V)
    y = (alpha_m(V))/(alpha_m(V) + beta_m(V));
end
