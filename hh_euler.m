%--------------------------------------------------------------------------
% Hannah Lee
% AMATH 383 Spring 2020
% University of Washington
% hh_euler aims to confirm the calculations of hh_model by solving the 
% model numerically using the forward Euler approximations instead of the
% built-in ODE45 solver utilized in hh_model.
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

% Maximum conductances (1952 original HH model constants)
g_K = 36;
g_Na = 120;
g_leak = 0.3;

c_m = 1;
V_init = 0;

%---------------------
% Initializing vectors
%---------------------

T = 200;            
dt = .01;
t = 0:dt:T;         
steps = length(t);  

V = zeros(1, steps);
n = zeros(1, steps);
m = zeros(1, steps);
h = zeros(1, steps);
I_ext = zeros(1, steps);

V(1) = V_init;
n(1) = 0.3;
m(1) = 0.1;
h(1) = 0.6;

%--------------------
% Input current array       
%--------------------

I_ext(1:steps) = 10;

%--------------------------------------------
% Solving HH model using forward Euler method
%--------------------------------------------

for k = 1:steps - 1
    V(k + 1) = V(k) + dt*((1/c_m)*(g_Na*(m(k)^3)*h(k)*(E_Na - V(k))...
               + g_K*(n(k)^4)*(E_K - V(k))...
               + g_leak*(E_leak - V(k))+ I_ext(k)));
    m(k + 1) = m(k) + dt*(alpha_m(V(k))*(1 - m(k)) - beta_m(V(k))*m(k));
    h(k + 1) = h(k) + dt*(alpha_h(V(k))*(1 - h(k)) - beta_h(V(k))*h(k));
    n(k + 1) = n(k) + dt*(alpha_n(V(k))*(1 - n(k)) - beta_n(V(k))*n(k));

end

%-----------------------------------------------------------------
% Plotting membrane potential, gating variables, and input current
%-----------------------------------------------------------------

% Membrane Potential alone
figure
plot(t, V)
xlabel('Time (ms)');
ylabel('Membrane Potential');
title('Voltage vs. Time');

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
plot(t, I_ext);
xlabel('Time (ms)');
ylabel('Input Current (mV)');
title('Input Current vs. Time');

% Subplots of 3 gating variables and input current
figure
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
plot(t, I_ext);
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
plot(t, I_ext);
title('Input Current');






