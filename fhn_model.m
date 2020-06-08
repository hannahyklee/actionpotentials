%--------------------------------------------------------------------------
% Hannah Lee
% AMATH 383 Spring 2020
% University of Washington
% fhn_model numerically integrates and plots solutions of V(t) for the 
% Fitzhugh-Nagumo model for neuron action potentials. fhn_model also
% generates phase-plane portraits so that analysis of the solutions can be
% done.
%--------------------------------------------------------------------------

clear all; close all; clc;

%-------------------
% Defining Constants
%-------------------

a = 0.08;
b = 0.7;
c = 0.8;
I = 0.5;

%-----------------------------------------------------
% Solving FHN model using MATLAB builtin function ODE45
%-----------------------------------------------------

fhn_eqns = @(t, vars)[
    vars(1) - (vars(1).^3)./3 - vars(2) + I;
    a.*(vars(1) + b - c.*vars(2))];
 
tspan = [0, 200];
vars_0 = [I; 0]; % initial values for V, w

[t, vals] = ode45(fhn_eqns, tspan, vars_0);

V = vals(:, 1);     % F(V, w)
w = vals(:, 2);     % G(V, w)


%-----------------------------------------------------------------
% Plotting membrane potential, gating variables, and input current
%-----------------------------------------------------------------

% Membrane Potential alone
figure
plot(t, V, 'r')
xlabel('Time');
ylabel('Membrane Potential');
title('Voltage vs. Time');
savefig('fhn_membrane')

% Recovery variable w
figure
plot(t,w, 'b');
xlabel('Time');
title('w vs. Time');

%-----------------
% Phase plane plot
%-----------------

diff_eqs = @(t,Y)[ Y(1) - (Y(1).^3)./3 - Y(2) + I;
    a.*(Y(1) + b - c.*Y(2))];

[x, y] = meshgrid(-2:0.5:2, -2:0.5:2);

u = zeros(size(x));
v = zeros(size(x));

t = 0;
for i = 1:numel(x)
    Yprime = diff_eqs(t, [x(i); y(i)]);
    u(i) = Yprime(1);
    v(i) = Yprime(2);
end

figure
quiver(x,y,u,v); 
xlabel('V');
ylabel('w');
title('Phase Plane for the FitzHugh-Nagumo Model');
axis tight equal;

% plotting nullclines:

% nullcline 1: F(V, w) = 0, dV/dt = 0:
n_V = @(x)(x.*(1 - x.^2/3) + I);

% nullcline 2: G(V, w) = 0, dw/dt = 0:
n_w = @(x)((x + b)./c);

hold 'on';
tspan2 = (-2:.01:2);
plot(tspan2, n_V(tspan2), 'k');

hold 'on';
tspan3 = (-2:.01:1);
plot(tspan3, n_w(tspan3), 'k');

% finding fixed point:

diff = @(V)(n_V(V) - n_w(V));

Vstar = fzero(diff, 0);

% finding eigenvalues of Jacobian:

eigen1 = -0.5 * (Vstar^2 + a*c - 1) ...
        + 0.5 * sqrt((Vstar^2 + a*c - 1)^2 - 4*(-a*c + a + a*c*Vstar^2));
 
eigen2 = -0.5 * (Vstar^2 + a*c - 1) ...
        - 0.5 * sqrt((Vstar^2 + a*c - 1)^2 - 4*(-a*c + a + a*c*Vstar^2));
    
det = eigen1 * eigen2;

trace = eigen1 + eigen2;




