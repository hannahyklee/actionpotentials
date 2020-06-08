%--------------------------------------------------------------------------
% Hannah Lee
% AMATH 383 Spring 2020
% University of Washington
% alpha_m is a function that calculates the alpha value for gating variable
% m with input V.
%--------------------------------------------------------------------------

function a_m = alpha_m(V)
    a_m = (0.1.*(25 - V))./(exp((25 - V)./10) - 1);
end