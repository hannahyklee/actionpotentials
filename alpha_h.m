%--------------------------------------------------------------------------
% Hannah Lee
% AMATH 383 Spring 2020
% University of Washington
% alpha_h is a function that calculates the alpha value for gating variable
% h with input V.
%--------------------------------------------------------------------------

function a_h = alpha_h(V)
    a_h = 0.07.*exp(-V./20);
end