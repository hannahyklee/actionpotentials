%--------------------------------------------------------------------------
% Hannah Lee
% AMATH 383 Spring 2020
% University of Washington
% beta_m is a function that calculates the beta value for gating variable
% m with input V.
%--------------------------------------------------------------------------

function b_m = beta_m(V)
    b_m = 4*exp(-V./18);
end
