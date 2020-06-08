%--------------------------------------------------------------------------
% Hannah Lee
% AMATH 383 Spring 2020
% University of Washington
% beta_n is a function that calculates the beta value for gating variable
% n with input V.
%--------------------------------------------------------------------------

function b_n = beta_n(V)
    b_n = (0.125.*exp(-V./80));
end
