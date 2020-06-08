%--------------------------------------------------------------------------
% Hannah Lee
% AMATH 383 Spring 2020
% University of Washington
% beta_h is a function that calculates the beta value for gating variable
% h with input V.
%--------------------------------------------------------------------------

function b_h = beta_h(V)
    b_h = 1./(exp((30 - V)./10) + 1);
end 