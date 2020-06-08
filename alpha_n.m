%--------------------------------------------------------------------------
% Hannah Lee
% AMATH 383 Spring 2020
% University of Washington
% alpha_n is a function that calculates the alpha value for gating variable
% n with input V.
%--------------------------------------------------------------------------

function a_n = alpha_n(V)
    a_n = (0.01.*(10-V))./(exp((10-V)./10) - 1);
end