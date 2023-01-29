% RasSoln: This function outputs the phasor for given aquifer flow parameter inputs using the fully-confined and leaky-confined analytical models developed by:
% Rasmussen, T. C., Haborak, K. G., & Young, M. H. (2003). Estimating aquifer hydraulic properties using sinusoidal pumping at the Savannah River site, South Carolina
% USA. Hydrogeology Journal, 11(4), 466?482. https://doi.org/10.1007/s10040-003-0255-7

% Inputs:
%   test_list (num_tests x 4) - Matrix for each individual oscillatory flow experiment
%       Column 1 - Oscillation period [s]
%       Column 2 - Angular frequency [1/s]
%       Column 3 - Maximum pumping rate [m^3/s]
%       Column 4 - Radial distance [m]
%   s - Vector of natural log transformed parameter inputs
%       (2 x 1) for fully-confined model [T, S]
%       (3 x 1) for leaky-confined model [T, S, L]
%   soln - String describing which analytical solution to use 
%       'confined' or 'leaky'

% Outputs:
%   y_mod (2*num_obs x 1) - Vector of simulated phasor values separated into real and imaginary components. Odd elements are real component and even elements are imaginary component. 

% Code developed by Jeremy Patterson
% Created Dec 2018; Updated May 2022

function [y_mod] = RasSoln(test_list, s, soln)

s = exp(s);
D = s(1) / s(2);

omega = test_list(:,2);
Q_max = test_list(:,3);

r = test_list(:,4);
num_obs = numel(r);

phasor_mod = zeros(num_obs, 1);
for j = 1 : num_obs

    if strcmp(soln, 'leaky') == 1
        B_sq = s(1) / s(3);
        arg = sqrt(((1i * r(j)^2 * omega(j)) / D) +...
                    (r(j)^2 / B_sq));

    elseif strcmp(soln, 'confined') == 1
        arg = sqrt((1i * r(j)^2 * omega(j)) / D);
    end

    phasor_mod(j) =  Q_max(j) / (2 * pi * s(1)) * besselk(0, arg);
end
y_mod = zeros(2*num_obs,1);
y_mod(1:2:end-1) = real(phasor_mod);
y_mod(2:2:end) = imag(phasor_mod);
