function [J] = jacob(p, delta, test_list, soln)

% jacob: This function computes the Jacobian (i.e., parameter sensitivity) matrix for T and S (confined analysis) or T, S, and L (leaky analysis) using the analytical solutions 
% developed by Rasmussen et al., (2013). The Jacobian is then used in linearized error propagation to quantify the uncertainty of estimated aquifer parameters. 
%
% Inputs:
%   p - (num_param x 1) is the vector of aquifer parameters
%   delta - (num_param x 1) is the small incremental value by which parameters are changed to determine parameter sensitivity
%   test_list - (num_obs x 4) is the list of all pumping tests being analyzed. The columns of the matrix are:
%               [Period AngularFrequency(omega) Q_max r]
%   soln - Is a string that tells which analytical solution is being used. Should be 'confined' or 'leaky'
%
% Outputs:  
%   J - (numdata x numparam) Jacobian (i.e., parameter sensitivity) matrix 

% Code by Jeremy Patterson
% Created Dec 2020; Updated May 2022

% Modeled phasors with unchanged parameters
[coeffs_base] = RasSoln(test_list, p, soln);
    
for i = 1 : numel(p)
    pj = p; 
    pj(i) = p(i) + delta(i);
    
    % Modeled phasors with perturbed parameters
    [coeffs_mod] = RasSoln(test_list, pj, soln); 
    
    J(:,i) = (coeffs_mod - coeffs_base) ./ delta(i);   
end


