function [y_hat, varargout] = MultiFreqPEST(test_list, y_obs, cov_y, s_init, lambda, delta, idx, soln)

num_obs = numel(idx);

for i = 1 : num_obs
    data_cov{i} = reshape(cov_y(idx(i),:), 2, 2);
    data_unc(i,:) = 1.96 * sqrt(diag(data_cov{i}));
end
R = blkdiag(data_cov{1:end});
R_inv = inv(R);
data_err = sqrt(mean(diag(R)));

y = zeros(2*num_obs, 1);
y(1:2:end-1) = y_obs(idx,1); 
y(2:2:end) = y_obs(idx,2);

% LM Inversion and Modeled Data
[pest, ~, out_flag] = Lev_Marq(test_list(idx,1:4), s_init, y, R_inv, lambda, delta, soln);
if out_flag == 0
    s_hat = Lev_Marq(test_list(idx,1:4), pest, y, R_inv, lambda, delta, soln);
else
    s_hat = pest;
end

% Parameter Estimate Uncertainty
J = jacob(s_hat, delta, test_list(idx,1:4), soln);
s_cov = inv(J' * R_inv * J);
D_unc = 1.96 * sqrt([1 -1] * s_cov(1:2,1:2) * [1; -1]);
s_unc = [1.96 * sqrt(diag(s_cov)); D_unc];

y_hat = RasSoln(test_list(idx,1:4), s_hat, soln);

obj_func = (1/2) * ((y-y_hat)' * R_inv * (y-y_hat));
rmse = sqrt(mean((y-y_hat).^2));

if strcmp(soln, 'confined') == 1
    varargout = {rmse, data_err, data_unc, [s_hat(1); s_hat(2); s_hat(1)-s_hat(2)], s_unc};
elseif strcmp(soln, 'leaky') == 1
    varargout = {rmse, data_err, data_unc, [s_hat(1); s_hat(2); s_hat(3); s_hat(1)-s_hat(2)], s_unc};
else
    error('Pick a valid analytical solution')
end