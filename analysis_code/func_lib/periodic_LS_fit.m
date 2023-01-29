function [obj_func, varargout] = periodic_LS_fit(times,sig,periods)

nreqout = 1;

n = numel(sig);
num_per = numel(periods);
num_times = size(times,1);

if numel(sig) ~= num_times
    error('Size of signal is incorrect')
end

G = zeros(num_times,2*num_per);
Gcol = 1;
for p = 1:1:num_per
    G(:,[Gcol Gcol+1]) = [cos(2*pi/periods(p).*times), -sin(2*pi/periods(p).*times)];
    Gcol = Gcol + 2;
end

c_opt = (G'*G)\(G'*sig);

sig_sq = ((sig - G*c_opt)'*(sig - G*c_opt)) / n; %Mean squared error
obj_func = inv(G'*G) .* sig_sq; %Data covariance matrix

phasors_opt = zeros(num_per,1);
c_loc = 1;
for p = 1:1:num_per
    phasors_opt(p) = c_opt(c_loc) + 1i*c_opt(c_loc+1);
    c_loc = c_loc + 2;
end

if nargout > nreqout
    varargout{1} = phasors_opt;
    varargout{2} = sig_sq;
end
