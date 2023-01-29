function [condrelz_field, varargout] = simple_condrelz(x,data_loc, data_val, GCF_func, varargin)

%simple_condrelz: Kriging equations used to generate a conditional
%realization of a parameter field given data at specified locations.
%Utilizes the "function estimation" approach.
%Note: As a global kriging / conditional realization method, this code is
%for relatively small-scale problems only (1000's of positions or fewer).
%Other methods such as spectral approaches may be used for larger problem
%sizes. Additionally, neighborhood kriging may be employed with care.
%[condrelz_field,[times]] =
%simple_condrelz(x,data_loc,data_val,GCF_func,[drift_func])
%where:
%   -condrelz_field is the conditional realization of the property values
%   at the locations x
%   -times (optional) is a 2 x 1 vector that returns 1) the amount of time
%   required to compute the functional coefficients for the kriging
%   function estimate form, and 2) the time to krig all datapoints using
%   these values
%   -x is the list of locations (unk x dim matrix) where the property needs
%   to be conditionally simulated
%   -data_loc is the list of locations (n x dim matrix) where the
%   property's value has been measured (assumed error-free measurements)
%   -data_val is the list of measured property values at the data_loc
%   locations (n x 1 vector)
%   -GCF_func is the generalized covariance function or covariance function
%   which gives the expected covariance for two points at a specified
%   distance. This function must accept matrix arguments for distances.
%   -drift_func (optional) is a function which takes an (n x dim) list of
%   points and returns an (n x p) list of drift function values for each
%   point. If drift_func is not supplied a constant unknown mean is
%   assumed.

drift_func = [];

nreqin = 4;
if nargin > nreqin
    drift_func = varargin{1};
end

num_unk = size(x,1);
dim = size(x,2);
num_data = size(data_loc,1);

%Create an unconditional realization at the measurement points and the
%unknown locations
full_loclist = [x; data_loc];
fulllist_dist = dimdist(full_loclist);
fulllist_covar = GCF_func(fulllist_dist);
Qpf = chol(fulllist_covar);
uncond_relz = Qpf'*randn((num_unk+num_data),1);
uncond_unk = uncond_relz(1:num_unk);
uncond_data = uncond_relz((num_unk+1):(num_unk+num_data));

%Form matrix to find weights (xis and beta)
tstart = clock;
data_dist = dimdist(data_loc);
data_covar = GCF_func(data_dist);
if ~isempty(drift_func)
    X = drift_func(data_loc);
else
    X = ones(num_data,1);
end
p = size(X,2);
A = [0.5*data_covar X; X' zeros(p,p)];
b = [data_val - uncond_data; zeros(p,1)];
xibeta_vec = A\b;
times(1) = etime(clock,tstart);

%Perform kriging for all points

tstart = clock;
condrelz_field = zeros(num_unk,1);
for i = 1:num_unk
    x2data_dist = dimdist(x(i,:),data_loc);
    x2data_covar = GCF_func(x2data_dist);
    if ~isempty(drift_func)
        xdrift = drift_func(x(i,:));
    else
        xdrift = 1;
    end
    condrelz_field(i) = [0.5*x2data_covar xdrift]*xibeta_vec + uncond_unk(i);
end
times(2) = etime(clock,tstart);

if nargout > 1
    varargout{1} = times;
end