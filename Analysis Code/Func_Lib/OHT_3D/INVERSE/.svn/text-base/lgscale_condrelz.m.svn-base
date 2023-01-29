function [condrelz_field, varargout] = lgscale_condrelz(x,x_griddim,data_loc, data_val, GCF_func, varargin)

%lgscale_condrelz: Kriging equations used to generate a conditional
%realization of a parameter field given data at specified locations.
%Utilizes the "function estimation" approach. Code by Michael Cardiff.
%This code generates an unconditional realization using Toeplitz tricks on
%a fine, equispaced grid which is then interpolated to the locations of the
%measurements and unknowns. The conditional realization is thus somewhat
%APPROXIMATE. If the unknown locations (x) are on a regular equispaced
%grid, then these may be used. Otherwise, an equispaced grid is
%automatically generated that covers the region spanned by the unknowns and
%measurements.
%[condrelz_field,[times]] =
%lgscale_condrelz(x,x_griddim,data_loc,data_val,GCF_func,[drift_func])
%where:
%   -condrelz_field is the conditional realization of the property values
%   at the locations x
%   -times (optional) is a 2 x 1 vector that returns 1) the amount of time
%   required to compute the functional coefficients for the kriging
%   function estimate form, and 2) the time to krig all datapoints using
%   these values
%   -x is the list of locations (unk x dim matrix) where the property needs
%   to be conditionally simulated
%   -x_griddim: If x is an equispaced grid, x_griddim is a (dim x 1) vector
%   that represents the number of points along each dimension. If x is not
%   an equispaced grid, set x_griddim = [];
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

nreqin = 5;
if nargin > nreqin
    drift_func = varargin{1};
end

%Maximum number of points in the equi-spaced grids. This can be adjusted
%up or down based on memory availability. Finer grids result in more
%accurate conditional realizations.
equigrid_maxsize = 1e6;

num_unk = size(x,1);
dim = size(x,2);
num_data = size(data_loc,1);

%Create an unconditional realization at the measurement points and the
%unknown locations
if ~isempty(x_griddim)
    %If the unknown locations are an equispaced grid that span the
    %measurement locations, then use this to generate an unconditional
    %realization using toeplitz tricks and then interpolate to the data 
    %locations.
    
    dist_firstrow = dimdist(x,x(1,:));
    
    cov_firstrow = GCF_func(dist_firstrow);
    uncond_unk = toepmat_vector_math(cov_firstrow, 'r', [],dim, x_griddim);
    %Only need one conditional realization, toepmat_vector_math returns 2
    uncond_unk = uncond_unk(:,1); 
    if dim == 1
        uncond_data = interp1(x,uncond_unk,data_loc);
    elseif dim == 2
        xg1 = reshape(x(:,1),x_griddim);
        xg2 = reshape(x(:,2),x_griddim);
        ug = reshape(uncond_unk,x_griddim);
        uncond_data = interp2(xg1,xg2,ug,data_loc(:,1),data_loc(:,2));
    elseif dim == 3
        xg1 = reshape(x(:,1),x_griddim);
        xg2 = reshape(x(:,2),x_griddim);
        xg3 = reshape(x(:,3),x_griddim);
        ug = reshape(uncond_unk,x_griddim);
        uncond_data = interp3(xg1,xg2,xg3,ug,data_loc(:,1),data_loc(:,2),data_loc(:,3));
    end
else
    %In here - For the case when the unknown locations are NOT an
    %equispaced grid, then create an unconditional realization on a
    %GENERATED equispaced grid with the largest possible number of
    %elements, then interpolate to both unknown locations and data
    %locations
    points_perdim = floor(equigrid_maxsize^(1/dim));
    cgrid_griddim = points_perdim*ones(1,dim);
    all_locs = [x; data_loc];
    inp_list = '';
    out_list = '';
    bounds = zeros(dim,2);
    for i = 1:1:dim
        bounds(i,1) = min(all_locs(:,i));
        bounds(i,2) = max(all_locs(:,i));
    end
    %TO GO HERE - use outputs to calculate first row of distance matrix,
    %first row of covariance matrix. Call toepmat_vector_math, and finally
    %interpolate to all locations using appropriate interpolation routines,
    %as above.
    if dim == 1
        xg1 = (linspace(bounds(1,1),bounds(1,2),points_perdim))';
        cgrid_coords = xg1;
    elseif dim == 2
        [xg1,xg2] = meshgrid(linspace(bounds(1,1),bounds(1,2),points_perdim),...
            linspace(bounds(2,1),bounds(2,2),points_perdim));
        cgrid_coords = [reshape(xg1,numel(xg1),1) reshape(xg2,numel(xg2),1)];
    elseif dim == 3
        [xg1,xg2,xg3] = meshgrid(linspace(bounds(1,1),bounds(1,2),points_perdim),...
            linspace(bounds(2,1),bounds(2,2),points_perdim),...
            linspace(bounds(3,1),bounds(3,2),points_perdim));
        cgrid_coords = [reshape(xg1,numel(xg1),1) reshape(xg2,numel(xg2),1) ...
            reshape(xg3,numel(xg3),1)];
    end
    
    dist_firstrow = dimdist(cgrid_coords,cgrid_coords(1,:));
%     dist_firstrow = zeros(num_unk,1);
%     for i = 1:1:dim
%         dist_firstrow = dist_firstrow + (cgrid_coords(1,i) - x(:,i)).^2;
%     end
%     dist_firstrow = dist_firstrow.^.5;
     
    cov_firstrow = GCF_func(dist_firstrow);
    
    uncond_grid = toepmat_vector_math(cov_firstrow, 'r', [],dim, cgrid_griddim);
    %Only need one conditional realization, toepmat_vector_math returns 2
    uncond_grid = uncond_grid(:,1); 

    if dim == 1
        uncond_unk = interp1(xg1,uncond_grid,x);
        uncond_data = interp1(xg1,uncond_grid,data_loc);
    elseif dim == 2
        ug = reshape(uncond_grid,cgrid_griddim);
        uncond_unk = interp2(xg1,xg2,ug,x(:,1),x(:,2));
        uncond_data = interp2(xg1,xg2,ug,data_loc(:,1),data_loc(:,2));
    elseif dim == 3
        ug = reshape(uncond_grid,cgrid_griddim);
        uncond_unk = interp3(xg1,xg2,xg3,ug,x(:,1),x(:,2),x(:,3));
        uncond_data = interp3(xg1,xg2,xg3,ug,data_loc(:,1),data_loc(:,2),data_loc(:,3));
    end
end

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