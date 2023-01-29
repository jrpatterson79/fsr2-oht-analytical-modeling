function [Qproduct] = covar_product_aper(Q_aper_row,invec,numx,numy,varargin)

%covar_product_aper: Function for computing product of a vector with the
%covariance matrix for the model parameter fracture aperture. Code assumes that
%the covariance matrix Q has the structure Q = [Q_aper], where Q_aper is determined
%by a spatial covariance function, i.e. the variogram. 
%This code utilizes the toeplitz approach to matrix-vector products, thus
%only the first rows of the Q matrices need to be supplied
%
%[Qproduct] = covar_product_aper(Q_aper_row,invec,numx,numy,{numz})
%
%Where:
%   -Qproduct is the product of the Q_aper matrix with invec 
%   -Q_aper_row (numcells x 1) is the first row of the covariance matrix
%   for aperture values, where numcells is the number of cells in the model.
%   -invec (numcells x 1) is the vector being multiplied by Q_aper
%   -numx is the number of cells in the x direction
%   -numy is the number of cells in the y direction
%   -numz (optional) is the number of cells in the z direction
%
% Code by Michael Cardiff
% 10/2014, Last Updated: 1/2015 

num_reqin = 4;
dim = 2;

if nargin > num_reqin
    numz = varargin{1};
    dim = 3;
end

numaper = numel(Q_aper_row);

if dim == 2
    Qproduct = toepmat_vector_math(Q_aper_row,'*',invec,2,[numy numx]);
else
    Qproduct = toepmat_vector_math(Q_aper_row,'*',invec,3,[numy numx numz]);
end
   