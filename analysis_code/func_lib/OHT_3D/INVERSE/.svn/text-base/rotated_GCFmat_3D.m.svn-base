function GCF_mat = rotated_GCFmat_3D(dirdistmat,GCF_func,GCF_params,angles)

%rotated_GCFmat_3D: Function which takes as input the x,y, and z distances
%between locations in a field and calculates the generalized covariance
%values (or variance values), given an *anisotropic* variogram with
%arbitrary principle axes.
%
%Syntax:
%   [GCF_mat] = rotated_GCFmat_3D(dirdistmat,GCF_func,GCF_params,angles)
%   where:
%       -GCF_mat is the output generalized covariance or variance values
%       for the points
%       -dirdistmat is a (num_pts1 x num_pts2 x dim) matrix representing
%       the distances between two setst of points
%       -GCF_func is the anisotropic GCF function that calculates GCF
%       values in the principle coordinate system. For example, GCF_func =
%       @(distmat) var*(((distmat(:,:,1)./2).^2 + (distmat(:,:,2)./5).^2 +
%       (distmat(:,:,3)./9).^2).^.5) could be used as a linear variogram
%       with different correlation lengths along each principle direction.
%       -angles is a (3 x 1) vector giving the rotation angles that relate
%       the points coordinate system to the principle coordinate system.



dims = size(dirdistmat,3);
num_pts1 = size(dirdistmat,1);
num_pts2 = size(dirdistmat,2);

Rx = [1 0 0; 0 cosd(angles(1)) -sind(angles(1)); 0 sind(angles(1)) cosd(angles(1))];
Ry = [cosd(angles(2)) 0 sind(angles(2)); 0 1 0; -sind(angles(2)) 0 cosd(angles(2))];
Rz = [cosd(angles(3)) -sind(angles(3)) 0; sind(angles(3)) cosd(angles(3)) 0; 0 0 1];

%Reshape to vectors to perform rotation mathematics, then back into the
%matrices.

dirdist_list = zeros(num_pts1*num_pts2,dims);
for i = 1:1:dims
    dirdist_list(:,i) = reshape(dirdistmat(:,:,i),num_pts1*num_pts2,1);
end

rot_dirdist_list = (Rx*Ry*Rz*dirdist_list')';

rot_dirdistmat = zeros(num_pts1,num_pts2,dims);
for i = 1:1:dims
    rot_dirdistmat(:,:,i) = reshape(rot_dirdist_list(:,i),num_pts1,num_pts2);
end

GCF_mat = GCF_func(rot_dirdistmat,GCF_params);


