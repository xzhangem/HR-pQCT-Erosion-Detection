function [c_x, c_y, c_z, hlf_l, hlf_w, hlf_h] = BoundBoXLoc(cluster_coe)
%% Parameter Setting
% Input:
%   cluter_coe: (num, 3) matrix of the coefficients of one erosion
% Output:
%   (c_x, c_y, c_z) is the coefficient of the centroid;
%   (hlf_l, hlf_w, hlf_h) is the half size of each side of bounding box

x_max = max(cluster_coe(:,1));  x_min = min(cluster_coe(:,1));
y_max = max(cluster_coe(:,2));  y_min = min(cluster_coe(:,2));
z_max = max(cluster_coe(:,3));  z_min = min(cluster_coe(:,3));

c_x = (x_max + x_min) / 2; hlf_l = (x_max - x_min) / 2;
c_y = (y_max + y_min) / 2; hlf_w = (y_max - y_min) / 2;
c_z = (z_max + z_min) / 2; hlf_h = (z_max - z_min) / 2;


end