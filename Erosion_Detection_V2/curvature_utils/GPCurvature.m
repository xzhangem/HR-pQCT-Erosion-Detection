function [result] = GPCurvature(vol, C0_sigma, C0_filtsize)

%% Parameter setting
% Input:
%   vol: input binary segmentation volume
%   C0_sigma: sigma for Gaussian filter
%   C0_filtersize: filter size for Gaussian filter

% Output:
%   result: raw erosion detection result with K>0 & -1<H<0 at give scale
%   principle_matrix: at the possible erosion part two principe curvatures
%       are stored as k_min and k_max

contour = bwperim(vol);
contour = single(contour);

[vx, vy, vz] = size(vol);

D = bwdist(contour);
for i = 1 : vx
    for j = 1 : vy
        for k = 1 : vz
            if vol(i,j,k) == 0
                D(i,j,k) = -D(i,j,k);
            end
        end
    end
end

D = imgaussfilt3(D, C0_sigma, 'FilterSize', C0_filtsize);

[g_x, g_y, g_z] = gradient(D);


norm_g = (g_x.^2 + g_y.^2 + g_z.^2).^0.5;
n_x = -g_x ./ norm_g; n_y = -g_y ./ norm_g; n_z = -g_z ./ norm_g;

norm_vector(:,:,:,1) = n_x .* contour; norm_vector(:,:,:,2) = n_y .* contour; 
norm_vector(:,:,:,3) = n_z .* contour;


I = eye(3);
P_matrix = zeros(vx, vy, vz, 3, 3); % for contour

for i = 1 : vx
    for j = 1 : vy
        for k = 1 : vz
            if contour(i, j, k) == 1
                norm_v = squeeze(norm_vector(i, j, k, :));
                P_matrix(i, j, k, :, :) = I - norm_v * norm_v';
            end
        end
    end
end

%--------------------------------------------------------------------------

[g_xx, g_xy, g_xz] = gradient(g_x);
[g_yx, g_yy, g_yz] = gradient(g_y);
[g_zx, g_zy, g_zz] = gradient(g_z);

g_xy = 0.5 * (g_xy + g_yx);
g_yz = 0.5 * (g_yz + g_zy);
g_zx = 0.5 * (g_xz + g_zx);

Hessian_matrix(:,:,:,1,1) = g_xx .* contour; Hessian_matrix(:,:,:,1,2) = g_xy .* contour; Hessian_matrix(:,:,:,1,3) = g_zx .* contour;
Hessian_matrix(:,:,:,2,1) = g_xy .* contour; Hessian_matrix(:,:,:,2,2) = g_yy .* contour; Hessian_matrix(:,:,:,2,3) = g_yz .* contour;
Hessian_matrix(:,:,:,3,1) = g_zx .* contour; Hessian_matrix(:,:,:,3,2) = g_yz .* contour; Hessian_matrix(:,:,:,3,3) = g_zz .* contour;

fprintf('Hessian Finished. \n');

G_matrix = zeros(vx, vy, vz, 3, 3);
Curvature_matrix = zeros(vx, vy, vz, 2); % Gaussian curvature for 1st dim, mean curvature for 2nd dim. 

for i = 1 : vx
    for j = 1 : vy
        for k = 1 : vz
            if contour(i, j, k) == 1
                G_submatrix = - squeeze(P_matrix(i,j,k,:,:)) * squeeze(Hessian_matrix(i,j,k,:,:)) ...
                    * squeeze(P_matrix(i,j,k,:,:)) / norm_g(i,j,k);
                G_matrix(i,j,k,:,:) = G_submatrix;
                %fprintf('[%d, %d, %d] G_matrix finished. \n', i, j, k);
                T = trace(G_submatrix);
                F = norm(G_submatrix, 'fro');
                %fprintf('%d, %d \n', T, F);
                gauss_curvature = (T^2 - F^2) / 2;
                mean_curvature = T / 2;
                Curvature_matrix(i,j,k,1) = gauss_curvature;    Curvature_matrix(i,j,k,2) = mean_curvature; 
            end
        end
    end
end

result = (Curvature_matrix(:,:,:,1) > 0.0) .* (Curvature_matrix(:,:,:,2) > -1) .* (Curvature_matrix(:,:,:,2) < 0);
result = single(result);
result = result .* contour;

%Curvature_matrix(:,:,:,1) = Curvature_matrix(:,:,:,1) .* result;
%Curvature_matrix(:,:,:,2) = Curvature_matrix(:,:,:,2) .* result;
end