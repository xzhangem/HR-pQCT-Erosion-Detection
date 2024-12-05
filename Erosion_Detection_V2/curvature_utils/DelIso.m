function [de_erosion, de_curvature_tensor] = DelIso(erosion, curvature_tensor, iso_r)
%% Parameter Setting
% Input: 
%   erosion: the raw erosion detection result
%   curvature_tensor: the curvature-scale 5D-tensor indicating curvature
%       information of the erosion
%   iso_r: radium for determining isolated points
% Output:
%   de_erosion: the erosion after deleting isoloated parts
%   de_curvature_tensor: the curvature_tensor after deleting isolated parts

bound_len = 25;
[vx, vy, vz, len, ~] = size(curvature_tensor);
vol_mask = zeros(vx, vy, vz);
vol_mask(bound_len:vx-bound_len, bound_len:vy-bound_len, bound_len:vz-bound_len) = 1;
de_erosion = erosion .* vol_mask;

tensor_mask = zeros(vx, vy, vz, len, 4);
tensor_mask(bound_len:vx-bound_len, bound_len:vy-bound_len, bound_len:vz-bound_len,:,:) = 1;
de_curvature_tensor = curvature_tensor .* tensor_mask;

for i = 1 + iso_r : vx - iso_r
    for j = 1 + iso_r : vy - iso_r
        for k = 1 + iso_r : vz - iso_r
            if de_erosion(i,j,k) == 1
                submat = de_erosion(i-iso_r:i+iso_r, j-iso_r:j+iso_r, k-iso_r:k+iso_r);
                s = sum(sum(sum(submat)));
                if s < 5
                    de_erosion(i,j,k) = 0;
                    de_curvature_tensor(i,j,k,:,:) = 0;
                end
            end
        end
    end
end

end