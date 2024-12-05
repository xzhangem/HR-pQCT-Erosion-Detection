clc;
clear all;
close all;

erosion = niftiread('../Weibull_refine/3926_M25478_MCP3_bl_o_Weibull_20_refine.nii');
n = 9;
sigma =  8 / sqrt(2); %(sqrt(2))^n;
filtersize = FilterSize(sigma);
[erosion, principle_curvature] = GPCurvature(erosion, sigma, filtersize);
[vx, vy, vz] = size(erosion);
bound_len = 30;
vol_mask  = zeros(vx, vy, vz);
vol_mask(bound_len:vx-bound_len, bound_len:vy-bound_len, bound_len:vz-bound_len) = 1;
erosion = erosion .* vol_mask;

[label_erosion, erosion_num] = SpatialCluster(erosion, 7);

erosion_coe_cell = ClusterIndx(label_erosion, erosion_num);

num_per_erosion = zeros(erosion_num, 1);
for i = 1 : erosion_num
    a = size(erosion_coe_cell{i});
    num_per_erosion(i, 1) = a(1);
end

boundingbox_vol = zeros(vx, vy, vz);

for i = 1 : erosion_num
    [c_x, c_y, c_z, hlf_l, hlf_w, hlf_h] = BoundBoXLoc(erosion_coe_cell{i});
    x_min = floor(c_x - hlf_l);     x_max = ceil(c_x + hlf_l);
    y_min = floor(c_y - hlf_w);     y_max = ceil(c_y + hlf_w);
    z_min = floor(c_z - hlf_h);     z_max = ceil(c_z + hlf_h);
    
    boundingbox_vol(x_min:x_max, y_min:y_max, z_min:z_max) = i;
end

ellipsoid_vol = zeros(vx, vy, vz);

for i = 1 : erosion_num
    if num_per_erosion(i, 1) > 60
        erosion_part = erosion_coe_cell{i};
        erosion_part = erosion_part';
        [E, c] = lowner(erosion_part, 0.001);
        c = c';
        eig_value = eig(E);
        curvature_ratio = max(eig_value) / min(eig_value);
        fprintf('%d Erosion with %d points: %f.\n', i, num_per_erosion(i,1), curvature_ratio);
        if curvature_ratio < 10
            ell_v = ((4 * pi) / 3) / (eig_value(1) * eig_value(2) * eig_value(3));
            fprintf('%d Erosion Volume: %f. \n', i, ell_v);
            for p = 1 : vx
                for q = 1 : vy
                    for r = 1 : vz
                        vec = [p, q, r] - c;
                        if vec * E * vec' <= 1
                            ellipsoid_vol(p,q,r) = i;
                        end
                    end
                end
            end
        end
    end
end

niftiwrite(ellipsoid_vol, '3926_3_s5_erosion.nii');

function w = FilterSize(s)
w = 2 * ceil(2 * s) + 1;
end