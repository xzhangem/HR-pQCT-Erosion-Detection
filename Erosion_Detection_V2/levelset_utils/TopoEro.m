function [region_vol, ellipsoid_vol] = TopoEro(new_vol, ratio_thresh)
addpath('../curvature_utils/');
erosion_min = 1;
erosion_max = 100000;
[poss_erosion, num] = bwlabeln(new_vol, 26);
stats = regionprops(poss_erosion);
area = cat(1, stats.Area);
%niftiwrite(poss_erosion, 'new_vol.nii');
%[sort_area, ind] = sort(area, 'descend');
erosion_coe_cell = ClusterIndx(poss_erosion, num);
num_per_erosion = zeros(num, 1);

[vx, vy, vz] = size(new_vol);

for i = 1 : num
    a = size(erosion_coe_cell{i});
    num_per_erosion(i, 1) = a(1);
end

ellipsoid_vol = zeros(vx, vy, vz);

for i = 1 : num
    if (num_per_erosion(i, 1) > erosion_min) && (num_per_erosion(i, 1) < erosion_max)  
        erosion_part = erosion_coe_cell{i};
        erosion_part = erosion_part';
        [E, c] = lowner(erosion_part, 0.001);
        E_isnan = isnan(E); E_isinf = isinf(E);
        if (sum(sum(E_isnan)) == 0) && (sum(sum(E_isinf)) == 0)
            c = c';
            eig_value = eig(E);
            curvature_ratio = max(eig_value) / min(eig_value);
            fprintf('%d Erosion with %d points: %f.\n', i, num_per_erosion(i,1), curvature_ratio);
            if (curvature_ratio < ratio_thresh) && (curvature_ratio > 0.0)
                ell_v = ((4 * pi) / 3) / (eig_value(1) * eig_value(2) * eig_value(3));
                fprintf('%d Erosion Volume: %f. \n', i, ell_v);
                for p = 1 : vx
                    for q = 1 : vy
                        for r = 1 : vz
                            vec = [p, q, r] - c;
                            if vec * E * vec' <= 1
                                ellipsoid_vol(p, q, r) = i;
                            end
                        end
                    end
                end
            end
        end
    end
end
%niftiwrite(ellipsoid_vol, 'ellip_vol.nii');
val_vec = unique(ellipsoid_vol);
val_vec = val_vec(2:end);
region_vol = zeros(vx, vy, vz);
for i = 1 : length(val_vec)
    region_vol = region_vol + val_vec(i) * single(poss_erosion == val_vec(i));
end
end