function ellipsoid_vol = EllipsoidShape(erosion, radium, erosion_min, ratio_thresh)

% Input:
%   erosion: the raw erosion volume only with voxel-wise information
%   radium: radium for spatial clustering
%   erosion_min: threshold for erasing extremely small erosions
%   ratio_thresh: threshold ratio for curvature threshold
% Output:
%   ellipsoid_vol: detected erosion in form of ellipsoid

[vx, vy, vz] = size(erosion);

[label_erosion, erosion_num] = SpatialCluster(erosion, radium);
erosion_coe_cell = ClusterIndx(label_erosion, erosion_num);

num_per_erosion = zeros(erosion_num, 1);

for i = 1 : erosion_num
    a = size(erosion_coe_cell{i});
    num_per_erosion(i, 1) = a(1);
end

ellipsoid_vol = zeros(vx, vy, vz);

for i = 1 : erosion_num
    if (num_per_erosion(i, 1) > erosion_min)  
        erosion_part = erosion_coe_cell{i};
        erosion_part = erosion_part';
        [E, c] = lowner(erosion_part, 0.001);
        E_isnan = isnan(E); E_isinf = isinf(E);
        if (sum(sum(E_isnan)) == 0) && (sum(sum(E_isinf)) == 0)
            c = c';
            eig_value = eig(E);
            curvature_ratio = max(eig_value) / min(eig_value);
            fprintf('%d Erosion with %d points: %f.\n', i, num_per_erosion(i,1), curvature_ratio);
            if (curvature_ratio < ratio_thresh) && (curvature_ratio > 0) 
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