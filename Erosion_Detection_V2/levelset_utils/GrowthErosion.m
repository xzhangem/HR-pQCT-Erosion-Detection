function GrowthErosion(save_dir, vol, new_vol, num_thresh, ratio_thresh)

if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

timestep = 1; % 1
mu = 0.2 / timestep; % weight of the distance regularization term
iter = 10; % iteration of the global distribution levelset
lambda = 1; % weight of the weighted length term
alpha = 1; % weight of the weighted area/dilation term
epsilon = 0.2; % width of Dirac Delta function and Heaviside function


[erosion_lab, num] = bwlabeln(new_vol, 26);
[vx, vy, vz] = size(erosion_lab);
for l = 1 : num
    erosion_vol = single(erosion_lab == l);
    phi = ErosionGrow(erosion_vol, vol, mu, lambda, alpha, epsilon, timestep, 5);
    ero_vol = single(phi > 0);
    erosion_num = sum(ero_vol(:));
    tag = 1;
    erosion_coe = zeros(3, erosion_num);
    for i = 1 : vx
        for j = 1 : vy
            for k = 1 : vz
                if ero_vol(i,j,k) == 1
                    erosion_coe(1, tag) = i; erosion_coe(2, tag) = j; erosion_coe(3, tag) = k;
                    tag = tag + 1;
                end
            end
        end
    end
    [E, c] = lowner(erosion_coe, 0.001);
    E_isnan = isnan(E); E_isinf = isinf(E);
    if (sum(sum(E_isnan)) == 0) && (sum(sum(E_isinf)) == 0)
        c = c';
        eig_value = eig(E);
        curvature_ratio = max(eig_value) / min(eig_value);
        fprintf('Erosion with %d points: %f.\n', erosion_num, curvature_ratio);
        if (curvature_ratio < ratio_thresh) && (curvature_ratio > 0.0) && (erosion_num > num_thresh)
            ero_name = [save_dir, 'ero_', num2str(l), '.nii'];
            poss_erosion = zeros(vx+2, vy+2, vz+2);
            poss_erosion(2:vx+1, 2:vy+1, 2:vz+1) = ero_vol;
            niftiwrite(poss_erosion, ero_name);
        end
    end
end


end