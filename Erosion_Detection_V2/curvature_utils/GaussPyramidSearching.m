function [start_scale, raw_erosion, curvature_scale_tensor] = GaussPyramidSearching(vol, ratio, init_sigma)

%% Parameter Setting
% Input:
%   vol: input binary segmentation volume
%   th_ratio: threshold ratio determining the last scale factor
%   initial_sigma: initial sigma for searching
% Output:
%   start_scale: the smallest significant scale for erosion recognition
%   raw_erosion: binary volume for indicating erosion parts
%   curvature_scale_tensor: size of (vx, vy, vz, scale_len, 4), the last 
%       dimension is (k_min, k_max, start_scale, end_scale) 
%       (only for raw erosion parts)

[vx, vy, vz] = size(vol);
contour = single(bwperim(vol));
contour_sum = sum(sum(sum(contour)));
resolution_thresh = contour_sum * ratio;
%init_fs = FilterSize(init_sigma);


n = 0;
tag = true;
sqrt_2 = sqrt(2);
while (tag)
    sigma = init_sigma * sqrt_2^n;
    filtersize = FilterSize(sigma);
    [erosion, ~] = GPCurvature(vol, sigma, filtersize);
    n = n + 1;
    fprintf('Scale: %d-----ratio: %f. \n', sigma, sum(sum(sum(erosion))) / contour_sum);
    if sum(sum(sum(erosion))) < resolution_thresh
        tag = false;
    end
end

sigma_tail = round(sigma);
sigma_head = round(sigma / sqrt_2);
start_scale = sigma_head;


sigma_op = sigma_tail;
[head_Gauss_indx, head_principle_mat] = GPCurvature(vol, sigma_head, FilterSize(sigma_head));

scale_interval = 1;
scale_len = sigma_op - sigma_head + 1;
indx_list = cell(scale_len, 1);
principle_list = cell(scale_len, 1);
indx_list{1} = head_Gauss_indx;
principle_list{1} = head_principle_mat;
scale_indx = 2;

for i = sigma_head+1 : scale_interval : sigma_op
    [indx_mat, principle_mat] = GPCurvature(vol, i, FilterSize(i));
    indx_list{scale_indx} = indx_mat;
    principle_list{scale_indx} = principle_mat;
    scale_indx = scale_indx + 1;
end

indx_sum_matrix = zeros(vx, vy, vz);

for i = 1 : scale_len
    indx_sum_matrix = indx_sum_matrix + indx_list{i};
end

raw_erosion = single(indx_sum_matrix > 0);

level_volume = zeros(vx, vy, vz, 2); % 3rd dimension is (start_scale, end_scale)

for i = 1 : vx
    for j = 1 : vy
        for k = 1 : vz
            for l = scale_len : -scale_interval : 1
                bool_mat = indx_list{l};
                if bool_mat(i, j, k) == 1 && level_volume(i, j, k, 2) == 0
                     level_volume(i, j, k, 2) = sigma_head + l - 1;
                     break;
                end
            end
            for m = 1 : scale_interval : scale_len
                bool_mat = indx_list{m};
                if bool_mat(i, j, k) == 1 && level_volume(i, j, k, 1) == 0
                    level_volume(i, j, k, 1) = sigma_head + m - 1;
                    break;
                end
            end
        end
    end
end


curvature_scale_tensor = zeros(vx, vy, vz, scale_len, 4); % (k_min, k_max, start_scale, end_scale)

for i = 1 : vx
    for j = 1 : vy
        for k = 1 : vz
            for l = 1 : scale_len
                curvature_scale_tensor(i,j,k,l,1) = principle_list{l}(i,j,k,1);
                curvature_scale_tensor(i,j,k,l,2) = principle_list{l}(i,j,k,2);
                curvature_scale_tensor(i,j,k,l,3) = level_volume(i,j,k,1);
                curvature_scale_tensor(i,j,k,l,4) = level_volume(i,j,k,2);
            end
        end
    end
end

end



function w = FilterSize(s)
w = 2 * ceil(2 * s) + 1;
end


function sigma = BiSection(vol, head, tail, lower_val, upper_val)

if (tail - head < 1)
    sigma = (head + tail) / 2;
else
    mid = (head + tail) / 2;
    [erosion, ~] = GPCurvature(vol, mid, FilterSize(mid));
    ero_num = sum(sum(sum(erosion)));
    if (ero_num >= lower_val) && (ero_num <= upper_val)
        sigma = mid;
    elseif ero_num > upper_val
        sigma = BiSection(vol, mid, tail, lower_val, upper_val);
    else
        sigma = BiSection(vol, head, mid, lower_val, upper_val);
    end
end

end


