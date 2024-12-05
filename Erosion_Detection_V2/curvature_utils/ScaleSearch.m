function [scale_interval, volume_set] = ScaleSearch(vol, init_sigma, thresh_ratio, N)

%% Parameter Setting
% Input:
%   vol: input binary segmentation volume
%   init_sigma: initial sigma for searching 
%   thresh_ratio: threshold ratio determining the last scale factor
%   N: maximum trace-back iteration 

% Output:
%   scale_interval: plausible scale set;
%   volume_set: binary volume set on the scale set 

[vx, vy, vz] = size(vol);
contour = single(bwperim(vol));
contour_sum = sum(sum(sum(contour)));
resolution_thresh = contour_sum * thresh_ratio;

n = 0;
tag = true;
sqrt_2 = sqrt(2);

while (tag)
    sigma = init_sigma * sqrt_2^n;
    filtersize = FilterSize(sigma);
    [erosion, ~] = GPCurvature(vol, sigma, filtersize);
    
    
    if n == 0
        erosion_set = erosion;
    else
        erosion_set = cat(4, erosion_set, erosion);
    end
    
    fprintf('Scale: %d ------ ratio: %f \n', sigma, sum(sum(sum(erosion))) / contour_sum);
    if sum(sum(sum(erosion))) > resolution_thresh
        n = n + 1;
    else
        tag = false;
    end    
end

indx_max = n - 1;
sigma_max = init_sigma * sqrt_2^(indx_max);
back_tag = true;
t = 0;
N = min(N, n-1);
while (t < N & back_tag)
    now_indx = indx_max - t;
    now_erosion = erosion_set(:, :, :, now_indx+1);
    pre_erosion = erosion_set(:, :, :, now_indx);
    now_sum = sum(sum(sum(now_erosion)));
    pre_sum = sum(sum(sum(pre_erosion)));
    if pre_sum > now_sum
        t = t + 1;
    else
        back_tag = false;
    end
end

inverse_t_interval = -t : 0;

if t >= N
    inverse_t_interval = inverse_t_interval(end - N + 1: end);
end

scale_interval = init_sigma * sqrt_2.^(indx_max + inverse_t_interval);

volume_set = erosion_set(:, :, :, indx_max + inverse_t_interval);
end


function w = FilterSize(s)
w = 2 * ceil(2 * s) + 1;
end