function detect_vol = ErosionDetection(refine_seg, input_sigma, input_radium)
addpath('levelset_utils');
addpath('curvature_utils');

%% Parameter Setting
% Erosion Detection Parameter Part
sigma = input_sigma; % scale space parameter, default 4
clust_radium = round(input_radium); % cluster parameter, default 5
%iso_r = 1; % radium for determining isolated points need to be erased
erosion_min = 120; % threshold for erasing extremely small erosions
ratio_thresh = 11; % threshold ratio for curvature threshold

sqrt_2 = sqrt(2);
filtersize = FilterSize(sigma);

[vx, vy, vz] = size(refine_seg);

[erosion] = GPCurvature(refine_seg, sigma, filtersize);
erosion = VolumeTail(erosion, 5);
detect_ellp = EllipsoidShape(erosion, clust_radium, erosion_min, ratio_thresh);
detect_ellp = single(detect_ellp > 0);
[detect_vol, ~] = bwlabeln(detect_ellp, 26);
detect_vol = single(detect_vol);

end

function w = FilterSize(s)
w = 2 * ceil(2 * s) + 1;
end


function t_vol = VolumeTail(vol, bound_len)
[vx, vy, vz] = size(vol);
vol_mask  = zeros(vx, vy, vz);
vol_mask(bound_len:vx-bound_len, bound_len:vy-bound_len, bound_len:vz-bound_len) = 1;
t_vol = vol .* vol_mask;
end

