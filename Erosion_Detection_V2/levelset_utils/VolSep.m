function [up_crop, down_crop, nozero_maxh] = VolSep(input_vol)
%addpath('../levelset_utils');

[vx, vy, vz] = size(input_vol);
%origin_vol = input_vol(2 : vx - 1, 2 : vy - 1, 2 : vz - 1);
origin_vol = input_vol;

vol = origin_vol;
[vx, vy, vz] = size(vol);

%for i = 1 : vz
%    vol(:, :, i) = single(imfill(vol(:, :, i)));
%end
sep_tag = 0;
sep_check = 1;

while (sep_tag == 0)
    %sep_check
    %sep_tag
    if sep_check == 1
        test_vol = vol;
    else
        test_se = strel('cuboid', [2*sep_check-1, 2*sep_check-1, 2*sep_check-1]);
        test_vol = single(imopen(vol, test_se));
    end
    [test_label, ~] = bwlabeln(test_vol, 26);
    test_stats = regionprops(test_label);
    test_area = cat(1, test_stats.Area);
    [~, test_ind] = sort(test_area, 'descend');
    %length(test_ind)
    if length(test_ind) == 1
        sep_check = sep_check + 1;
    else
        test_bone_1 = single(test_label == test_ind(1));
        test_bone_2 = single(test_label == test_ind(2));
        z_sum_1 = sum(sum(test_bone_1));
        z_sum_2 = sum(sum(test_bone_2));
        zerosum_1 = single(z_sum_1 == 0);
        zerosum_2 = single(z_sum_2 == 0);
        if (sum(zerosum_1(:)) == 0) || (sum(zerosum_2(:)) == 0)
            sep_check = sep_check + 1;
        else
            [~, ~, test_z_1] = MaskCent(test_bone_1);
            [~, ~, test_z_2] = MaskCent(test_bone_2);
            if test_z_1 > test_z_2
                bone_up = test_bone_1;
                bone_down = test_bone_2;
            else
                bone_up = test_bone_2;
                bone_down = test_bone_1;
            end
            sep_tag = 1;
        end
    end
end


%{
se = strel('cuboid', [17, 17, 17]);
vol = single(imopen(vol, se));
size(vol)
[con_lab, num] = bwlabeln(vol, 26);
stats = regionprops(con_lab);
area = cat(1, stats.Area);
[sort_area, ind] = sort(area, 'descend');

bone_1 = single(con_lab == ind(1));
bone_2 = single(con_lab == ind(2));

[~, ~, z_1] = MaskCent(bone_1);
[~, ~, z_2] = MaskCent(bone_2);

if z_1 > z_2
    bone_up = bone_1;
    bone_down = bone_2;
else
    bone_up = bone_2;
    bone_down = bone_1;
end
%}

max_h = zeros(vx, vy);
min_h = zeros(vx, vy);

for i = 1 : vx
    for j =  1 : vy
        max_h(i, j) = max(squeeze(bone_down(i, j, :)) .* single((1 : vz)'));
        upz_vec = squeeze(bone_up(i, j, :)) .* single((1 : vz)');
        if sum(upz_vec) == 0
            min_h(i, j) = vz;
        else
            upz_vec(find(upz_vec == 0)) = [];
            min_h(i, j) = min(upz_vec);
        end
    end
end

max_mask = single(max_h > 0);
contour = single(bwperim(max_mask));
[D, idx] = bwdist(contour);

%indx_x = zeros(vx, vy);
%indx_y = zeros(vx, vy);

[row, col] = ind2sub(size(max_h), idx);

nozero_maxh = max_h;

for i = 1 : vx
    for j = 1 : vy
        if max_h(i, j) == 0
            x_ind = row(i, j);
            y_ind = col(i, j);
            nozero_maxh(i, j) = max_h(x_ind, y_ind);
        end
    end
end

res_h = min_h - max_h;
min_gap = ceil(min(res_h(:)) / 2);
nozero_maxh = nozero_maxh + min(2, min_gap);

down_crop = zeros(vx, vy, vz);
for i = 1 : vx
    for j = 1 : vy
        for k = 1 : nozero_maxh(i,j)
            down_crop(i, j, k) = origin_vol(i, j, k);
        end
    end
end

up_crop = zeros(vx, vy, vz);
for i = 1 : vx
    for j = 1 : vy
        for k = (nozero_maxh(i,j)+1) : vz
            up_crop(i, j, k) = origin_vol(i, j, k);
        end
    end
end
%niftiwrite(up_crop, 'up_crop.nii');

%{
res_h = min_h - max_h;
min_gap = ceil(min(res_h(:)) / 2);

niftiwrite(max_h, 'max_h.nii');

mean_sep = (min_h + max_h) / 2;
for i = 1 : vx
    for j =  1 : vy
        if (min_h(i, j) == vz) && (max_h(i, j) ~= 0)
            mean_sep(i, j) = max(max_h(:)) + min_gap;
        elseif (max_h(i, j) == 0) && (min_h(i, j) ~= vz)
            temp_vec = min_h;
            temp_vec(find(temp_vec == 0)) = [];
            mean_sep(i, j) = min(temp_vec(:)) - min_gap;
        else
            mean_sep(i, j) = mean_sep(i, j);
        end
    end
end

down_crop = zeros(vx, vy, vz);
for i = 1 : vx
    for j = 1 : vy
        for k = 1 : floor(mean_sep(i, j))
            down_crop(i, j, k) = origin_vol(i, j, k);
        end
    end
end

up_crop = zeros(vx, vy, vz);
for i = 1 : vx
    for j = 1 : vy
        for k = ceil(mean_sep(i, j)) : vz
            up_crop(i, j, k) = origin_vol(i, j, k);
        end
    end
end
%}
end

function [cx, cy, cz] = MaskCent(mask)
[vx, vy, vz] = size(mask);
[X, Y, Z] = meshgrid(1:vx, 1:vy, 1:vz);
X = permute(X, [2, 1, 3]);
Y = permute(Y, [2, 1, 3]);
Z = permute(Z, [2, 1, 3]);
mask_X = X .* mask;
mask_Y = Y .* mask;
mask_Z = Z .* mask;
cx = sum(mask_X(:)) / sum(mask(:));
cy = sum(mask_Y(:)) / sum(mask(:));
cz = sum(mask_Z(:)) / sum(mask(:));
end