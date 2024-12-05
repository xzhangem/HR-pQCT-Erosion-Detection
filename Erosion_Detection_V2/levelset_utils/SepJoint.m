function sep_vol = SepJoint(vol)
[vx, vy, vz] = size(vol);
CONNECT_BONE = 0;

[con_lab, ~] = bwlabeln(vol, 26);
stats = regionprops(con_lab);
area = cat(1, stats.Area);
[~, ind] = sort(area, 'descend');

if length(ind) < 2
    CONNECT_BONE = 1;
else
    part_bone_1 = single(con_lab == ind(1));
    part_bone_2 = single(con_lab == ind(2));
    z_sum_1 = sum(sum(part_bone_1));
    z_sum_2 = sum(sum(part_bone_2));
    zerosum_1 = single(z_sum_1 == 0);
    zerosum_2 = single(z_sum_2 == 0);
    if (sum(zerosum_1(:)) == 0) || (sum(zerosum_2(:)) == 0)
        CONNECT_BONE = 1;
    end
end

if CONNECT_BONE == 0
    [~, ~, part_z_1] = MaskCent(part_bone_1);
    [~, ~, part_z_2] = MaskCent(part_bone_2);
    if part_z_1 > part_z_2
        bone_up = part_bone_1;
    else
        bone_up = part_bone_2;
    end
    up_crop = bone_up;
else
    fill_vol = IsoDiffusionFill(vol, 2, 1); %6
    [~, ~, init_hmap] = VolSep(fill_vol);
    h_map = HeightLevel(init_hmap, fill_vol, 1, 8, 0.1, 20000);
    up_crop = zeros(vx, vy, vz);
    for i = 1 : vx
        for j = 1 : vy
            for k = round(h_map(i,j)) : vz
                up_crop(i, j, k) = vol(i, j, k);
            end
        end
    end
end

sep_vol = up_crop;
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

