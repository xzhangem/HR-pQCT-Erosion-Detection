clc;
clear all;
close all;

vol = niftiread('cvh.nii');
[vx, vy, vz] = size(vol);
mask = single(vol > 0);

voxel_num = sum(sum(sum(mask)));
erosion_set = unique(vol);
erosion_set = erosion_set(2:end);
erosion_set = sort(erosion_set);

erosion_num = length(erosion_set);
for j = 1 : erosion_num
    contour_part = single(vol == erosion_set(j));
    contour_num = sum(sum(sum(contour_part)));
    %x_vector = zeros(contour_num, 1); y_vector = zeros(contour_num, 1);
    %z_vector = zeros(contour_num, 1);
    cood_vector = zeros(contour_num, 3);
    count = 1;
    for p = 1 : vx
        for q = 1 : vy
            for r = 1 : vz
                if contour_part(p, q, r) == 1
                    %x_vector(count) = p; y_vector(count) = q; 
                    %z_vector(count) = r; count = count + 1;
                    cood_vector(count, 1) = p; cood_vector(count, 2) = q;
                    cood_vector(count, 3) = r; count = count + 1;
                end
            end
        end
    end
    [K, V] = convhull(cood_vector);
    conv_erosion = zeros(vx, vy, vz);
    convhull_point = cood_vector(K);
    for l = 1 : length(convhull_point)
        conv_erosion(convhull_point(l,1), convhull_point(l,2), convhull_point(l,3)) = 1;
    end
end


