clc;
clear all;
close all;

detect_file = '../detect_erosion_thr_10/';
seg_file = '../Gamma_refine/';

seg_dir = dir(fullfile(seg_file, '*.nii'));
seg_num = length(seg_dir);

for i = 1 : seg_num
    name = seg_dir(i).name;
    split_name = strsplit(name, '.');
    seg_name = [seg_file, name];
    detect_name = [detect_file, split_name{1}, '_ero.nii'];
    seg_vol = niftiread(seg_name);
    ellp_vol = niftiread(detect_name);
    contour = single(bwperim(seg_vol));
    [vx, vy, vz] = size(seg_vol);
    erosion_set = unique(ellp_vol); 
    erosion_set = sort(erosion_set);
    convhull_vol = zeros(vx, vy, vz);
    if (length(erosion_set) ~= 1)
        erosion_set = erosion_set(2:end);
        erosion_num = length(erosion_set);
        for j = 1 : erosion_num
            region_part = single(ellp_vol == erosion_set(j));
            contour_part = region_part .* contour;
            contour_num = sum(sum(sum(contour_part)));
            x_vector = zeros(contour_num, 1); y_vector = zeros(contour_num, 1);
            z_vector = zeros(contour_num, 1);
            count = 1;
            for p = 1 : vx
                for q = 1 : vy
                    for r = 1 : vz
                        if contour_part(p, q, r) == 1
                            x_vector(count) = p; y_vector(count) = q; 
                            z_vector(count) = r; count = count + 1;
                        end
                    end
                end
            end
            K = convhull(x_vector, y_vector, z_vector);
            conv_erosion = zeros(vx, vy, vz);
            for l = 1 : length(K)
                conv_erosion(K(l,1), K(l,2), K(l,3)) = 1;
            end
            %conv_erosion = imfill(conv_erosion, 'holes');
            conv_erosion = conv_erosion * erosion_set(j);
            convhull_vol = convhull_vol + conv_erosion;
        end
    end
end

