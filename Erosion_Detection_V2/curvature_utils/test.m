clc;
clear all;
close all;

th_ratio = 0.008; % threshold ratio determining the last scale factor for significant erosion searching
init_sigma = 0.5; % initial sigma for multiscale searching
iso_r = 1; % radium for determining isolated points need to be erased
adj_radium = 20; % radium for neighborhood searching

num_th = 250;

refine_seg = niftiread('../Gamma_refine/3243_M25616_MCP3_bl_o_Gamma.nii');

[start_scale, raw_erosion, curvature_scale_tensor] = GaussPyramidSearching(refine_seg, th_ratio, init_sigma);
[raw_erosion, curvature_scale_tensor] = DelIso(raw_erosion, curvature_scale_tensor, iso_r);
[label_erosion, erosion_num] = SpatialCluster(raw_erosion, adj_radium); 
erosion_coe_cell = ClusterIndx(label_erosion, erosion_num);

num_per_erosion = zeros(erosion_num, 1);
for i = 1 : erosion_num
    a = size(erosion_coe_cell{i});
    num_per_erosion(i, 1) = a(1);
end

%num_per_erosion = single(num_per_erosion > num_th) .* num_per_erosion;

curvature_ratio_cell = cell(erosion_num, 1);
curvature_length_cell = cell(erosion_num, 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : erosion_num
    if num_per_erosion(i,1) > 0
        curvature_ratio_cell{i} = zeros(num_per_erosion(i,1), 1);
        curvature_length_cell{i} = zeros(num_per_erosion(i,1), 1);
    end
end

cluster_counter = ones(erosion_num, 1);

for i = 1 : erosion_num
    if num_per_erosion(i,1) > 0
        length = num_per_erosion(i, 1);
        coe_vector = erosion_coe_cell{i};
        for j = 1 : length
            x = coe_vector(j,1); y = coe_vector(j,2); z = coe_vector(j,3);
            if curvature_scale_tensor(x, y, z, 1, 3) == start_scale
                curvature_ratio_cell{i}(j, 1) = curvature_scale_tensor(x,y,z,1,1) / curvature_scale_tensor(x,y,z,1,2);
                curvature_length_cell{i}(j, 1) = (curvature_scale_tensor(x,y,z,1,1)^2 + ...
                    curvature_scale_tensor(x,y,z,1,2)^2)^0.5;
            end
        end
        curvature_ratio_cell{i}(find(curvature_ratio_cell{i} == 0)) = [];
        curvature_length_cell{i}(find(curvature_length_cell{i} == 0)) = [];
    end
end

