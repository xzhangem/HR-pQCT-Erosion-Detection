clc;
clear all;
close all;

detect_file = '../detect_erosion_thr_25/';
erosion_file = '../erosion_label/';

IOU_threshold = 0.2;

erosion_dir = dir(fullfile(erosion_file, '*.nii'));
erosion_num = length(erosion_dir);

ROC_matrix = zeros(erosion_num, 4); % TP, FN, FP, TN

for i = 1 : erosion_num
    name = erosion_dir(i).name;
    split_name = strsplit(name, '_e.');
    erosion_name = [erosion_file, name];
    detect_name = [detect_file, split_name{1}, '_o_Gamma_ero.nii'];
    erosion_vol = niftiread(erosion_name);
    detect_vol = niftiread(detect_name);
    
    erosion_set = unique(erosion_vol);
    erosion_set = erosion_set(2:end);
    
    detect_set = unique(detect_vol);
    if length(detect_set) == 1
        ROC_matrix(i, 1) = 0; ROC_matrix(i, 2) = length(erosion_set);
        ROC_matrix(i, 3) = 0; ROC_matrix(i, 4) = 0;
    else
        detect_set = detect_set(2:end);
        e_num = length(erosion_set); d_num = length(detect_set);
        IOU_matrix = zeros(e_num, d_num);
        for j = 1 : e_num
            erosion_mask = single(erosion_vol == erosion_set(j));
            for k = 1 : d_num
                detect_mask = single(detect_vol == detect_set(k));
                cap_mask = single(erosion_mask & detect_mask); 
                cup_mask = single(erosion_mask | detect_mask);
                IOU = sum(sum(sum(cap_mask))) / sum(sum(sum(cup_mask)));
                IOU_matrix(j, k) = IOU;
            end
        end
        fprintf('%d th volume. \n', i);
        
        for k = 1 : d_num
            col_max = max(IOU_matrix(:,k));
            IOU_col_mask = single(IOU_matrix(:,k) == col_max);
            IOU_matrix(:, k) = IOU_matrix(:, k) .* IOU_col_mask;
        end
        IOU_mask = single(IOU_matrix > IOU_threshold);
        
        
        for k = 1 : e_num
            row_max = max(IOU_matrix(k,:));
            IOU_row_mask = single(IOU_matrix(k,:) == row_max);
            IOU_matrix(k, :) = IOU_matrix(k, :) .* IOU_row_mask;
        end
        
        IOU_matrix
        
        ROC_matrix(i, 1) = sum(sum(sum(IOU_mask))); % TP
        ROC_matrix(i, 2) = e_num - ROC_matrix(i, 1); % FN
        ROC_matrix(i, 3) = sum(single(sum(IOU_mask) == 0)); % FP
    end
end

total_vec = sum(ROC_matrix);
total_TP = total_vec(1); total_FN = total_vec(2); total_FP = total_vec(3);
precision = total_TP / (total_TP + total_FP);
recall = total_TP / (total_TP + total_FN);
fprintf('precision: %.3f; recall: %.3f.\n', precision, recall)