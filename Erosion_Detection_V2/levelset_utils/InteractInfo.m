function [interact_vol, interact_total_num, interact_ratio, interact_slice_max, maxlength] = InteractInfo(fill_vol, dect_vol)
[vx, vy, vz] = size(fill_vol);

label_vec = unique(dect_vol);
label_vec = label_vec(2 : end);

interact_vol = zeros(vx, vy, vz);
fill_curve = single(bwperim(fill_vol));

for i = 1 : length(label_vec)
    dect_region = single(dect_vol == label_vec(i));
    dect_curve = single(bwperim(dect_region));
    interact_vol = interact_vol + label_vec(i) * single((dect_curve + fill_curve) > 1);
end

%niftiwrite(interact_vol, 'interact_vol.nii');

interact_total_num = zeros(length(label_vec), 1);
for i = 1 : length(label_vec)
    interact_total_num(i) = sum(sum(sum(single(interact_vol == label_vec(i)))));
end

interact_ratio = zeros(length(label_vec), 1);
for i = 1 : length(label_vec)
    interact_ratio(i) = interact_total_num(i) / sum(sum(sum(single(dect_vol == label_vec(i)))));
end

interact_slice_max = zeros(length(label_vec), 1);
maxlength = zeros(length(label_vec), 1);
for k = 1 : vz
    t_slice = squeeze(interact_vol(:, :, k));
    for i = 1 : length(label_vec)
        p_sum = sum(sum(single(t_slice == label_vec(i))));
        p_length = CurveLength(single(t_slice == label_vec(i)));
        if interact_slice_max(i) < p_sum
            interact_slice_max(i) = p_sum;
        end
        if maxlength(i) < p_length
            maxlength(i) = p_length;
        end
    end
end

end