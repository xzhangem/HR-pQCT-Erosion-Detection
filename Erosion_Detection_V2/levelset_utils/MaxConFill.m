function fill_up_vol = MaxConFill(up_vol)
[vx, vy, vz] = size(up_vol);
% se = strel('sphere', 17); % 9? 11
% full_upcrop = single(imclose(up_vol, se));
full_upcrop = IsoDiffusionFill(up_vol, 26, 1);
for i = 1 : vz
    im = squeeze(full_upcrop(:, :, i));
    im = imfill(im);
    full_upcrop(:, :, i) = single(im);
end
up_res = full_upcrop - up_vol;

contour_se = strel('cuboid', [9,9,9]); %9
up_res = single(imopen(up_res, contour_se));
[con_lab, ~] = bwlabeln(up_res, 26);
%niftiwrite(con_lab, 'new_up_conlab.nii');
stats = regionprops(con_lab);
area = cat(1, stats.Area);
[~, ind] = sort(area, 'descend');

if length(ind) < 1
    part = zeros(vx, vy, vz);
else
    part = single(con_lab == ind(1));
end

fill_up_vol = single((part + up_vol) > 0);
fill_up_vol = single(imclose(fill_up_vol, contour_se));
for i = 1 : vz
    im = squeeze(fill_up_vol(:, :, i));
    im = imfill(im);
    fill_up_vol(:, :, i) = single(im);
end

end