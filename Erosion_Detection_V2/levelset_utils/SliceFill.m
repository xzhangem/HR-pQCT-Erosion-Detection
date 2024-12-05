function [up_vol, new_vol] = SliceFill(vol)
%addpath('../levelset_utils');

[vx, vy, vz] = size(vol);
up_vol = single(vol);
%se = strel('cuboid', [7, 7, 7]);
%up_vol = single(imclose(vol, se));
for i = 1 : vz
    im = squeeze(up_vol(:, :, i));
    im = imfill(im);
    up_vol(:, :, i) = single(im);
end

xyfill_vol = up_vol;

for i = 1 : vx
    im = squeeze(up_vol(i, :, :));
    im = imfill(im);
    up_vol(i, :, :) = single(im);
end

for i = 1 : vy
    im = squeeze(up_vol(:, i, :));
    im = imfill(im);
    up_vol(:, i, :) = single(im);
end

for i = 1 : vz
    im = squeeze(up_vol(:, :, i));
    im = imfill(im);
    up_vol(:, :, i) = single(im);
end

new_vol = single(up_vol - xyfill_vol);

end