function up_crop = UpHeight(vol)
[~, ~, init_h] = VolSep(vol);

[vx, vy, vz] = size(vol);
c_vol = vol(2 : vx-1, 2 : vy-1, 2 : vz-1);
c_vol = MaxConRegion(c_vol);
c_vol = MaxConFill(c_vol);
c_vol = MaxConFill(c_vol);

h_map = HeightLevel(init_h, c_vol, 1, 10, 0.1, 20000);

[vx, vy, vz] = size(c_vol);

up_crop = zeros(vx, vy, vz);

for i = 1 : vx
    for j = 1 : vy
        for k = round(h_map(i,j)) : vz
            up_crop(i, j, k) = c_vol(i, j, k);
        end
    end
end

up_crop = MaxConRegion(up_crop);

end