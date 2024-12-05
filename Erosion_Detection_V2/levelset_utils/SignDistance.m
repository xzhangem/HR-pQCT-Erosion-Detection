function distant_vol = SignDistance(process_vol)

[vx, vy, vz] = size(process_vol);
contour = single(bwperim(process_vol));
distant_vol = bwdist(contour);
for i = 1 : vx
    for j = 1 : vy
        for k = 1 : vz
            if process_vol(i,j,k) == 0
                distant_vol(i,j,k) = -distant_vol(i,j,k);
            end
        end
    end
end

end