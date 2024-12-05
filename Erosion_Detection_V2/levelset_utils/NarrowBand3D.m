function [nb_indx_x, nb_indx_y, nb_indx_z] = NarrowBand3D(phi, r)
l = 2 * r + 1;
[vx, vy, vz] = size(phi);
%fprintf('[%d, %d, %d]\n', vx, vy, vz);
zero_indx = zeros(vx, vy, vz);
nrwbnd_indx = zeros(vx, vy, vz);
zero_indx(2:vx-1, 2:vy-1, 2:vz-1) = (phi(1:vx-2, 2:vy-1, 2:vz-1) ...
.* phi(3:vx, 2:vy-1, 2:vz-1) < 0) + (phi(2:vx-1, 1:vy-2, 2:vz-1) .* phi(2:vx-1, 3:vy, 2:vz-1) < 0) ...
+ (phi(2:vx-1, 2:vy-1, 1:vz-2) .* phi(2:vx-1, 2:vy-1, 3:vz) < 0);

for i = 2 : vx-1
    for j = 2 : vy-1
        for k = 2 : vz-1
            if (zero_indx(i, j, k) ~= 0)
                nrwbnd_indx(i-r:i+r, j-r:j+r, k-r:k+r) = ones(l,l,l);
            end
        end
    end
end

nb_indx_x = [];
nb_indx_y = [];
nb_indx_z = [];


for i = 1 : vx
    for j = 1 : vy
        for k = 1 : vz
            if (nrwbnd_indx(i, j, k) ~= 0)
                nb_indx_x = [nb_indx_x; i];
                nb_indx_y = [nb_indx_y; j];
                nb_indx_z = [nb_indx_z; k];
            end
        end
    end
end

end