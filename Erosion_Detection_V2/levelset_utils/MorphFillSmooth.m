function pro = MorphFillSmooth(vol, kernel, iter, kernel_type)
%% Parameter Setting
timestep = 5;
mu = 0.2 / timestep; % Weight of the distance regulartization term
inner_iter = 1;
outer_iter = 1; % 15 1st
lambda = 0.4; % Weight of the length term
alpha = 1; % Wetight of the data fitting term
epsilon = 2; % Width of the Dirac Delta function and Heaviside function

[vx, vy, vz] = size(vol);

if kernel_type == "cuboid"
    se = strel('cuboid', [kernel, kernel, kernel]);
else
    se = strel('sphere', floor(kernel / 2));
end

%se = strel('cuboid', [kernel, kernel, kernel]);
process_vol = imclose(vol, se);
process_vol = imfill(process_vol, 'holes');

for sl = 1 : vz
    tem = process_vol(:, :, sl);
    tem = imfill(tem, 'holes');
    process_vol(:, :, sl) = tem;
end

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
phi = distant_vol;

if iter == 0
    pro = process_vol;
else
    phi = LSEdgeSmooth(phi, process_vol, mu, lambda, alpha, epsilon, timestep, iter);

    pro = single(phi > 0);
    for sl = 1 : vz
        tem = pro(:,:,sl);
        tem = imfill(tem, 'holes');
        pro(:,:,sl) = tem;
    end
end


end 