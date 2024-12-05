function result_vol = IsoDiffusionFill(vol, timestep, iter_num)



pad_v = ceil(timestep) + 2;
vol = padarray(vol, [pad_v,pad_v,pad_v], 'symmetric', 'both');
[vx, vy, vz] = size(vol);
for iter = 1 : iter_num
    isosurface = SignDistance(vol);
    isosurface = NeumannBoundCond(isosurface);
    [g_x, g_y, g_z] = gradient(isosurface);
    norm_mat = (g_x.^2 + g_y.^2 + g_z.^2).^0.5;
    isosurface = isosurface + timestep * norm_mat;
    vol = single(isosurface > 0);
    for sl = 1 : vz
        tem = vol(:, :, sl);
        tem = imfill(tem, 'holes');
        vol(:, :, sl) = tem;
    end   
    
    isosurface = SignDistance(vol);
    [g_x, g_y, g_z] = gradient(isosurface);
    norm_mat = (g_x.^2 + g_y.^2 + g_z.^2).^0.5;
    isosurface = isosurface - (timestep * 1) * norm_mat;
    vol = single(isosurface > 0);
    for sl = 1 : vz
        tem = vol(:, :, sl);
        tem = imfill(tem, 'holes');
        vol(:, :, sl) = tem;
    end    
       
end


result_vol = vol(1+pad_v : vx-pad_v, 1+pad_v : vy-pad_v, 1+pad_v : vz-pad_v);

end

function g = NeumannBoundCond(f)
[nx, ny, nz] = size(f);
g = f;
g([1 nx], [1 ny], [1 nz]) = g([3 nx-2], [3, ny-2], [3, nz-2]);
g([1 nx], 2:end-1, 2:end-1) = g([3 nx-2], 2:end-1, 2:end-1);
g(2:end-1, [1 ny], 2:end-1) = g(2:end-1, [3 ny-2], 2:end-1);
g(2:end-1, 2:end-1, [1 nz]) = g(2:end-1, 2:end-1, [3 nz-2]);
end
