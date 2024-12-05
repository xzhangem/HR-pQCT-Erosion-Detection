function h_map = HeightLevel(init_h, seg, mu, lambda, timestep, iter_num)
% \mu \int_{\Omega} h(x,y) I(x, y, h(x,y)) dx dy + \lambda \int_{\Omega}
% |\nabla h(x,y)| dx dy
% \lambda should be larger to ensure the height should be at similar level

h_map = init_h;
[vx, vy, vz] = size(seg);
[~, ~, gz] = gradient(seg);
for it = 1 : iter_num
    h_map = NBCond(h_map);
    for i = 1 : vx
        for j = 1 : vy
            if h_map(i, j) < 1
                h_map(i, j) = 1;
            elseif h_map(i, j) > vz
                h_map(i, j) = vz;
            else
                h_map(i, j) = h_map(i, j);
            end
        end
    end
    [gh_x, gh_y] = gradient(h_map);
    s = sqrt(gh_x.^2 + gh_y.^2);
    small_num = 1e-9;
    N_ghx = gh_x ./ (s + small_num);
    N_ghy = gh_y ./ (s + small_num);
    curvature = Div(N_ghx, N_ghy);    
    h_map_floor = floor(h_map);
    h_map_ceil = ceil(h_map);
    gz_h = zeros(vx, vy);
    seg_h = zeros(vx, vy);
    for i = 1 : vx
        for j = 1 : vy
            gz_h(i, j) = (h_map(i, j) - h_map_floor(i, j)) * gz(i, j, h_map_floor(i,j)) + (h_map_ceil(i,j) - h_map(i, j)) * ...
                gz(i, j, h_map_ceil(i, j));
            seg_h(i, j) = (h_map(i, j) - h_map_floor(i, j)) * seg(i, j, h_map_floor(i,j)) + (h_map_ceil(i,j) - h_map(i, j)) * ...
                seg(i, j, h_map_ceil(i, j));
        end
    end
    h_map = h_map - timestep * (mu * (h_map .* gz_h + seg_h) - lambda * curvature);
end
for i = 1 : vx
   for j = 1 : vy
        if h_map(i, j) < 1
            h_map(i, j) = 1;
        elseif h_map(i, j) > vz
            h_map(i, j) = vz;
        else
            h_map(i, j) = h_map(i, j);
        end
   end
end
end

function f = Div(x, y)
[x_x, ~] = gradient(x);
[~, y_y] = gradient(y);
f = x_x + y_y;
end

function g = NBCond(f)
[nx, ny] = size(f);
g = f;
g([1 nx], [1 ny]) = g([3 nx-2], [3, ny-2]);
g([1 nx], 2:end-1) = g([3 nx-2], 2:end-1);
g(2:end-1, [1 ny]) = g(2:end-1, [3 ny-2]);
g(2:end-1, 2:end-1) = g(2:end-1, 2:end-1);
end
