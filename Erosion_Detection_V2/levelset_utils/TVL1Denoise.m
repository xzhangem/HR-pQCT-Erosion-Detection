function result = TVL1Denoise(vol, iter, lambda, epsilon, timestep)
% Input: 
%   vol: 3D bone volume;
%   iter: iteration number for optimization;
%   lambda: coefficient of data term;
%   epsilon: small number for numerial issues;
%   timestep: step length for iteration;
% Output:
%   result: denoised 3D bone volume;

[vx, vy, vz] = size(vol);
u = vol;
vol = NeumannBoundCond(vol);
[vol_x, vol_y, vol_z] = gradient(vol);

for it = 1 : iter
    fprintf('Volume smoothing with TVL1 (%d of %d in total).\n', it, iter);
    [u_x, u_y, u_z] = gradient(u);
    s = sqrt(u_x.^2 + u_y.^2 + u_z.^2);
    N_x = u_x ./ (s + epsilon);
    N_y = u_y ./ (s + epsilon);
    N_z = u_z ./ (s + epsilon);
    m_curvature = Div(N_x, N_y, N_z);
    sigma = int3D((u - vol).^2);
    adp_lambda = int3D(s - vol_x .* N_x - vol_y .* N_y - vol_z .* N_z);
    adp_lambda = -adp_lambda / (sigma + epsilon);
    if abs(adp_lambda) < lambda
        lam = adp_lambda;
    else
        lam = lambda;
    end
    u = u + timestep * (m_curvature - lam * (u - vol)); 
end
result = u;
end

function f = Div(x, y, z)
[x_x, ~, ~] = gradient(x);
[~, y_y, ~] = gradient(y);
[~, ~, z_z] = gradient(z);
f = x_x + y_y + z_z;
end


function f = int3D(x)
f = sum(sum(sum(x)));
end


function g = NeumannBoundCond(f)
[nx, ny, nz] = size(f);
g = f;
g([1 nx], [1 ny], [1 nz]) = g([3 nx-2], [3, ny-2], [3, nz-2]);
g([1 nx], 2:end-1, 2:end-1) = g([3 nx-2], 2:end-1, 2:end-1);
g(2:end-1, [1 ny], 2:end-1) = g(2:end-1, [3 ny-2], 2:end-1);
g(2:end-1, 2:end-1, [1 nz]) = g(2:end-1, 2:end-1, [3 nz-2]);
end