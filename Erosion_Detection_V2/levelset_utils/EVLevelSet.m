function phi = EVLevelSet(phi_0, vol, g, mu, lambda, alpha, beta_v, gamma_v, epsilon, timestep, iter)
% Input:
%   phi_0: initial level set function in this stage
%   vol: double volume for updating external enengy
%   g: 3D surface indicator function
%   mu: weight of the distance regularization term
%   lambda: weight of the weighted length term
%   alpha: weight of the weighted area term
%   beta: weight of the weighted internal Gumbel Distribution
%   gamma: weight of the weighted external Gaussian
%   epsilon: width of Dirac Delta function and Heaviside function
%   timestep: time step
%   iter: number of iteration 
%   bool_global: bool value. 1 for Global lvlset, 0 for local evolution.

% Output:
%   phi: the updated level set function

% This script models the bone region as the Extreme Value distribution and use MLE
% to estimate the location and scale. Brent iterative updating is used for 
% solving the equation. Scale parameter 'sc' and location
% parameter 'lc' 

vol = double(vol);


phi = double(phi_0);
[vx, vy, vz] = gradient(g);
[lx, ly, lz] = size(vol);

for i = 1 : iter
    phi = NeumannBoundCond(phi);
    [phi_x, phi_y, phi_z] = gradient(phi);
    s = sqrt(phi_x.^2 + phi_y.^2 + phi_z.^2);
    small_number = 1e-7;
    Nx = phi_x ./ (s + small_number);
    Ny = phi_y ./ (s + small_number);
    Nz = phi_z ./ (s + small_number);
    curvature = Div(Nx, Ny, Nz);
    dist_regu = DistReg(phi);
    dirac_phi = SmoothDirac(phi, epsilon);
    heavi_phi = SmoothHeavi(phi, epsilon);
    area_term = dirac_phi .* g;
    edge_term = dirac_phi .* (vx .* Nx + vy .* Ny + vz .* Nz) + dirac_phi .* g .* curvature;
    c_ex = (int3D(vol .* (1 - heavi_phi))) / (int3D(1 - heavi_phi));
    sigma_ex = (int3D((vol - c_ex).^2 .* (1 - heavi_phi))) / (int3D(1 - heavi_phi));
    ext_term = ((vol - c_ex).^2 / (2 * (sigma_ex + small_number)) + 0.5 * log(sigma_ex));
    
    % Estimate location shift in the area heavi_phi > 0
    
    fore_mask = single(heavi_phi > 0);
    param = EVParam(vol, fore_mask);
    lc = param(1); sc = param(2); 
    fprintf('location & scale: %f, %f. \n', lc, sc);
    
    %M = 100;
    %phi_mask = single((phi >= -M) & (phi <= M));
    
    %inner_log_vol = LogVol(1 + sh * (vol - lc) / sc);
    inter_term = log(sc) - (vol - lc) / sc + exp((vol - lc) / sc);
    phi = phi + timestep * (mu * dist_regu + lambda * edge_term - alpha * area_term - beta_v * inter_term + ...
        gamma_v * ext_term);
    
end

end


function f = Div(x, y, z)
[x_x, ~, ~] = gradient(x);
[~, y_y, ~] = gradient(y);
[~, ~, z_z] = gradient(z);
f = x_x + y_y + z_z;
end

function f = DistReg(phi)
[phi_x, phi_y, phi_z] = gradient(phi);
s = sqrt(phi_x.^2 + phi_y.^2 + phi_z.^2);
a = (s>=0) & (s<=1);
b = (s>1);
ps = a .* sin(2 * pi * s) / (2 * pi) + b .* (s-1);
dps=((ps~=0) .* ps + (ps==0)) ./ ((s~=0) .* s + (s==0));
f = Div(dps .* phi_x - phi_x, dps .* phi_y - phi_y, dps .* phi_z - phi_z) + 4 * del2(phi);
end

function f = SmoothDirac(phi, epsilon)
f = (1 / (2 * epsilon)) * (1 + cos(pi * phi / epsilon));
b = (phi <= epsilon) & (phi >= -epsilon);
f = f .* b;
end

function f = SmoothHeavi(phi, epsilon)
f = (1 / 2) * (1 + phi / epsilon + (1 / pi) * sin(pi * phi / epsilon));
b = (phi <= epsilon) & (phi >= -epsilon);
a = (phi > epsilon);
c = (phi < -epsilon);
f = f .* b + 1. * a + 0. * c;
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

function beta = BetaEst(k)
syms x;
beta = vpasolve(k * (gamma(x + 1))^2 == gamma(1 + 2 * x));
beta = 1 / beta;
end

function g = LogVol(x)
[vx, vy, vz] = size(x);
g = zeros(vx, vy, vz);

for i = 1 : vx
    for j = 1 : vy
        for k = 1 : vz
            if x(i, j, k) <= 0 
                g(i, j, k) = 0;
            else
                g(i, j, k) = log(x(i, j, k));
            end
        end
    end
end

end

function para = EVParam(vol, mask)
[vx, vy, vz] = size(vol);
num = sum(sum(sum(mask)));
vol_vec = zeros(num, 1);
tag = 1;
for i = 1 : vx
    for j = 1 : vy
        for k = 1 : vz
            if mask(i, j, k) == 1
                vol_vec(tag) = vol(i,j,k);
                tag = tag + 1;
            end
        end
    end
end
para = gevfit(vol_vec);
end