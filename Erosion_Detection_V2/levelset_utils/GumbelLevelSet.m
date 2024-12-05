function phi = GumbelLevelSet(phi_0, vol, g, mu, lambda, alpha, beta, gamma, epsilon, timestep, iter, bool_global)
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

% This script models the bone region as the Gumbel distribution and use MLE
% to estimate the location and scale. Brent iterative updating is used for 
% solving the equation. 

vol = double(vol);

phi = double(phi_0);
[vx, vy, vz] = gradient(g);
[lx, ly, lz] = size(vol);

for k = 1 : iter
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
    %ext_term = dirac_phi .* ((vol - c_ex).^2 / (2 * (sigma_ex + small_number)) + 0.5 * log(sigma_ex));
    ext_term = (vol - c_ex).^2 / (2 * (sigma_ex + small_number)) + 0.5 * log(sigma_ex);

    inter_mean = int3D(vol .* heavi_phi) / int3D(heavi_phi);
    [init_mu, init_sigma] = GumbelInit(vol, heavi_phi);
    %{
    x = sym('x');
    F = x - inter_mean + int3D(vol .* exp(-vol / x) .* heavi_phi) / ... 
    int3D(exp(-vol / x) .* heavi_phi);
    func = matlabFunction(F, 'Vars', x);
    sol_sigma = fsolve(func, init_sigma);
    sol_mu = -sol_sigma * log(int3D(exp(-vol / sol_sigma) .* heavi_phi) / ...
        int3D(heavi_phi));
    %}
    sol_mu = init_mu;
    sol_sigma = init_sigma;
    %inter_term = dirac_phi .* (exp((-vol + sol_mu) / sol_mu) + log(sol_sigma) + (vol - sol_mu) / sol_sigma);

    inter_term = exp((-vol + sol_mu) / sol_mu) + log(sol_sigma) + (vol - sol_mu) / sol_sigma;

    
    
    phi = phi + timestep * (mu * dist_regu + lambda * edge_term - alpha * area_term - beta * inter_term + ...
    gamma * ext_term);
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


%{
function f = SmoothDirac(phi, epsilon)
f = (epsilon / phi) ./ (epsilon^2. + phi.^2);
end

function f = SmoothHeavi(phi, epsilon)
f = 0.5 * (1 + (2/pi) * atan(phi / epsilon));
end
%}

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

function [mu, sigma] = GumbelInit(vol, heavi_phi)
euler_const = 0.57721;
mean = int3D(vol .* heavi_phi) / int3D(heavi_phi);
var = int3D((vol - mean).^2 .* heavi_phi) / int3D(heavi_phi);
sigma = sqrt(6 * var) / pi;
mu = mean - sigma * euler_const;
end
