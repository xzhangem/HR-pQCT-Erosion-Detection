function phi = LSThreshold(phi_0, vol, g, mu, lambda, alpha, beta, gamma, epsilon, timestep, iter)

% Input: 
%   phi_0: inital level set function in this stage
%   vol: double volume for updating mean for external energy
%   g: 3D surface indicator function
%   mu: weight of the distance regularization term
%   lambda: weight of the weighted length term 
%   alpha: weight of the weighted area term
%   beta: weight of the weighted internal Gaussian
%   gamma: weight of the weighted external Gaussian
%   epsilon: width of Dirac Delta function and Heaviside function
%   timestep: time step
%   iter: number of iterations

%Output:
%   phi: the updated level set function

phi = phi_0;
[vx, vy, vz] = gradient(g);
for k = 1 : iter
    phi = NeumannBoundCond(phi);
    [phi_x, phi_y, phi_z] = gradient(phi);
    s = sqrt(phi_x.^2 + phi_y.^2 + phi_z.^2);
    small_number = 1e-9;
    Nx = phi_x ./ (s + small_number);
    Ny = phi_y ./ (s + small_number);
    Nz = phi_z ./ (s + small_number);
    curvature = Div(Nx, Ny, Nz);
    dist_regu = DistReg(phi);
    dirac_phi = SmoothDirac(phi, epsilon);
    heavi_phi = SmoothHeavi(phi, epsilon);
    area_term = dirac_phi .* g;
    %edge_term = dirac_phi .* (vx .* Nx + vy .* Ny + vz .* Nz) + dirac_phi .* g .* curvature;
    %edge_term = curvature;
    edge_term = dirac_phi .* curvature;
    %area_term = dirac_phi;
    c_in = (int3D(vol .* heavi_phi)) / (int3D(heavi_phi));
    c_ex = (int3D(vol .* (1 - heavi_phi))) / (int3D(1 - heavi_phi));
    %sigma_in = (int3D((vol - c_in).^2 .* heavi_phi)) / (int3D(heavi_phi));
    %sigma_ex = (int3D((vol - c_ex).^2 .* (1 - heavi_phi))) / (int3D(1 - heavi_phi));
    inter_term = (vol - c_in).^2;
    ext_term = (vol - c_ex).^2;
    %sigma_term = log(sigma_in / sigma_ex);
    %bound_term = BoundConstrain(phi);
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

function f = BoundConstrain(phi)
a = (phi < 0);
b = (phi > 1);
c = (phi >= 0) & (phi <= 1);
f = -2. * a + 2. * b + 0. * c;  
end
