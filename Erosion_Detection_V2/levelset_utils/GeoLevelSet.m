function phi = GeoLevelSet(phi_0, g, mu, lambda, alpha, epsilon, timestep, iter)

% Input: 
%   phi_0: inital level set function in this stage
%   g: 3D surface indicator function
%   mu: weight of the distance regularization term
%   lambda: weight of the weighted length term 
%   alpha: weight of the weighted area term
%   epsilon: width of Dirac Delta function and Heaviside function
%   timestep: time step
%   iter: number of iterations

%Output:
%   phi: the updated level set function

region = 0;
phi = phi_0;
unit_vol = ones(size(phi));
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
    %heavi_phi = SmoothHeavi(phi, epsilon);
    area_term = dirac_phi .* g;
    edge_term = dirac_phi .* (vx .* Nx + vy .* Ny + vz .* Nz) + dirac_phi .* g .* curvature;

    phi = phi + timestep * (mu * dist_regu + lambda * edge_term - alpha * area_term);
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
%{
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
%}

%%{
function f = SmoothDirac(phi, epsilon)
f = (epsilon / pi) ./ (epsilon^2. + phi.^2);
end

function f = SmoothHeavi(phi, epsilon)
f = 0.5 * (1 + (2/pi) * atan(phi / epsilon));
end
%%}

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
